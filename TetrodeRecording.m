classdef TetrodeRecording < handle
	properties
		Files = []
		Path = []
		Part = [1, 1]
		Notes
		FrequencyParameters
		ChannelMap
		SelectedChannels
		Spikes
		DigitalEvents
	end

	properties (Transient)
		Amplifier
		BoardDigIn
	end

	%----------------------------------------------------
	%		Methods
	%----------------------------------------------------
	methods
		function obj = TetrodeRecording()
		end

		function Preview(obj, varargin)
			p = inputParser;
			addParameter(p, 'Channels', [], @isnumeric);
			addParameter(p, 'ChunkSize', 10, @isnumeric);
			addParameter(p, 'HideResults', false, @islogical);
			parse(p, varargin{:});
			channels 		= p.Results.Channels;
			chunkSize 		= p.Results.ChunkSize;
			hideResults 	= p.Results.HideResults;

			obj.ReadFiles(chunkSize);
			if isempty(channels)
				channels = obj.MapChannelID([obj.Amplifier.Channels.NativeOrder]);
			end
			obj.SpikeDetect(channels, 'NumSigmas', 4, 'WaveformWindow', [-0.5, 0.5]);
			obj.SpikeSort(channels, 'ClusterMethod', 'kmeans', 'FeatureMethod', 'PCA', 'Dimension', 3);
			if ~hideResults
				obj.PlotAllChannels();
			end
		end

		function expName = GetExpName(obj)
			expName = strsplit(obj.Path, '\');
			expName = expName{end - 1};
		end

		function SelectFiles(obj)
			[obj.Files, obj.Path, ~] = uigetfile('*.rhd', 'Select an RHD2000 Data File', 'MultiSelect', 'on');
		end

		function ReadFiles(obj, varargin)
			p = inputParser;
			addOptional(p, 'ChunkSize', 2, @isnumeric);
			addParameter(p, 'Chunks', 'first');
			addParameter(p, 'SubstractMean', true, @islogical);
			addParameter(p, 'RemoveTransient', true, @islogical);
			addParameter(p, 'SpikeDetect', false, @islogical);
			addParameter(p, 'DigitalDetect', false, @islogical);
			parse(p, varargin{:});
			chunkSize = p.Results.ChunkSize;
			chunks = p.Results.Chunks;
			substractMean = p.Results.SubstractMean;
			removeTransient = p.Results.RemoveTransient;
			spikeDetect = p.Results.SpikeDetect;
			digitalDetect = p.Results.DigitalDetect;

			files = obj.Files;
			numChunks = ceil(length(files)/chunkSize);
			if ischar(chunks)
				switch lower(chunks)
					case 'all'
						chunks = 1:numChunks;
					case 'remaining'
						chunks = 2:numChunks;
					case 'first'
						chunks = 1;
				end
			end

			for iChunk = chunks
				obj.ClearCache();
				TetrodeRecording.TTS(['Processing chunk ', num2str(iChunk), '/', num2str(numChunks), ':\n']);
				obj.ReadRHD(obj.Files((iChunk - 1)*chunkSize + 1:min(iChunk*chunkSize, length(obj.Files))))
				obj.MapChannels();
                if substractMean
                    obj.SubstractMean();
                end
                if removeTransient
					obj.RemoveTransient();
				end
				if spikeDetect
					obj.SpikeDetect([obj.Spikes.Channel], 'Append', true);
					obj.GetDigitalData('Append', true);
				else
					obj.GetDigitalData('Append', false);
				end
			end
			if digitalDetect
				obj.GetDigitalEvents(true);
			end
		end

		function ReadRHD(obj, files)
			TetrodeRecording.TTS(['	Loading data:\n'])
			for iFile = 1:length(files)
				filename = [obj.Path, files{iFile}];

				if exist(filename, 'file') ~= 2
					error(['File ', filename, ' not found.']);
				end

				fid = fopen(filename, 'r');
				s = dir(filename);
				filesize = s.bytes;

				% Check 'magic number' at beginning of file to make sure this is an Intan Technologies RHD2000 data file.
				magic_number = fread(fid, 1, 'uint32');
				if magic_number ~= hex2dec('c6912702')
					error('Unrecognized file type.');
				end

				% Read version number.
				data_file_main_version_number = fread(fid, 1, 'int16');
				data_file_secondary_version_number = fread(fid, 1, 'int16');

				if (data_file_main_version_number == 1)
					num_samples_per_data_block = 60;
				else
					num_samples_per_data_block = 128;
				end

				% Read information of sampling rate and amplifier frequency settings.
				sample_rate = fread(fid, 1, 'single');
				dsp_enabled = fread(fid, 1, 'int16');
				actual_dsp_cutoff_frequency = fread(fid, 1, 'single');
				actual_lower_bandwidth = fread(fid, 1, 'single');
				actual_upper_bandwidth = fread(fid, 1, 'single');

				desired_dsp_cutoff_frequency = fread(fid, 1, 'single');
				desired_lower_bandwidth = fread(fid, 1, 'single');
				desired_upper_bandwidth = fread(fid, 1, 'single');

				% This tells us if a software 50/60 Hz notch filter was enabled during
				% the data acquisition.
				notch_filter_mode = fread(fid, 1, 'int16');
				notch_filter_frequency = 0;
				if (notch_filter_mode == 1)
					notch_filter_frequency = 50;
				elseif (notch_filter_mode == 2)
					notch_filter_frequency = 60;
				end

				desired_impedance_test_frequency = fread(fid, 1, 'single');
				actual_impedance_test_frequency = fread(fid, 1, 'single');

				% Place notes in data strucure
				notes = struct( ...
					'note1', TetrodeRecording.ReadQString(fid), ...
					'note2', TetrodeRecording.ReadQString(fid), ...
					'note3', TetrodeRecording.ReadQString(fid) );
			
				% If data file is from GUI v1.1 or later, see if temperature sensor data
				% was saved.
				num_temp_sensor_channels = 0;
				if ((data_file_main_version_number == 1 && data_file_secondary_version_number >= 1) ...
					|| (data_file_main_version_number > 1))
					num_temp_sensor_channels = fread(fid, 1, 'int16');
				end

				% If data file is from GUI v1.3 or later, load eval board mode.
				eval_board_mode = 0;
				if ((data_file_main_version_number == 1 && data_file_secondary_version_number >= 3) ...
					|| (data_file_main_version_number > 1))
					eval_board_mode = fread(fid, 1, 'int16');
				end

				% If data file is from v2.0 or later (Intan Recording Controller),
				% load name of digital reference channel.
				if (data_file_main_version_number > 1)
					objTemp(iFile).ReferenceChannel = TetrodeRecording.ReadQString(fid);
				end

				% Place frequency-related information in data structure.
				frequency_parameters = struct( ...
					'AmplifierSampleRate', sample_rate, ...
					'AuxInputDampleRate', sample_rate / 4, ...
					'SupplyVoltageSampleRate', sample_rate / num_samples_per_data_block, ...
					'BoardADCSampleRate', sample_rate, ...
					'BoardDigInSampleRate', sample_rate, ...
					'DesiredDSPCutoffFrequency', desired_dsp_cutoff_frequency, ...
					'ActualDspCutoffFrequency', actual_dsp_cutoff_frequency, ...
					'DSPEnabled', dsp_enabled, ...
					'DesiredLowerBandwidth', desired_lower_bandwidth, ...
					'ActualLowerBandwidth', actual_lower_bandwidth, ...
					'DesiredUpperBandwidth', desired_upper_bandwidth, ...
					'ActualUpperBandwidth', actual_upper_bandwidth, ...
					'NotchFilterFrequency', notch_filter_frequency, ...
					'DesiredImpedanceTestFrequency', desired_impedance_test_frequency, ...
					'ActualImpedanceTestFrequency', actual_impedance_test_frequency );

				% Define data structure for spike trigger settings.
				spike_trigger_struct = struct( ...
					'VoltageTriggerMode', {}, ...
					'VoltageThreshold', {}, ...
					'DigitalTriggerChannel', {}, ...
					'DigitalEdgePolarity', {} );

				% Define data structure for data channels.
				channel_struct = struct( ...
					'NativeChannelName', {}, ...
					'CustomChannelName', {}, ...
					'NativeOrder', {}, ...
					'CustomOrder', {}, ...
					'BoardStream', {}, ...
					'ChipChannel', {}, ...
					'PortName', {}, ...
					'PortPrefix', {}, ...
					'PortNumber', {}, ...
					'ElectrodeImpedanceMagnitude', {}, ...
					'ElectrodeImpedancePhase', {} );

				new_channel = struct(channel_struct);

				% Create structure arrays for each type of data channel.
				objTemp(iFile).Amplifier.Channels = struct(channel_struct);
				objTemp(iFile).BoardDigIn.Channels = struct(channel_struct);

				amplifier_index = 1;
				aux_input_index = 1;
				supply_voltage_index = 1;
				board_adc_index = 1;
				board_dig_in_index = 1;
				board_dig_out_index = 1;

				% Read signal summary from data file header.
				number_of_signal_groups = fread(fid, 1, 'int16');

				for signal_group = 1:number_of_signal_groups
					signal_group_name = TetrodeRecording.ReadQString(fid);
					signal_group_prefix = TetrodeRecording.ReadQString(fid);
					signal_group_enabled = fread(fid, 1, 'int16');
					signal_group_num_channels = fread(fid, 1, 'int16');
					signal_group_num_amp_channels = fread(fid, 1, 'int16');

					if (signal_group_num_channels > 0 && signal_group_enabled > 0)
						new_channel(1).PortName = signal_group_name;
						new_channel(1).PortPrefix = signal_group_prefix;
						new_channel(1).PortNumber = signal_group;
						for signal_channel = 1:signal_group_num_channels
							new_channel(1).NativeChannelName = TetrodeRecording.ReadQString(fid);
							new_channel(1).CustomChannelName = TetrodeRecording.ReadQString(fid);
							new_channel(1).NativeOrder = fread(fid, 1, 'int16');
							new_channel(1).CustomOrder = fread(fid, 1, 'int16');
							signal_type = fread(fid, 1, 'int16');
							channel_enabled = fread(fid, 1, 'int16');
							new_channel(1).ChipChannel = fread(fid, 1, 'int16');
							new_channel(1).BoardStream = fread(fid, 1, 'int16');
							fread(fid, 4, 'int16');
							% new_trigger_channel(1).VoltageTriggerMode = fread(fid, 1, 'int16');
							% new_trigger_channel(1).VoltageThreshold = fread(fid, 1, 'int16');
							% new_trigger_channel(1).DigitalTriggerChannel = fread(fid, 1, 'int16');
							% new_trigger_channel(1).DigitalEdgePolarity = fread(fid, 1, 'int16');
							new_channel(1).ElectrodeImpedanceMagnitude = fread(fid, 1, 'single');
							new_channel(1).ElectrodeImpedancePhase = fread(fid, 1, 'single');
							
							if (channel_enabled)
								switch (signal_type)
									case 0
										objTemp(iFile).Amplifier.Channels(amplifier_index) = new_channel;
										amplifier_index = amplifier_index + 1;
									case 1
										aux_input_index = aux_input_index + 1;
									case 2
										supply_voltage_index = supply_voltage_index + 1;
									case 3
										board_adc_index = board_adc_index + 1;
									case 4
										objTemp(iFile).BoardDigIn.Channels(board_dig_in_index) = new_channel;
										board_dig_in_index = board_dig_in_index + 1;
									case 5
										board_dig_out_index = board_dig_out_index + 1;
									otherwise
										error('Unknown channel type');
								end
							end
							
						end
					end
				end

				% Summarize contents of data file.
				num_amplifier_channels = amplifier_index - 1;
				num_aux_input_channels = aux_input_index - 1;
				num_supply_voltage_channels = supply_voltage_index - 1;
				num_board_adc_channels = board_adc_index - 1;
				num_board_dig_in_channels = board_dig_in_index - 1;
				num_board_dig_out_channels = board_dig_out_index - 1;

				% Determine how many samples the data file contains.

				% Each data block contains num_samples_per_data_block amplifier samples.
				bytes_per_block = num_samples_per_data_block * 4;  % timestamp data
				bytes_per_block = bytes_per_block + num_samples_per_data_block * 2 * num_amplifier_channels;
				% Auxiliary inputs are sampled 4x slower than amplifiers
				bytes_per_block = bytes_per_block + (num_samples_per_data_block / 4) * 2 * num_aux_input_channels;
				% Supply voltage is sampled once per data block
				bytes_per_block = bytes_per_block + 1 * 2 * num_supply_voltage_channels;
				% Board analog inputs are sampled at same rate as amplifiers
				bytes_per_block = bytes_per_block + num_samples_per_data_block * 2 * num_board_adc_channels;
				% Board digital inputs are sampled at same rate as amplifiers
				if (num_board_dig_in_channels > 0)
					bytes_per_block = bytes_per_block + num_samples_per_data_block * 2;
				end
				% Board digital outputs are sampled at same rate as amplifiers
				if (num_board_dig_out_channels > 0)
					bytes_per_block = bytes_per_block + num_samples_per_data_block * 2;
				end
				% Temp sensor is sampled once per data block
				if (num_temp_sensor_channels > 0)
					bytes_per_block = bytes_per_block + 1 * 2 * num_temp_sensor_channels; 
				end

				% How many data blocks remain in this file?
				data_present = 0;
				bytes_remaining = filesize - ftell(fid);
				if (bytes_remaining > 0)
					data_present = 1;
				end

				num_data_blocks = bytes_remaining / bytes_per_block;

				num_amplifier_samples = num_samples_per_data_block * num_data_blocks;
				num_aux_input_samples = (num_samples_per_data_block / 4) * num_data_blocks;
				num_supply_voltage_samples = 1 * num_data_blocks;
				num_board_adc_samples = num_samples_per_data_block * num_data_blocks;
				num_board_dig_in_samples = num_samples_per_data_block * num_data_blocks;
				num_board_dig_out_samples = num_samples_per_data_block * num_data_blocks;

				if (data_present)
					% Pre-allocate memory for data.
					objTemp(iFile).Amplifier.Timestamps = zeros(1, num_amplifier_samples);
					objTemp(iFile).Amplifier.Data = zeros(num_amplifier_channels, num_amplifier_samples);
					objTemp(iFile).BoardDigIn.Data = zeros(num_board_dig_in_channels, num_board_dig_in_samples);
					board_dig_in_raw = zeros(1, num_board_dig_in_samples);
					board_dig_out_raw = zeros(1, num_board_dig_out_samples);

					% Read sampled data from file.
					tic, TetrodeRecording.TTS(['		Reading file ', num2str(iFile), '/', num2str(length(files)), ' (''', files{iFile},''')...']);

					amplifier_index = 1;
					aux_input_index = 1;
					supply_voltage_index = 1;
					board_adc_index = 1;
					board_dig_in_index = 1;
					board_dig_out_index = 1;

					print_increment = 10;
					percent_done = print_increment;
					for i=1:num_data_blocks
						% In version 1.2, we moved from saving timestamps as unsigned
						% integeters to signed integers to accomidate negative (adjusted)
						% timestamps for pretrigger data.
						if ((data_file_main_version_number == 1 && data_file_secondary_version_number >= 2) ...
						|| (data_file_main_version_number > 1))
							objTemp(iFile).Amplifier.Timestamps(amplifier_index:(amplifier_index + num_samples_per_data_block - 1)) = fread(fid, num_samples_per_data_block, 'int32');
						else
							objTemp(iFile).Amplifier.Timestamps(amplifier_index:(amplifier_index + num_samples_per_data_block - 1)) = fread(fid, num_samples_per_data_block, 'uint32');
						end
						if (num_amplifier_channels > 0)
							objTemp(iFile).Amplifier.Data(:, amplifier_index:(amplifier_index + num_samples_per_data_block - 1)) = fread(fid, [num_samples_per_data_block, num_amplifier_channels], 'uint16')';
						end
						if (num_aux_input_channels > 0)
							fread(fid, [(num_samples_per_data_block / 4), num_aux_input_channels], 'uint16')';
						end
						if (num_supply_voltage_channels > 0)
							fread(fid, [1, num_supply_voltage_channels], 'uint16')';
						end
						if (num_temp_sensor_channels > 0)
							fread(fid, [1, num_temp_sensor_channels], 'int16')';
						end
						if (num_board_adc_channels > 0)
							fread(fid, [num_samples_per_data_block, num_board_adc_channels], 'uint16')';
						end
						if (num_board_dig_in_channels > 0)
							board_dig_in_raw(board_dig_in_index:(board_dig_in_index + num_samples_per_data_block - 1)) = fread(fid, num_samples_per_data_block, 'uint16');
						end
						if (num_board_dig_out_channels > 0)
							% board_dig_out_raw(board_dig_out_index:(board_dig_out_index + num_samples_per_data_block - 1)) = fread(fid, num_samples_per_data_block, 'uint16');
							fread(fid, num_samples_per_data_block, 'uint16');
						end

						amplifier_index = amplifier_index + num_samples_per_data_block;
						aux_input_index = aux_input_index + (num_samples_per_data_block / 4);
						supply_voltage_index = supply_voltage_index + 1;
						board_adc_index = board_adc_index + num_samples_per_data_block;
						board_dig_in_index = board_dig_in_index + num_samples_per_data_block;
						board_dig_out_index = board_dig_out_index + num_samples_per_data_block;

						fraction_done = 100 * (i / num_data_blocks);
						if (fraction_done >= percent_done)
							% TetrodeRecording.TTS('\t%d%% done...\n', percent_done);
							percent_done = percent_done + print_increment;
						end
					end

					% Make sure we have read exactly the right amount of data.
					bytes_remaining = filesize - ftell(fid);
					if (bytes_remaining ~= 0)
						error('Error: End of file not reached.');
					end
				end

				% Close data file.
				fclose(fid);

				if (data_present)
					% Extract digital input channels to separate variables.
					for i=1:num_board_dig_in_channels
						mask = 2^(objTemp(iFile).BoardDigIn.Channels(i).NativeOrder) * ones(size(board_dig_in_raw));
						objTemp(iFile).BoardDigIn.Data(i, :) = (bitand(board_dig_in_raw, mask) > 0);
					end

					% Scale voltage levels appropriately.
					objTemp(iFile).Amplifier.Data = 0.195 * (objTemp(iFile).Amplifier.Data - 32768); % units = microvolts

					% Scale time steps (units = seconds).
					objTemp(iFile).Amplifier.Timestamps = objTemp(iFile).Amplifier.Timestamps / sample_rate;
					objTemp(iFile).BoardDigIn.Timestamps = objTemp(iFile).Amplifier.Timestamps;
				end
				TetrodeRecording.TTS(['Done(', num2str(toc, '%.2f'), ' seconds).\n'])
			end

			% Count number of samples in all files
			obj.Amplifier.NumSamples = 0;
			obj.BoardDigIn.NumSamples = 0;
			for iFile = 1:length(files)
				obj.Amplifier.NumSamples = obj.Amplifier.NumSamples + size(objTemp(iFile).Amplifier.Timestamps, 2);
				obj.BoardDigIn.NumSamples = obj.BoardDigIn.NumSamples + size(objTemp(iFile).BoardDigIn.Timestamps, 2);
			end

			% Combine files
			tic, TetrodeRecording.TTS('	Concatenating data...');

			obj.Notes = notes;
			obj.FrequencyParameters = frequency_parameters;
			obj.Amplifier.Channels = objTemp(1).Amplifier.Channels;
			obj.BoardDigIn.Channels = objTemp(1).BoardDigIn.Channels;

			obj.Amplifier.Timestamps = zeros(1, obj.Amplifier.NumSamples);
			obj.Amplifier.Data = zeros(length(obj.Amplifier.Channels), obj.Amplifier.NumSamples);
			obj.BoardDigIn.Timestamps = zeros(1, obj.BoardDigIn.NumSamples);
			obj.BoardDigIn.Data = zeros(length(obj.BoardDigIn.Channels), obj.BoardDigIn.NumSamples);

			iSample.Amplifier = 0;
			iSample.BoardDigIn = 0;
			for iFile = 1:length(files)
				obj.Amplifier.Timestamps(1, iSample.Amplifier + 1:iSample.Amplifier + size(objTemp(iFile).Amplifier.Timestamps, 2)) = objTemp(iFile).Amplifier.Timestamps;
				obj.Amplifier.Data(:, iSample.Amplifier + 1:iSample.Amplifier + size(objTemp(iFile).Amplifier.Timestamps, 2)) = objTemp(iFile).Amplifier.Data;
				obj.BoardDigIn.Timestamps(1, iSample.BoardDigIn + 1:iSample.BoardDigIn + size(objTemp(iFile).BoardDigIn.Timestamps, 2)) = objTemp(iFile).BoardDigIn.Timestamps;
				obj.BoardDigIn.Data(:, iSample.BoardDigIn + 1:iSample.BoardDigIn + size(objTemp(iFile).BoardDigIn.Timestamps, 2)) = objTemp(iFile).BoardDigIn.Data;

				iSample.Amplifier = iSample.Amplifier + size(objTemp(iFile).Amplifier.Timestamps, 2);
				iSample.BoardDigIn = iSample.BoardDigIn + size(objTemp(iFile).BoardDigIn.Timestamps, 2);
			end
			TetrodeRecording.TTS(['Done(', num2str(toc, '%.2f'), ' seconds).\n'])
		end

		% Used to preview a small portion of loaded data. Will remove used data from workspace.
		function TrimData(obj, numSamples)
			obj.BoardDigIn.NumSamples = numSamples;
			obj.BoardDigIn.Timestamps = obj.BoardDigIn.Timestamps(1:numSamples);
			obj.BoardDigIn.Data = obj.BoardDigIn.Data(:, 1:numSamples);
			obj.Amplifier.NumSamples = numSamples;
			obj.Amplifier.Timestamps = obj.Amplifier.Timestamps(1:numSamples);
			obj.Amplifier.Data = obj.Amplifier.Data(:, 1:numSamples);
		end

		% Amplifier data arranged in tetrode order (i.e., first four elements in array is first tetrode)
		function MapChannels(obj, EIBMap, headstageType)
			tic, TetrodeRecording.TTS('	Remapping amplifier channels in tetrode order...');

			if nargin < 3
				headstageType = 'intan';
			end
			if nargin < 2
				EIBMap = [15, 13, 11, 9, 7, 5, 3, 1, 2, 4, 6, 8, 10, 12, 14, 16];
				EIBMap = [EIBMap, EIBMap + 16];
				[~, EIBMap] = sort(EIBMap);
			end

			switch headstageType
				case 'intan'
					HSMap = [25:32, 1:8, 24:-1:9];
				case 'open ephys'
					HSMap = [26, 25, 27, 28, 29, 30, 31, 1, 32, 3, 2, 5, 4, 7, 6, 8, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9];
				otherwise
					error('Unrecognized headstage type.'); 
			end

			tetrodeMap = HSMap(EIBMap);
			recordedChannels = [obj.Amplifier.Channels.NativeOrder] + 1;
			
			obj.Amplifier.DataMapped = NaN(32, size(obj.Amplifier.Data, 2));
			iChn = 0;
			for targetChannel = tetrodeMap
				iChn = iChn + 1;
				sourceChannel = find(recordedChannels == targetChannel, 1);
				if ~isempty(sourceChannel)
					obj.Amplifier.DataMapped(iChn, :) = obj.Amplifier.Data(sourceChannel, :);
				end
			end
			obj.ChannelMap.ElectrodeInterfaceBoard = EIBMap;
			obj.ChannelMap.Headstage = HSMap;
			obj.ChannelMap.Tetrode = tetrodeMap;
			obj.Amplifier.Data = obj.Amplifier.DataMapped; 
			obj.Amplifier = rmfield(obj.Amplifier, 'DataMapped');
			TetrodeRecording.TTS(['Done(', num2str(toc, '%.2f'), ' seconds).\n'])
		end

		% Convert zero-based raw channel id to remapped tetrode id
		function tetrodeChannel = MapChannelID(obj, rawChannel)
			tetrodeChannel = [];
			for iChannel = rawChannel
				tetrodeChannel = [tetrodeChannel, find(obj.ChannelMap.Tetrode == (iChannel + 1))];
			end
		end

		% Substract by 32 chn mean
		function SubstractMean(obj)
			dataMean = nanmean(obj.Amplifier.Data, 1);
			obj.Amplifier.Data = obj.Amplifier.Data - dataMean;
		end

		% Remove lick/press-transient by setting signal to zero
		function RemoveTransient(obj, digChannel, dilate, dilateSize)
			if nargin < 2
				digChannel = [1, 2]; % Default lick channel is 2 on intan board. This was used for Daisy1. Also do this for lever press (Chn 1)
			end
			if nargin < 3
				dilate = true;
			end
			if nargin < 4
				dilateSize = round(30*obj.FrequencyParameters.AmplifierSampleRate/1000); % By default, lick digital event is extended 15 ms to the left and right.
			end

			% Verify if digital input sample rate is identical to amplifier sample rate
			if obj.FrequencyParameters.AmplifierSampleRate ~= obj.FrequencyParameters.BoardDigInSampleRate
				error('Board digital input has a different sampling rate from amplifier. Aborted.');
			end

			tic, TetrodeRecording.TTS('	Removing lick/touch-related transients. This might take a while...');
			if dilate
				mask = false(1, obj.Amplifier.NumSamples);
				for thisChannel = digChannel
					thisChannel = find([obj.BoardDigIn.Channels.NativeOrder] == thisChannel);
					mask = mask | logical(movmean(obj.BoardDigIn.Data(thisChannel, :), dilateSize));
				end
			else
				for thisChannel = digChannel
					thisChannel = find([obj.BoardDigIn.Channels.NativeOrder] == thisChannel);
					mask = logical(obj.BoardDigIn.Data(thisChannel, :));
				end
			end

			obj.Amplifier.Data(:, mask) = 0;
			TetrodeRecording.TTS(['Done(', num2str(toc, '%.2f'), ' seconds).\n'])
		end

		% Expand waveform window, fill unavailable data with NaN
		function varargout = GetWaveforms(obj, channel, waveformWindow, varargin)
			p = inputParser;
			addRequired(p, 'Channel', @isnumeric);
			addRequired(p, 'WaveformWindow', @(x) isnumeric(x) && length(x) == 2);
			addOptional(p, 'Index', NaN, @isnumeric);
			addParameter(p, 'IndexType', 'SampleIndex', @ischar);
			parse(p, channel, waveformWindow, varargin{:});
			channel = p.Results.Channel;
			waveformWindow = p.Results.WaveformWindow;
			index = p.Results.Index;
			indexType = p.Results.IndexType;

			if length(index) == 1 && isnan(index)
				index = obj.Spikes(channel).SampleIndex;
				indexType = 'SampleIndex';
			end

			switch indexType
				case 'SampleIndex'
					sampleIndex = index;
					% Only interpolate if sample index in non-integer
					if sum(rem(sampleIndex, 1) == 0) == length(sampleIndex)
						timestamps = obj.Amplifier.Timestamps(sampleIndex);
					else
						timestamps = interp1(1:length(obj.Amplifier.Timestamps), obj.Amplifier.Timestamps, sampleIndex, 'linear');
					end
				case 'Timestamps'
					% Always interpolate if input index in timestamps. This will always be slower.
					timestamps = index;
					sampleIndex = interp1(obj.Amplifier.Timestamps, 1:length(obj.Amplifier.Timestamps), timestamps, 'linear');
				otherwise
					error(['Unrecognized index type: ''', indexType, ''', must be ''SampleIndex'' or ''Timestamps''.'])
			end

			sampleRate = obj.FrequencyParameters.AmplifierSampleRate/1000;
			waveforms = NaN(length(sampleIndex), 1 + round(sampleRate*(waveformWindow(2) - waveformWindow(1))));
			t = waveformWindow(1):1/sampleRate:waveformWindow(2);
			i = sampleRate*t;
			for iWaveform = 1:length(sampleIndex)
				iQuery = sampleIndex(iWaveform) + i;
				if (sum(rem(iQuery, 1) == 0) == length(iQuery)) && min(iQuery) > 0 && max(iQuery) <= size(obj.Amplifier.Data, 2)
					waveforms(iWaveform, :) = obj.Amplifier.Data(channel, iQuery);
				else
					waveforms(iWaveform, :) = interp1(1:size(obj.Amplifier.Data, 2), obj.Amplifier.Data(channel, :), iQuery, 'pchip', NaN);
				end
			end

			% Output
			varargout = {waveforms, t, timestamps, sampleIndex};
		end

		% Detect spike-like waveforms by simple (or advanced) thresholding
		function SpikeDetect(obj, channels, varargin)
			p = inputParser;
			addRequired(p, 'Channels', @isnumeric);
			addParameter(p, 'NumSigmas', 4, @isnumeric);
			addParameter(p, 'Direction', 'negative', @ischar);
			addParameter(p, 'WaveformWindow', [-1.25, 1.25], @isnumeric);
			addParameter(p, 'Append', false, @islogical);
			parse(p, channels, varargin{:});
			channels = p.Results.Channels;
			append = p.Results.Append;

			TetrodeRecording.TTS('	Detecting spikes:\n');
			sampleRate = obj.FrequencyParameters.AmplifierSampleRate/1000;

			for iChannel = channels
				% If append mode
				if ~append
					numSigmas = p.Results.NumSigmas;
					directionMode = p.Results.Direction;
					waveformWindow = p.Results.WaveformWindow;
				else
					numSigmas = obj.Spikes(iChannel).Threshold.NumSigmas;
					directionMode = obj.Spikes(iChannel).Threshold.Direction;
					waveformWindow = obj.Spikes(iChannel).WaveformWindow;
				end
				% Auto-threshold for spikes
				% median(abs(x))/0.6745 is a better estimation of noise amplitude than std() when there are spikes -- Quiroga, 2004
				threshold = numSigmas*nanmedian(abs(obj.Amplifier.Data(iChannel, :)))/0.6745;
				switch lower(directionMode)
					case 'negative'
						direction = -1;
					case 'positive'
						direction = 1;
					case 'auto'
						direction = sign(median(obj.Amplifier.Data(iChannel, abs(obj.Amplifier.Data(iChannel, :)) > 1.5*threshold))); % Check if spikes are positive or negative
					otherwise
						error(['Unrecognized spike detection mode ''', directionMode, '''.'])
				end
				tic, TetrodeRecording.TTS(['		Channel ', num2str(iChannel), ' (', num2str(char(952)), ' = ', num2str(numSigmas), num2str(char(963)), ' = ', num2str(direction*threshold), ')...']);
				
				% Find spikes
				[~, sampleIndex] = findpeaks(direction*obj.Amplifier.Data(iChannel, :), 'MinPeakHeight', threshold, 'MinPeakProminence', threshold);

				% Extract waveforms
				[waveforms, t] = obj.GetWaveforms(iChannel, waveformWindow, sampleIndex, 'IndexType', 'SampleIndex');

				% Align waveforms to peak
				numWaveforms = size(waveforms, 1);
				i = sampleRate*t;
				[~, maxIndex] = max(direction*waveforms, [], 2);
				alignmentShift = i(maxIndex);
				[waveforms, t, timestamps, sampleIndex] = obj.GetWaveforms(iChannel, waveformWindow, sampleIndex + alignmentShift, 'IndexType', 'SampleIndex');

				% Store data
				if ~append
					obj.Spikes(iChannel).Channel = iChannel;

					obj.Spikes(iChannel).SampleIndex = sampleIndex;
					obj.Spikes(iChannel).Timestamps = timestamps;
					obj.Spikes(iChannel).Waveforms = waveforms;

					obj.Spikes(iChannel).WaveformTimestamps = t;
					obj.Spikes(iChannel).WaveformWindow = waveformWindow;
					obj.Spikes(iChannel).Threshold.NumSigmas = numSigmas;
					obj.Spikes(iChannel).Threshold.Threshold = direction*threshold;
					obj.Spikes(iChannel).Threshold.Direction = directionMode;
				else
					obj.Spikes(iChannel).SampleIndex = [obj.Spikes(iChannel).SampleIndex, sampleIndex + length(obj.DigitalEvents.Timestamps)];
					obj.Spikes(iChannel).Timestamps = [obj.Spikes(iChannel).Timestamps, timestamps];
					obj.Spikes(iChannel).Waveforms = [obj.Spikes(iChannel).Waveforms; waveforms];
					obj.Spikes(iChannel).Threshold.Threshold = [obj.Spikes(iChannel).Threshold.Threshold, direction*threshold];
				end

				TetrodeRecording.TTS(['Done(', num2str(numWaveforms), ' waveforms, ', num2str(toc, '%.2f'), ' seconds).\n'])
			end
		end

		function GetDigitalData(obj, varargin)
			p = inputParser;
			addParameter(p, 'ChannelCue', 4, @isnumeric);
			addParameter(p, 'ChannelPress', 2, @isnumeric);
			addParameter(p, 'ChannelLick', 3, @isnumeric);
			addParameter(p, 'ChannelReward', 5, @isnumeric);
			addParameter(p, 'Append', false, @islogical);
			parse(p, varargin{:});
			channelCue 		= p.Results.ChannelCue;
			channelPress 	= p.Results.ChannelPress;
			channelLick 	= p.Results.ChannelLick;
			channelReward 	= p.Results.ChannelReward;
			append 			= p.Results.Append;

			if append
				obj.DigitalEvents.Data = [obj.DigitalEvents.Data, obj.BoardDigIn.Data([channelCue, channelPress, channelLick, channelReward], :)];
				obj.DigitalEvents.Timestamps = [obj.DigitalEvents.Timestamps, obj.BoardDigIn.Timestamps];
			else
				obj.DigitalEvents.Data = obj.BoardDigIn.Data([channelCue, channelPress, channelLick, channelReward], :);
				obj.DigitalEvents.Timestamps = obj.BoardDigIn.Timestamps;
			end
		end

		function GetDigitalEvents(obj, clearCache)
			if nargin < 2
				clearCache = false;
			end

			tic, TetrodeRecording.TTS('	Extracting digital events...');
			% Get digital signals
			[obj.DigitalEvents.CueOn, obj.DigitalEvents.CueOff] 		= TetrodeRecording.FindEdges(obj.DigitalEvents.Data(1, :), obj.DigitalEvents.Timestamps);
			[obj.DigitalEvents.PressOn, obj.DigitalEvents.PressOff] 	= TetrodeRecording.FindEdges(obj.DigitalEvents.Data(2, :), obj.DigitalEvents.Timestamps);
			[obj.DigitalEvents.LickOn, obj.DigitalEvents.LickOff] 		= TetrodeRecording.FindEdges(obj.DigitalEvents.Data(3, :), obj.DigitalEvents.Timestamps);
			[obj.DigitalEvents.RewardOn, obj.DigitalEvents.RewardOff] 	= TetrodeRecording.FindEdges(obj.DigitalEvents.Data(4, :), obj.DigitalEvents.Timestamps);

			if clearCache
				obj.DigitalEvents.Timestamps = [];
				obj.DigitalEvents.Data = [];
			end

			TetrodeRecording.TTS(['Done(', num2str(toc, '%.2f'), ' seconds).\n'])
		end

		% This compresses data by ~ 20 times
		function ClearCache(obj)
			obj.Amplifier = [];
			obj.BoardDigIn = [];
			mem = memory();
			TetrodeRecording.TTS(['Cached data cleared. System memory: ', num2str(round(mem.MemUsedMATLAB/1024^2)), ' MB used (', num2str(round(mem.MemAvailableAllArrays/1024^2)), ' MB available).\n']);
		end

		% SpikeSort: PCA & Cluster
		function SpikeSort(obj, channels, varargin)
			p = inputParser;
			addRequired(p, 'Channels', @isnumeric);
			addParameter(p, 'Dimension', 10, @isnumeric);
			addParameter(p, 'FeatureMethod', 'WaveletTransform', @ischar);
			addParameter(p, 'WaveformWindow', [], @isnumeric);
			addParameter(p, 'ClusterMethod', 'gaussian', @ischar);
			addParameter(p, 'NumClusters', 3, @isnumeric);
			parse(p, channels, varargin{:});
			channels = p.Results.Channels;
			dimension = p.Results.Dimension;
			featureMethod = p.Results.FeatureMethod;
			waveformWindow = p.Results.WaveformWindow;
			clusterMethod = p.Results.ClusterMethod;
			numClusters = p.Results.NumClusters;

			obj.RemoveNaNs(channels);
			obj.FeatureExtract(channels, 'WaveformWindow', waveformWindow, 'Method', featureMethod, 'Dimension', dimension);
			obj.Cluster(channels, 'Method', clusterMethod, 'NumClusters', numClusters);
		end

		function RemoveNaNs(obj, channels)
			for iChannel = channels
				if isempty(obj.Spikes(iChannel).Waveforms)
					continue
				end
				iWaveformToDiscard = sum(isnan(obj.Spikes(iChannel).Waveforms), 2) > 0;
				obj.Spikes(iChannel).Waveforms(iWaveformToDiscard, :) = [];
				obj.Spikes(iChannel).Timestamps(iWaveformToDiscard) = [];
				obj.Spikes(iChannel).SampleIndex(iWaveformToDiscard) = [];
			end
		end

		function FeatureExtract(obj, channels, varargin)
			p = inputParser;
			addRequired(p, 'Channels', @isnumeric);
			addParameter(p, 'WaveformWindow', [], @isnumeric);
			addParameter(p, 'Method', 'WaveletTransform', @ischar);
			addParameter(p, 'WaveDecLevel', 4, @isnumeric);
			addParameter(p, 'Dimension', 10, @isnumeric);
			parse(p, channels, varargin{:});
			channels = p.Results.Channels;
			waveformWindow = p.Results.WaveformWindow;
			method = p.Results.Method;
			waveDecLevel = p.Results.WaveDecLevel;
			dimension = p.Results.Dimension;

			switch lower(method)
				case 'wavelettransform'
					methodDisplayName = 'wavelet transform';
				case 'pca'
					methodDisplayName = 'PCA';
				otherwise
					error('Unrecognized feature extraction method.')
			end

			tic, TetrodeRecording.TTS(['	Extracting waveform features (', methodDisplayName, '):\n']);

			for iChannel = channels
				if isempty(obj.Spikes(iChannel).Waveforms)
					continue
				end

				tic, TetrodeRecording.TTS(['		Channel ', num2str(iChannel), '...']);

				obj.Spikes(iChannel).Feature.Method = method;

				if isempty(waveformWindow)
					thisWaveformWindow = obj.Spikes(iChannel).WaveformWindow;
				else
					thisWaveformWindow = waveformWindow;
				end
					
				switch lower(method)
					case 'wavelettransform'
						obj.Spikes(iChannel).Feature.Parameters = struct('Dimension', dimension, 'WaveformWindow', thisWaveformWindow, 'Level', waveDecLevel);
						[obj.Spikes(iChannel).Feature.Coeff, obj.Spikes(iChannel).Feature.Stats] = obj.WaveletTransform(iChannel, 'Level', waveDecLevel, 'WaveformWindow', thisWaveformWindow, 'Dimension', dimension);
					case 'pca'
						obj.Spikes(iChannel).Feature.Parameters = struct('Dimension', dimension, 'WaveformWindow', thisWaveformWindow);
						[obj.Spikes(iChannel).Feature.Coeff, obj.Spikes(iChannel).Feature.Stats] = obj.PCA(iChannel, 'WaveformWindow', thisWaveformWindow, 'Dimension', dimension);
				end

				TetrodeRecording.TTS(['Done(', num2str(toc, '%.2f'), ' seconds).\n'])
			end
		end

		function [coeff, stats] = PCA(obj, channel, varargin)
			p = inputParser;
			addRequired(p, 'Channel', @isnumeric);
			addParameter(p, 'WaveformWindow', [], @isnumeric);
			addParameter(p, 'Dimension', 10, @isnumeric);
			parse(p, channel, varargin{:});
			channel = p.Results.Channel;
			waveformWindow = p.Results.WaveformWindow;
			dimension = p.Results.Dimension;

			if isempty(waveformWindow)
				inWindow = true(size(obj.Spikes(channel).WaveformTimestamps));
			else
				inWindow = obj.Spikes(channel).WaveformTimestamps >= waveformWindow(1) & obj.Spikes(channel).WaveformTimestamps <= waveformWindow(2);
			end
			[stats, coeff] = pca(obj.Spikes(channel).Waveforms(:, inWindow));
			coeff = coeff(:, 1:dimension);
			stats = struct('Basis', stats(:, 1:dimension));
		end

		function [coeff, stats] = WaveletTransform(obj, channel, varargin)
			p = inputParser;
			addRequired(p, 'Channel', @isnumeric);
			addParameter(p, 'Level', 4, @isnumeric);
			addParameter(p, 'WaveformWindow', [], @isnumeric);
			addParameter(p, 'Dimension', 10, @isnumeric);
			parse(p, channel, varargin{:});
			channel = p.Results.Channel;
			level = p.Results.Level;
			waveformWindow = p.Results.WaveformWindow;
			dimension = p.Results.Dimension;

			if isempty(waveformWindow)
				inWindow = true(size(obj.Spikes(channel).WaveformTimestamps));
			else
				inWindow = obj.Spikes(channel).WaveformTimestamps >= waveformWindow(1) & obj.Spikes(channel).WaveformTimestamps <= waveformWindow(2);
			end

			% Wavelet transform
			numWaveforms = length(obj.Spikes(channel).Timestamps);
			waveformLength = sum(inWindow);
			coeff = zeros(numWaveforms, waveformLength);
			for iWaveform = 1:numWaveforms
				[c, ~] = wavedec(obj.Spikes(channel).Waveforms(iWaveform, inWindow), level, 'haar');	% Haar wavelet decomposition
				coeff(iWaveform, :) = c(1:waveformLength);
			end

			% Select informative coefficients via KS test. 10 coefficients with largest deviation from normality are selected.
			ksstat = zeros(1, waveformLength);
			for iCoeff = 1:waveformLength
				thisCoeff = coeff(:, iCoeff);
				% Discard outliers (3 sigma)
				thisCoeff = thisCoeff(thisCoeff > mean(thisCoeff) - 3*std(thisCoeff) & thisCoeff < mean(thisCoeff) + 3*std(thisCoeff));

				if (sum(thisCoeff == 0) ~= length(thisCoeff)) && length(thisCoeff) > dimension
					[~, ~, ksstat(iCoeff)] = kstest((thisCoeff - mean(thisCoeff))./std(thisCoeff)); % Check normality for each coefficient
				else
					ksstat(iCoeff) = 0;
				end
			end
			% Pick 10 least 'normal' coefficients (largest ksstat)
			[~, I] = sort(ksstat, 'descend');
			I = I(1:min(dimension, waveformLength));
			coeff = coeff(:, I);
			stats = struct('ksstat', ksstat(I));
		end

		function Cluster(obj, channels, varargin)
			p = inputParser;
			addRequired(p, 'Channels', @isnumeric);
			addParameter(p, 'Method', 'kmeans', @ischar);
			addParameter(p, 'Dimension', [], @isnumeric);
			addParameter(p, 'Clusters', [], @isnumeric);
			addParameter(p, 'NumClusters', [], @isnumeric);
			parse(p, channels, varargin{:});
			channels = p.Results.Channels;
			method = p.Results.Method;
			dimension = p.Results.Dimension;
			clusters = p.Results.Clusters;
			numClusters = p.Results.NumClusters;

			switch lower(method)
				case 'kmeans'
					methodDisplayName = 'k-means';
				case 'gaussian'
					methodDisplayName = 'Gaussian mixture model';
				case 'spc'
					methodDisplayName = 'superparamagnetic';
				otherwise
					error('Unrecognized clustering method. Must be ''kmeans'', ''gaussian'', or ''SPC''.')			
			end			
			tic, TetrodeRecording.TTS(['	Clustering (', methodDisplayName, '):\n']);
			for iChannel = channels
				numWaveforms = size(obj.Spikes(iChannel).Waveforms, 1);
				if isempty(clusters)
					selected = true(1, numWaveforms);
				else
					selected = ismember(obj.Spikes(iChannel).Cluster.Classes, clusters);
				end

				if isempty(dimension) || dimension > size(obj.Spikes(iChannel).Feature.Coeff, 2)
					thisDimension = size(obj.Spikes(iChannel).Feature.Coeff, 2);
				else
					thisDimension = dimension;
				end
				obj.Spikes(iChannel).Cluster.Method = method;
				if isempty(obj.Spikes(iChannel).Waveforms)
					continue
				end
				tic, TetrodeRecording.TTS(['		Channel ', num2str(iChannel), '...']);
				switch lower(method)
					case 'kmeans'
						[classesSelected, obj.Spikes(iChannel).Cluster.Stats] = obj.KMeans(iChannel, 'NumClusters', numClusters, 'Dimension', thisDimension, 'SelectedWaveforms', selected);
					case 'gaussian'
						gm = fitgmdist(obj.Spikes(iChannel).Feature.Coeff(selected, 1:thisDimension), numClusters, 'RegularizationValue', 0.001);
						classesSelected = cluster(gm, obj.Spikes(iChannel).Feature.Coeff(selected, 1:thisDimension));
					case 'spc'
						[classesSelected, obj.Spikes(iChannel).Cluster.Stats] = obj.SPC(iChannel, 'Dimension', thisDimension, 'SelectedWaveforms', selected);
				end

				if nnz(selected) ~= numWaveforms
					classesUntouched = obj.Spikes(iChannel).Cluster.Classes(~selected);
					numClassesTotal = length(unique(classesUntouched)) + length(unique(classesSelected));
					map = setdiff(1:numClassesTotal, classesUntouched);
					classesSelectedCopy = classesSelected;
					for iClass = unique(classesSelectedCopy)'
						classesSelected(classesSelectedCopy == iClass) = map(iClass);
					end
				end
				obj.Spikes(iChannel).Cluster.Classes(selected) = classesSelected;

				TetrodeRecording.TTS(['Done(', num2str(toc, '%.2f'), ' seconds).\n'])
			end
		end

		% kmeans clustering
		function varargout = KMeans(obj, channel, varargin)
			p = inputParser;
			addRequired(p, 'Channel', @isnumeric);
			addParameter(p, 'Dimension', [], @isnumeric);
			addParameter(p, 'NumClusters', [], @isnumeric);
			addParameter(p, 'MaxNumClusters', 7, @isnumeric);
			addParameter(p, 'SelectedWaveforms', [], @(x) (islogical(x) | isnumeric(x)));
			parse(p, channel, varargin{:});
			channel 			= p.Results.Channel;
			dimension 			= p.Results.Dimension;
			numClusters			= p.Results.NumClusters;
			maxNumClusters 		= p.Results.MaxNumClusters;
			selectedWaveforms	= p.Results.SelectedWaveforms;

			feature = obj.Spikes(channel).Feature.Coeff;

			if isempty(selectedWaveforms)
				selectedWaveforms = true(1, size(feature, 1));
			end
			if isempty(dimension) || dimension > size(feature, 2)
				dimension = size(feature, 2);
			end

			feature = feature(selectedWaveforms, 1:dimension);

			warning('off')

			if isempty(numClusters)
				stats = evalclusters(feature, 'kmeans', 'CalinskiHarabasz', 'KList', 2:maxNumClusters);
				numClusters = stats.OptimalK;
			else
				stats = [];
			end
			classes = kmeans(feature, numClusters, 'Replicates', 10);

			warning('on')
			varargout = {classes, stats};
		end

		% Superparamagnetic clustering
		function varargout = SPC(obj, channel, varargin)
			p = inputParser;
			addRequired(p, 'Channel', @isnumeric);
			addParameter(p, 'Dimension', [], @isnumeric);
			addParameter(p, 'MinTemp', 0, @isnumeric);
			addParameter(p, 'MaxTemp', 0.251, @isnumeric);
			addParameter(p, 'TempStep', 0.01, @isnumeric);
			addParameter(p, 'NumSteps', 3, @isnumeric);
			addParameter(p, 'SWCycles', 100, @isnumeric);
			addParameter(p, 'KNearestNeighbours', 11, @isnumeric);
			addParameter(p, 'MinClusterSize', NaN, @isnumeric);
			addParameter(p, 'MinClusterSizeRatio', 0.0015, @isnumeric);
			addParameter(p, 'MaxNumClusters', 13, @isnumeric);
			addParameter(p, 'MaxNumWaveforms', 40000, @isnumeric); % If too many spikes use template matching for extra spikes
			addParameter(p, 'SelectedWaveforms', [], @(x) (islogical(x) | isnumeric(x)));
			parse(p, channel, varargin{:});
			channel 			= p.Results.Channel;
			dimension 			= p.Results.Dimension;
			minTemp 			= p.Results.MinTemp;
			maxTemp 			= p.Results.MaxTemp;
			tempStep 			= p.Results.TempStep;
			numSteps 			= p.Results.NumSteps;
			minClusterSize 		= p.Results.MinClusterSize;
			minClusterSizeRatio = p.Results.MinClusterSizeRatio;
			maxNumClusters 		= p.Results.MaxNumClusters;
			maxNumWaveforms		= p.Results.MaxNumWaveforms;
			selectedWaveforms	= p.Results.SelectedWaveforms;

			feature = obj.Spikes(channel).Feature.Coeff;
			if isempty(selectedWaveforms)
				selectedWaveforms = true(1, size(feature, 1));
			end
			feature = feature(selectedWaveforms, :);
			waveforms = obj.Spikes(channel).Waveforms(selectedWaveforms, :);
			numWaveforms = size(feature, 1);

			if numWaveforms > maxNumWaveforms
				templateMatching = true;
			else
				templateMatching = false;
			end

			if templateMatching
				indicesSPC = randperm(numWaveforms, maxNumWaveforms);
				indicesTemplateMatching = setdiff(1:numWaveforms, indicesSPC);
				featureSPC = feature(indicesSPC, :);
				waveformsSPC = waveforms(indicesSPC, :);
				waveformsTemplateMatching = waveforms(indicesTemplateMatching, :);
			else
				featureSPC = feature;
			end

			if isnan(minClusterSize)
				minClusterSize = round(size(featureSPC, 1)*minClusterSizeRatio);
			end

			if isempty(dimension) || dimension > size(feature, 2)
				dimension = size(feature, 2);
			end

			fileIn = 'temp_in';
			fileOut = 'temp_out';
			save(fileIn, 'featureSPC', '-ascii');

			classesSPC = [];
			clu = [];
			tree = [];

			temps = minTemp:tempStep:maxTemp;
			for iTemp = 1:numSteps:length(temps)
				thisMinTemp = temps(iTemp);
				thisMaxTemp = temps(min(length(temps), iTemp + numSteps - 1));
				fid = fopen(sprintf('%s.run', fileOut), 'wt');
				fprintf(fid, 'NumberOfPoints: %s\n', num2str(size(featureSPC, 1)));
				fprintf(fid, 'DataFile: %s\n', fileIn);
				fprintf(fid, 'OutFile: %s\n', fileOut);
				fprintf(fid, 'Dimensions: %s\n', num2str(dimension));
				fprintf(fid, 'MinTemp: %s\n', num2str(thisMinTemp));
				fprintf(fid, 'MaxTemp: %s\n', num2str(thisMaxTemp));
				fprintf(fid, 'TempStep: %s\n', num2str(tempStep));
				fprintf(fid, 'SWCycles: %s\n', num2str(p.Results.SWCycles));
				fprintf(fid, 'KNearestNeighbours: %s\n', num2str(p.Results.KNearestNeighbours));
				fprintf(fid, 'MSTree|\n');
				fprintf(fid, 'DirectedGrowth|\n');
				fprintf(fid, 'SaveSuscept|\n');
				fprintf(fid, 'WriteLables|\n');
				fprintf(fid, 'WriteCorFile~\n');
				fclose(fid);

				[status,result] = dos(sprintf('"%s" %s.run', which('cluster_64.exe'), fileOut));

				clu = [clu; load([fileOut, '.dg_01.lab'])];
				tree = [tree; load([fileOut, '.dg_01'])];

				delete([fileOut, '.dg_01.lab'])
				delete([fileOut, '.dg_01'])
				delete([fileOut, '.run']);
				delete([fileOut, '*.mag']);
				delete([fileOut, '*.edges']);
				delete([fileOut, '*.param']);
				delete([fileOut, '*.knn']);

				% Find proper temperature (when the 4 biggest clusters no longer significantly increase in size)
				if size(tree, 1) > 1
					dSizeCluster = diff(tree(:, 5:8));
					stable = [false; sum(dSizeCluster <= minClusterSize, 2) == 4];
					if nnz(stable) > 0
						classesSPC = clu(find(stable, 1), 3:end);
						break
					end
				end
			end

			delete(fileIn);

			for iCluster = 1:max(classesSPC)
				if nnz(classesSPC == iCluster) < minClusterSize
					classesSPC(classesSPC == iCluster) = 0;
				end
			end

			classesSPC = classesSPC + 1;
			classesSPC = classesSPC';

			if templateMatching
				% Build templates
				templates = zeros(length(unique(classesSPC)), size(waveformsSPC, 2));
				for iCluster = 1:max(classesSPC)
					templates(iCluster, :) = mean(waveformsSPC(classesSPC == iCluster, :), 1);
				end

				% Match to template via corr
				nearestIndex = knnsearch(waveformsSPC, waveformsTemplateMatching, 'k', 1);
				classesTemplateMatching = classesSPC(nearestIndex);

				classes = zeros(numWaveforms, 1);
				classes(indicesSPC) = classesSPC;
				classes(indicesTemplateMatching) = classesTemplateMatching;
			else
				classes = classesSPC;
			end

			% Output
			stats = struct('Classes', clu, 'ClusterSizes', tree, 'MinClusterSize', minClusterSize, 'BestTempIndex', find(stable, 1), 'IndicesSPC', indicesSPC, 'IndicesTemplateMatching', indicesTemplateMatching);
			varargout = {classes, stats};			
		end

		function ClusterMerge(obj, channel, mergeList, varargin)
			p = inputParser;
			addRequired(p, 'Channel', @isnumeric);
			addRequired(p, 'MergeList', @iscell);
			parse(p, channel, mergeList, varargin{:});
			channel = p.Results.Channel;
			mergeList = p.Results.MergeList;

			newCluster = NaN(size(obj.Spikes(channel).Cluster.Classes));
			for iNewCluster = 1:length(mergeList)
				newCluster(ismember(obj.Spikes(channel).Cluster.Classes, mergeList{iNewCluster})) = iNewCluster;
			end
			obj.Spikes(channel).Cluster.Classes = newCluster;
		end

		function ClusterRemove(obj, channel, discardList, varargin)
			p = inputParser;
			addRequired(p, 'Channel', @isnumeric);
			addRequired(p, 'DiscardList', @isnumeric);
			parse(p, channel, discardList, varargin{:});
			channel = p.Results.Channel;
			discardList = p.Results.DiscardList;

			% Discard unwanted clusters
			iWaveformToDiscard = ismember(obj.Spikes(channel).Cluster.Classes, discardList);
			obj.Spikes(channel).Waveforms(iWaveformToDiscard, :) = [];
			obj.Spikes(channel).Timestamps(iWaveformToDiscard) = [];
			obj.Spikes(channel).SampleIndex(iWaveformToDiscard) = [];
			obj.Spikes(channel).Feature.Coeff(iWaveformToDiscard, :) = [];
			obj.Spikes(channel).Cluster.Classes(iWaveformToDiscard) = [];

			% Renumber clusters from 1 to numClusters
			obj.ClusterMerge(channel, num2cell(unique(obj.Spikes(channel).Cluster.Classes)));
		end

		function ClusterDecimate(obj, channel, cluster, varargin)
			p = inputParser;
			addRequired(p, 'Channel', @isnumeric);
			addRequired(p, 'Cluster', @isnumeric);
			addOptional(p, 'LuckyNumber', 10, @isnumeric);
			parse(p, channel, cluster, varargin{:});
			channel = p.Results.Channel;
			cluster = p.Results.Cluster;
			luckyNumber = p.Results.LuckyNumber;

			inCluster = find(obj.Spikes(channel).Cluster.Classes == cluster);
			notInLuck = inCluster(~ismember(1:length(inCluster), 1:luckyNumber:length(inCluster)));

			obj.Spikes(channel).SampleIndex(notInLuck) 			= [];
			obj.Spikes(channel).Timestamps(notInLuck) 			= [];
			obj.Spikes(channel).Waveforms(notInLuck, :) 		= [];
			obj.Spikes(channel).Feature.Coeff(notInLuck, :)		= [];
			obj.Spikes(channel).Cluster.Classes(notInLuck) 		= [];
		end

		% Plot raw channel
		function PlotChannelRaw(obj, channel, xRange, yRange, digInChannels)
			if nargin < 3
				xRange = NaN;
			end
			if nargin < 4
				yRange = [-110, 110];
			end
			if nargin < 5
				digInChannels = [4, 2, 5, 3]; % CUE, LEVER, REWARD, LICK
			end
			digInScaled = obj.BoardDigIn.Data(digInChannels, :).*[100; 500; 200; 350]; % CUE, LEVER, REWARD, LICK

			obj.Amplifier.DataMean = nanmean(obj.Amplifier.Data, 1);

			figure('Name', ['Channel ', num2str(channel)], 'Units', 'normalized', 'Position', [0, 0, 1, 0.6])

			subplot(2, 1, 1)
			hold on
			plot(obj.Amplifier.Timestamps, obj.Amplifier.Data(channel, :));
			plot(obj.BoardDigIn.Timestamps, digInScaled, 'r:', 'LineWidth', 2);
			hold off
			if ~isnan(xRange)
				xlim(xRange)
			end
			ylim(yRange)
			xlabel('Time (s)')
			ylabel('\muV')
			title('Raw')

			subplot(2, 1, 2)
			hold on
			plot(obj.Amplifier.Timestamps, obj.Amplifier.Data(channel, :) - obj.Amplifier.DataMean);
			plot(obj.BoardDigIn.Timestamps, digInScaled, 'r:', 'LineWidth', 2);
			hold off
			if ~isnan(xRange)
				xlim(xRange)
			end
			ylim(yRange)
			xlabel('Time (s)')
			ylabel('\muV')
			title('Substracted by 32-Chn mean')
		end

		function PlotWaveforms(obj, channel, varargin)
			p = inputParser;
			addRequired(p, 'Channel', @isnumeric);
			addOptional(p, 'NumWaveforms', 50, @isnumeric);
			addParameter(p, 'Clusters', [], @isnumeric);
			addParameter(p, 'SelectedCluster', [], @isnumeric);
			addParameter(p, 'ReferenceCluster', [], @isnumeric);
			addParameter(p, 'YLim', [], @isnumeric);
			addParameter(p, 'FrameRate', 120, @isnumeric);
			addParameter(p, 'MaxShown', 200, @isnumeric);
			addParameter(p, 'WaveformWindow', [], @isnumeric);
			addParameter(p, 'PlotMean', false, @islogical);
			addParameter(p, 'Ax', []);
			parse(p, channel, varargin{:});
			channel 			= p.Results.Channel;
			numWaveforms 		= p.Results.NumWaveforms;
			clusters 			= p.Results.Clusters;
			selectedCluster 	= p.Results.SelectedCluster;
			referenceCluster 	= p.Results.ReferenceCluster;
			yRange 				= p.Results.YLim;
			frameRate 			= p.Results.FrameRate;
			maxShown 			= p.Results.MaxShown;
			waveformWindow 		= p.Results.WaveformWindow;
			plotMean 			= p.Results.PlotMean;
			ax 					= p.Results.Ax;

			waveforms = obj.Spikes(channel).Waveforms;
			numWaveformsTotal = size(waveforms, 1);

			if isempty(yRange)
				yRange = [min(waveforms(:)), max(waveforms(:))];
			end

			if isempty(waveformWindow)
				waveformWindow = obj.Spikes(channel).WaveformWindow;
			end
			t = obj.Spikes(channel).WaveformTimestamps;
			score = obj.Spikes(channel).Feature.Coeff;
			if isempty(clusters)
				clusters = nonzeros(unique(obj.Spikes(channel).Cluster.Classes));
			end

			if isfield(obj.Spikes(channel), 'Cluster')
				clusterID = obj.Spikes(channel).Cluster.Classes;
				toKeep = ismember(clusterID, clusters);
				clusterID = clusterID(toKeep);
				waveforms = waveforms(toKeep, :);
				score = score(toKeep, :);
				numWaveformsTotal = size(waveforms, 1);
			else
				clusterID = ones(numWaveformsTotal, 1);
			end

			if isempty(ax) || nnz(ishandle(ax)) == 0
				hFigure = figure('Units', 'Normalized', 'OuterPosition', [0.7, 0, 0.3, 1]);
				hAxes1 = subplot(2, 1, 1);
				hAxes2 = subplot(2, 1, 2);
			else
				hAxes1 = ax(1);
				hAxes2 = ax(2);
				if ishandle(hAxes1)
					hFigure = hAxes1.Parent;
				else
					hFigure = hAxes2.Parent;
				end
			end

			if ishandle(hAxes1)
				axes(hAxes1)
				hold(hAxes1, 'on')
				hLegends = [];
				for iCluster = unique(nonzeros(clusterID))'
					inCluster = clusterID == iCluster;
					if sum(inCluster) == 0
						continue
					end
					percentage = 100*sum(inCluster)/size(obj.Spikes(channel).Waveforms, 1);
					percentage = [num2str(percentage, '%.1f'), '%'];
					count = sum(inCluster);
					if count < 1000
						count = num2str(count);
					else
						count = [num2str(count/1000, '%.1f'), 'k'];
					end

					if isempty(selectedCluster)
						[thisColor, ~] = TetrodeRecording.GetColorAndStyle(iCluster);
					else
						if iCluster == selectedCluster
							[thisColor, ~] = TetrodeRecording.GetColorAndStyle(iCluster);
							if thisColor == 'k'
								thisColor = 'r';
							end
						else
							thisColor = 'k';
						end
					end
					dispName = ['Cluster ', num2str(iCluster), ' (', percentage, ' | ', count, ')'];
					h = scatter3(hAxes1, score(inCluster, 1), score(inCluster, 2), score(inCluster, 3), 1, thisColor, 'DisplayName', dispName);
					if isempty(selectedCluster) || (~isempty(selectedCluster) && iCluster == selectedCluster)
						hLegends = [hLegends, h];
						if iCluster == selectedCluster
							title(hAxes2, dispName)
						end
					end
				end
				hold(hAxes1, 'off')
				xlabel(hAxes1, '1st Coefficient')
				ylabel(hAxes1, '2nd Coefficient')
				zlabel(hAxes1, '3rd Coefficient')
				if isempty(selectedCluster)
					title(hAxes1, 'Feature space')
				end
				legend(hLegends, 'Location', 'Best')
			end

			if ishandle(hAxes2)
				axes(hAxes2)
				hold(hAxes2, 'on')
				hLegends = [];
				hAxes2.UserData.hWaveforms = [];
				hAxes2.UserData.iWaveform = 0;
				hLegends = [hLegends, line(hAxes2, 'XData', [obj.Spikes(channel).WaveformWindow(1), 0], 'YData', repmat(mean(obj.Spikes(channel).Threshold.Threshold), [1, 2]), 'Color', 'k', 'LineWidth', 3, 'DisplayName', ['Threshold (', num2str(mean(obj.Spikes(channel).Threshold.Threshold)), ' \muV)'])];
				legend(hLegends, 'AutoUpdate', 'off', 'Location', 'Best')
				xlim(hAxes2, waveformWindow)
				ylim(hAxes2, yRange)

				if frameRate > 0
					hTimer = timer(...
						'ExecutionMode', 'FixedSpacing',...
					 	'Period', round((1/frameRate)*1000)/1000,...
					 	'TimerFcn', {@TetrodeRecording.OnPlotChannelRefresh, hAxes2, t, waveforms, numWaveforms, numWaveformsTotal, clusterID}...
					 	);
					hAxes2.UserData.hTimer = hTimer;
					hFigure.KeyPressFcn = {@TetrodeRecording.OnKeyPress, hTimer};
					hFigure.CloseRequestFcn = {@TetrodeRecording.OnFigureClosed, hTimer};
					start(hTimer);
				else
					xlabel(hAxes2, 'Time (ms)');
					ylabel(hAxes2, 'Voltage (\muV)');
					if isempty(selectedCluster)
						title(hAxes2, 'Waveforms');
					end
					for iCluster = unique(nonzeros(clusterID))'
						% Skip this cluster in some modes
						if (~isempty(selectedCluster) && ~ismember(iCluster, [selectedCluster, referenceCluster]))
							continue
						end

						thisWaveforms = waveforms(clusterID==iCluster, :);
						
						if iCluster == referenceCluster
							thisColor = 'k';
							thisStyle = '--';
						else
							[thisColor, thisStyle] = TetrodeRecording.GetColorAndStyle(iCluster);
						end
						if plotMean
							thisMean = mean(thisWaveforms, 1);
							thisStd = std(thisWaveforms, 0, 1);
							thisUp = prctile(thisWaveforms, 95, 1);
							thisDown = prctile(thisWaveforms, 5, 1);
							line(hAxes2, t, thisMean, 'LineStyle', thisStyle, 'Color', thisColor);
							patch(hAxes2, [t, t(end:-1:1)], [thisDown, thisUp(end:-1:1)], thisColor,...
								'FaceAlpha', 0.15, 'EdgeColor', 'none');
							patch(hAxes2, [t, t(end:-1:1)], [thisMean - thisStd, thisMean(end:-1:1) + thisStd(end:-1:1)], thisColor,...
								'FaceAlpha', 0.4, 'EdgeColor', 'none');
							line(hAxes2, t, [thisMean - thisStd; thisMean + thisStd], 'LineStyle', '--', 'Color', thisColor);
							line(hAxes2, t, [thisDown; thisUp], 'LineStyle', ':', 'Color', thisColor);
						else
							if size(thisWaveforms, 1) > maxShown
								percentShown = max(1, round(100*(maxShown/size(thisWaveforms, 1))));
							else
								percentShown = 100;
							end
							thisWaveforms = thisWaveforms(1:ceil(100/percentShown):end, :);
							line(hAxes2, t, thisWaveforms, 'LineStyle', thisStyle, 'Color', thisColor);
						end
					end
				end
			end
		end

		function Raster(obj, channels, varargin)
			p = inputParser;
			addRequired(p, 'Channels', @isnumeric);
			addOptional(p, 'Reference', 'CueOn', @ischar);
			addOptional(p, 'Event', 'PressOn', @ischar);
			addOptional(p, 'Exclude', 'LickOn', @ischar);
			addParameter(p, 'AlignTo', 'Event', @ischar);
			addParameter(p, 'Sort', true, @islogical);
			addParameter(p, 'ExtendedWindow', [-2, 2], @isnumeric);
			addParameter(p, 'Clusters', [], @isnumeric);
			addParameter(p, 'XLim', [], @isnumeric);
			addParameter(p, 'Ax', []);
			parse(p, channels, varargin{:});
			channels = p.Results.Channels;
			reference = p.Results.Reference;
			event = p.Results.Event;
			exclude = p.Results.Exclude;
			alignTo = p.Results.AlignTo;
			doSort = p.Results.Sort;
			extendedWindow = p.Results.ExtendedWindow;
			clusters = p.Results.Clusters;
			xRange = p.Results.XLim;
			ax = p.Results.Ax;

			if ~isempty(reference)
				referenceDisplayName = reference;
				reference = obj.DigitalEvents.(reference);
			else
				reference = [];
			end
			if ~isempty(event)
				eventDisplayName = event;
				event = obj.DigitalEvents.(event);
			else
				event = [];
			end
			if ~isempty(exclude)
				exclude = obj.DigitalEvents.(exclude);
			else
				exclude = [];
			end

			% Recalculate timestamps relative to cue on
			% Find the first lever press (true first movement: before any lick/press has occured since cue on)

			% Get spikes between two reference and first event
			[reference, event, ~, ~] = TetrodeRecording.FindFirstInTrial(reference, event, exclude);

			% Bin spikes into trials
			for iChannel = channels
				[spikes, trials] = obj.GetSpikesByTrial(iChannel, reference, event, extendedWindow, clusters);
				% Sort trials by time to movement
				if doSort
					[~, I] = sort(reference - event);
					trialsSorted = changem(trials, 1:length(unique(I)), I);	% !!! changem requires mapping toolbox.
				else
					trialsSorted = trials;
				end

				if isempty(ax)
					hFigure = figure('Units', 'Normalized', 'OuterPosition', [0, 0, 0.75, 1]);
					hAxes = axes(hFigure);
				else
					hAxes = ax;
					axes(hAxes);
				end

				if strcmpi(alignTo, 'Event') || strcmpi(alignTo, 'Movement')
					spikesRelative = spikes - event(trials);
					eventRelative = reference(trials) - event(trials);
					eventDisplayNameRelative = referenceDisplayName; 
					referenceDisplayNameRelative = eventDisplayName; 
				elseif strcmpi(alignTo, 'Reference') || strcmpi(alignTo, 'Cue')
					spikesRelative = spikes - reference(trials);
					eventRelative = event(trials) - reference(trials);
					eventDisplayNameRelative = eventDisplayName; 
					referenceDisplayNameRelative = referenceDisplayName; 
				end
				hold on
				plot(hAxes, spikesRelative, trialsSorted, '.',...
					'MarkerSize', 5,...
					'MarkerEdgeColor', 'k',...
					'MarkerFaceColor', 'k',...
					'LineWidth', 1.5,...
					'DisplayName', 'Spike'...
				)
				plot(hAxes, eventRelative, trialsSorted, '.',...
					'MarkerSize', 10,...
					'MarkerEdgeColor', 'r',...
					'MarkerFaceColor', 'r',...
					'LineWidth', 1.5,...
					'DisplayName', eventDisplayNameRelative...
				)
				title(hAxes, 'Spike raster')
				xlabel(hAxes, ['Time relative to ', referenceDisplayNameRelative, ' (s)'])
				ylabel(hAxes, 'Trial')
				legend(hAxes, 'Location', 'Best');
				if ~isempty(xRange)
					xlim(hAxes, xRange);
				end
				hold off
			end
		end

		function PETH(obj, channels, varargin)
			p = inputParser;
			addRequired(p, 'Channels', @isnumeric);
			addOptional(p, 'Reference', 'CueOn', @ischar);
			addOptional(p, 'Event', 'PressOn', @ischar);
			addOptional(p, 'Exclude', 'LickOn', @ischar);
			addParameter(p, 'Clusters', [], @isnumeric);
			addParameter(p, 'MinTrialLength', 0, @isnumeric);
			addParameter(p, 'Bins', 5, @isnumeric);
			addParameter(p, 'BinMethod', 'percentile', @ischar);
			addParameter(p, 'SpikeRateWindow', 100, @isnumeric);
			addParameter(p, 'ExtendedWindow', [0, 0], @isnumeric);
			addParameter(p, 'Ax', []);
			parse(p, channels, varargin{:});
			channels = p.Results.Channels;
			reference = p.Results.Reference;
			event = p.Results.Event;
			exclude = p.Results.Exclude;
			clusters = p.Results.Clusters;
			minTrialLength = p.Results.MinTrialLength;
			nBins = p.Results.Bins;
			binMethod = p.Results.BinMethod;
			spikeRateWindow = p.Results.SpikeRateWindow;
			extendedWindow = p.Results.ExtendedWindow;
			ax = p.Results.Ax;

			if ~isempty(reference)
				referenceDisplayName = reference;
				reference = obj.DigitalEvents.(reference);
			else
				reference = [];
			end
			if ~isempty(event)
				eventDisplayName = event;
				event = obj.DigitalEvents.(event);
			else
				event = [];
			end
			if ~isempty(exclude)
				exclude = obj.DigitalEvents.(exclude);
			else
				exclude = [];
			end

			[reference, event, ~, ~] = TetrodeRecording.FindFirstInTrial(reference, event, exclude);

			% Bin trials according to trial length (t(event) - t(reference))
			trialLength = event - reference;
			switch binMethod
				case 'percentile'
					edges = prctile(trialLength(trialLength > minTrialLength), 0:(100/nBins):100);
				case 'equal'
					edges = linspace(max(minTrialLength, min(trialLength)), max(trialLength), nBins + 1);
				otherwise
					error('Unrecognized bin method.')
			end
			[NTrials, ~, bins] = histcounts(trialLength, edges);

			for iChannel = channels
				[spikes, trials] = obj.GetSpikesByTrial(iChannel, reference, event, extendedWindow, clusters);

				if isempty(ax)
					hAxes = axes(figure());
				else
					hAxes = ax;
					axes(hAxes);
				end
				hold on
				for iBin = 1:nBins
					inBin = ismember(trials, find(bins == iBin));
					cutOff = -edges(iBin);
					spikesRelative = spikes(inBin) - event(trials(inBin)); % Spike times relative to event
					spikesRelative = spikesRelative(spikesRelative > (cutOff + extendedWindow(1)) & spikesRelative < extendedWindow(2));
					% Windows for estimating spike rate
					thisEdges = (cutOff + extendedWindow(1)):(spikeRateWindow/1000):extendedWindow(2);
					if length(thisEdges) < 3
						continue
					end
					% if thisEdges(end) < 0
					% 	thisEdges(end + 1) = 0;
					% end
					thisSpikeRate = histcounts(spikesRelative, thisEdges);
					thisSpikeRate = (1000*thisSpikeRate/spikeRateWindow)/NTrials(iBin);
					thisCenters = (thisEdges(1:end - 1) + thisEdges(2:end))/2;

					[thisColor, thisStyle] = TetrodeRecording.GetColorAndStyle(iBin);
					plot(hAxes, thisCenters, thisSpikeRate, [thisColor, thisStyle],...
						'DisplayName', ['[', num2str(-cutOff, 2), ' s, ', num2str(edges(iBin + 1), 2), ' s] (', num2str(NTrials(iBin)), ' trials)'],...
						'LineWidth', 2.5);
				end
				xlabel(hAxes, ['Time relative to ', eventDisplayName, ' (s)']);
				ylabel(hAxes, 'Mean firing rate (Hz)')
				legend(hAxes, 'Location', 'Best');
				title(hAxes, 'Peri-event time histogram');
				hold off
			end
		end

		function PlotChannel(obj, channel, varargin)
			p = inputParser;
			addRequired(p, 'Channel', @isnumeric);
			addParameter(p, 'Reference', 'CueOn', @ischar);
			addParameter(p, 'Event', 'PressOn', @ischar);
			addParameter(p, 'Exclude', 'LickOn', @ischar);
			addParameter(p, 'Clusters', [], @isnumeric);
			addParameter(p, 'WaveformWindow', [], @isnumeric);
			addParameter(p, 'ExtendedWindow', [0, 0], @isnumeric);
			addParameter(p, 'MinTrialLength', 0, @isnumeric);
			addParameter(p, 'Bins', 4, @isnumeric);
			addParameter(p, 'BinMethod', 'percentile', @ischar);
			addParameter(p, 'SpikeRateWindow', 100, @isnumeric);
			addParameter(p, 'RasterXLim', [], @isnumeric);
			addParameter(p, 'WaveformYLim', [-200, 200], @isnumeric);
			addParameter(p, 'FontSize', 8, @isnumeric);
			addParameter(p, 'PrintMode', false, @islogical);
			addParameter(p, 'FrameRate', 0, @isnumeric);
			addParameter(p, 'Fig', []);
			parse(p, channel, varargin{:});
			iChannel 		= p.Results.Channel;
			reference 		= p.Results.Reference;
			event 			= p.Results.Event;
			exclude 		= p.Results.Exclude;
			clusters 		= p.Results.Clusters;
			waveformWindow 	= p.Results.WaveformWindow;
			extendedWindow 	= p.Results.ExtendedWindow;
			minTrialLength 	= p.Results.MinTrialLength;
			bins 			= p.Results.Bins;
			binMethod 		= p.Results.BinMethod;
			spikeRateWindow = p.Results.SpikeRateWindow;
			rasterXLim 		= p.Results.RasterXLim;
			waveformYLim	= p.Results.WaveformYLim;
			fontSize 		= p.Results.FontSize;
			printMode 		= p.Results.PrintMode;
			frameRate 		= p.Results.FrameRate;
			hFigure 		= p.Results.Fig;

			allChannels = [obj.Spikes.Channel];
			if isempty(iChannel)
				iChannel = allChannels(1);
			end

			xRatio 	= 0.4;
			yRatio 	= 0.6;
			xMargin = 0.04;
			yMargin = 0.06;
			fMargin = 0.04;
			buttonWidth = 0.04;
			buttonHeight = 0.03;
			buttonSpacing = 0.0075;
			hUp 	= (1 - 2*fMargin)*yRatio - 2*yMargin;
			hUpHalf = (hUp - 2*yMargin)/2;
			hDown 	= (1 - 2*fMargin)*(1 - yRatio) - 2*yMargin;
			w 		= 1 - 2*fMargin - 2*xMargin;
			wLeft	= (1 - 2*fMargin)*xRatio - 2*xMargin;
			wRight	= (1 - 2*fMargin)*(1 - xRatio) - 2*xMargin;

			expName = obj.GetExpName();				
			displayName = [expName, ' (Channel ', num2str(iChannel), ')'];

			if ~isempty(hFigure)
				clf(hFigure)
			else
				hFigure	= figure('Units', 'Normalized', 'Position', [0, 0, 1, 1], 'GraphicsSmoothing', 'on');
				hFigure.UserData.PlotMean = true;
				hFigure.UserData.ReferenceCluster = [];
			end
			hFigure.Name = displayName;

			if printMode
				fontSize = 14;
				hFigure.Units = 'pixels';
				hFigure.InnerPosition = [0, 0, 1332, 999];
			end
			set(hFigure, 'DefaultAxesFontSize', fontSize);

			hWaveform 	= subplot('Position', [fMargin + xMargin, fMargin + 5*yMargin + hDown + hUpHalf, wLeft, hUpHalf]);
			hPCA 		= subplot('Position', [fMargin + xMargin, fMargin + 3*yMargin + hDown, wLeft, hUpHalf]);
			hRaster		= subplot('Position', [fMargin + xMargin + wLeft + 2*xMargin, fMargin + 3*yMargin + hDown, wRight, hUp]);
			hPETH 		= subplot('Position', [fMargin + xMargin, fMargin + yMargin, w, hDown]);

			if isempty(clusters)
				hFigure.UserData.SelectedClusters = unique(obj.Spikes(iChannel).Cluster.Classes);
			else
				hFigure.UserData.SelectedClusters = clusters;
			end

			hTitle = suptitle(displayName);

			hButtonPlot = uicontrol(...
				'Style', 'pushbutton',...
				'String', 'Plot ...',...
				'Callback', {@obj.GUIPlotClusters, iChannel, p, hFigure, hWaveform, hPCA, hRaster, hPETH},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Position', [xMargin, 1 - yMargin, buttonWidth, min(yMargin, buttonHeight)]);
			hButtonPlot.Position(1) = hButtonPlot.Position(1) + hButtonPlot.Position(3);
			hPrev = hButtonPlot;

			hButtonMerge = uicontrol(...
				'Style', 'pushbutton',...
				'String', 'Merge ...',...
				'Callback', {@obj.GUIMergeClusters, iChannel, p, hFigure, hWaveform, hPCA, hRaster, hPETH},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Position', hPrev.Position);
			hButtonMerge.Position(1) = hButtonMerge.Position(1) + hButtonMerge.Position(3) + buttonSpacing;
			hPrev = hButtonMerge;

			hButtonRemove = uicontrol(...
				'Style', 'pushbutton',...
				'String', 'Remove ...',...
				'Callback', {@obj.GUIRemoveClusters, iChannel, p, hFigure, hWaveform, hPCA, hRaster, hPETH},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Position', hPrev.Position);
			hButtonRemove.Position(1) = hButtonRemove.Position(1) + hButtonRemove.Position(3) + buttonSpacing;
			hPrev = hButtonRemove;

			hButtonDecimate = uicontrol(...
				'Style', 'pushbutton',...
				'String', 'Decimate ...',...
				'Callback', {@obj.GUIDecimateClusters, iChannel, p, hFigure, hWaveform, hPCA, hRaster, hPETH},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Position', hPrev.Position);
			hButtonDecimate.Position(1) = hButtonDecimate.Position(1) + hButtonDecimate.Position(3) + buttonSpacing;
			hPrev = hButtonDecimate;

			hButtonRecluster = uicontrol(...
				'Style', 'pushbutton',...
				'String', 'Recluster ...',...
				'Callback', {@obj.GUIRecluster, iChannel, p, hFigure, hWaveform, hPCA, hRaster, hPETH},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Position', hPrev.Position);
			hButtonRecluster.Position(1) = hButtonRecluster.Position(1) + hButtonRecluster.Position(3) + buttonSpacing;
			hPrev = hButtonRecluster;

			hButtonPlotMean = uicontrol(...
				'Style', 'togglebutton',...
				'String', 'Mean',...
				'Callback', {@obj.GUIPlotMean, iChannel, p, hFigure, hWaveform, hPCA, hRaster, hPETH},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Value', hFigure.UserData.PlotMean,...
				'Position', [buttonSpacing, hWaveform.Position(2) + hWaveform.Position(4) - buttonHeight, buttonWidth, buttonHeight]);
			hPrev = hButtonPlotMean;

			hButtonSelRef = uicontrol(...
				'Style', 'pushbutton',...
				'String', 'Reference ...',...
				'Callback', {@obj.GUISelRef, iChannel, p, hFigure, hWaveform, hPCA, hRaster, hPETH},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Position', hPrev.Position);
			hButtonSelRef.Position(2) = hButtonSelRef.Position(2) - hButtonSelRef.Position(4) - buttonHeight;
			hPrev = hButtonSelRef;

			hButtonPlotAllClusters = uicontrol(...
				'Style', 'pushbutton',...
				'String', 'Expand ...',...
				'Callback', {@obj.GUIPlotAllClusters, iChannel, hFigure},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Position', hPrev.Position);
			hButtonPlotAllClusters.Position(2) = hButtonPlotAllClusters.Position(2) - hButtonPlotAllClusters.Position(4) - buttonHeight;
			hPrev = hButtonPlotAllClusters;

			hButtonSavePlot = uicontrol(...
				'Style', 'pushbutton',...
				'String', 'Clipboard',...
				'Callback', {@obj.GUISavePlot, hFigure},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Position', hPrev.Position);
			hButtonSavePlot.Position(2) = hButtonSavePlot.Position(2) - hButtonSavePlot.Position(4) - buttonHeight;
			hPrev = hButtonSavePlot;

			nextChn = find(allChannels == iChannel) + 1;
			nextChn = mod(nextChn, length(allChannels));
			if nextChn == 0
				nextChn = length(allChannels);
			end
			nextChn = allChannels(nextChn);

			prevChn = find(allChannels == iChannel) - 1;
			prevChn = mod(prevChn, length(allChannels));
			if prevChn == 0
				prevChn = length(allChannels);
			end
			prevChn = allChannels(prevChn);

			hButtonNextChn = uicontrol(hFigure,...
				'Style', 'pushbutton',...
				'String', 'Next Chn',...
				'Callback', {@(~, ~, iChannel, hFigure) obj.PlotChannel(iChannel, 'Fig', hFigure, 'PrintMode', printMode), nextChn, hFigure},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Position', [1 - xMargin - 2*buttonWidth, 1 - yMargin, buttonWidth, min(yMargin, buttonHeight)]);
			hPrev = hButtonNextChn;

			hAxesTextCurChn = axes('Position', hPrev.Position, 'Visible', 'off');
			hTextCurChn = text(hAxesTextCurChn, 0.5, 0.5,...
				['Chn ', num2str(find(allChannels == iChannel)), '/', num2str(length(allChannels))],...
				'HorizontalAlignment', 'center');
			hAxesTextCurChn.Position(1) = hAxesTextCurChn.Position(1) - hAxesTextCurChn.Position(3) - buttonSpacing;
			hPrev = hAxesTextCurChn;

			hButtonDeleteChn = uicontrol(hFigure,...
				'Style', 'pushbutton',...
				'String', 'Delete Chn',...
				'Callback', {@obj.GUIDeleteChannel, iChannel, nextChn, hFigure},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Position', hPrev.Position);
			hButtonDeleteChn.Position(2) = hButtonDeleteChn.Position(2) - hButtonDeleteChn.Position(4);

			hButtonPrevChn = uicontrol(hFigure,...
				'Style', 'pushbutton',...
				'String', 'Prev Chn',...
				'Callback', {@(~, ~, iChannel, hFigure) obj.PlotChannel(iChannel, 'Fig', hFigure, 'PrintMode', printMode), prevChn, hFigure},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Position', hPrev.Position);
			hButtonPrevChn.Position(1) = hButtonPrevChn.Position(1) - hButtonPrevChn.Position(3) - buttonSpacing;
			hPrev = hButtonPrevChn;

			obj.GUIBusy(hFigure, true);
			obj.ReplotChannel(iChannel, p, hFigure, hWaveform, hPCA, hRaster, hPETH);
			obj.GUIBusy(hFigure, false);
		end

		function PlotAllClusters(obj, channel, varargin)
			p = inputParser;
			addRequired(p, 'Channel', @isnumeric);
			addParameter(p, 'Clusters', [], @isnumeric);
			addParameter(p, 'ReferenceCluster', [], @isnumeric);
			addParameter(p, 'WaveformWindow', [], @isnumeric);
			addParameter(p, 'WaveformYLim', [-200, 200], @isnumeric);
			addParameter(p, 'FontSize', 8, @isnumeric);
			addParameter(p, 'FrameRate', 0, @isnumeric);
			addParameter(p, 'PlotMean', true, @islogical);
			addParameter(p, 'Fig', []);
			addParameter(p, 'FigMain', []);
			parse(p, channel, varargin{:});
			iChannel 			= p.Results.Channel;
			clusters 			= p.Results.Clusters;
			referenceCluster 	= p.Results.ReferenceCluster;
			waveformWindow 		= p.Results.WaveformWindow;
			waveformYLim		= p.Results.WaveformYLim;
			fontSize 			= p.Results.FontSize;
			frameRate 			= p.Results.FrameRate;
			plotMean 			= p.Results.PlotMean;
			hFigure 			= p.Results.Fig;
			hFigureMain 		= p.Results.FigMain;

			if isempty(clusters)
				clusters = unique(obj.Spikes(channel).Cluster.Classes);
			end
			numClusters = length(clusters);
			numCols = 4;

			xPadding = 0.03*numCols/numClusters;
			yPadding = 0.03;
			wAx = (1 - (numClusters + 1)*xPadding)/numClusters;
			hAx = (1 - 3*yPadding)/2;
			ySlider = 0.025;

			if ~isempty(hFigureMain) && isfield(hFigureMain.UserData, 'PlotMean')
				plotMean = hFigureMain.UserData.PlotMean;
			end

			if isempty(hFigure)
				hFigure	= figure('Units', 'Normalized', 'Position', [0, 0.1, 1, 0.7], 'GraphicsSmoothing', 'off');
				hPanel = uipanel('Parent', hFigure, 'Position', [0, ySlider, numClusters/numCols, 1 - ySlider]);
				if numClusters > numCols
					hSlider = uicontrol('Style', 'slider', 'Units', 'Normalized', 'Position', [0, 0, 1, ySlider],...
						'Min', 0, 'Max', numClusters/numCols - 1, 'SliderStep', (numClusters/numCols - 1)./[numClusters, numClusters/numCols],...
						'Callback', {@(hSlider, ~) set(hPanel, 'Position', [-hSlider.Value, ySlider, numClusters/numCols, 1 - ySlider])});
				end
				hWaveform = gobjects(1, numClusters);
				hFeature = gobjects(1, numClusters);
				for iAx = 1:numClusters
					hWaveform(iAx) = axes(hPanel, 'Units', 'Normalized', 'OuterPosition', [iAx*xPadding + (iAx - 1)*wAx, 2*yPadding + hAx, wAx, hAx]);
					hFeature(iAx) = axes(hPanel, 'Units', 'Normalized', 'OuterPosition', [iAx*xPadding + (iAx - 1)*wAx, yPadding, wAx, hAx]);
					obj.PlotWaveforms(channel, 'Clusters', clusters, 'SelectedCluster', iAx, 'ReferenceCluster', referenceCluster,...
						'WaveformWindow', waveformWindow, 'YLim', waveformYLim, 'FrameRate', frameRate, 'PlotMean', plotMean,...
						'Ax', [hFeature(iAx), hWaveform(iAx)]);
				end
			end
		end

		function ReplotChannel(obj, iChannel, p, hFigure, hWaveform, hPCA, hRaster, hPETH)
			reference 		= p.Results.Reference;
			event 			= p.Results.Event;
			exclude 		= p.Results.Exclude;
			waveformWindow 	= p.Results.WaveformWindow;
			extendedWindow 	= p.Results.ExtendedWindow;
			minTrialLength 	= p.Results.MinTrialLength;
			bins 			= p.Results.Bins;
			binMethod 		= p.Results.BinMethod;
			spikeRateWindow = p.Results.SpikeRateWindow;
			rasterXLim 		= p.Results.RasterXLim;
			waveformYLim	= p.Results.WaveformYLim;
			frameRate 		= p.Results.FrameRate;

			clusters = hFigure.UserData.SelectedClusters;
			plotMean = hFigure.UserData.PlotMean;
			referenceCluster = hFigure.UserData.ReferenceCluster;

			% Stop refresh if frameRate > 0
			if isfield(hWaveform.UserData, 'hTimer')
				if isvalid(hWaveform.UserData.hTimer)
					stop(hWaveform.UserData.hTimer);
					delete(hWaveform.UserData.hTimer);
				end
			end

			% Clear axes
			cla(hPCA)
			cla(hWaveform)
			cla(hRaster)
			cla(hPETH)

			if length(clusters) == 1
				selectedCluster = clusters;
				clusters = [];
				if selectedCluster == referenceCluster
					referenceCluster = [];
				end
			else
				selectedCluster = [];
				referenceCluster = [];
			end

			% Replot newly selected clusters
			obj.PlotWaveforms(iChannel, 'Clusters', clusters, 'WaveformWindow', waveformWindow,...
				'YLim', waveformYLim, 'FrameRate', frameRate, 'PlotMean', plotMean,...
				'ReferenceCluster', referenceCluster, 'SelectedCluster', selectedCluster,...
				'Ax', [hPCA, hWaveform]);
			obj.Raster(iChannel, reference, event, exclude, 'Clusters', clusters,...
				'AlignTo', 'Event', 'ExtendedWindow', extendedWindow, 'XLim', rasterXLim,...
				'Ax', hRaster, 'Sort', true);
			obj.PETH(iChannel, reference, event, exclude, 'Clusters', clusters,...
				'MinTrialLength', minTrialLength, 'Bins', bins, 'BinMethod', binMethod,...
				'SpikeRateWindow', spikeRateWindow, 'ExtendedWindow', extendedWindow,...
				'Ax', hPETH);
		end

		function GUIPlotClusters(obj, hButton, evnt, iChannel, p, hFigure, hWaveform, hPCA, hRaster, hPETH)
			obj.GUIBusy(hFigure, true);
			liststr = cellfun(@num2str, num2cell(unique(obj.Spikes(iChannel).Cluster.Classes)), 'UniformOutput', false);

			[clusters, ok] = listdlg(...
				'PromptString', 'Plot clusters:',...
				'SelectionMode', 'multiple',...
				'OKString', 'Plot',...
				'ListString', liststr,...
				'InitialValue', []);

			if (ok && ~isempty(clusters))
				clusters = cellfun(@str2num, liststr(clusters));
				% Replot clusters
				hFigure.UserData.SelectedClusters = clusters;
				obj.ReplotChannel(iChannel, p, hFigure, hWaveform, hPCA, hRaster, hPETH);
			end
			obj.GUIBusy(hFigure, false);
		end

		function GUIRemoveClusters(obj, hButton, evnt, iChannel, p, hFigure, hWaveform, hPCA, hRaster, hPETH)
			obj.GUIBusy(hFigure, true);
			liststr = cellfun(@num2str, num2cell(unique(obj.Spikes(iChannel).Cluster.Classes)), 'UniformOutput', false);

			[clusters, ok] = listdlg(...
				'PromptString', 'Remove clusters:',...
				'SelectionMode', 'multiple',...
				'OKString', 'Remove',...
				'ListString', liststr,...
				'InitialValue', []);

			if (ok && ~isempty(clusters))
				clusters = cellfun(@str2num, liststr(clusters));
				answer = questdlg(...
					['Permanently remove selected clusters (', mat2str(clusters), ')?'],...
					'Remove Cluster(s)',...
					'Remove', 'Cancel',...
					'Cancel');
				if strcmpi(answer, 'Remove')
					% Remove selected clusters
					obj.ClusterRemove(iChannel, clusters);
					hFigure.UserData.SelectedClusters = unique(obj.Spikes(iChannel).Cluster.Classes);

					% Replot clusters
					obj.ReplotChannel(iChannel, p, hFigure, hWaveform, hPCA, hRaster, hPETH);
				end
			end
			obj.GUIBusy(hFigure, false);
		end

		function GUIDecimateClusters(obj, hButton, evnt, iChannel, p, hFigure, hWaveform, hPCA, hRaster, hPETH)
			obj.GUIBusy(hFigure, true);
			liststr = cellfun(@num2str, num2cell(unique(obj.Spikes(iChannel).Cluster.Classes)), 'UniformOutput', false);

			[clusters, ok] = listdlg(...
				'PromptString', 'Decimate clusters:',...
				'SelectionMode', 'multiple',...
				'OKString', 'Decimate',...
				'ListString', liststr,...
				'InitialValue', []);

			if (ok && ~isempty(clusters))
				clusters = cellfun(@str2num, liststr(clusters));
				answer = questdlg(...
					['Permanently decimate selected clusters (', mat2str(clusters), ')?'],...
					'Decimate Cluster(s)',...
					'Decimate', 'Cancel',...
					'Cancel');
				if strcmpi(answer, 'Decimate')
					% Decimate selected clusters
					for iCluster = clusters
						obj.ClusterDecimate(iChannel, iCluster);
					end

					% Replot clusters
					obj.ReplotChannel(iChannel, p, hFigure, hWaveform, hPCA, hRaster, hPETH);
				end
			end
			obj.GUIBusy(hFigure, false);
		end

		function GUIMergeClusters(obj, hButton, evnt, iChannel, p, hFigure, hWaveform, hPCA, hRaster, hPETH)
			obj.GUIBusy(hFigure, true);
			liststr = cellfun(@num2str, num2cell(unique(obj.Spikes(iChannel).Cluster.Classes)), 'UniformOutput', false);

			[clusters, ok] = listdlg(...
				'PromptString', 'Merge clusters:',...
				'SelectionMode', 'multiple',...
				'OKString', 'Merge',...
				'ListString', liststr,...
				'InitialValue', []);

			if (ok && ~isempty(clusters))
				clusters = cellfun(@str2num, liststr(clusters));
				answer = questdlg(...
					['Merge selected clusters (', mat2str(clusters), ')?'],...
					'Merge Clusters',...
					'Merge', 'Cancel',...
					'Cancel');
				if strcmpi(answer, 'Merge')
					% Merge selected clusters
					allClusters = unique(obj.Spikes(iChannel).Cluster.Classes);
					mergeList = [{clusters}; num2cell(allClusters(~ismember(allClusters, clusters)))];
					obj.ClusterMerge(iChannel, mergeList);
					hFigure.UserData.SelectedClusters = unique(obj.Spikes(iChannel).Cluster.Classes);

					% Replot clusters
					obj.ReplotChannel(iChannel, p, hFigure, hWaveform, hPCA, hRaster, hPETH);
				end
			end
			obj.GUIBusy(hFigure, false);
		end

		function GUIRecluster(obj, hButton, evnt, iChannel, p, hFigure, hWaveform, hPCA, hRaster, hPETH)
			obj.GUIBusy(hFigure, true);
			liststr = cellfun(@num2str, num2cell(unique(obj.Spikes(iChannel).Cluster.Classes)), 'UniformOutput', false);

			[clusters, ok] = listdlg(...
				'PromptString', 'Recluster clusters:',...
				'SelectionMode', 'multiple',...
				'OKString', 'Recluster',...
				'ListString', liststr,...
				'InitialValue', []);

			if (ok && ~isempty(clusters))
				clusters = cellfun(@str2num, liststr(clusters));
				answer = questdlg(...
					['Recluster selected clusters (', mat2str(clusters), ')?'],...
					'Recluster Clusters',...
					'Recluster', 'Cancel',...
					'Cancel');
				if strcmpi(answer, 'Recluster')
					% Recluster selected clusters
					clusterMethod = obj.Spikes(iChannel).Cluster.Method;
					obj.Cluster(iChannel, 'Clusters', clusters, 'Method', clusterMethod);
					hFigure.UserData.SelectedClusters = unique(obj.Spikes(iChannel).Cluster.Classes);

					% Replot clusters
					obj.ReplotChannel(iChannel, p, hFigure, hWaveform, hPCA, hRaster, hPETH);
				end
			end
			obj.GUIBusy(hFigure, false);
		end

		function GUISelRef(obj, hButton, evnt, iChannel, p, hFigure, hWaveform, hPCA, hRaster, hPETH)
			obj.GUIBusy(hFigure, true);
			liststr = cellfun(@num2str, num2cell(unique(obj.Spikes(iChannel).Cluster.Classes)), 'UniformOutput', false);

			[referenceCluster, ok] = listdlg(...
				'PromptString', 'Select cluster as reference:',...
				'SelectionMode', 'single',...
				'OKString', 'Reference',...
				'CancelString', 'No Reference',...
				'ListString', liststr,...
				'InitialValue', []);

			if ok
				hFigure.UserData.ReferenceCluster = cellfun(@str2num, liststr(referenceCluster));
			else
				hFigure.UserData.ReferenceCluster = [];
			end
			obj.ReplotChannel(iChannel, p, hFigure, hWaveform, hPCA, hRaster, hPETH);

			obj.GUIBusy(hFigure, false);
		end

		function GUIPlotMean(obj, hButton, evnt, iChannel, p, hFigure, hWaveform, hPCA, hRaster, hPETH)
			obj.GUIBusy(hFigure, true);
			if logical(hButton.Value) ~= hFigure.UserData.PlotMean
				hFigure.UserData.PlotMean = logical(hButton.Value);
				obj.ReplotChannel(iChannel, p, hFigure, hWaveform, hPCA, hRaster, hPETH);
			end
			obj.GUIBusy(hFigure, false);
		end

		function GUIPlotAllClusters(obj, hButton, evnt, iChannel, hFigure)
			obj.GUIBusy(hFigure, true);
			liststr = cellfun(@num2str, num2cell(unique(obj.Spikes(iChannel).Cluster.Classes)), 'UniformOutput', false);

			[referenceCluster, ok] = listdlg(...
				'PromptString', 'Select cluster as reference:',...
				'SelectionMode', 'single',...
				'OKString', 'Reference',...
				'ListString', liststr,...
				'InitialValue', []);

			if ok
				referenceCluster = cellfun(@str2num, liststr(referenceCluster));
				% Expand and plot all clusters
				obj.PlotAllClusters(iChannel, 'Clusters', hFigure.UserData.SelectedClusters, 'ReferenceCluster', referenceCluster, 'FigMain', hFigure);
			end
			obj.GUIBusy(hFigure, false);
		end

		function GUISavePlot(obj, hButton, evnt, hFigure)
			hButtons = findobj(hFigure.Children, 'Type', 'uicontrol');
			set(hButtons, 'Visible', 'off');
			print(hFigure, '-clipboard', '-dbitmap')
			set(hButtons, 'Visible', 'on');
		end

		function GUIDeleteChannel(obj, hButton, evnt, iChannel, nextChn, hFigure)
			obj.GUIBusy(hFigure, true);
			answer = questdlg(...
				['Permanently delete current channel (', num2str(iChannel), ')?'],...
				'Delete channel',...
				'Delete', 'Cancel',...
				'Cancel');
			if strcmpi(answer, 'Delete')
				% Delete channel
				for field = fieldnames(obj.Spikes)'
					obj.Spikes(iChannel).(field{1}) = [];
				end				
				if nextChn == iChannel
					close(hFigure)
				else
					obj.PlotChannel(nextChn, 'Fig', hFigure)
				end
			end
			obj.GUIBusy(hFigure, false);
		end

		function GUIBusy(obj, hFigure, busy)
			hButtons = findobj(hFigure.Children, 'Type', 'uicontrol');

			if busy
				set(hButtons, 'Enable', 'off')
			else
				set(hButtons, 'Enable', 'on')
			end
		end

		function PlotAllChannels(obj, varargin)
			p = inputParser;
			addParameter(p, 'Channels', [], @isnumeric);
			addParameter(p, 'PercentShown', 5, @isnumeric); % What percentage of waveforms are plotted (0 - 100)
			addParameter(p, 'MaxShown', 200, @isnumeric);
			addParameter(p, 'Fontsize', 8, @isnumeric);
			addParameter(p, 'PlotMethod', 'all', @ischar); % 'all', 'mean'
			addParameter(p, 'YLim', [-400, 400], @isnumeric); % 'all', 'mean'
			parse(p, varargin{:});
			channels 		= p.Results.Channels;
			percentShown 	= p.Results.PercentShown;
			maxShown 		= p.Results.MaxShown;
			fontSize 		= p.Results.Fontsize;
			plotMethod 		= p.Results.PlotMethod;
			yRange 			= p.Results.YLim;

			if isempty(channels)
				channels = [obj.Spikes.Channel];
			end

			expName = obj.GetExpName();

			hFigure	= figure('Units', 'Normalized', 'Position', [0, 0, 1, 1], 'Name', expName, 'DefaultAxesFontSize', fontSize,...
				'GraphicsSmoothing', 'off');
			hFigure.UserData.SelectedChannels = false(32, 1);
			hAxes = gobjects(1, 35);

			for iChannel = channels
				hAxes(iChannel)	= subplot(5, 7, iChannel);
				hAxes(iChannel).UserData.Channel = iChannel;
				hAxes(iChannel).ButtonDownFcn = @obj.OnAxesClicked;
				xlabel(hAxes(iChannel), 'Time (ms)');
				ylabel(hAxes(iChannel), 'Voltage (\muV)');
				title(hAxes(iChannel), ['Channel ', num2str(iChannel)]);
				clusterID = obj.Spikes(iChannel).Cluster.Classes;
				for iCluster = unique(nonzeros(clusterID))'
					thisWaveforms = obj.Spikes(iChannel).Waveforms(clusterID==iCluster, :);
					switch lower(plotMethod)
						case 'all'
							thisWaveforms = thisWaveforms(1:ceil(100/percentShown):end, :);
							if size(thisWaveforms, 1) > maxShown
								thisWaveforms = thisWaveforms(randperm(size(thisWaveforms, 1), maxShown), :);
							end
						case 'mean'
							thisWaveforms = mean(thisWaveforms, 1);
					end
					[thisColor, thisStyle] = TetrodeRecording.GetColorAndStyle(iCluster);
					line(hAxes(iChannel), obj.Spikes(iChannel).WaveformTimestamps, thisWaveforms, 'LineStyle', thisStyle, 'Color', thisColor);
				end
				ylim(hAxes(iChannel), yRange);
				drawnow
			end

			suptitle(expName);
		end

		function OnAxesClicked(obj, hAxes, evnt)
			channel = hAxes.UserData.Channel;
			button = evnt.Button; % 1, 2, 3 (LMB, MMB, RMB click)
			hFigure = hAxes.Parent;

			switch button
				case 1
					hFigure.UserData.SelectedChannels(channel) = true;
					hAxes.Box 						= 'on';
					hAxes.LineWidth 				= 2;
					hAxes.XColor 					= 'r';
					hAxes.YColor 					= 'r';
					hAxes.Title.Color 				= 'r';
					hAxes.TitleFontSizeMultiplier 	= 1.6;
					selectedChannels = transpose(find(hFigure.UserData.SelectedChannels));
					obj.SelectedChannels = selectedChannels;
					disp(['Selected channels: ', mat2str(selectedChannels)])
					drawnow
				case 3
					hFigure.UserData.SelectedChannels(channel) = false;
					hAxes.Box 						= 'off';
					hAxes.LineWidth 				= 0.5;
					hAxes.XColor 					= 'k';
					hAxes.YColor 					= 'k';
					hAxes.Title.Color 				= 'k';
					hAxes.TitleFontSizeMultiplier 	= 1.1;
					selectedChannels = transpose(find(hFigure.UserData.SelectedChannels));
					obj.SelectedChannels = selectedChannels;
					disp(['Selected channels: ', mat2str(selectedChannels)])
					drawnow
			end
		end

		% Sort spikes and digital events into trial structure
		function [spikes, trials] = GetSpikesByTrial(obj, channel, reference, event, extendedWindow, clusters)
			if nargin < 5
				extendedWindow = [-1, 0];
			end

			if nargin < 6
				clusters = [];
			end

			edges = [reference + extendedWindow(1); event + extendedWindow(2)];
			edges = edges(:);
			if isempty(edges)
				spikes = [];
				trials = [];
			else
				if ~isempty(clusters)
					spikes = obj.Spikes(channel).Timestamps(ismember(obj.Spikes(channel).Cluster.Classes, clusters));
				else
					spikes = obj.Spikes(channel).Timestamps;
				end
				[~, ~, bins] = histcounts(spikes, edges);
				oddBins = rem(bins, 2) ~= 0;	% Spikes in odd bins occur between reference and event, should keep these spikes
				spikes = spikes(oddBins);
				trials = (bins(oddBins) + 1)/2;
			end
		end
	end

	methods (Static)
		function previewObj = BatchPreview()
			previewObj = TetrodeRecording();
			dirs = uipickfiles('Prompt', 'Select (multiple) folders...');
			dirs = dirs(isfolder(dirs));
			for iDir = 1:length(dirs)
				files = dir([dirs{iDir}, '\*.rhd']);
				stepSize = round(length(files)/6);
				files = {files(stepSize:stepSize:5*stepSize).name};
				previewObj(iDir) = TetrodeRecording();
				previewObj(iDir).Path = [dirs{iDir}, '\'];
				previewObj(iDir).Files = files;
				previewObj(iDir).Preview('HideResults', true);
			end
			for iDir = 1:length(dirs)
				previewObj(iDir).PlotAllChannels();
				previewObj(iDir).ClearCache();
			end
			TetrodeRecording.RandomWords();
		end

		function BatchProcess(previewObj, varargin)
			p = inputParser;
			addRequired(p, 'Obj', @(x) isa(x, 'TetrodeRecording'));
			addParameter(p, 'ChunkSize', 10, @isnumeric);
			addParameter(p, 'NumSigmas', 4, @isnumeric);
			addParameter(p, 'WaveformWindow', [-1, 1], @isnumeric);
			addParameter(p, 'FeatureMethod', 'WaveletTransform', @ischar);
			addParameter(p, 'ClusterMethod', 'kmeans', @ischar);
			addParameter(p, 'Dimension', 10, @isnumeric);
			addParameter(p, 'Prefix', 'tr_', @ischar);
			parse(p, previewObj, varargin{:});
			previewObj 		= p.Results.Obj;
			chunkSize 		= p.Results.ChunkSize;
			numSigmas 		= p.Results.NumSigmas;
			waveformWindow 	= p.Results.WaveformWindow;
			featureMethod 	= p.Results.FeatureMethod;
			clusterMethod 	= p.Results.ClusterMethod;
			dimension 		= p.Results.Dimension;
			prefix 			= p.Results.Prefix;

			selectedChannels = {previewObj.SelectedChannels};
			allPaths = {previewObj.Path};
			for iDir = 1:length(selectedChannels)
				channels = selectedChannels{iDir};
				if ~isempty(channels)
					TetrodeRecording.ProcessFolder(allPaths{iDir}, chunkSize, channels, numSigmas, waveformWindow, featureMethod, clusterMethod, dimension, prefix);
				end
			end
			TetrodeRecording.RandomWords();
		end

		function ProcessFolder(thisPath, chunkSize, channels, numSigmas, waveformWindow, featureMethod, clusterMethod, dimension, prefix)
			tr = TetrodeRecording();
			tr.Path = thisPath;
			files = dir([tr.Path, '*.rhd']);
			tr.Files = {files.name};
			tr.ReadFiles(chunkSize, 'DigitalDetect', true);
			tr.SpikeDetect(channels, 'NumSigmas', numSigmas, 'WaveformWindow', waveformWindow);
			tr.ReadFiles(chunkSize, 'Chunks', 'remaining', 'SpikeDetect', true, 'DigitalDetect', true);
			tr.SpikeSort(channels, 'FeatureMethod', featureMethod, 'ClusterMethod', clusterMethod, 'Dimension', dimension);
			tr.ClearCache();
			TetrodeRecording.BatchSave(tr, 'Prefix', prefix, 'DiscardData', false, 'MaxChannels', 5);
		end

		function tr = BatchLoad()
			files = uipickfiles('Prompt', 'Select .mat files containing TetrodeRecording objects to load...', 'Type', {'*.mat', 'MAT-files'});
			for iFile = 1:length(files)
				S(iFile) = load(files{iFile}, 'tr');
			end
			% Merge multi-part files
			for iFile = 1:length(files)
				part = S(iFile).tr.Part;
				if nnz(part > 1) == 2
					partOne = cellfun(@(tr) strcmpi(tr.Path, S(iFile).tr.Path) && tr.Part(1) == 1, {S.tr});
					channels = [S(iFile).tr.Spikes.Channel];
					S(partOne).tr.Spikes(channels) = S(iFile).tr.Spikes(channels);
					S(partOne).tr.Part(2) = S(partOne).tr.Part(2) - 1;
					S(iFile).tr.Part = [];
				end
			end

			partOne = cellfun(@(tr) ~isempty(tr.Part), {S.tr});

			tr = [S(partOne).tr];
		end

		function BatchSave(TR, varargin)
			p = inputParser;
			addParameter(p, 'Prefix', '', @ischar);
			addParameter(p, 'DiscardData', false, @islogical);
			addParameter(p, 'MaxChannels', [], @isnumeric);
			parse(p, varargin{:});
			prefix = p.Results.Prefix;
			discardData = p.Results.DiscardData;
			maxChannels = p.Results.MaxChannels;

			for iTr = 1:length(TR)
				tr = TR(iTr);
				expName = tr.GetExpName();
				if discardData
					tr.Spikes = [];
					tr.DigitalEvents = [];
				end
				if ~isfolder([tr.Path, '..\SpikeSort'])
					mkdir([tr.Path, '..\SpikeSort'])
				end
				partition = false;
				if ~isempty(maxChannels)
					allChannels = [tr.Spikes.Channel];
					numParts = ceil(length(allChannels)/maxChannels);
					if numParts > 1
						partition = true;
					end
				end

				if partition
					spikes = tr.Spikes;
					try
						for iPart = 1:numParts
							tr.Part = [iPart, numParts];
							for iChannel = 1:length(tr.Spikes)
								for field = fieldnames(tr.Spikes)'
									tr.Spikes(iChannel).(field{1}) = [];
								end
							end
							for iChannel = allChannels((iPart - 1)*maxChannels + 1:min(length(allChannels), iPart*maxChannels))
								tr.Spikes(iChannel) = spikes(iChannel);
							end
							file = [tr.Path, '..\SpikeSort\', prefix, expName, '(', num2str(iPart), ')', '.mat'];
							save(file, 'tr');
						end
					catch ME
						tr.Spikes = spikes;
						rethrow(ME)
					end
					tr.Spikes = spikes;
				else
					file = [tr.Path, '..\SpikeSort\', prefix, expName, '.mat'];
					save(file, 'tr');					
				end
			end
		end

		function [thisColor, thisStyle] = GetColorAndStyle(iClass)
			colors = 'rgbcmyk';
			styles = {'-', '--', ':', '-.'};

			iStyle = ceil(iClass/length(colors));
			iColor = mod(iClass, length(colors));
			if iColor == 0
				iColor = length(colors);
			end
			thisColor = colors(iColor);
			thisStyle = styles{iStyle};
		end

		function a = ReadQString(fid)
			% Read Qt style QString.  The first 32-bit unsigned number indicates the length of the string (in bytes).  If this number equals 0xFFFFFFFF, the string is null.
			a = '';
			length = fread(fid, 1, 'uint32');
			if length == hex2num('ffffffff')
				return;
			end
			% convert length from bytes to 16-bit Unicode words
			length = length / 2;

			for i=1:length
				a(i) = fread(fid, 1, 'uint16');
			end
		end

		function [eventOn, eventOff] = FindEdges(event, t)
			grad = diff([0, event]);
			eventOn = grad == 1;
			eventOff = grad == -1;

			if nargin >= 2
				eventOn = t(eventOn);
				eventOff = t(eventOff);
			else
				eventOn = find(eventOn);
				eventOff = find(eventOff);
			end
		end

		% Find first event after reference
		function [reference, event, iReference, iEvent] = FindFirstInTrial(reference, event, eventExclude)
			edges = [reference(1:end - 1), max(event(end), reference(end))];
			[~, ~, bins] = histcounts(event, edges);
			event = event(bins ~= 0);
			bins = nonzeros(bins);
			[iReference, iEvent] = unique(bins, 'first');
			reference = reference(iReference);
			event = event(iEvent);

			if ~isempty(eventExclude)
				% Filter out trials where mouse licked before pressing
				edges = [reference; event];
				edges = edges(:);
				[~, ~, bins] = histcounts(eventExclude, edges);
				oddBins = rem(bins, 2) ~= 0;	% Licks in odd bins occur between cue on and lever press, should exclude these trials
				toRemove = (unique(bins(oddBins)) + 1)/2;
				event(toRemove) = [];
				reference(toRemove) = [];
				iReference(toRemove) = [];
				iEvent(toRemove) = [];
			end
		end

		function OnPlotChannelRefresh(~, ~, hAxes, t, waveforms, numWaveforms, numWaveformsTotal, clusterID)
			hAxes.UserData.iWaveform = hAxes.UserData.iWaveform + 1;
			iWaveform = hAxes.UserData.iWaveform;

			if iWaveform > numWaveformsTotal
				stop(hAxes.UserData.hTimer);
				delete(hAxes.UserData.hTimer);
				return
			else
				if iWaveform == 1
					xlabel(hAxes, 'Time (ms)');
					ylabel(hAxes, 'Voltage (\muV)');
					title(hAxes, 'Waveforms');
				end
				if length(hAxes.UserData.hWaveforms) < numWaveforms
					[thisColor, thisStyle] = TetrodeRecording.GetColorAndStyle(clusterID(iWaveform));
					hAxes.UserData.hWaveforms = [hAxes.UserData.hWaveforms, line(hAxes, 'XData', t, 'YData', waveforms(iWaveform, :), 'LineStyle', thisStyle, 'Color', thisColor, 'DisplayName', ['Waveform (Cluster ', num2str(clusterID(iWaveform)), ')'])];
				else
					iHandle = mod(iWaveform, numWaveforms);
					if iHandle == 0
						iHandle = numWaveforms;
					end
					hThisWaveform = hAxes.UserData.hWaveforms(iHandle);
					hThisWaveform.YData = waveforms(iWaveform, :);
					[hThisWaveform.Color, hThisWaveform.LineStyle] = TetrodeRecording.GetColorAndStyle(clusterID(iWaveform));
				end
				drawnow
			end	
		end

		function OnKeyPress(~, evnt, hTimer)
			if isvalid(hTimer)
				if strcmpi(evnt.Key, 'space')
					switch lower(hTimer.Running)
						case 'on'
							stop(hTimer);
						case 'off'
							start(hTimer);
					end
				end
			end
		end

		function OnFigureClosed(~, ~, hTimer)
			if isvalid(hTimer)
				stop(hTimer);
				delete(hTimer);
			end
			delete(gcf);
		end

		function TTS(txt, speak)
			fprintf(1, txt);
			if nargin < 2
				speak = false;
			end
			if speak
				txt = strsplit(txt, '\');
				txt = txt{1};
				tts(txt);
			end
		end

		function RandomWords()
			words = {...
				'49 times, we fought that beast. It had a chicken head with duck feet, with a woman''s face too.',...
				'Every day I worry all day about what''s waiting in the bushes of love. Something''s waiting in the bushes for us, something''s waiting in the bushes of love.',...
				'Now run, run, run, jump! I can be a backpack while you run.',...
				'Rocking, rocking and rolling. Down to the beach I''m strolling. But the seagulls poke at my head, not fun! I said "Seagulls... mmgh! Stop it now!"',...
				'Some day when you are older you could get hit by a boulder. While you''re lying there screaming come help me please, the seagulls come poke your knees.',...
				'Even though it looks like it''s the future. It''s really a long long time ago when there were knights and they got into fights using sabres of light!"',...
				'Twenty nights in the ice is a long time when there''s hostiles on the hill. It''s not about what they want, you just gotta walk your walk.'...
			};
			TetrodeRecording.TTS([words{randi(length(words))}, '\n'], true)
		end
	end
end
