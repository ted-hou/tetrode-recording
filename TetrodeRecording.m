classdef TetrodeRecording < handle
	properties
		Files
		Path
		Notes
		FrequencyParameters
		ChannelMap
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
			obj.SelectFiles();
		end

		function ReadFiles(obj, varargin)
			p = inputParser;
			addOptional(p, 'ChunkSize', 2, @isnumeric);
			addParameter(p, 'Chunks', 'first');
			addParameter(p, 'SubstractMean', true, @islogical);
			addParameter(p, 'SpikeDetect', false, @islogical);
			addParameter(p, 'DigitalDetect', false, @islogical);
			parse(p, varargin{:});
			chunkSize = p.Results.ChunkSize;
			chunks = p.Results.Chunks;
			substractMean = p.Results.SubstractMean;
			spikeDetect = p.Results.SpikeDetect;
			digitalDetect = p.Results.DigitalDetect;

			files = obj.Files;
			numChunks = ceil(length(files)/chunkSize);
			if ischar(chunks)
				if strcmp(chunks, 'all')
					chunks = 1:numChunks;
				elseif strcmp(chunks, 'remaining')
					chunks = 2:numChunks;
				elseif strcmp(chunks, 'first')
					chunks = 1;
				end
			end

			for iChunk = chunks
				obj.ClearCache();
				TetrodeRecording.TTS(['Loading chunk ', num2str(iChunk), '/', num2str(numChunks), '...\n']);
				obj.ReadRHD(obj.Files((iChunk - 1)*chunkSize + 1:min(iChunk*chunkSize, length(obj.Files))))
				obj.MapChannels();
                if substractMean
                    obj.SubstractMean();
                end
				obj.RemoveTransient();
				if spikeDetect
					for iChannel = [obj.Spikes.Channel]
						waveformWindow = obj.Spikes(iChannel).WaveformWindow;
						waveformWindowExtended = obj.Spikes(iChannel).WaveformWindowExtended;
						threshold = obj.Spikes(iChannel).Threshold;
						obj.SpikeDetect(...
							iChannel, threshold.Trigger,...
							'ExitThreshold', threshold.Exit,...
							'ExitWindow', threshold.ExitWindow,...
							'ExclusionThreshold', threshold.Exclusion,...
							'WaveformWindow', waveformWindow,...
							'waveformWindowExtended', waveformWindowExtended,...
							'Append', true);
					end
					obj.GetDigitalData('Append', true);
				else
					obj.GetDigitalData('Append', false);
				end
			end
			if digitalDetect
				obj.GetDigitalEvents(true);
			end
		end

		function SelectFiles(obj)
			[obj.Files, obj.Path, ~] = uigetfile('*.rhd', 'Select an RHD2000 Data File', 'MultiSelect', 'on');
		end

		function ReadRHD(obj, files)
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
					TetrodeRecording.TTS(['	Reading data from file ', num2str(iFile), '/', num2str(length(files)), ' (''', files{iFile},''')...\n']);

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
			end

			% Count number of samples in all files
			obj.Amplifier.NumSamples = 0;
			obj.BoardDigIn.NumSamples = 0;
			for iFile = 1:length(files)
				obj.Amplifier.NumSamples = obj.Amplifier.NumSamples + size(objTemp(iFile).Amplifier.Timestamps, 2);
				obj.BoardDigIn.NumSamples = obj.BoardDigIn.NumSamples + size(objTemp(iFile).BoardDigIn.Timestamps, 2);
			end

			% Combine files
			tic, TetrodeRecording.TTS('All files loaded. Concatenating data...');

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
			TetrodeRecording.TTS(['	Done(', num2str(toc, 2), ' seconds).\n'])
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
			tic, TetrodeRecording.TTS('Remapping amplifier channels in tetrode order...');

			if nargin < 3
				headstageType = 'intan';
			end
			if nargin < 2
				EIBMap = [15, 13, 11, 9, 7, 5, 3, 1, 2, 4, 6, 8, 10, 12, 14, 16];
				EIBMap = [EIBMap, EIBMap + 16];
				[~, EIBMap] = sort(EIBMap);
			end

			if strcmp(headstageType, 'intan')
				HSMap = [25:32, 1:8, 24:-1:9];
			elseif strcmp(headstageType, 'open ephys')
				HSMap = [26, 25, 27, 28, 29, 30, 31, 1, 32, 3, 2, 5, 4, 7, 6, 8, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9];
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
			TetrodeRecording.TTS(['Done(', num2str(toc, 2), ' seconds).\n'])
		end

		% Convert zero-based raw channel id to remapped tetrode id
		function tetrodeChannel = MapChannelID(obj, rawChannel)
			tetrodeChannel = find(obj.ChannelMap.Tetrode == (rawChannel + 1));
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
				dilate = false;
			end
			if nargin < 4
				dilateSize = 12; % By default, lick digital event is extended 2 ms to the left and right.
			end

			% Verify if digital input sample rate is identical to amplifier sample rate
			if obj.FrequencyParameters.AmplifierSampleRate ~= obj.FrequencyParameters.BoardDigInSampleRate
				error('Board digital input has a different sampling rate from amplifier. Aborted.');
			end

			tic, TetrodeRecording.TTS('Removing lick/touch-related transients. This might take a while...');
			if dilate
				mask = false(1, obj.Amplifier.NumSamples);
				dilateSE = ones(round(obj.FrequencyParameters.AmplifierSampleRate*2*dilateSize/1000) + 1);
				for thisChannel = digChannel
					thisChannel = find([obj.BoardDigIn.Channels.NativeOrder] == thisChannel);
					mask = mask | logical(imdilate(obj.BoardDigIn.Data(thisChannel, :), dilateSE));
				end
			else
				for thisChannel = digChannel
					thisChannel = find([obj.BoardDigIn.Channels.NativeOrder] == thisChannel);
					mask = logical(obj.BoardDigIn.Data(thisChannel, :));
				end
			end
			
			obj.Amplifier.Data(:, mask) = 0;
			TetrodeRecording.TTS(['Done(', num2str(toc, 2), ' seconds).\n'])
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

			if strcmp(indexType, 'SampleIndex')
				sampleIndex = index;
				% Only interpolate if sample index in non-integer
				if sum(rem(sampleIndex, 1) == 0) == length(sampleIndex)
					timestamps = obj.Amplifier.Timestamps(sampleIndex);
				else
					timestamps = interp1(1:length(obj.Amplifier.Timestamps), obj.Amplifier.Timestamps, sampleIndex, 'linear');
				end
			elseif strcmp(indexType, 'Timestamps')
				% Always interpolate if input index in timestamps. This will always be slower.
				timestamps = index;
				sampleIndex = interp1(obj.Amplifier.Timestamps, 1:length(obj.Amplifier.Timestamps), timestamps, 'linear');
			else
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
		function SpikeDetect(obj, channels, triggerThreshold, varargin)
			p = inputParser;
			addRequired(p, 'Channels', @isnumeric);
			addRequired(p, 'TriggerThreshold', @isnumeric);
			addParameter(p, 'ExitThreshold', NaN, @isnumeric);
			addParameter(p, 'ExitWindow', [0, NaN], @(x) isnumeric(x) && length(x) == 2);
			addParameter(p, 'ExclusionThreshold', NaN, @isnumeric);
			addParameter(p, 'WaveformWindow', [-0.5, 0.5], @isnumeric);
			addParameter(p, 'WaveformWindowExtended', [-1, 2], @isnumeric);
			addParameter(p, 'Append', false, @islogical);
			parse(p, channels, triggerThreshold, varargin{:});
			channels = p.Results.Channels;
			triggerThreshold = p.Results.TriggerThreshold;
			exitThreshold = p.Results.ExitThreshold;
			exitWindow = p.Results.ExitWindow;
			exclusionThreshold = p.Results.ExclusionThreshold;
			waveformWindow = p.Results.WaveformWindow;
			waveformWindowExtended = p.Results.WaveformWindowExtended;
			append = p.Results.Append;

			tic, TetrodeRecording.TTS('Detecting spikes...');
			sampleRate = obj.FrequencyParameters.AmplifierSampleRate/1000;
			if isnan(exitWindow(2))
				exitWindow(2) = waveformWindow(2);
			end

			for iChannel = channels
				% Find triggerThreshold-crossing events
				[~, sampleIndex] = findpeaks(sign(triggerThreshold)*obj.Amplifier.Data(iChannel, :), 'MinPeakHeight', abs(triggerThreshold));

				% Extract waveforms
				[waveforms, t] = obj.GetWaveforms(iChannel, waveformWindow, sampleIndex, 'IndexType', 'SampleIndex');
				if ~isnan(exitThreshold)
					exitWaveforms = obj.GetWaveforms(iChannel, exitWindow, sampleIndex, 'IndexType', 'SampleIndex');
				end

				% Align waveforms to peak
				numWaveforms = size(waveforms, 1);
				i = sampleRate*t;

				% Figure out via consensus if max amplitude is positive or negative.
				[~, maxIndex] = max(abs(waveforms), [], 2);
				peak = zeros(numWaveforms, 1);
				for iWaveform = 1:numWaveforms
					peak(iWaveform) = waveforms(iWaveform, maxIndex(iWaveform));
				end
				if sum(peak > 0) >= sum(peak < 0)
					[~, maxIndex] = max(waveforms, [], 2);
				else
					[~, maxIndex] = max(-waveforms, [], 2);
				end
				alignmentShift = i(maxIndex);

				[waveforms, t] = obj.GetWaveforms(iChannel, waveformWindow, sampleIndex + alignmentShift, 'IndexType', 'SampleIndex');
				[waveformsExtended, tExtended, timestamps, sampleIndex] = obj.GetWaveforms(iChannel, waveformWindowExtended, sampleIndex + alignmentShift, 'IndexType', 'SampleIndex');

				% Filter waveforms by optional criteria
				includeWaveform = true(numWaveforms, 1);
				for iWaveform = 1:numWaveforms
					% (Optional) Only keep waveforms that cross a second exitThreshold after reaching the first triggerThreshold
					if ~isnan(exitThreshold)
						if sum(sign(exitWaveforms(iWaveform, :) - exitThreshold) == sign(exitThreshold)) == 0
							includeWaveform(iWaveform) = false;
							continue
						end
					end
					% (Optional) Filter out waveforms with max amplitude exceeding exclusionThreshold
					if ~isnan(exclusionThreshold)
						if sum(sign(waveforms(iWaveform, :) - exclusionThreshold) == sign(exclusionThreshold)) > 0
							includeWaveform(iWaveform) = false;
							continue
						end
					end
				end

				% Remove invalid waveforms
				waveforms = waveforms(includeWaveform, :);
				waveformsExtended = waveformsExtended(includeWaveform, :);
				sampleIndex = sampleIndex(includeWaveform);
				timestamps = timestamps(includeWaveform);

				% Store data
				if ~append
					obj.Spikes(iChannel).Channel = iChannel;

					obj.Spikes(iChannel).SampleIndex = sampleIndex;
					obj.Spikes(iChannel).Timestamps = timestamps;
					obj.Spikes(iChannel).Waveforms = waveforms;
					obj.Spikes(iChannel).WaveformsExtended = waveformsExtended;

					obj.Spikes(iChannel).WaveformTimestamps = t;
					obj.Spikes(iChannel).WaveformTimestampsExtended = tExtended;
					obj.Spikes(iChannel).WaveformWindow = waveformWindow;
					obj.Spikes(iChannel).WaveformWindowExtended = waveformWindowExtended;
					obj.Spikes(iChannel).Threshold.Trigger = triggerThreshold;
					obj.Spikes(iChannel).Threshold.Exit = exitThreshold;
					obj.Spikes(iChannel).Threshold.ExitWindow = exitWindow;
					obj.Spikes(iChannel).Threshold.Exclusion = exclusionThreshold;
					obj.Spikes(iChannel).AlignmentShift = alignmentShift;
				else
					obj.Spikes(iChannel).SampleIndex = [obj.Spikes(iChannel).SampleIndex, sampleIndex + length(obj.DigitalEvents.Timestamps)];
					obj.Spikes(iChannel).Timestamps = [obj.Spikes(iChannel).Timestamps, timestamps];
					obj.Spikes(iChannel).Waveforms = [obj.Spikes(iChannel).Waveforms; waveforms];
					obj.Spikes(iChannel).WaveformsExtended = [obj.Spikes(iChannel).WaveformsExtended; waveformsExtended];
				end

				TetrodeRecording.TTS(['Done(', num2str(toc, 2), ' seconds).\n'])
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

			tic, TetrodeRecording.TTS('Extracting digital events...');
			% Get digital signals
			[obj.DigitalEvents.CueOn, obj.DigitalEvents.CueOff] 		= TetrodeRecording.FindEdges(obj.DigitalEvents.Data(1, :), obj.DigitalEvents.Timestamps);
			[obj.DigitalEvents.PressOn, obj.DigitalEvents.PressOff] 	= TetrodeRecording.FindEdges(obj.DigitalEvents.Data(2, :), obj.DigitalEvents.Timestamps);
			[obj.DigitalEvents.LickOn, obj.DigitalEvents.LickOff] 		= TetrodeRecording.FindEdges(obj.DigitalEvents.Data(3, :), obj.DigitalEvents.Timestamps);
			[obj.DigitalEvents.RewardOn, obj.DigitalEvents.RewardOff] 	= TetrodeRecording.FindEdges(obj.DigitalEvents.Data(4, :), obj.DigitalEvents.Timestamps);

			if clearCache
				obj.DigitalEvents.Timestamps = [];
				obj.DigitalEvents.Data = [];
			end

			TetrodeRecording.TTS(['Done(', num2str(toc, 2), ' seconds).\n'])
		end

		% This compresses data by ~ 20 times
		function ClearCache(obj)
			obj.Amplifier = [];
			obj.BoardDigIn = [];
			mem = memory();
			TetrodeRecording.TTS(['Cached data cleared. System memory: ', num2str(round(mem.MemUsedMATLAB/1024^2)), ' MB used (', num2str(round(mem.MemAvailableAllArrays/1024^2)), ' MB available).\n']);
		end

		% SpikeSort: PCA & Cluster
		function SpikeSort(obj, channels, k, varargin)
			p = inputParser;
			addRequired(p, 'Channels', @isnumeric);
			addRequired(p, 'K', @isnumeric);
			addParameter(p, 'Method', 'kmeans', @ischar);
			addParameter(p, 'Clusters', [], @isnumeric);
			addParameter(p, 'WaveformWindow', [], @isnumeric);
			addParameter(p, 'HideResults', false, @islogical);
			parse(p, channels, k, varargin{:});
			channels = p.Results.Channels;
			k = p.Results.K;
			method = p.Results.Method;
			clusters = p.Results.Clusters;
			waveformWindow = p.Results.WaveformWindow;
			hideResults = p.Results.HideResults;

			for iChannel = channels
                if isfield(obj.Spikes(iChannel), 'Cluster')
                    if isempty(clusters)
                        thisClusters = unique(obj.Spikes(iChannel).Cluster);
                    else
                        thisClusters = clusters;
                    end
                else
                    thisClusters = clusters;
                end
                if isempty(waveformWindow)
                	thisWaveformWindow = obj.Spikes(iChannel).WaveformWindow;
                else
                	thisWaveformWindow = waveformWindow;
                end
				obj.RemoveNaNs(channels);
				obj.PCA(channels, thisClusters, thisWaveformWindow);
				obj.Cluster(channels, k, method, thisClusters);
				if ~hideResults
					obj.PlotChannel(iChannel, 'WaveformWindow', thisWaveformWindow);
				end
			end
		end

		function RemoveNaNs(obj, channels)
			for iChannel = channels
				iWaveformToDiscard = sum(isnan(obj.Spikes(iChannel).Waveforms), 2) > 0;
				obj.Spikes(iChannel).Waveforms(iWaveformToDiscard, :) = [];
				obj.Spikes(iChannel).WaveformsExtended(iWaveformToDiscard, :) = [];
				obj.Spikes(iChannel).Timestamps(iWaveformToDiscard) = [];
				obj.Spikes(iChannel).SampleIndex(iWaveformToDiscard) = [];
			end
		end

		function PCA(obj, channels, clusters, waveformWindow)
			tic, TetrodeRecording.TTS('Principal component analysis...');
			for iChannel = channels
                if ~isempty(clusters)
                    inCluster = ismember(obj.Spikes(iChannel).Cluster, clusters);
                else
                    inCluster = true(size(obj.Spikes(iChannel).Timestamps));
                end
                inWindow = obj.Spikes(iChannel).WaveformTimestamps >= waveformWindow(1) & obj.Spikes(iChannel).WaveformTimestamps <= waveformWindow(2);
				obj.Spikes(iChannel).PCA = [];
				[obj.Spikes(iChannel).PCA.Coeff, obj.Spikes(iChannel).PCA.Score(inCluster, :)] = pca(obj.Spikes(iChannel).Waveforms(inCluster, inWindow));
			end
			TetrodeRecording.TTS(['Done(', num2str(toc, 2), ' seconds).\n'])
		end

		function Cluster(obj, channels, k, method, clusters)
			tic, TetrodeRecording.TTS('Clustering...');
			for iChannel = channels
                if isfield(obj.Spikes(iChannel), 'Cluster')
                    if isempty(obj.Spikes(iChannel).Cluster)
                        inCluster = true(size(obj.Spikes(iChannel).Timestamps));
                    else
                        inCluster = ismember(obj.Spikes(iChannel).Cluster, clusters);
                    end
                else
                    inCluster = true(size(obj.Spikes(iChannel).Timestamps));
                end
				
				obj.Spikes(iChannel).Cluster = [];
				if strcmp(method, 'kmeans')
					obj.Spikes(iChannel).Cluster(inCluster) = kmeans(obj.Spikes(iChannel).PCA.Score(inCluster, :), k);
				elseif strcmp(method, 'gaussian')
					gm = fitgmdist(obj.Spikes(iChannel).PCA.Score(inCluster, :), k);
					obj.Spikes(iChannel).Cluster(inCluster) = cluster(gm, obj.Spikes(iChannel).PCA.Score(inCluster, :));
				else
					error('Unrecognized clustering method. Must be ''kmeans'' or ''gaussian''.')			
				end
			end
			TetrodeRecording.TTS(['Done(', num2str(toc, 2), ' seconds).\n'], toc >= 8)
		end

		function ClusterMerge(obj, channel, mergeList, varargin)
			p = inputParser;
			addRequired(p, 'Channel', @isnumeric);
			addRequired(p, 'MergeList', @iscell);
			addParameter(p, 'HideResults', false, @islogical);
			parse(p, channel, mergeList, varargin{:});
			channel = p.Results.Channel;
			mergeList = p.Results.MergeList;
			hideResults = p.Results.HideResults;

			newCluster = NaN(size(obj.Spikes(channel).Cluster));
			for iNewCluster = 1:length(mergeList)
				newCluster(ismember(obj.Spikes(channel).Cluster, mergeList{iNewCluster})) = iNewCluster;
			end
			obj.Spikes(channel).Cluster = newCluster;

			if ~hideResults
				obj.PlotChannel(channel);
			end
		end

		function ClusterRemove(obj, channel, discardList, varargin)
			p = inputParser;
			addRequired(p, 'Channel', @isnumeric);
			addRequired(p, 'DiscardList', @isnumeric);
			addParameter(p, 'HideResults', false, @islogical);
			parse(p, channel, discardList, varargin{:});
			channel = p.Results.Channel;
			discardList = p.Results.DiscardList;
			hideResults = p.Results.HideResults;

			% Discard unwanted clusters
			iWaveformToDiscard = ismember(obj.Spikes(channel).Cluster, discardList);
			obj.Spikes(channel).Waveforms(iWaveformToDiscard, :) = [];
			obj.Spikes(channel).WaveformsExtended(iWaveformToDiscard, :) = [];
			obj.Spikes(channel).Timestamps(iWaveformToDiscard) = [];
			obj.Spikes(channel).SampleIndex(iWaveformToDiscard) = [];
			obj.Spikes(channel).PCA.Score(iWaveformToDiscard, :) = [];
			obj.Spikes(channel).Cluster(iWaveformToDiscard) = [];

			% Renumber clusters from 1 to numClusters
			obj.ClusterMerge(channel, num2cell(unique(obj.Spikes(channel).Cluster)), 'HideResults', hideResults);
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

		% Plot spikes, run after ThresholdWaveforms
		function PlotChannel(obj, channel, varargin)
			p = inputParser;
			addRequired(p, 'Channel', @isnumeric);
			addOptional(p, 'NumWaveforms', 30, @isnumeric);
			addParameter(p, 'Clusters', [], @isnumeric);
			addParameter(p, 'YLim', [-500, 500], @isnumeric);
			addParameter(p, 'FrameRate', 30, @isnumeric);
			addParameter(p, 'WaveformWindow', [], @isnumeric);
			parse(p, channel, varargin{:});
			channel = p.Results.Channel;
			numWaveforms = p.Results.NumWaveforms;
			clusters = p.Results.Clusters;
			yRange = p.Results.YLim;
			frameRate = p.Results.FrameRate;
			waveformWindow = p.Results.WaveformWindow;

			waveforms = obj.Spikes(channel).WaveformsExtended;
			numWaveformsTotal = size(waveforms, 1);
			if isempty(waveformWindow)
				waveformWindow = obj.Spikes(channel).WaveformWindow;
			end
			waveformWindowExtended = obj.Spikes(channel).WaveformWindowExtended;
			t = obj.Spikes(channel).WaveformTimestampsExtended;
			score = obj.Spikes(channel).PCA.Score;
			if isempty(clusters)
				clusters = nonzeros(unique(obj.Spikes(channel).Cluster));
			end

			if isfield(obj.Spikes(channel), 'Cluster')
				clusterID = obj.Spikes(channel).Cluster;
				toKeep = ismember(clusterID, clusters);
				clusterID = clusterID(toKeep);
				waveforms = waveforms(toKeep, :);
				score = score(toKeep, :);
				numWaveformsTotal = size(waveforms, 1);
			else
				clusterID = ones(numWaveformsTotal, 1);
			end

			hFigure = figure('Units', 'Normalized', 'OuterPosition', [0.7, 0, 0.3, 1]);
			colors = 'rgbcmyk';		
			subplot(2, 1, 1)
			hold on
			for iCluster = unique(nonzeros(clusterID))'
				inCluster = clusterID == iCluster;
				percentage = round(100*sum(inCluster)/size(obj.Spikes(channel).Waveforms, 1));
				scatter3(score(inCluster, 1), score(inCluster, 2), score(inCluster, 3), 1, colors(iCluster), 'DisplayName', ['Cluster ', num2str(iCluster), ' (', num2str(percentage), '%)'])
			end
			hold off
			axis equal
			xlabel('1st Principal Component')
			ylabel('2nd Principal Component')
			zlabel('3rd Principal Component')
			title('Spike Sorting')
			legend()
			hAxes = subplot(2, 1, 2);
			hold on
			hLegends = [];
			hAxes.UserData.hWaveforms = [];
			hAxes.UserData.iWaveform = 0;
			if sum(waveformWindowExtended == waveformWindow) < 2
				hLegends = [hLegends, line(repmat(waveformWindow(1), [1, 2]), yRange, 'Color', 'r', 'DisplayName', ['SpikeSort Window (', num2str(waveformWindow(1)), ' to ', num2str(waveformWindow(2)), ' ms)'])];
				line(repmat(waveformWindow(2), [1, 2]), yRange, 'Color', 'r')
			end
			hLegends = [hLegends, line([obj.Spikes(channel).WaveformWindow(1), 0], repmat(obj.Spikes(channel).Threshold.Trigger, [1, 2]), 'Color', 'k', 'LineWidth', 3, 'DisplayName', ['Trigger Threshold (', num2str(obj.Spikes(channel).Threshold.Trigger), ' \muV)'])];
			if ~isnan(obj.Spikes(channel).Threshold.Exit)
				hLegends = [hLegends, line(obj.Spikes(channel).Threshold.ExitWindow, repmat(obj.Spikes(channel).Threshold.Exit, [1, 2]), 'Color', 'b', 'LineWidth', 3, 'DisplayName', ['Exit Threshold (', num2str(obj.Spikes(channel).Threshold.Exit), ' \muV)'])];
			end
			if ~isnan(obj.Spikes(channel).Threshold.Exclusion)
				hLegends = [hLegends, line(obj.Spikes(channel).WaveformWindow, repmat(obj.Spikes(channel).Threshold.Exclusion, [1, 2]), 'Color', 'm', 'LineWidth', 3, 'DisplayName', ['Exclusion Threshold (', num2str(obj.Spikes(channel).Threshold.Exclusion), ' \muV)'])];
			end
			legend(hLegends, 'AutoUpdate', 'off')
			ylim(yRange)

			hTimer = timer(...
				'ExecutionMode', 'FixedSpacing',...
			 	'Period', round((1/frameRate)*1000)/1000,...
			 	'TimerFcn', {@TetrodeRecording.OnPlotChannelRefresh, hAxes, t, waveforms, numWaveforms, numWaveformsTotal, colors, clusterID}...
			 	);
			hAxes.UserData.hTimer = hTimer;
			hFigure.KeyPressFcn = {@TetrodeRecording.OnKeyPress, hTimer};
			hFigure.CloseRequestFcn = {@TetrodeRecording.OnFigureClosed, hTimer};
			start(hTimer);
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
			if strcmp(binMethod, 'percentile')
				edges = prctile(trialLength(trialLength > minTrialLength), 0:(100/nBins):100);
			elseif strcmp(binMethod, 'equal')
				edges = linspace(max(minTrialLength, min(trialLength)), max(trialLength), nBins + 1);
			end
			[NTrials, ~, bins] = histcounts(trialLength, edges);

			colors = 'rgbcmyk';
			styles = {'-', '--', ':', '-.'};

			for iChannel = channels
				[spikes, trials] = obj.GetSpikesByTrial(iChannel, reference, event, [0, 0], clusters);

				if isempty(ax)
					hAxes = axes(figure());
				else
					hAxes = ax;
				end
				hold on
				for iBin = 1:nBins
					inBin = ismember(trials, find(bins == iBin));
					cutOff = -edges(iBin);
					spikesRelative = spikes(inBin) - event(trials(inBin)); % Spike times relative to event
					spikesRelative = spikesRelative(spikesRelative > cutOff & spikesRelative < 0);
					% Windows for estimating spike rate
					thisEdges = cutOff:(spikeRateWindow/1000):0;
					if length(thisEdges) < 3
						continue
					end
					% if thisEdges(end) < 0
					% 	thisEdges(end + 1) = 0;
					% end
					thisSpikeRate = histcounts(spikesRelative, thisEdges);
					thisSpikeRate = (1000*thisSpikeRate/spikeRateWindow)/NTrials(iBin);
					thisCenters = (thisEdges(1:end - 1) + thisEdges(2:end))/2;

					iStyle = ceil(iBin/length(colors));
					iColor = mod(iBin, length(colors));
					if iColor == 0
						iColor = 7;
					end
					thisColor = colors(iColor);
					thisStyle = styles{iStyle};
					plot(hAxes, thisCenters, thisSpikeRate, [thisColor, thisStyle],...
						'DisplayName', ['[', num2str(-cutOff, 2), ' s, ', num2str(edges(iBin + 1), 2), ' s] (', num2str(NTrials(iBin)), ' trials)'],...
						'LineWidth', 2.5);
				end
				xlabel(hAxes, ['Time relative to ', eventDisplayName, ' (s)']);
				ylabel(hAxes, 'Mean firing rate (Hz)')
				legend(hAxes, 'Location', 'southwest');
				hold off
			end
		end

		% Sort spikes and digital events into trial structure
		function Raster(obj, channels, varargin)
			p = inputParser;
			addRequired(p, 'Channels', @isnumeric);
			addOptional(p, 'Reference', 'CueOn', @ischar);
			addOptional(p, 'Event', 'PressOn', @ischar);
			addOptional(p, 'Exclude', 'LickOn', @ischar);
			addParameter(p, 'ExtendedWindow', [-2, 2], @isnumeric);
			addParameter(p, 'Clusters', [], @isnumeric);
			parse(p, channels, varargin{:});
			channels = p.Results.Channels;
			reference = p.Results.Reference;
			event = p.Results.Event;
			exclude = p.Results.Exclude;
			extendedWindow = p.Results.ExtendedWindow;
			clusters = p.Results.Clusters;

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
				[~, I] = sort(reference - event);
				trialsSorted = changem(trials, 1:length(unique(I)), I);	% !!! changem requires mapping toolbox.

				figure('Units', 'Normalized', 'OuterPosition', [0, 0, 0.75, 1]);
				hAxes = subplot(1, 2, 1);
				hold on
				plot(hAxes, spikes - reference(trials), trialsSorted, '.',...
					'MarkerSize', 5,...
					'MarkerEdgeColor', 'k',...
					'MarkerFaceColor', 'k',...
					'LineWidth', 1.5,...
					'DisplayName', 'Spike'...
				)
				plot(hAxes, event(trials) - reference(trials), trialsSorted, '.',...
					'MarkerSize', 10,...
					'MarkerEdgeColor', 'b',...
					'MarkerFaceColor', 'b',...
					'LineWidth', 1.5,...
					'DisplayName', eventDisplayName...
				)
				expName = strsplit(obj.Path, '\');
				expName = expName{end - 1};
				title([expName, ' - Channel ', num2str(iChannel)], 'Interpreter', 'none')
				xlabel(['Time relative to ', referenceDisplayName, ' (s)'])
				ylabel('Trial')
				legend('Location', 'northeast');
				hold off

				hAxes = subplot(1, 2, 2);
				hold on
				plot(hAxes, spikes - event(trials), trialsSorted, '.',...
					'MarkerSize', 5,...
					'MarkerEdgeColor', [.1, .1, .1],...
					'MarkerFaceColor', [.1, .1, .1],...
					'LineWidth', 1.5,...
					'DisplayName', 'Spike'...
				)
				plot(hAxes, reference(trials) - event(trials), trialsSorted, '.',...
					'MarkerSize', 10,...
					'MarkerEdgeColor', 'r',...
					'MarkerFaceColor', 'r',...
					'LineWidth', 1.5,...
					'DisplayName', referenceDisplayName...
				)
				xlabel(['Time relative to ', eventDisplayName,' (s)'])
				ylabel('Trial')
				legend('Location', 'northwest');
				hold off
			end
		end

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
					spikes = obj.Spikes(channel).Timestamps(ismember(obj.Spikes(channel).Cluster, clusters));
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

		function OnPlotChannelRefresh(~, ~, hAxes, t, waveforms, numWaveforms, numWaveformsTotal, colors, clusterID)
			hAxes.UserData.iWaveform = hAxes.UserData.iWaveform + 1;
			iWaveform = hAxes.UserData.iWaveform;

			if iWaveform > numWaveformsTotal
				stop(hAxes.UserData.hTimer);
				delete(hAxes.UserData.hTimer);
				return
			else
				if length(hAxes.UserData.hWaveforms) >= numWaveforms
					delete(hAxes.UserData.hWaveforms(1));
					hAxes.UserData.hWaveforms = hAxes.UserData.hWaveforms(2:end);
				end
				hAxes.UserData.hWaveforms = [hAxes.UserData.hWaveforms, plot(hAxes, t, waveforms(iWaveform, :), 'LineStyle', '-', 'Color', colors(clusterID(iWaveform)), 'DisplayName', ['Waveform (Cluster ', num2str(clusterID(iWaveform)), ')'])];
				if iWaveform == 1
					xlabel(hAxes, 'Time (ms)');
					ylabel(hAxes, 'Voltage (\muV');
					title(hAxes, 'Waveforms');
				end
			end	
		end

		function OnKeyPress(~, evnt, hTimer)
			if isvalid(hTimer)
				if strcmp(evnt.Key, 'space')
					if strcmp(hTimer.Running, 'on')
						stop(hTimer);
					elseif strcmp(hTimer.Running, 'off')
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
	end
end
