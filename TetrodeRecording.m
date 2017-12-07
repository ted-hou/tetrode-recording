classdef TetrodeRecording < handle
	properties
		Files
		Path
		Notes
		FrequencyParameters
		Amplifier
		BoardDigIn
		Spikes
	end

	%----------------------------------------------------
	%		Methods
	%----------------------------------------------------
	methods
		function obj = TetrodeRecording(chunkSize)
			obj.SelectFiles();

			if nargin < 2
				chunkSize = input('Chunk size: ');
			end
			obj.ReadRHD(obj.Files(1:min(chunkSize, length(obj.Files))))
			obj.MapChannels();
			obj.RemoveTransient();
		end

		function SelectFiles(obj)
			[obj.Files, obj.Path, ~] = uigetfile('*.rhd', 'Select an RHD2000 Data File', 'MultiSelect', 'on');
		end

		function ReadRHD(obj, files)
			% if nargin < 2
			% 	[files, path, ~] = uigetfile('*.rhd', 'Select an RHD2000 Data File', 'MultiSelect', 'on');
			% end

			for iFile = 1:length(files)
				filename = [obj.Path, files{iFile}];

				if exist(filename, 'file') ~= 2
					error(['File ', filename, ' not found.']);
					return
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
					fprintf(1, ['Reading data from file ', num2str(iFile), '/', num2str(length(files)), ' (''', files{iFile},''')...\n']);

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
							% fprintf(1, '\t%d%% done...\n', percent_done);
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
			fprintf(1, 'All files loaded. Concatenating data...');

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
			fprintf(1, 'Done!\n');
		end

		% Amplifier data arranged in tetrode order (i.e., first four elements in array is first tetrode)
		function MapChannels(obj, EIBMap, headstageType)
			if isfield(obj.Amplifier, 'ChannelMap')
				warning('Amplifier data already remapped in tetrode order. I refuse to do it again.');
				return
			end

			fprintf(1, 'Remapping amplifier channels in tetrode order...');

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
			
			obj.Amplifier.DataMapped = nan(32, size(obj.Amplifier.Data, 2));
			iChn = 0;
			for targetChannel = tetrodeMap
				iChn = iChn + 1;
				sourceChannel = find(recordedChannels == targetChannel, 1);
				if ~isempty(sourceChannel)
					obj.Amplifier.DataMapped(iChn, :) = obj.Amplifier.Data(sourceChannel, :);
				end
			end
			obj.Amplifier.ChannelMap.ElectrodeInterfaceBoard = EIBMap;
			obj.Amplifier.ChannelMap.Headstage = HSMap;
			obj.Amplifier.ChannelMap.Tetrode = tetrodeMap;
			obj.Amplifier.Data = obj.Amplifier.DataMapped; 
			obj.Amplifier = rmfield(obj.Amplifier, 'DataMapped');
			fprintf(1, 'Done!\n');
		end

		% Convert zero-based raw channel id to remapped tetrode id
		function tetrodeChannel = GetMappedChannelID(obj, rawChannel)
			tetrodeChannel = find(obj.Amplifier.ChannelMap.Tetrode == (rawChannel + 1));
		end

		function PlotChannel(obj, channel, xRange, yRange, digInChannels)
			if nargin < 4
				yRange = [-110, 110];
			end
			if nargin < 5
				digInChannels = [4, 2, 5, 3]; % CUE, LEVER, REWARD, LICK
			end
			digInScaled = obj.BoardDigIn.Data(digInChannels, :).*[100; 500; 200; 350]; % CUE, LEVER, REWARD, LICK

			obj.Amplifier.DataMean = nanmean(obj.Amplifier.Data, 1);

			tetrodeColors = {'c-', 'm-', 'y-', 'b-', 'c--', 'm--', 'y--', 'b--'};

			figure('Units', 'normalized', 'Position', [0, 0, 1, 0.6])

			subplot(2, 1, 1)
			hold on
			plot(obj.Amplifier.Timestamps, obj.Amplifier.Data(channel, :));
			plot(obj.BoardDigIn.Timestamps, digInScaled, 'r:', 'LineWidth', 2);
			hold off
			xlim(xRange)
			ylim(yRange)
			xlabel('Time (s)')
			ylabel('\muV')
			title('Raw')

			subplot(2, 1, 2)
			hold on
			plot(obj.Amplifier.Timestamps, obj.Amplifier.Data(channel, :) - obj.Amplifier.DataMean);
			plot(obj.BoardDigIn.Timestamps, digInScaled, 'r:', 'LineWidth', 2);
			hold off
			xlim(xRange)
			ylim(yRange)
			xlabel('Time (s)')
			ylabel('\muV')
			title('Substracted by 32-Chn mean')
		end

		% Remove lick/press-transient by setting signal to zero
		function RemoveTransient(obj, digChannel, dilate)
			if nargin < 2
				digChannel = [1, 2]; % Default lick channel is 2 on intan board. This was used for Daisy1. Also do this for lever press (Chn 1)
			end
			if nargin < 3
				dilate = 12; % By default, lick digital event is extended 2 ms to the left and right.
			end

			% Verify if digital input sample rate is identical to amplifier sample rate
			if obj.FrequencyParameters.AmplifierSampleRate ~= obj.FrequencyParameters.BoardDigInSampleRate
				error('Board digital input has a different sampling rate from amplifier. Aborted.');
				return
			end

			mask = logical(zeros(1, obj.Amplifier.NumSamples));
			dilateSE = ones(round(obj.FrequencyParameters.AmplifierSampleRate*2*dilate/1000) + 1);
			for thisChannel = digChannel
				thisChannel = find([obj.BoardDigIn.Channels.NativeOrder] == thisChannel);
				mask = mask | logical(imdilate(obj.BoardDigIn.Data(thisChannel, :), dilateSE));
			end
			obj.Amplifier.Data(:, mask) = 0;
		end

		% Threshold and plot spikes
		function ThresholdAndPlotSpikes(obj, channel, threshold, numSpikes, spikeWindow)
			if nargin < 4
				numSpikes = 30;
			end
			if nargin < 5
				spikeWindow = [0.5, 1];
			end
			[waveforms, t, spikeID] = ThresholdSpikes(obj, channel, threshold, spikeWindow);
			figure();
			hold on
			for iSpike = 1:numSpikes
				plot(t, waveforms(iSpike, :), 'k');
				xlabel('Time (ms)');
				ylabel('Voltage (\muV');
				title('Spikes');
				pause(0.01);
			end
			hold off			
		end

		function [waveforms, t, spikeID] = ThresholdSpikes(obj, channel, threshold, spikeWindow)
			if nargin < 4
				spikeWindow = [0.5, 1];
			end
			% Find threshold-crossing events
			spikeID = find([false, diff(obj.Amplifier.Data(channel, :) > threshold) == sign(threshold)]);
			numSpikes = length(spikeID);

			sampleRate = obj.FrequencyParameters.AmplifierSampleRate;
			t = -spikeWindow(1):1000/sampleRate:spikeWindow(2);
			waveforms = nan(numSpikes, 1 + sampleRate*sum(spikeWindow)/1000);
			for iSpike = 1:numSpikes
				startIndex = spikeID(iSpike) - round(sampleRate*spikeWindow(1)/1000);
				endIndex = spikeID(iSpike) + round(sampleRate*spikeWindow(2)/1000);
				if startIndex >= 1 && endIndex <= obj.Amplifier.NumSamples
					waveforms(iSpike, :) = obj.Amplifier.Data(channel, startIndex:endIndex);
				end
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
	end
end
