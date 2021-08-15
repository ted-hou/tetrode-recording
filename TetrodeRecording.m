classdef TetrodeRecording < handle
	properties
		System = ''
		Files = []
		Path = []
		Part = [1, 1]
		Notes
		FrequencyParameters
		ChannelMap
		SelectedChannels
		Spikes
		DigitalEvents
		AnalogIn
		StartTime
	end

	properties (Transient)
		Amplifier
		BoardDigIn
		BoardADC
		NEV
		NSx
	end

	%----------------------------------------------------
	%		Methods
	%----------------------------------------------------
	methods
		function obj = TetrodeRecording(data, timestamps, sampleRate)
			if nargin > 0
				obj.Amplifier.Data = data;
				obj.Amplifier.Timestamps = timestamps;
				obj.FrequencyParameters.AmplifierSampleRate = sampleRate;
				obj.Files = 'mouse_yyyymmdd';
				obj.Path = 'Z:\mouse_yyyymmdd\';
			end
		end

		function Preview(obj, varargin)
			p = inputParser;
			addParameter(p, 'Duration', [300, 360], @isnumeric); % start:stop in seconds
			addParameter(p, 'Channels', [], @isnumeric);
			addParameter(p, 'ChunkSize', 10, @isnumeric);
			addParameter(p, 'HideResults', false, @islogical);
			addParameter(p, 'Rig', 1, @isnumeric);
			parse(p, varargin{:});
			duration 		= p.Results.Duration;
			channels 		= p.Results.Channels;
			chunkSize 		= p.Results.ChunkSize;
			hideResults 	= p.Results.HideResults;
			rig 			= p.Results.Rig;

			if isempty(obj.Path)
				obj.SelectFiles();
			end
			obj.ReadFiles(chunkSize, 'Rig', rig, 'Duration', duration, 'NumSigmas', 2.5, 'NumSigmasReturn', 1.25, 'NumSigmasReject', 40, 'WaveformWindow', [-0.35, 0.35]);
			obj.SpikeSort(channels, 'ClusterMethod', 'kmeans', 'FeatureMethod', 'PCA', 'Dimension', 3);
			if ~hideResults
				obj.PlotAllChannels();
			end
		end

		function expName = GetExpName(obj)
			expName = strsplit(obj.Path, '\');
			expName = expName{end - 1};
		end

		function startTime = GetStartTime(obj)
			if isempty(obj.StartTime)
				if isempty(obj.NEV)
					NEV = openNEV([obj.Path, obj.Files{1}], 'nosave', 'nomat');
					startTime = NEV.MetaTags.DateTimeRaw;
				else
					startTime = obj.NEV.MetaTags.DateTimeRaw;
				end
				startTime(7) = startTime(7) + startTime(8)/1000;
				startTime = startTime([1,2,4,5,6,7]);
				startTime = datetime(startTime, 'TimeZone', 'UTC');
				obj.StartTime = startTime;
			else
				startTime = obj.StartTime;
			end
		end

		function SelectFiles(obj)
			[obj.Files, obj.Path, filterindex] = uigetfile({'*.nev; *.ns*; *.ccf', 'Blackrock Files (*.nev, *.ns*, *.ccf)'; '*.rhd', 'Intan RHD2000 Files (*.rhd)'}, 'Select files', 'MultiSelect', 'on');
			switch filterindex
				case 1
					obj.System = 'Blackrock';
				case 2
					obj.System = 'Intan';
			end
				
		end

		function ReadFiles(obj, varargin)
			p = inputParser;
			addOptional(p, 'ChunkSize', 2, @isnumeric);
			addParameter(p, 'Duration', [], @isnumeric); % [0, 60] to read 0 - 60 seconds of data. For blackrock only
			addParameter(p, 'SubstractMean', false, @islogical);
			addParameter(p, 'RemoveTransient', false, @islogical);
			addParameter(p, 'DigitalChannels', 'auto', @(x) iscell(x) || ischar(x)); % 'auto', or custom, e.g. {'Cue', 4; 'Press', 2; 'Lick', 3; 'Reward', 5}
			addParameter(p, 'AnalogChannels', 'auto', @(x) iscell(x) || ischar(x)); % 'auto', or custom, e.g. {'AccX', 1; 'AccY', 2; 'AccZ', 3;}
			addParameter(p, 'Channels', [], @isnumeric);
			addParameter(p, 'ChannelsOnRig', NaN, @isnumeric);
			addParameter(p, 'NumSigmas', 4, @isnumeric);
			addParameter(p, 'NumSigmasReturn', [], @isnumeric);
			addParameter(p, 'NumSigmasReject', [], @isnumeric);
			addParameter(p, 'Direction', 'negative', @ischar);
			addParameter(p, 'WaveformWindow', [-1.25, 1.25], @isnumeric);
			addParameter(p, 'Rig', 1, @isnumeric);
			addParameter(p, 'DetectSpikes', true, @islogical);
			addParameter(p, 'DetectEvents', true, @islogical);
			parse(p, varargin{:});
			chunkSize 		= p.Results.ChunkSize;
			duration 		= p.Results.Duration;
			substractMean 	= p.Results.SubstractMean;
			removeTransient = p.Results.RemoveTransient;
			digitalChannels = p.Results.DigitalChannels;
			analogChannels 	= p.Results.AnalogChannels;
			channels 		= p.Results.Channels;
			channelsOnRig	= p.Results.ChannelsOnRig;
			numSigmas 		= p.Results.NumSigmas;
			numSigmasReturn = p.Results.NumSigmasReturn;
			numSigmasReject = p.Results.NumSigmasReject;
			direction 		= p.Results.Direction;
			waveformWindow 	= p.Results.WaveformWindow;
			rig 			= p.Results.Rig;
			detectSpikes	= p.Results.DetectSpikes;
			detectEvents	= p.Results.DetectEvents;

			files = obj.Files;
			numChunks = ceil(length(files)/chunkSize);

			switch lower(obj.System)
				case 'intan'
					for iChunk = 1:numChunks
						TetrodeRecording.TTS(['Processing chunk ', num2str(iChunk), '/', num2str(numChunks), ':\n']);
						obj.ReadIntan(obj.Files((iChunk - 1)*chunkSize + 1:min(iChunk*chunkSize, length(obj.Files))))
						obj.GenerateChannelMap('HeadstageType', 'Intan');
						if substractMean
							obj.SubstractMean();
						end
						if removeTransient
							obj.RemoveTransient();
						end
						if isempty(channels)
							channels = obj.MapChannel_RawToTetrode([obj.Amplifier.Channels.NativeOrder] + 1);
						end
						obj.SpikeDetect(channels, 'NumSigmas', numSigmas, 'NumSigmasReturn', numSigmasReturn, 'NumSigmasReject', numSigmasReject, 'WaveformWindow', waveformWindow, 'Direction', direction, 'Append', iChunk > 1);
						obj.GetDigitalData('ChannelLabels', digitalChannels, 'Append', iChunk > 1);
						obj.GetAnalogData('ChannelLabels', analogChannels, 'Append', iChunk > 1);
						obj.ClearCache();
					end
					obj.GetDigitalEvents(true);
				case 'blackrock'
					if ~detectEvents
						digitalChannels = {};
					else
						if rig == 1
							digitalChannels = {'Cue', 0; 'Reward', 1; 'Lick', 2; 'Press', 3; 'Stim', 4};
						elseif rig == 2
							digitalChannels = {'Cue', 15; 'Reward', 14; 'Lick', 13; 'Press', 12; 'Stim', 11};
						end
					end
					% First read
					if isnan(channelsOnRig)
						obj.ReadBlackrock('DigitalChannels', digitalChannels, 'Duration', duration);
						obj.ChannelMap.Rig1 = find(ismember([obj.NSx.ElectrodesInfo.ElectrodeID], 1:32));
						obj.ChannelMap.Rig2 = find(ismember([obj.NSx.ElectrodesInfo.ElectrodeID], 65:96));
						if rig == 1
							channelsOnRig = obj.ChannelMap.Rig1;
						elseif rig == 2
							channelsOnRig = obj.ChannelMap.Rig2;
						end
						obj.Amplifier.Data = obj.Amplifier.Data(channelsOnRig, :);
						if detectSpikes
							obj.SpikeDetect(1:size(obj.Amplifier.Data, 1), 'NumSigmas', numSigmas, 'NumSigmasReturn', numSigmasReturn, 'NumSigmasReject', numSigmasReject, 'WaveformWindow', waveformWindow, 'Direction', direction, 'Append', false);
						end
					% Second read
					else
						channelsToRead = channelsOnRig(channels);
						obj.ReadBlackrock('Channels', channelsToRead, 'DigitalChannels', digitalChannels, 'Duration', duration);
						if detectSpikes
							obj.SpikeDetect(1:size(obj.Amplifier.Data, 1), 'NumSigmas', numSigmas, 'NumSigmasReturn', numSigmasReturn, 'NumSigmasReject', numSigmasReject, 'WaveformWindow', waveformWindow, 'Direction', direction, 'Append', false);
						end
					end
					if detectSpikes
						obj.ClearCache();
					end
			end
			obj.GetStartTime();
		end

		function ReadBlackrock(obj, varargin)
			p = inputParser;
			addParameter(p, 'Channels', [], @isnumeric);
			addParameter(p, 'DigitalChannels', {'Lever', 1; 'Lick', 2; 'Cue', 4; 'Reward', 7}, @iscell); % {'Lever', 1; 'Lick', 2; 'Cue', 4; 'Reward', 7}
			addParameter(p, 'Duration', [], @isnumeric); % [0, 60] seconds
			parse(p, varargin{:});
			channels 		= p.Results.Channels;
			digitalChannels = p.Results.DigitalChannels;
			duration 		= p.Results.Duration;

			if ~isempty(duration) && length(duration) == 2
				duration = ['t:', num2str(duration(1)), ':', num2str(duration(2))];
			else
				duration = [];
			end

			tic, TetrodeRecording.TTS(['	Loading data (', obj.GetExpName(), ')...'])
			% Load NEV
			isNEV = contains(obj.Files, '.nev'); 
			if nnz(isNEV) == 1
				filename = [obj.Path, obj.Files{isNEV}];
				if isempty(duration)
					obj.NEV = openNEV(filename, 'nosave', 'nomat');
				else
					obj.NEV = openNEV(filename, 'nosave', 'nomat', duration);
				end
			elseif nnz(isNEV) > 1
				warning('More than one NEV file specified so none will be loaded. Digital events might be missing.')
			else
				warning('No NEV file loaded. Digital events might be missing.')
			end

			% Load NSx
			isNSx = contains(obj.Files, '.ns'); 
			if nnz(isNSx) == 1
				filename = [obj.Path, obj.Files{isNSx}];
				if isempty(channels)
					if isempty(duration)
						obj.NSx = openNSx(filename);
					else
						obj.NSx = openNSx(filename, 'duration', duration, 'sec');
					end
				else
					if isempty(duration)
						obj.NSx = openNSx(filename, 'channels', channels);
					else
						obj.NSx = openNSx(filename, 'channels', channels, 'duration', duration, 'sec');
					end
				end
			elseif nnz(isNSx) > 1
				warning('More than one NSx file specified so none will be loaded.')
			else
				warning('No NSx file specified.')
			end
			TetrodeRecording.TTS(['Done(', num2str(toc, '%.2f'), ' seconds).\n'])

			% Process amplifier data
			sampleRate = obj.NSx.MetaTags.SamplingFreq;
			obj.FrequencyParameters.AmplifierSampleRate = sampleRate;
			numSamples = obj.NSx.MetaTags.DataPoints(end);

			obj.Amplifier.NumSamples = numSamples;
			obj.Amplifier.Timestamps = (0:(numSamples - 1))/sampleRate;
			if iscell(obj.NSx.Data)
				obj.Amplifier.Data = [obj.NSx.Data{end}];
				obj.FrequencyParameters.SysInitDelay.NumSamples = length(obj.NSx.Data{1});
				obj.FrequencyParameters.SysInitDelay.Duration = length(obj.NSx.Data{1})/sampleRate;
				obj.FrequencyParameters.SysInitDelay.DataTrimmed = true;
				warning(['For some readon data was truncated into ', num2str(length(obj.NSx.Data)), ' chunk (', mat2str(cellfun(@(x) length(x)/sampleRate, obj.NSx.Data)),'). Timestamps: ', mat2str(obj.NSx.MetaTags.Timestamp),'. If its only ywo parts and the first chunk is a few seconds long or less, then its discarded data recorded by rig 1 before rig 2 starts recording. If it''s more than a few seconds long or there are more than one data chunk, then you''re screwed because data acquisition was paused and then manually resumed - possibly due to disk being full.']);
			else
				obj.FrequencyParameters.SysInitDelay.NumSamples = NaN;
				obj.FrequencyParameters.SysInitDelay.Duration = NaN;
				obj.FrequencyParameters.SysInitDelay.DataTrimmed = false;
				obj.Amplifier.Data = obj.NSx.Data;
			end
			obj.NSx.Data = [];

			% Parse digital events
			if ~isempty(digitalChannels)
				digitalData = flip(dec2bin(obj.NEV.Data.SerialDigitalIO.UnparsedData), 2);
				digitalTimestamps = obj.NEV.Data.SerialDigitalIO.TimeStampSec;

				for iChannel = 1:size(digitalChannels, 1)
					iBit = digitalChannels{iChannel, 2} + 1;
					channelName = digitalChannels{iChannel, 1};
					[obj.DigitalEvents.([channelName, 'On']), obj.DigitalEvents.([channelName, 'Off'])] = obj.FindEdges(transpose(digitalData(:, iBit)), digitalTimestamps);
				end
			end
		end

		function ReadDigitalEvents(obj)
			% Determine digital channels
			rig = TetrodeRecording.GetRig(obj.Path);
			if rig == 1
				digitalChannels = {'Cue', 0; 'Reward', 1; 'Lick', 2; 'Press', 3; 'Stim', 4};
			elseif rig == 2
				digitalChannels = {'Cue', 15; 'Reward', 14; 'Lick', 13; 'Press', 12; 'Stim', 11};
			end

			% Load NEV
			isNEV = contains(obj.Files, '.nev'); 
			if nnz(isNEV) == 1
				filename = [obj.Path, obj.Files{isNEV}];
				obj.NEV = openNEV(filename, 'nosave', 'nomat');
			elseif nnz(isNEV) > 1
				warning('More than one NEV file specified so none will be loaded. Digital events might be missing.')
			else
				warning('No NEV file loaded. Digital events might be missing.')
			end

			% Parse digital events
			digitalData = flip(dec2bin(obj.NEV.Data.SerialDigitalIO.UnparsedData), 2);
			digitalTimestamps = obj.NEV.Data.SerialDigitalIO.TimeStampSec;

			for iChannel = 1:size(digitalChannels, 1)
				iBit = digitalChannels{iChannel, 2} + 1;
				channelName = digitalChannels{iChannel, 1};
				[obj.DigitalEvents.([channelName, 'On']), obj.DigitalEvents.([channelName, 'Off'])] = obj.FindEdges(transpose(digitalData(:, iBit)), digitalTimestamps);
			end
			obj.NEV = [];
		end

		function ReadIntan(obj, files)
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
				objTemp(iFile).BoardADC.Channels = struct(channel_struct);

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
										objTemp(iFile).BoardADC.Channels(board_adc_index) = new_channel;
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
					objTemp(iFile).BoardADC.Data = zeros(num_board_adc_channels, num_board_adc_samples);
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
					if (eval_board_mode == 1)
						objTemp(iFile).BoardADC.Data = 152.59e-6 * (objTemp(iFile).BoardADC.Data - 32768); % units = volts
					elseif (eval_board_mode == 13) % Intan Recording Controller
						objTemp(iFile).BoardADC.Data = 312.5e-6 * (objTemp(iFile).BoardADC.Data - 32768); % units = volts    
					else
						objTemp(iFile).BoardADC.Data = 50.354e-6 * objTemp(iFile).BoardADC.Data; % units = volts
					end
					% Scale time steps (units = seconds).
					objTemp(iFile).Amplifier.Timestamps = objTemp(iFile).Amplifier.Timestamps / sample_rate;
					objTemp(iFile).BoardDigIn.Timestamps = objTemp(iFile).Amplifier.Timestamps;
					objTemp(iFile).BoardADC.Timestamps = objTemp(iFile).Amplifier.Timestamps;
				end
				TetrodeRecording.TTS(['Done(', num2str(toc, '%.2f'), ' seconds).\n'])
			end

			% Count number of samples in all files
			obj.Amplifier.NumSamples = 0;
			obj.BoardDigIn.NumSamples = 0;
			obj.BoardADC.NumSamples = 0;
			for iFile = 1:length(files)
				obj.Amplifier.NumSamples = obj.Amplifier.NumSamples + size(objTemp(iFile).Amplifier.Timestamps, 2);
				obj.BoardDigIn.NumSamples = obj.BoardDigIn.NumSamples + size(objTemp(iFile).BoardDigIn.Timestamps, 2);
				obj.BoardADC.NumSamples = obj.BoardADC.NumSamples + size(objTemp(iFile).BoardADC.Timestamps, 2);
			end

			% Combine files
			tic, TetrodeRecording.TTS('	Concatenating data...');

			obj.Notes = notes;
			obj.FrequencyParameters = frequency_parameters;
			obj.Amplifier.Channels = objTemp(1).Amplifier.Channels;
			obj.BoardDigIn.Channels = objTemp(1).BoardDigIn.Channels;
			obj.BoardADC.Channels = objTemp(1).BoardADC.Channels;

			obj.Amplifier.Timestamps = zeros(1, obj.Amplifier.NumSamples);
			obj.Amplifier.Data = zeros(length(obj.Amplifier.Channels), obj.Amplifier.NumSamples);
			obj.BoardDigIn.Timestamps = zeros(1, obj.BoardDigIn.NumSamples);
			obj.BoardDigIn.Data = zeros(length(obj.BoardDigIn.Channels), obj.BoardDigIn.NumSamples);
			obj.BoardADC.Timestamps = zeros(1, obj.BoardADC.NumSamples);
			obj.BoardADC.Data = zeros(length(obj.BoardADC.Channels), obj.BoardADC.NumSamples);

			iSample.Amplifier = 0;
			iSample.BoardDigIn = 0;
			iSample.BoardADC = 0;
			for iFile = 1:length(files)
				obj.Amplifier.Timestamps(1, iSample.Amplifier + 1:iSample.Amplifier + size(objTemp(iFile).Amplifier.Timestamps, 2)) = objTemp(iFile).Amplifier.Timestamps;
				obj.Amplifier.Data(:, iSample.Amplifier + 1:iSample.Amplifier + size(objTemp(iFile).Amplifier.Timestamps, 2)) = objTemp(iFile).Amplifier.Data;
				obj.BoardDigIn.Timestamps(1, iSample.BoardDigIn + 1:iSample.BoardDigIn + size(objTemp(iFile).BoardDigIn.Timestamps, 2)) = objTemp(iFile).BoardDigIn.Timestamps;
				obj.BoardDigIn.Data(:, iSample.BoardDigIn + 1:iSample.BoardDigIn + size(objTemp(iFile).BoardDigIn.Timestamps, 2)) = objTemp(iFile).BoardDigIn.Data;
				obj.BoardADC.Timestamps(1, iSample.BoardADC + 1:iSample.BoardADC + size(objTemp(iFile).BoardADC.Timestamps, 2)) = objTemp(iFile).BoardADC.Timestamps;
				obj.BoardADC.Data(:, iSample.BoardADC + 1:iSample.BoardADC + size(objTemp(iFile).BoardADC.Timestamps, 2)) = objTemp(iFile).BoardADC.Data;

				iSample.Amplifier = iSample.Amplifier + size(objTemp(iFile).Amplifier.Timestamps, 2);
				iSample.BoardDigIn = iSample.BoardDigIn + size(objTemp(iFile).BoardDigIn.Timestamps, 2);
				iSample.BoardADC = iSample.BoardADC + size(objTemp(iFile).BoardADC.Timestamps, 2);
			end
			TetrodeRecording.TTS(['Done(', num2str(toc, '%.2f'), ' seconds).\n'])
		end

		% Used to preview a small portion of loaded data. Will remove used data from workspace.
		function TrimData(obj, numSamples)
			if ~isempty(obj.BoardDigIn)
				obj.BoardDigIn.NumSamples = numSamples;
				obj.BoardDigIn.Timestamps = obj.BoardDigIn.Timestamps(1:numSamples);
				obj.BoardDigIn.Data = obj.BoardDigIn.Data(:, 1:numSamples);
			end
			obj.Amplifier.NumSamples = numSamples;
			obj.Amplifier.Timestamps = obj.Amplifier.Timestamps(1:numSamples);
			obj.Amplifier.Data = obj.Amplifier.Data(:, 1:numSamples);
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

			tic, TetrodeRecording.TTS('	Removing lick/touch-related transients...');
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
			channel 		= p.Results.Channel;
			waveformWindow 	= p.Results.WaveformWindow;
			index 			= p.Results.Index;
			indexType 		= p.Results.IndexType;

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
			t = [flip(0:-1/sampleRate:waveformWindow(1)), 1/sampleRate:1/sampleRate:waveformWindow(2)];
			waveforms = NaN(length(sampleIndex), length(t));
			i = sampleRate*t;
			for iWaveform = 1:length(sampleIndex)
				iQuery = sampleIndex(iWaveform) + i;
				if (sum(rem(iQuery, 1) == 0) == length(iQuery)) && min(iQuery) > 0 && max(iQuery) <= size(obj.Amplifier.Data, 2)
					% waveforms(iWaveform, :) = obj.Amplifier.Data(obj.MapChannel_TetrodeToRecorded(channel), iQuery);
					waveforms(iWaveform, :) = obj.Amplifier.Data(channel, iQuery);
				else
					% waveforms(iWaveform, :) = interp1(1:size(obj.Amplifier.Data, 2), double(obj.Amplifier.Data(obj.MapChannel_TetrodeToRecorded(channel), :)), iQuery, 'pchip', NaN);
					waveforms(iWaveform, :) = interp1(1:size(obj.Amplifier.Data, 2), double(obj.Amplifier.Data(channel, :)), iQuery, 'pchip', NaN);
				end
			end

			% Output
			varargout = {waveforms, t, timestamps, sampleIndex};
		end

		% Delete waveforms by sampleIndex/simpleIndex (and additionally by cluster if specified)
		function DeleteWaveforms(obj, channel, index, varargin)
			p = inputParser;
			addRequired(p, 'Channel', @isnumeric);
			addRequired(p, 'Index', @isnumeric);
			addParameter(p, 'IndexType', 'Default', @ischar);
			addParameter(p, 'Clusters', [], @isnumeric);
			parse(p, channel, index, varargin{:});
			channel 	= p.Results.Channel;
			index 		= p.Results.Index;
			indexType	= p.Results.IndexType;
			clusters	= p.Results.Clusters;

			switch lower(indexType)
				case 'default'
					selected = find(ismember(1:length(obj.Spikes(channel).SampleIndex), index));
				case 'sampleindex'
					selected = find(ismember(obj.Spikes(channel).SampleIndex, index));
				case 'threshold'
					selected = find(obj.Spikes(channel).SampleIndex <= index);
			end
			if ~isempty(clusters)
				classes = obj.Spikes(channel).Cluster.Classes(selected);
				selected = selected(ismember(classes, clusters));
			end

			% Discard unwanted clusters
			obj.Spikes(channel).Waveforms(selected, :) = [];
			obj.Spikes(channel).Timestamps(selected) = [];
			obj.Spikes(channel).SampleIndex(selected) = [];
			obj.Spikes(channel).Feature.Coeff(selected, :) = [];
			obj.Spikes(channel).Cluster.Classes(selected) = [];

			% Remove empty clusters and consolidate cluster number
			if length(unique(obj.Spikes(channel).Cluster.Classes)) < max(obj.Spikes(channel).Cluster.Classes)
				obj.ClusterMerge(channel, num2cell(unique(obj.Spikes(channel).Cluster.Classes)));
			end
		end

		% Detect spike by simple thresholding
		function SpikeDetect(obj, channels, varargin)
			p = inputParser;
			addRequired(p, 'Channels', @isnumeric);
			addParameter(p, 'NumSigmas', 4, @isnumeric); % Spike detection threshold = (this*sigma*direction). Sigma is estimated noise standard deviation.
			addParameter(p, 'NumSigmasReturn', [], @isnumeric); % [] to disable. Waveform must return to this*sigma*direction after crossing threshold, helps remove noisy periods with non-zero baseline.
			addParameter(p, 'NumSigmasReject', [], @isnumeric); % [] to disable. Reject huge waveforms that exceed this many sigmas in either direction
			addParameter(p, 'Direction', 'negative', @ischar);
			addParameter(p, 'WaveformWindow', [-0.35, 0.35], @isnumeric);
			addParameter(p, 'Append', false, @islogical);
			parse(p, channels, varargin{:});
			channels 		= p.Results.Channels;
			numSigmas 		= p.Results.NumSigmas;
			numSigmasReturn = p.Results.NumSigmasReturn;
			numSigmasReject = p.Results.NumSigmasReject;
			directionMode 	= p.Results.Direction;
			waveformWindow 	= p.Results.WaveformWindow;
			append 			= p.Results.Append;

			TetrodeRecording.TTS('	Detecting spikes:\n');
			sampleRate = obj.FrequencyParameters.AmplifierSampleRate/1000;

			for iChannel = channels
				% Auto-threshold for spikes
				% median(abs(x))/0.6745 is a better estimation of noise amplitude than std() when there are spikes -- Quiroga, 2004
				% shoud be median(abx(x - median(x))) but for ephys usually median(x) = 0;
				sigma = nanmedian(abs(obj.Amplifier.Data(iChannel, :)))/0.6745;
				threshold = numSigmas*sigma;
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
				[~, sampleIndex] = findpeaks(double(direction*obj.Amplifier.Data(iChannel, :)), 'MinPeakHeight', threshold, 'MinPeakProminence', threshold);

				% Extract waveforms
				[waveforms, t] = obj.GetWaveforms(iChannel, waveformWindow, sampleIndex, 'IndexType', 'SampleIndex');

				% Align waveforms to peak
				numWaveforms = size(waveforms, 1);
				i = sampleRate*t;
				[~, maxIndex] = max(direction*waveforms, [], 2);
				alignmentShift = i(maxIndex);
				[waveforms, t, timestamps, sampleIndex] = obj.GetWaveforms(iChannel, waveformWindow, sampleIndex + alignmentShift, 'IndexType', 'SampleIndex');

				% Reject waveforms that do not return to a certain level after crossing threshold
				if ~isempty(numSigmasReturn)
					if direction > 0
						selected = min(waveforms(:, t > 0), [], 2) <= numSigmasReturn*sigma;
					else
						selected = max(waveforms(:, t > 0), [], 2) >= -numSigmasReturn*sigma;
					end
					waveforms = waveforms(selected, :);
					timestamps = timestamps(selected);
					sampleIndex = sampleIndex(selected);
				end

				% Reject waveforms that exceed a threshold
				if ~isempty(numSigmasReject)
					selected = max(abs(waveforms), [], 2) < abs(numSigmasReject*sigma);
					waveforms = waveforms(selected, :);
					timestamps = timestamps(selected);
					sampleIndex = sampleIndex(selected);
				end

				% Double data so we can do divisions and stuff
				waveforms = double(waveforms);

				% Store data
				if ~append
					obj.Spikes(iChannel).Channel = iChannel;

					obj.Spikes(iChannel).SampleIndex = sampleIndex;
					obj.Spikes(iChannel).Timestamps = timestamps;
					obj.Spikes(iChannel).Waveforms = waveforms;

					obj.Spikes(iChannel).WaveformTimestamps = t;
					obj.Spikes(iChannel).WaveformWindow = waveformWindow;
					obj.Spikes(iChannel).Threshold.NumSigmas = numSigmas;
					obj.Spikes(iChannel).Threshold.NumSigmasReturn = numSigmasReturn;
					obj.Spikes(iChannel).Threshold.NumSigmasReject = numSigmasReject;
					obj.Spikes(iChannel).Threshold.Threshold = direction*threshold;
					obj.Spikes(iChannel).Threshold.ThresholdReturn = direction*numSigmasReturn*sigma;
					obj.Spikes(iChannel).Threshold.ThresholdReject = [abs(numSigmasReject*sigma); -abs(numSigmasReject*sigma)];
					obj.Spikes(iChannel).Threshold.Direction = directionMode;
				else
					obj.Spikes(iChannel).SampleIndex = [obj.Spikes(iChannel).SampleIndex, sampleIndex + length(obj.DigitalEvents.Timestamps)];
					obj.Spikes(iChannel).Timestamps = [obj.Spikes(iChannel).Timestamps, timestamps];
					obj.Spikes(iChannel).Waveforms = [obj.Spikes(iChannel).Waveforms; waveforms];
					obj.Spikes(iChannel).Threshold.Threshold = [obj.Spikes(iChannel).Threshold.Threshold, direction*threshold];
					obj.Spikes(iChannel).Threshold.ThresholdReturn = [obj.Spikes(iChannel).Threshold.ThresholdReturn, direction*numSigmasReturn*sigma];
					obj.Spikes(iChannel).Threshold.ThresholdReject = [obj.Spikes(iChannel).Threshold.ThresholdReject, [abs(numSigmasReject*sigma); -abs(numSigmasReject*sigma)]];
				end

				TetrodeRecording.TTS(['Done(', num2str(numWaveforms), ' waveforms, ', num2str(toc, '%.2f'), ' seconds).\n'])
			end
		end

		function GetAnalogData(obj, varargin)
			p = inputParser;
			addParameter(p, 'ChannelLabels', {}, @(x) iscell(x) || ischar(x)); %e.g. {'AccelerometerX', 4; 'AccelerometerY', 2; 'AccelerometerZ', 3}
			addParameter(p, 'Append', false, @islogical);
			parse(p, varargin{:});
			channelLabels 	= p.Results.ChannelLabels;
			append 			= p.Results.Append;

			if ischar(channelLabels) && strcmpi(channelLabels, 'auto')
				channelLabels = transpose([{obj.BoardADC.Channels.CustomChannelName}; num2cell(1:length(obj.BoardADC.Channels))]);
			end

			if isempty(channelLabels)
				obj.AnalogIn.ChannelNames = {};
				obj.AnalogIn.Data = [];
				obj.AnalogIn.Timestamps = [];
				return
			end

			channels = cell2mat(channelLabels(:, 2));
			channelNames = channelLabels(:, 1);

			if append
				obj.AnalogIn.Data = [obj.AnalogIn.Data, obj.BoardADC.Data(channels, :)];
				obj.AnalogIn.Timestamps = [obj.AnalogIn.Timestamps, obj.BoardADC.Timestamps];
			else
				obj.AnalogIn.ChannelNames = channelNames;
				obj.AnalogIn.Data = obj.BoardADC.Data(channels, :);
				obj.AnalogIn.Timestamps = obj.BoardADC.Timestamps;
			end
		end

		function GetDigitalData(obj, varargin)
			p = inputParser;
			addParameter(p, 'ChannelLabels', {}, @(x) iscell(x) || ischar(x)); % {'Cue', 4; 'Press', 2; 'Lick', 3; 'Reward', 5}
			addParameter(p, 'Append', false, @islogical);
			parse(p, varargin{:});
			channelLabels 	= p.Results.ChannelLabels;
			append 			= p.Results.Append;

			if ischar(channelLabels) && strcmpi(channelLabels, 'auto')
				channelLabels = transpose([{obj.BoardDigIn.Channels.CustomChannelName}; num2cell(1:length(obj.BoardDigIn.Channels))]);
			end

			if isempty(channelLabels)
				obj.DigitalEvents.ChannelNames = {};
				obj.DigitalEvents.Data = [];
				obj.DigitalEvents.Timestamps = [];
				return
			end

			channels = cell2mat(channelLabels(:, 2));
			channelNames = channelLabels(:, 1);

			if append
				obj.DigitalEvents.Data = [obj.DigitalEvents.Data, obj.BoardDigIn.Data(channels, :)];
				obj.DigitalEvents.Timestamps = [obj.DigitalEvents.Timestamps, obj.BoardDigIn.Timestamps];
			else
				obj.DigitalEvents.ChannelNames = channelNames;
				obj.DigitalEvents.Data = obj.BoardDigIn.Data(channels, :);
				obj.DigitalEvents.Timestamps = obj.BoardDigIn.Timestamps;
			end
		end

		function GetDigitalEvents(obj, clearCache)
			if nargin < 2
				clearCache = false;
			end

			% Get digital signals
			channelNames = obj.DigitalEvents.ChannelNames;
			if isempty(channelNames)
				return
			end

			tic, TetrodeRecording.TTS('	Extracting digital events...');

			for iChannel = 1:length(channelNames)
				[obj.DigitalEvents.([channelNames{iChannel}, 'On']), obj.DigitalEvents.([channelNames{iChannel}, 'Off'])] = TetrodeRecording.FindEdges(obj.DigitalEvents.Data(iChannel, :), obj.DigitalEvents.Timestamps);
			end
			
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
			obj.BoardADC = [];
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

			if isempty(channels)
				channels = [obj.Spikes.Channel];
			end

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
			channel 		= p.Results.Channel;
			level 			= p.Results.Level;
			waveformWindow 	= p.Results.WaveformWindow;
			dimension 		= p.Results.Dimension;

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

				if numWaveforms == 0
					continue
				end

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
			addParameter(p, 'MaxNumClusters', 5, @isnumeric);
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

		function PlotWaveforms(obj, channel, varargin)
			p = inputParser;
			addRequired(p, 'Channel', @isnumeric);
			addOptional(p, 'NumWaveforms', 50, @isnumeric);
			addParameter(p, 'Clusters', [], @isnumeric);
			addParameter(p, 'SelectedCluster', [], @isnumeric);
			addParameter(p, 'ReferenceCluster', [], @isnumeric);
			addParameter(p, 'YLim', [], @(x) isnumeric(x) || ischar(x));
			addParameter(p, 'FrameRate', 0, @isnumeric);
			addParameter(p, 'MaxShown', 500, @isnumeric);
			addParameter(p, 'WaveformWindow', [], @isnumeric);
			addParameter(p, 'SelectedSampleIndex', [], @isnumeric);
			addParameter(p, 'FeatureAxisLim', [], @isnumeric);
			addParameter(p, 'Polygon', [], @isnumeric);
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
			selectedSampleIndex	= p.Results.SelectedSampleIndex;
			featureAxisLim		= p.Results.FeatureAxisLim;
			polygon				= p.Results.Polygon;
			plotMean 			= p.Results.PlotMean;
			ax 					= p.Results.Ax;

			if isempty(waveformWindow)
				waveformWindow = obj.Spikes(channel).WaveformWindow;
			end
			t 			= obj.Spikes(channel).WaveformTimestamps;
			sampleIndex	= obj.Spikes(channel).SampleIndex;
			waveforms 	= obj.Spikes(channel).Waveforms;
			classes 	= obj.Spikes(channel).Cluster.Classes;
			feature 	= obj.Spikes(channel).Feature.Coeff;
			threshold 	= obj.Spikes(channel).Threshold.Threshold;

			% If blackrock, convert data to 'uV'
			if strcmpi(obj.System, 'Blackrock')
				waveforms = waveforms/4;
				threshold = double(threshold)/4;
			end

			if ~isempty(clusters)
				selected 	= ismember(classes, clusters);
				sampleIndex = sampleIndex(selected);
				waveforms 	= waveforms(selected, :);
				classes 	= classes(selected);
				feature 	= feature(selected, :);
			end

			if ~isempty(selectedSampleIndex)
				selected 	= ismember(sampleIndex, selectedSampleIndex);
				sampleIndex = sampleIndex(selected);
				waveforms 	= waveforms(selected, :);
				classes 	= classes(selected);
				feature 	= feature(selected, :);
			end

			numWaveformsTotal = size(waveforms, 1);

			if isempty(yRange)
				yRange = [min(waveforms(:)), max(waveforms(:))];
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
				for iCluster = unique(nonzeros(classes))'
					inCluster = classes == iCluster;
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
					h = scatter3(hAxes1, feature(inCluster, 1), feature(inCluster, 2), feature(inCluster, 3), 1, thisColor, 'DisplayName', dispName);
					if isempty(selectedCluster) || (~isempty(selectedCluster) && iCluster == selectedCluster)
						hLegends = [hLegends, h];
						if iCluster == selectedCluster
							title(hAxes2, dispName)
						end
					end
				end
				hold(hAxes1, 'off')

				if isempty(featureAxisLim)
					axis(hAxes1, 'auto')
				else
					xlim(hAxes1, featureAxisLim(1, :));
					ylim(hAxes1, featureAxisLim(2, :));
					zlim(hAxes1, featureAxisLim(3, :));
				end
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
				hLegends = [hLegends, line(hAxes2, 'XData', [obj.Spikes(channel).WaveformWindow(1), 0], 'YData', repmat(mean(threshold), [1, 2]), 'Color', 'k', 'LineWidth', 3, 'DisplayName', ['Threshold (', num2str(mean(threshold)), ' \muV)'])];
				legend(hLegends, 'AutoUpdate', 'off', 'Location', 'Best')
				xlim(hAxes2, waveformWindow)
				ylim(hAxes2, yRange)

				if frameRate > 0
					hTimer = timer(...
						'ExecutionMode', 'FixedSpacing',...
					 	'Period', round((1/frameRate)*1000)/1000,...
					 	'TimerFcn', {@TetrodeRecording.OnPlotChannelRefresh, hAxes2, t, waveforms, numWaveforms, numWaveformsTotal, classes}...
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
					for iCluster = unique(nonzeros(classes))'
						% Skip this cluster in some modes
						if (~isempty(selectedCluster) && ~ismember(iCluster, [selectedCluster, referenceCluster]))
							continue
						end

						thisWaveforms = waveforms(classes==iCluster, :);
						
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
							if isempty(maxShown) || maxShown == 0 || size(thisWaveforms, 1) <= maxShown
								percentShown = 100;
							else
								[~, I] = sort(sum(thisWaveforms, 2));
								thisWaveforms = thisWaveforms(I(1:ceil(size(thisWaveforms, 1)/maxShown):end), :);
							end
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
			addParameter(p, 'SelectedSampleIndex', [], @isnumeric);
			addParameter(p, 'Ax', []);
			parse(p, channels, varargin{:});
			channels 			= p.Results.Channels;
			reference 			= p.Results.Reference;
			event 				= p.Results.Event;
			exclude 			= p.Results.Exclude;
			alignTo 			= p.Results.AlignTo;
			doSort 				= p.Results.Sort;
			extendedWindow 		= p.Results.ExtendedWindow;
			clusters 			= p.Results.Clusters;
			xRange 				= p.Results.XLim;
			selectedSampleIndex	= p.Results.SelectedSampleIndex;
			ax 					= p.Results.Ax;

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
			reference = sort(reference);
			event = sort(event);
			[reference, event] = TetrodeRecording.FindFirstInTrial(reference, event, exclude);

			% Bin spikes into trials
			for iChannel = channels
				[sampleIndex, spikes, trials] = obj.GetSpikesByTrial(iChannel, 'Reference', reference, 'Event', event, 'Clusters', clusters,...
					'Window', extendedWindow, 'WindowReference', 'StartAndEnd');

				if ~isempty(selectedSampleIndex)
					selected 	= ismember(sampleIndex, selectedSampleIndex);
					sampleIndex = sampleIndex(selected);
					spikes 		= spikes(selected);
					trials 		= trials(selected);
				end

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
				legend(hAxes, 'Location', 'NorthWest');%, 'FontSize', 12);
				if ~isempty(xRange)
					xlim(hAxes, xRange);
				end
				hold off

				hAxes.UserData.PlotParams = p.Results;
			end
		end

		function [trainOn, trainOff] = GetStimTrainTimestamps(obj)

			cueOn = sort(obj.DigitalEvents.CueOn);
			pulseOn = sort(obj.DigitalEvents.StimOn);
			pulseOff = sort(obj.DigitalEvents.StimOff);

			% Get the start and end timestamps of a stim train.
			[cueOnStimTrials, trainOn] = TetrodeRecording.FindFirstInTrial(cueOn, pulseOn);
			[~, trainOff] = TetrodeRecording.FindLastInTrial(cueOn, pulseOff);
		end

		function StimData = BacthProcessStimData(obj, ptr, channels, varargin)
			
		end

		% TODO: Clean up this messy poop code
		function StimData = ReadStimData(obj, ptr, channels, varargin)
			p = inputParser;
			addParameter(p, 'Window', [-20, 20], @isnumeric); % how many milliseconds before and after stim on to read.
			addParameter(p, 'MaxTrains', [], @isnumeric);
			parse(p, varargin{:});
			readWindow = p.Results.Window * 0.001;
			maxTrains = p.Results.MaxTrains;

			[trainOn, trainOff] = obj.GetStimTrainTimestamps();
			stimOn = obj.DigitalEvents.StimOn;
			stimOff = obj.DigitalEvents.StimOff;
			if (~isempty(maxTrains))
				trainOn = trainOn(1:maxTrains);
				trainOff = trainOff(1:maxTrains);
				stimOn = stimOn(stimOn <= trainOff(end));
				stimOff = stimOff(stimOff <= trainOff(end));
			end
				
			rig = TetrodeRecording.GetRig(obj.Path);
			if (rig == 1)
				channelMap = ptr.ChannelMap.Rig1;
			else
				channelMap = ptr.ChannelMap.Rig2;
			end

			% Convert to file channel
			% Terminology:
			%   - openNSx: channels, 1 : n, however many channels were recorded, empty channels are removed with no empty spacing. 2 rigs and 10 disabled channels: read channel 1:(64-10);
			%   - channelMap = ptr.ChannelMap(ptr.SelectedChannels): index: channel index in SpikeSort result. value: channel index for openNSx
			channelMap = channelMap(ptr.SelectedChannels);
			channels = channelMap(channels);

			sampleRate = obj.FrequencyParameters.AmplifierSampleRate;

			% First do a trial read to get channel mapping info
			filename = [obj.Path, obj.Files{end}];
			% NSx = openNSx(filename, 'read', 'duration', '1', 'sec');
			% ismember([NSx.ElectrodesInfo.ElectrodeID], rawChannels)

			% Max num of samples in a train for preallocating
			maxTrainLength = max(trainOff - trainOn) + diff(readWindow);
			maxTrainLength = ceil(sampleRate * maxTrainLength);

			StimData.Data = zeros(length(trainOn), length(channels), maxTrainLength);
			StimData.Timestamps = zeros(length(trainOn), maxTrainLength);

			% First read the first ten seconds to get sysInitDelay
			NSx = openNSx(filename, 'read', 'channels', channels, 'duration', [0, 10], 'sec');
			if iscell(NSx.Data)
				sysInitDelay = length(NSx.Data{1})/sampleRate;
				disp(['Truncated data. Substracting sysInitDelay ', num2str(sysInitDelay)]);
			else
				sysInitDelay = 0;
				disp(['ALL OKAY ', num2str(sysInitDelay)]);
			end

			for iTrain = 1:length(trainOn)
				startTime = tic();
				
				TetrodeRecording.TTS(['Reading train ', num2str(iTrain), '/', num2str(length(trainOn)), '...']);
				realReadWindow = sysInitDelay + readWindow + [trainOn(iTrain), trainOff(iTrain)];
				NSx = openNSx(filename, 'read', 'channels', channels, 'duration', realReadWindow, 'sec');
				numSamples = floor(NSx.MetaTags.DataPoints);
				firstSampleIndex = floor(NSx.MetaTags.Timestamp);

				if (iTrain == 1)
					electrodesInfo = NSx.ElectrodesInfo;
				end

				StimData.Data(iTrain, :, 1:numSamples) = NSx.Data;
                
				StimData.Timestamps(iTrain, 1:numSamples) = (firstSampleIndex : firstSampleIndex + numSamples - 1) / sampleRate;

				dur = toc(startTime);
				TetrodeRecording.TTS([num2str(dur*1000), ' ms.\n']);
			end

			StimData.StimOn = stimOn;
			StimData.StimOff = stimOff;
			StimData.TrainOn = trainOn;
			StimData.TrainOff = trainOff;
			StimData.Window = readWindow;
			StimData.ElectrodesInfo = electrodesInfo;
			StimData.Filename = filename;
			StimData.StartTime = obj.StartTime;
			StimData.SampleRate = sampleRate;
		end

		function RasterStim(obj, channels, varargin)
			p = inputParser;
			addRequired(p, 'Channels', @isnumeric);
			addParameter(p, 'ExtendedWindow', [-2, 2], @isnumeric);
			addParameter(p, 'Clusters', [], @isnumeric);
			addParameter(p, 'XLim', [], @isnumeric);
			addParameter(p, 'SelectedSampleIndex', [], @isnumeric);
			addParameter(p, 'Ax', []);
			parse(p, channels, varargin{:});
			channels 			= p.Results.Channels;
			extendedWindow 		= p.Results.ExtendedWindow;
			clusters 			= p.Results.Clusters;
			xRange 				= p.Results.XLim;
			selectedSampleIndex	= p.Results.SelectedSampleIndex;
			ax 					= p.Results.Ax;

			cueOn = sort(obj.DigitalEvents.CueOn);
			pulseOn = sort(obj.DigitalEvents.StimOn);
			pulseOff = sort(obj.DigitalEvents.StimOff);

			% Get spikes between first pulse start and last pulse end
			[cueOnStimTrials, trainOn] = TetrodeRecording.FindFirstInTrial(cueOn, pulseOn);
			[~, trainOff] = TetrodeRecording.FindLastInTrial(cueOn, pulseOff);

			% Bin spikes into trials
			for iChannel = channels
				[sampleIndex, spikes, trials] = obj.GetSpikesByTrial(iChannel, 'Reference', trainOn, 'Event', trainOff, 'Clusters', clusters,...
					'Window', extendedWindow, 'WindowReference', 'StartAndEnd');

				if ~isempty(selectedSampleIndex)
					selected 	= ismember(sampleIndex, selectedSampleIndex);
					sampleIndex = sampleIndex(selected);
					spikes 		= spikes(selected);
					trials 		= trials(selected);
				end

				if isempty(ax)
					hFigure = figure('Units', 'Normalized', 'OuterPosition', [0, 0, 0.75, 1]);
					hAxes = axes(hFigure);
				else
					hAxes = ax;
					axes(hAxes);
				end

				spikesRelative = spikes - trainOn(trials);

				edges = [trainOn, trainOff(end)];
				[~, ~, trialsPulseOn] = histcounts(pulseOn, edges);
				[~, ~, trialsPulseOff] = histcounts(pulseOff, edges);
                                
                
                trainOn1 = zeros(size(pulseOn));
                trainOn2 = zeros(size(pulseOff));
                trainOn1(trialsPulseOn > 0) = trainOn(trialsPulseOn(trialsPulseOn > 0));
                trainOn2(trialsPulseOff > 0) = trainOn(trialsPulseOff(trialsPulseOff > 0));
                
				pulseOnRelative = pulseOn - trainOn1;
				pulseOffRelative = pulseOff - trainOn2;

				hold on
				plot(hAxes, spikesRelative, trials, '.',...
					'MarkerSize', 5,...
					'MarkerEdgeColor', 'k',...
					'MarkerFaceColor', 'k',...
					'LineWidth', 1.5,...
					'DisplayName', 'Spike'...
				)
				plot(hAxes, pulseOnRelative, trialsPulseOn, '.',...
					'MarkerSize', 5,...
					'MarkerEdgeColor', 'b',...
					'MarkerFaceColor', 'b',...
					'LineWidth', 1.5,...
					'DisplayName', 'Pulse On'...
				)
				plot(hAxes, pulseOffRelative, trialsPulseOff, '.',...
					'MarkerSize', 5,...
					'MarkerEdgeColor', 'r',...
					'MarkerFaceColor', 'r',...
					'LineWidth', 1.5,...
					'DisplayName', 'Pulse Off'...
				)
				title(hAxes, 'Spike raster')
				xlabel(hAxes, ['Time relative to Stim On (s)'])
				ylabel(hAxes, 'Trial')
				legend(hAxes, 'Location', 'NorthWest');
				if ~isempty(xRange)
					xlim(hAxes, xRange);
				end
				hold off

				hAxes.UserData.PlotParams = p.Results;
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
			addParameter(p, 'SelectedSampleIndex', [], @isnumeric);
			addParameter(p, 'XLim', [], @isnumeric);
			addParameter(p, 'LineStyle', '-', @ischar);
			addParameter(p, 'Ax', []);
			parse(p, channels, varargin{:});
			channels 			= p.Results.Channels;
			reference 			= p.Results.Reference;
			event 				= p.Results.Event;
			exclude 			= p.Results.Exclude;
			clusters 			= p.Results.Clusters;
			minTrialLength 		= p.Results.MinTrialLength;
			nBins 				= p.Results.Bins;
			binMethod 			= p.Results.BinMethod;
			spikeRateWindow 	= p.Results.SpikeRateWindow;
			extendedWindow 		= p.Results.ExtendedWindow;
			selectedSampleIndex	= p.Results.SelectedSampleIndex;
			xRange 				= p.Results.XLim;			
			lineStyle			= p.Results.LineStyle;			
			ax 					= p.Results.Ax;

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

			reference = sort(reference);
			event = sort(event);
			[reference, event] = TetrodeRecording.FindFirstInTrial(reference, event, exclude);

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
				[sampleIndex, spikes, trials] = obj.GetSpikesByTrial(iChannel, 'Reference', reference, 'Event', event, 'Clusters', clusters,...
					'Window', extendedWindow, 'WindowReference', 'StartAndEnd');

				if ~isempty(selectedSampleIndex)
					selected = ismember(sampleIndex, selectedSampleIndex);
					sampleIndex = sampleIndex(selected);
					spikes 		= spikes(selected);
					trials 		= trials(selected);
				end

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

					[thisColor, thisStyle] = TetrodeRecording.GetColorAndStyle(iBin, 'Styles', {lineStyle});
					plot(hAxes, thisCenters, thisSpikeRate, [thisColor, thisStyle],...
						'DisplayName', ['[', num2str(-cutOff, 2), ' s, ', num2str(edges(iBin + 1), 2), ' s] (', num2str(NTrials(iBin)), ' trials)'],...
						'LineWidth', 1);
				end
				xlabel(hAxes, ['Time relative to ', eventDisplayName, ' (s)']);
				ylabel(hAxes, 'Mean firing rate (Spikes/s)')
				legend(hAxes, 'Location', 'Best');
				title(hAxes, 'Peri-event time histogram');
				if ~isempty(xRange)
					xlim(hAxes, xRange);
				end				
				hold off
			end
		end

		function varargout = PETHistCounts(obj, channel, varargin)
			p = inputParser;
			addRequired(p, 'Channel', @isnumeric);
			addParameter(p, 'Reference', 'CueOn', @ischar);
			addParameter(p, 'Event', 'PressOn', @ischar);
			addParameter(p, 'Exclude', 'LickOn', @ischar);
			addParameter(p, 'Clusters', [], @isnumeric);
			addParameter(p, 'TrialLength', 6, @isnumeric);
			addParameter(p, 'AllowedTrialLength', [0, Inf], @isnumeric);
			addParameter(p, 'SpikeRateWindow', 100, @isnumeric); % in ms
			addParameter(p, 'ExtendedWindow', 1, @isnumeric); % Extend window after event
			addParameter(p, 'MoveOnsetCorrection', [], @isnumeric); % Modify movement onset time for each trial, list should be as long as number of trials
			parse(p, channel, varargin{:});
			channel 			= p.Results.Channel;
			reference 			= p.Results.Reference;
			event 				= p.Results.Event;
			exclude 			= p.Results.Exclude;
			clusters 			= p.Results.Clusters;
			trialLength 		= p.Results.TrialLength;
			allowedTrialLength	= p.Results.AllowedTrialLength;
			spikeRateWindow 	= p.Results.SpikeRateWindow;
			extendedWindow 		= p.Results.ExtendedWindow;
			moveOnsetCorrection	= p.Results.MoveOnsetCorrection;

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

			reference = sort(reference);
			event = sort(event);
			[reference, event, ~, ~, toRemove] = TetrodeRecording.FindFirstInTrial(reference, event, exclude);

			if ~isempty(moveOnsetCorrection)
				moveOnsetCorrection(toRemove) = [];
				moveOnsetCorrection(moveOnsetCorrection>0 | moveOnsetCorrection<-4) = 0;
				event = event + moveOnsetCorrection;
			end

			sel = (event > trialLength) & ((event - reference) >= allowedTrialLength(1)) & ((event - reference) <= allowedTrialLength(2));
			reference = reference(sel);
			event = event(sel);

			if isempty(event)
				varargout = {[], [], 0};
				return
			end

			% Get spike timestamps
			spikes 	= obj.Spikes(channel).Timestamps;
			if ~isempty(clusters)
				spikes = spikes(ismember(obj.Spikes(channel).Cluster.Classes, clusters));
			end

			numTrials = length(event);
			centers = flip(extendedWindow:-spikeRateWindow/1000:-trialLength);
			centers = centers(2:end);
			% centers = centers(2:end) - spikeRateWindow/2000;
			spikeRate = NaN(numTrials, length(centers));

			for iTrial = 1:numTrials
				thisEdges = flip((event(iTrial) + extendedWindow):(-spikeRateWindow/1000):(event(iTrial) - trialLength));
				thisSpikeRate = histcounts(spikes, thisEdges);
				thisSpikeRate = 1000*thisSpikeRate/spikeRateWindow;
				spikeRate(iTrial, :) = thisSpikeRate;
			end

			meanSpikeRate = mean(spikeRate, 1);

			varargout = {meanSpikeRate, centers, numTrials, spikeRate, spikes, event, reference};
		end

		function varargout = PSTHistCounts(obj, channel, varargin)
			p = inputParser;
			addRequired(p, 'Channel', @isnumeric);
			addParameter(p, 'Clusters', [], @isnumeric);
			addParameter(p, 'ExtendedWindow', [-1, 1], @isnumeric); % [before trainOn, after trainOff]
			addParameter(p, 'SpikeRateWindow', 10, @isnumeric); % in ms
			parse(p, channel, varargin{:});
			channel 		= p.Results.Channel;
			clusters 		= p.Results.Clusters;
			extendedWindow 	= p.Results.ExtendedWindow;
			spikeRateWindow = p.Results.SpikeRateWindow;

			cueOn = sort(obj.DigitalEvents.CueOn);
			pulseOn = sort(obj.DigitalEvents.StimOn);
			pulseOff = sort(obj.DigitalEvents.StimOff);

			% Get spikes between first pulse start and last pulse end
			[cueOnStimTrials, trainOn] = TetrodeRecording.FindFirstInTrial(cueOn, pulseOn);
			[~, trainOff] = TetrodeRecording.FindLastInTrial(cueOn, pulseOff);

			if (length(trainOn) ~= length(trainOff))
				error('Something is wrong, num of trainOff events is different from trainOn events')
			end

			% For each Train, register the stimOn, stimOff events within
			trainEdges = reshape([trainOn(:)'; trainOff(:)'], [1, numel(trainOn) + numel(trainOff)]);
			[numPulsesInTrain, ~, pulseInTrain] = histcounts(pulseOn, trainEdges);
			pulseInTrain = (pulseInTrain + 1)/2; % Which train the pulse is in
			numPulsesInTrain = numPulsesInTrain(1:2:end); % How many pulses are in each train
			trainLength = trainOff - trainOn;
			trainLength = round(trainLength*100)/100; % How long trains are, rounded to nearest 0.01s

			% Convert session time into train time
			pulseOnRel = pulseOn - trainOn(pulseInTrain);
			pulseOffRel = pulseOff - trainOn(pulseInTrain);

			% Filter spikes by cluster
			spikes 	= obj.Spikes(channel).Timestamps;
			if ~isempty(clusters)
				spikes 	= spikes(ismember(obj.Spikes(channel).Cluster.Classes, clusters));
			end			

			% Classify trains by two numbers: length of train (rounded to 10 ms) and numPulsesInTrain
			[trainTypeID, ~, trainType] = unique(trainLength*1000 + numPulsesInTrain);

			PSTH = struct([]);
			% For each train type, calculate mean firing rate (PSTH), stimOn/Off times
			for iTrainType = 1:length(trainTypeID)
				ThisPSTH.TrainType = trainTypeID(iTrainType);
				firstTrainThisType = find(trainType == iTrainType, 1);
				ThisPSTH.TrainLength = trainLength(firstTrainThisType);
				ThisPSTH.NumPulses = numPulsesInTrain(firstTrainThisType);
				ThisPSTH.PulseOn = pulseOnRel(pulseInTrain == firstTrainThisType); % When do pulses come on for this type of stim train
				ThisPSTH.PulseOff = pulseOffRel(pulseInTrain == firstTrainThisType); % When do pulses go off for this type of stim train

				% Average spike rates for all trains in this train type
				trainsThisType = find(trainType == iTrainType);
				numTrainsThisType = length(trainsThisType);

				% 5 repetitions or less is too few, discard
				if (numTrainsThisType < 5)
					continue
				end

				ThisPSTH.NumTrains = numTrainsThisType;
				% Find binned spike rates from each train belonging to this type, then average them
				relEdgesThisType = extendedWindow(1):(spikeRateWindow/1000):(trainLength(firstTrainThisType) + extendedWindow(2));
				relCentersThisType = (relEdgesThisType(1:end - 1) + relEdgesThisType(2:end))/2;
				spikeRateThisType = NaN(length(trainsThisType), length(relCentersThisType));
				for iTrainRel = 1:length(trainsThisType)
					iTrain = trainsThisType(iTrainRel);
					thisEdges = (trainOn(iTrain) + extendedWindow(1)):(spikeRateWindow/1000):(trainOff(iTrain) + extendedWindow(2));
					thisSpikeRate = histcounts(spikes, thisEdges);
					thisSpikeRate = 1000*thisSpikeRate/spikeRateWindow;
					spikeRateThisType(iTrain, 1:length(thisSpikeRate)) = thisSpikeRate;
				end
				spikeRateThisType = nanmean(spikeRateThisType, 1);
				ThisPSTH.SpikeRate = spikeRateThisType;
				ThisPSTH.Timestamps = relEdgesThisType(1:end-1);

				if (isempty (PSTH))
					PSTH = ThisPSTH;
				else
					PSTH = [PSTH, ThisPSTH];
				end
			end

			varargout = {PSTH};
		end

		function PlotAllChannels(obj, varargin)
			p = inputParser;
			addParameter(p, 'Channels', [], @isnumeric);
			addParameter(p, 'PercentShown', 5, @isnumeric); % What percentage of waveforms are plotted (0 - 100)
			addParameter(p, 'MaxShown', 200, @isnumeric);
			addParameter(p, 'Fontsize', 8, @isnumeric);
			addParameter(p, 'PlotMethod', 'all', @ischar); % 'all', 'mean'
			addParameter(p, 'YLim', 'auto', @(x) isnumeric(x) || ischar(x)); % [-200, 200], or 'auto'
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
				if isempty(obj.Spikes(iChannel).Waveforms)
					continue
				end

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

		%% PlotUnit: Plot a single unit
		function varargout = PlotUnitSimple(obj, channel, unit, varargin)
			p = inputParser;
			addRequired(p, 'Channel', @isnumeric);
			addRequired(p, 'Unit', @isnumeric);
			addParameter(p, 'Reference', 'CueOn', @ischar);
			addParameter(p, 'Event', 'PressOn', @ischar);
			addParameter(p, 'Exclude', 'LickOn', @ischar);
			addParameter(p, 'Bins', 3, @isnumeric);
			addParameter(p, 'BinMethod', 'percentile', @ischar);
			addParameter(p, 'SpikeRateWindow', 100, @isnumeric);
			addParameter(p, 'RasterXLim', [-5, 0], @isnumeric);
			addParameter(p, 'RasterXLimStim', [-0.5, 0.5], @isnumeric);
			addParameter(p, 'PlotType', 'Raster', @ischar); % 'Raster', 'PETH'
			addParameter(p, 'AlignTo', 'Event', @ischar);
			addParameter(p, 'Position', [0, 0, 0.5, 1], @isnumeric);

			parse(p, channel, unit, varargin{:});
			iChannel 			= p.Results.Channel;
			iUnit 				= p.Results.Unit;
			reference 			= p.Results.Reference;
			event 				= p.Results.Event;
			exclude 			= p.Results.Exclude;
			bins 				= p.Results.Bins;
			binMethod 			= p.Results.BinMethod;
			spikeRateWindow 	= p.Results.SpikeRateWindow;
			rasterXLim 			= p.Results.RasterXLim;
			rasterXLimStim 		= p.Results.RasterXLimStim;

			hFigure = figure('Units', 'Normalized', 'Position', p.Results.Position, 'GraphicsSmoothing', 'on');
			hAx1 = subplot(2,1,1);
			hAx2 = subplot(2,1,2);

			switch lower(p.Results.PlotType)
				case 'raster'
					obj.Raster(iChannel, reference, event, exclude,...
						'Clusters', iUnit, 'XLim', rasterXLim, 'Ax', hAx1, 'ExtendedWindow', [-2, 2],...
						'AlignTo', p.Results.AlignTo);
				case 'peth'
					obj.PETH(iChannel, reference, event, exclude,...
						'Clusters', iUnit, 'XLim', rasterXLim, 'Ax', hAx1, 'ExtendedWindow', [-2, 2],...
						'Bins', bins, 'BinMethod', binMethod, 'SpikeRateWindow', spikeRateWindow);
			end

			obj.RasterStim(iChannel, 'Clusters', iUnit, 'Ax', hAx2, 'XLim', rasterXLimStim);

			expName = obj.GetExpName();				
			displayName = [expName, ' (Channel ', num2str(iChannel), ' Unit ', num2str(iUnit), ')'];
			hTitle = suptitle(displayName);

			title(hAx1, 'Lever Press')
			title(hAx2, 'Opto stim')

			varargout = {hFigure, hAx1, hAx2, hTitle};
		end
		

		% PlotChannel(iChannel, 'Reference', reference, 'Event', event, 'Exclude', exclude, 'Clusters', clusters, 'ReferenceCluster', referenceCluster, 'WaveformWindow', waveformWindow, 'ExtendedWindow', extendedWindow, 'MinTrialLength', minTrialLength, 'Bins', bins, 'BinMethod', binMethod, 'SpikeRateWindow', spikeRateWindow, 'RasterXLim', rasterXLim, 'WaveformYLim', waveformYLim, 'FontSize', fontSize, 'PrintMode', printMode, 'FrameRate', frameRate, 'Fig', hFigure)
		function varargout = PlotChannel(obj, channel, varargin)
			p = inputParser;
			addRequired(p, 'Channel', @isnumeric);
			addParameter(p, 'Reference', 'CueOn', @ischar);
			addParameter(p, 'Event', 'PressOn', @ischar);
			addParameter(p, 'Exclude', 'LickOn', @ischar);
			addParameter(p, 'Event2', '', @ischar);
			addParameter(p, 'Exclude2', '', @ischar);
			addParameter(p, 'Clusters', [], @isnumeric);
			addParameter(p, 'ReferenceCluster', [], @isnumeric);
			addParameter(p, 'WaveformWindow', [], @isnumeric);
			addParameter(p, 'ExtendedWindow', [0, 0], @isnumeric);
			addParameter(p, 'MinTrialLength', 0.5, @isnumeric);
			addParameter(p, 'Bins', 4, @isnumeric);
			addParameter(p, 'BinMethod', 'percentile', @ischar);
			addParameter(p, 'SpikeRateWindow', 100, @isnumeric);
			addParameter(p, 'RasterXLim', [], @isnumeric);
			addParameter(p, 'WaveformYLim', 'auto', @(x) isnumeric(x) || ischar(x)); % [-200, 200], 'auto', []
			addParameter(p, 'FontSize', 8, @isnumeric);
			addParameter(p, 'PrintMode', false, @islogical);
			addParameter(p, 'FrameRate', 0, @isnumeric);
			addParameter(p, 'PlotStim', false, @islogical);
			addParameter(p, 'Fig', []);
			parse(p, channel, varargin{:});
			iChannel 			= p.Results.Channel;
			reference 			= p.Results.Reference;
			event 				= p.Results.Event;
			exclude 			= p.Results.Exclude;
			event2 				= p.Results.Event2;
			exclude2 			= p.Results.Exclude2;
			clusters 			= p.Results.Clusters;
			referenceCluster	= p.Results.ReferenceCluster;
			waveformWindow 		= p.Results.WaveformWindow;
			extendedWindow 		= p.Results.ExtendedWindow;
			minTrialLength 		= p.Results.MinTrialLength;
			bins 				= p.Results.Bins;
			binMethod 			= p.Results.BinMethod;
			spikeRateWindow 	= p.Results.SpikeRateWindow;
			rasterXLim 			= p.Results.RasterXLim;
			waveformYLim		= p.Results.WaveformYLim;
			fontSize 			= p.Results.FontSize;
			printMode 			= p.Results.PrintMode;
			frameRate 			= p.Results.FrameRate;
			plotStim 			= p.Results.PlotStim;
			h.Figure 			= p.Results.Fig;

			allChannels = [obj.Spikes.Channel];
			if isempty(iChannel)
				iChannel = allChannels(1);
			end

			xRatio 	= 0.2;
			yRatio 	= 0.6;
			xMargin = 0.02;
			yMargin = 0.03;
			fMargin = 0.04;
			buttonWidth = 0.04;
			buttonHeight = 0.025;
			buttonXSpacing = 0.006;
			buttonYSpacing = 0.0075;
			hUp 	= (1 - 2*fMargin)*yRatio - 2*yMargin;
			hUpHalf = (hUp - 2*yMargin)/2;
			hDown 	= (1 - 2*fMargin)*(1 - yRatio) - 2*yMargin;
			hHalf 	= (1 - 2*fMargin - 5*yMargin)*0.5;
			w 		= 1 - 2*fMargin - 2*xMargin;
			wLeft	= (1 - 2*fMargin)*xRatio - 2*xMargin;
			wRight	= (1 - 2*fMargin)*(1 - xRatio) - 2*xMargin;

			expName = obj.GetExpName();				
			displayName = [expName, ' (Channel ', num2str(iChannel), ')'];

			if ~isempty(h.Figure)
				clf(h.Figure)
			else
				h.Figure = figure('Units', 'Normalized', 'Position', [0, 0, 1, 1], 'GraphicsSmoothing', 'on');
				h.Figure.UserData.PlotMean = true;
				h.Figure.UserData.ReferenceCluster = referenceCluster;
			end
			h.Figure.UserData.SelectedSampleIndex = []; % Selectively plot waveforms by sampleIndex. Reset everytime.
			h.Figure.Name = displayName;

			if printMode
				fontSize = 14;
				h.Figure.Units = 'pixels';
				% h.Figure.InnerPosition = [0, 0, 1332, 999];
				h.Figure.InnerPosition = [0, 0, 1776, 999];
			end
			set(h.Figure, 'DefaultAxesFontSize', fontSize);

			h.Waveform 	= subplot('Position', [fMargin + xMargin, fMargin + 5*yMargin + hDown + hUpHalf, wLeft, hUpHalf], 'Tag', 'Waveform');
			h.PCA 		= subplot('Position', [fMargin + xMargin, fMargin + 3*yMargin + hDown, wLeft, hUpHalf], 'Tag', 'PCA');
			h.Raster	= subplot('Position', [fMargin + xMargin + wLeft + 2*xMargin, fMargin + 3*yMargin + hHalf, wRight, hHalf], 'Tag', 'Raster');
			if (~isempty(event2) || plotStim)
				h.Raster2 = subplot('Position', [fMargin + xMargin + wLeft + 2*xMargin, fMargin + yMargin, wRight, hHalf], 'Tag', 'Raster2');
			else
				h.Raster2 = [];
			end
			h.PETH 		= subplot('Position', [fMargin + xMargin, fMargin + yMargin, wLeft, hDown], 'Tag', 'PETH');

			if isempty(clusters)
				h.Figure.UserData.SelectedClusters = unique(obj.Spikes(iChannel).Cluster.Classes);
			else
				h.Figure.UserData.SelectedClusters = clusters;
			end

			hTitle = suptitle(displayName);

			%----------------------------------------------------
			%		Buttons: 1st row (clusters)
			%----------------------------------------------------
			hAxesTextCluster = axes('Position', [buttonWidth, 1 - buttonHeight - buttonYSpacing, buttonWidth, min(yMargin, buttonHeight)], 'Visible', 'off');
			hTextCluster = text(hAxesTextCluster, 0.5, 0.5, 'Clusters:', 'Tag', 'HideWhenSaving',...
				'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
			hPrev = hAxesTextCluster;

			hButtonSelClusters = uicontrol(...
				'Style', 'pushbutton',...
				'String', 'Select ...',...
				'Callback', {@obj.GUIPlotClusters, iChannel, p, h},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Position', hPrev.Position);
			hPrev = hButtonSelClusters;
			hPrev.Position(1) = hPrev.Position(1) + hPrev.Position(3) + buttonXSpacing;

			hButtonSelRef = uicontrol(...
				'Style', 'pushbutton',...
				'String', 'Ref ...',...
				'Callback', {@obj.GUISelRef, iChannel, p, h},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Position', hPrev.Position);
			hPrev = hButtonSelRef;
			hPrev.Position(1) = hPrev.Position(1) + hPrev.Position(3) + buttonXSpacing;

			hButtonPlotAllClusters = uicontrol(...
				'Style', 'pushbutton',...
				'String', 'Expand ...',...
				'Callback', {@obj.GUIPlotAllClusters, iChannel, h},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Position', hPrev.Position);
			hPrev = hButtonPlotAllClusters;
			hPrev.Position(1) = hPrev.Position(1) + hPrev.Position(3) + buttonXSpacing;

			hButtonMerge = uicontrol(...
				'Style', 'pushbutton',...
				'String', 'Merge ...',...
				'Callback', {@obj.GUIMergeClusters, iChannel, p, h},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Position', hPrev.Position);
			hPrev = hButtonMerge;
			hPrev.Position(1) = hPrev.Position(1) + hPrev.Position(3) + buttonXSpacing;

			hButtonRemove = uicontrol(...
				'Style', 'pushbutton',...
				'String', 'Remove ...',...
				'Callback', {@obj.GUIRemoveClusters, iChannel, p, h},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Position', hPrev.Position);
			hPrev = hButtonRemove;
			hPrev.Position(1) = hPrev.Position(1) + hPrev.Position(3) + buttonXSpacing;

			hButtonDecimate = uicontrol(...
				'Style', 'pushbutton',...
				'String', 'Decimate ...',...
				'Callback', {@obj.GUIDecimateClusters, iChannel, p, h},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Position', hPrev.Position);
			hPrev = hButtonDecimate;
			hPrev.Position(1) = hPrev.Position(1) + hPrev.Position(3) + buttonXSpacing;

			hButtonRecluster = uicontrol(...
				'Style', 'pushbutton',...
				'String', 'Recluster ...',...
				'Callback', {@obj.GUIRecluster, iChannel, p, h},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Position', hPrev.Position);
			hPrev = hButtonRecluster;
			hPrev.Position(1) = hPrev.Position(1) + hPrev.Position(3) + buttonXSpacing;

			%----------------------------------------------------
			%		Buttons: 2nd row (spikes)
			%----------------------------------------------------
			hPrev = hAxesTextCluster;
			hAxesTextCluster = axes('Position', hPrev.Position, 'Visible', 'off');
			hTextCluster = text(hAxesTextCluster, 0.5, 0.5, 'Spikes:', 'Tag', 'HideWhenSaving',...
				'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
			hAxesTextCluster.Position(2) = hAxesTextCluster.Position(2) - hAxesTextCluster.Position(4) - buttonYSpacing;
			hPrev = hAxesTextCluster;

			hButtonSelSpikes = uicontrol(...
				'Style', 'pushbutton',...
				'String', 'Select ...',...
				'Callback', {@obj.GUISelectSpikes, iChannel, p, h},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Position', hPrev.Position);
			hPrev = hButtonSelSpikes;
			hPrev.Position(1) = hPrev.Position(1) + hPrev.Position(3) + buttonXSpacing;

			hButtonExcludeSpikes = uicontrol(...
				'Style', 'pushbutton',...
				'String', 'Exclude ...',...
				'Callback', {@obj.GUIExcludeSpikes, iChannel, p, h},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Position', hPrev.Position);
			hPrev = hButtonExcludeSpikes;
			hPrev.Position(1) = hPrev.Position(1) + hPrev.Position(3) + buttonXSpacing;

			hButtonSelReset = uicontrol(...
				'Style', 'pushbutton',...
				'String', 'Reset',...
				'Callback', {@obj.GUISelectAllSpikes, iChannel, p, h},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Position', hPrev.Position);
			hPrev = hButtonSelReset;
			hPrev.Position(1) = hPrev.Position(1) + hPrev.Position(3) + buttonXSpacing;

			hButtonDeleteSpikes = uicontrol(...
				'Style', 'pushbutton',...
				'String', 'Delete ...',...
				'Callback', {@obj.GUIDeleteSpikes, iChannel, p, h},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Position', hPrev.Position);
			hPrev = hButtonDeleteSpikes;
			hPrev.Position(1) = hPrev.Position(1) + hPrev.Position(3) + buttonXSpacing;

			hButtonSpikesToCluster = uicontrol(...
				'Style', 'pushbutton',...
				'String', 'To cluster ...',...
				'Callback', {@obj.GUISpikesToCluster, iChannel, p, h},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Position', hPrev.Position);
			hPrev = hButtonSpikesToCluster;
			hPrev.Position(1) = hPrev.Position(1) + hPrev.Position(3) + buttonXSpacing;

			%----------------------------------------------------
			%		Buttons: vertical (plot options)
			%----------------------------------------------------
			hButtonPlotMean = uicontrol(...
				'Style', 'togglebutton',...
				'String', 'Mean',...
				'Callback', {@obj.GUIPlotMean, iChannel, p, h},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Value', h.Figure.UserData.PlotMean,...
				'Position', [buttonXSpacing, h.Waveform.Position(2) + h.Waveform.Position(4) - buttonHeight, buttonWidth, buttonHeight]);
			hPrev = hButtonPlotMean;

			hButtonYLim = uicontrol(...
				'Style', 'pushbutton',...
				'String', 'YLim',...
				'Callback', {@(~, ~, ax) ylim(ax, 'auto'), h.Waveform},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Position', hPrev.Position);
			hPrev = hButtonYLim;
			hPrev.Position(2) = hPrev.Position(2) - hPrev.Position(4) - buttonYSpacing;

			hButtonSavePlot = uicontrol(...
				'Style', 'pushbutton',...
				'String', 'Clipboard',...
				'Callback', {@obj.GUISavePlot, h},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Position', hPrev.Position);
			hPrev = hButtonSavePlot;
			hPrev.Position(2) = hPrev.Position(2) - hPrev.Position(4) - buttonYSpacing;

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

			hButtonNextChn = uicontrol(h.Figure,...
				'Style', 'pushbutton',...
				'String', 'Next Chn',...
				'Callback', {@(~, ~) obj.PlotChannel(nextChn, 'Reference', reference, 'Event', event, 'Exclude', exclude, 'Event2', event2, 'Exclude2', exclude2, 'Clusters', clusters, 'ReferenceCluster', referenceCluster, 'WaveformWindow', waveformWindow, 'ExtendedWindow', extendedWindow, 'MinTrialLength', minTrialLength, 'Bins', bins, 'BinMethod', binMethod, 'SpikeRateWindow', spikeRateWindow, 'RasterXLim', rasterXLim, 'WaveformYLim', waveformYLim, 'FontSize', fontSize, 'PrintMode', printMode, 'FrameRate', frameRate, 'PlotStim', plotStim, 'Fig', h.Figure)},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Position', [1 - xMargin - 2*buttonWidth, 1 - yMargin, buttonWidth, min(yMargin, buttonHeight)]);
			hPrev = hButtonNextChn;

			hAxesTextCurChn = axes('Position', hPrev.Position, 'Visible', 'off');
			hTextCurChn = text(hAxesTextCurChn, 0.5, 0.5,...
				['Chn ', num2str(find(allChannels == iChannel)), '/', num2str(length(allChannels))],...
				'Tag', 'HideWhenSaving', 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
			hAxesTextCurChn.Position(1) = hAxesTextCurChn.Position(1) - hAxesTextCurChn.Position(3) - buttonXSpacing;
			hPrev = hAxesTextCurChn;

			hButtonDeleteChn = uicontrol(h.Figure,...
				'Style', 'pushbutton',...
				'String', 'Delete Chn',...
				'Callback', {@obj.GUIDeleteChannel, iChannel, nextChn, p, h},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Position', hPrev.Position);
			hButtonDeleteChn.Position(2) = hButtonDeleteChn.Position(2) - hButtonDeleteChn.Position(4);

			hButtonPrevChn = uicontrol(h.Figure,...
				'Style', 'pushbutton',...
				'String', 'Prev Chn',...
				'Callback', {@(~, ~) obj.PlotChannel(prevChn, 'Reference', reference, 'Event', event, 'Exclude', exclude, 'Event2', event2, 'Exclude2', exclude2, 'Clusters', clusters, 'ReferenceCluster', referenceCluster, 'WaveformWindow', waveformWindow, 'ExtendedWindow', extendedWindow, 'MinTrialLength', minTrialLength, 'Bins', bins, 'BinMethod', binMethod, 'SpikeRateWindow', spikeRateWindow, 'RasterXLim', rasterXLim, 'WaveformYLim', waveformYLim, 'FontSize', fontSize, 'PrintMode', printMode, 'FrameRate', frameRate, 'PlotStim', plotStim, 'Fig', h.Figure)},...
				'BusyAction', 'cancel',...
				'Units', 'Normalized',...
				'Position', hPrev.Position);
			hButtonPrevChn.Position(1) = hButtonPrevChn.Position(1) - hButtonPrevChn.Position(3) - buttonXSpacing;
			hPrev = hButtonPrevChn;

			obj.GUIBusy(h.Figure, true);
			obj.ReplotChannel(iChannel, p, h);

			if (~isempty(event2) || plotStim)
				numTrials = max(h.Raster.Children(1).YData);
				numTrials2 = max(h.Raster2.Children(1).YData);
				hRasterUp = (1 - 2*fMargin - 5*yMargin)*numTrials/(numTrials + numTrials2);
				hRasterDown = (1 - 2*fMargin - 5*yMargin)*numTrials2/(numTrials + numTrials2);
				h.Raster.Position(2) = fMargin + 3*yMargin + hRasterDown;
				h.Raster2.Position(2) = fMargin + yMargin;
				h.Raster.Position(4) = hRasterUp;
				h.Raster2.Position(4) = hRasterDown;
			else
				h.Raster.Position(2) = fMargin + yMargin;
				h.Raster.Position(4) = (1 - 2*fMargin - 2*yMargin); 
			end

			obj.GUIBusy(h.Figure, false);
			set(h.Figure, 'Visible', 'on');

			varargout = {h};			
		end

		function PlotAllClusters(obj, channel, varargin)
			p = inputParser;
			addRequired(p, 'Channel', @isnumeric);
			addParameter(p, 'Clusters', [], @isnumeric);
			addParameter(p, 'ReferenceCluster', [], @isnumeric);
			addParameter(p, 'WaveformWindow', [], @isnumeric);
			addParameter(p, 'WaveformYLim', [-200, 200], @(x) isnumeric(x) || ischar(x));
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

		function ReplotChannel(obj, iChannel, p, h)
			reference 		= p.Results.Reference;
			event 			= p.Results.Event;
			exclude 		= p.Results.Exclude;
			event2 			= p.Results.Event2;
			exclude2 		= p.Results.Exclude2;
			waveformWindow 	= p.Results.WaveformWindow;
			extendedWindow 	= p.Results.ExtendedWindow;
			minTrialLength 	= p.Results.MinTrialLength;
			bins 			= p.Results.Bins;
			binMethod 		= p.Results.BinMethod;
			spikeRateWindow = p.Results.SpikeRateWindow;
			rasterXLim 		= p.Results.RasterXLim;
			waveformYLim	= p.Results.WaveformYLim;
			frameRate 		= p.Results.FrameRate;
			plotStim 		= p.Results.PlotStim;

			clusters = h.Figure.UserData.SelectedClusters;
			plotMean = h.Figure.UserData.PlotMean;
			referenceCluster = h.Figure.UserData.ReferenceCluster;
			selectedSampleIndex = h.Figure.UserData.SelectedSampleIndex;

			% Stop refresh if frameRate > 0
			if isfield(h.Waveform.UserData, 'hTimer')
				if isvalid(h.Waveform.UserData.hTimer)
					stop(h.Waveform.UserData.hTimer);
					delete(h.Waveform.UserData.hTimer);
				end
			end

			% Clear axes
			hNames = fieldnames(h);
			for iHandle = 1:length(hNames)
				if isgraphics(h.(hNames{iHandle}), 'Axes')
					cla(h.(hNames{iHandle}))
				end
			end

			% Replot newly selected clusters
			if isgraphics(h.Raster, 'Axes')
				obj.Raster(iChannel, reference, event, exclude, 'Clusters', clusters,...
					'AlignTo', 'Event', 'ExtendedWindow', extendedWindow, 'XLim', rasterXLim,...
					'SelectedSampleIndex', selectedSampleIndex, 'Sort', true,...
					'Ax', h.Raster);
				obj.PETH(iChannel, reference, event, exclude, 'Clusters', clusters,...
					'MinTrialLength', minTrialLength, 'Bins', bins, 'BinMethod', binMethod,...
					'SpikeRateWindow', spikeRateWindow, 'ExtendedWindow', extendedWindow,...
					'SelectedSampleIndex', selectedSampleIndex,...
					'Ax', h.PETH);
			end

			if isgraphics(h.Raster2, 'Axes')
				if ~plotStim
					obj.Raster(iChannel, reference, event2, exclude2, 'Clusters', clusters,...
						'AlignTo', 'Event', 'ExtendedWindow', extendedWindow, 'XLim', rasterXLim,...
						'SelectedSampleIndex', selectedSampleIndex, 'Sort', true,...
						'Ax', h.Raster2);			
					obj.PETH(iChannel, reference, event2, exclude2, 'Clusters', clusters,...
						'MinTrialLength', minTrialLength, 'Bins', bins, 'BinMethod', binMethod,...
						'SpikeRateWindow', spikeRateWindow, 'ExtendedWindow', extendedWindow,...
						'SelectedSampleIndex', selectedSampleIndex,...
						'LineStyle', '--',...
						'Ax', h.PETH);
				else
					obj.RasterStim(iChannel, 'Clusters', clusters,...
						'ExtendedWindow', extendedWindow, 'XLim', [-1, 1],...
						'SelectedSampleIndex', selectedSampleIndex,...
						'Ax', h.Raster2);
				end
			end

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
			
			obj.PlotWaveforms(iChannel, 'Clusters', clusters, 'WaveformWindow', waveformWindow,...
				'YLim', waveformYLim, 'FrameRate', frameRate, 'PlotMean', plotMean,...
				'ReferenceCluster', referenceCluster, 'SelectedCluster', selectedCluster,...
				'SelectedSampleIndex', selectedSampleIndex,...
				'Ax', [h.PCA, h.Waveform]);
		end

		function GUIBusy(obj, hFigure, busy)
			hButtons = findobj(hFigure.Children, 'Type', 'uicontrol');

			if busy
				set(hButtons, 'Enable', 'off')
			else
				set(hButtons, 'Enable', 'on')
			end
		end

		function GUIPlotClusters(obj, hButton, evnt, iChannel, p, h)
			obj.GUIBusy(h.Figure, true);
			liststr = cellfun(@num2str, num2cell(unique(obj.Spikes(iChannel).Cluster.Classes)), 'UniformOutput', false);

			[clusters, ok] = listdlg(...
				'PromptString', 'Plot clusters:',...
				'SelectionMode', 'multiple',...
				'OKString', 'Plot',...
				'ListString', liststr,...
				'InitialValue', [],...
				'ListSize', [250, 150]);

			if (ok && ~isempty(clusters))
				clusters = cellfun(@str2num, liststr(clusters));
				% Replot clusters
				h.Figure.UserData.SelectedClusters = clusters;
				obj.ReplotChannel(iChannel, p, h);
			end
			obj.GUIBusy(h.Figure, false);
		end

		function GUIRemoveClusters(obj, hButton, evnt, iChannel, p, h)
			obj.GUIBusy(h.Figure, true);
			liststr = cellfun(@num2str, num2cell(unique(obj.Spikes(iChannel).Cluster.Classes)), 'UniformOutput', false);

			[clusters, ok] = listdlg(...
				'PromptString', 'Remove clusters:',...
				'SelectionMode', 'multiple',...
				'OKString', 'Remove',...
				'ListString', liststr,...
				'InitialValue', [],...
				'ListSize', [250, 150]);

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
					h.Figure.UserData.SelectedClusters = unique(obj.Spikes(iChannel).Cluster.Classes);
					h.Figure.UserData.SelectedSampleIndex = [];

					% Replot clusters
					obj.ReplotChannel(iChannel, p, h);
				end
			end
			obj.GUIBusy(h.Figure, false);
		end

		function GUIDecimateClusters(obj, hButton, evnt, iChannel, p, h)
			obj.GUIBusy(h.Figure, true);
			liststr = cellfun(@num2str, num2cell(unique(obj.Spikes(iChannel).Cluster.Classes)), 'UniformOutput', false);

			[clusters, ok] = listdlg(...
				'PromptString', 'Decimate clusters:',...
				'SelectionMode', 'multiple',...
				'OKString', 'Decimate',...
				'ListString', liststr,...
				'InitialValue', [],...
				'ListSize', [250, 150]);

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

					h.Figure.UserData.SelectedSampleIndex = [];

					% Replot clusters
					obj.ReplotChannel(iChannel, p, h);
				end
			end
			obj.GUIBusy(h.Figure, false);
		end

		function GUIMergeClusters(obj, hButton, evnt, iChannel, p, h)
			obj.GUIBusy(h.Figure, true);
			liststr = cellfun(@num2str, num2cell(unique(obj.Spikes(iChannel).Cluster.Classes)), 'UniformOutput', false);

			[clusters, ok] = listdlg(...
				'PromptString', 'Merge clusters:',...
				'SelectionMode', 'multiple',...
				'OKString', 'Merge',...
				'ListString', liststr,...
				'InitialValue', [],...
				'ListSize', [250, 150]);

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
					mergeList = [{clusters}, num2cell(allClusters(~ismember(allClusters, clusters)))];
					obj.ClusterMerge(iChannel, mergeList);
					h.Figure.UserData.SelectedClusters = unique(obj.Spikes(iChannel).Cluster.Classes);

					% Replot clusters
					obj.ReplotChannel(iChannel, p, h);
				end
			end
			obj.GUIBusy(h.Figure, false);
		end

		function GUIRecluster(obj, hButton, evnt, iChannel, p, h)
			obj.GUIBusy(h.Figure, true);
			liststr = cellfun(@num2str, num2cell(unique(obj.Spikes(iChannel).Cluster.Classes)), 'UniformOutput', false);

			[clusters, ok] = listdlg(...
				'PromptString', 'Recluster clusters:',...
				'SelectionMode', 'multiple',...
				'OKString', 'Recluster',...
				'ListString', liststr,...
				'InitialValue', [],...
				'ListSize', [250, 150]);

			if (ok && ~isempty(clusters))
				clusters = cellfun(@str2num, liststr(clusters));
				answer = questdlg(...
					['Recluster selected clusters (', mat2str(clusters), ')?'],...
					'Recluster Clusters',...
					'Recluster', 'Refeature & recluster', 'Cancel',...
					'Cancel');

				if strcmp(answer, 'Cancel')
					obj.GUIBusy(h.Figure, false);
					return					
				end

				% Re-extract features
				if strcmp(answer, 'Refeature & recluster')
					params = obj.Spikes(iChannel).Feature.Parameters;
					answerFeatureMethod = questdlg(...
						'Select feature extraction method',...
						'Refeature',...
						'PCA', 'WaveletTransform', 'Cancel',...
						'Cancel');

					if strcmp(answerFeatureMethod, 'Cancel')
						obj.GUIBusy(h.Figure, false);
						return					
					end
				end

				% Recluster
				if ismember(answer, {'Recluster', 'Refeature & recluster'})
					% Recluster selected clusters
					clusterMethod = obj.Spikes(iChannel).Cluster.Method;
					% obj.Cluster(iChannel, 'Clusters', clusters, 'Method', clusterMethod);
					answerClusterMethod = inputdlg({'Cluster method:', 'Number of clusters:'}, 'Recluster', 1, {clusterMethod, ''});
					if ismember(lower(answerClusterMethod{1}), {'kmeans', 'spc', 'gaussian'})
						clusterMethod = lower(answerClusterMethod{1});
					else
						warning(['Unrecognized clustering method. Using ''', clusterMethod, ''' instead.']);
					end
					numClusters = [];
					try
						numClusters = str2num(answerClusterMethod{2});
					catch
						numClusters = [];
					end
				end

				if strcmp(answer, 'Refeature & recluster')
					obj.FeatureExtract(iChannel, 'Method', answerFeatureMethod, 'Dimension', params.Dimension, 'WaveformWindow', params.WaveformWindow);
				end
				if ismember(answer, {'Recluster', 'Refeature & recluster'})
					obj.Cluster(iChannel, 'Clusters', clusters, 'Method', clusterMethod, 'NumClusters', numClusters);
					h.Figure.UserData.SelectedClusters = unique(obj.Spikes(iChannel).Cluster.Classes);
				end
				% Replot clusters
				obj.ReplotChannel(iChannel, p, h);
			end
			obj.GUIBusy(h.Figure, false);
		end

		function GUISelRef(obj, hButton, evnt, iChannel, p, h)
			obj.GUIBusy(h.Figure, true);
			liststr = cellfun(@num2str, num2cell(unique(obj.Spikes(iChannel).Cluster.Classes)), 'UniformOutput', false);

			[referenceCluster, ok] = listdlg(...
				'PromptString', 'Select cluster as reference:',...
				'SelectionMode', 'single',...
				'OKString', 'Reference',...
				'CancelString', 'No Reference',...
				'ListString', liststr,...
				'InitialValue', [],...
				'ListSize', [250, 150]);

			if ok
				h.Figure.UserData.ReferenceCluster = cellfun(@str2num, liststr(referenceCluster));
			else
				h.Figure.UserData.ReferenceCluster = [];
			end
			obj.ReplotChannel(iChannel, p, h);

			obj.GUIBusy(h.Figure, false);
		end

		function GUIPlotMean(obj, hButton, evnt, iChannel, p, h)
			obj.GUIBusy(h.Figure, true);
			if logical(hButton.Value) ~= h.Figure.UserData.PlotMean
				h.Figure.UserData.PlotMean = logical(hButton.Value);
				obj.ReplotChannel(iChannel, p, h);
			end
			obj.GUIBusy(h.Figure, false);
		end

		function GUIPlotAllClusters(obj, hButton, evnt, iChannel, h)
			obj.GUIBusy(h.Figure, true);
			liststr = cellfun(@num2str, num2cell(unique(obj.Spikes(iChannel).Cluster.Classes)), 'UniformOutput', false);

			[referenceCluster, ok] = listdlg(...
				'PromptString', 'Select cluster as reference:',...
				'SelectionMode', 'single',...
				'OKString', 'Reference',...
				'ListString', liststr,...
				'InitialValue', [],...
				'ListSize', [250, 150]);

			if ok
				referenceCluster = cellfun(@str2num, liststr(referenceCluster));
				% Expand and plot all clusters
				obj.PlotAllClusters(iChannel, 'Clusters', h.Figure.UserData.SelectedClusters, 'ReferenceCluster', referenceCluster, 'FigMain', h.Figure);
			end
			obj.GUIBusy(h.Figure, false);
		end

		function GUISavePlot(obj, hButton, evnt, h, varargin)
			p = inputParser;
			addParameter(p, 'Filename', '', @ischar); % '' to copy to clipboard
			addParameter(p, 'Reformat', 'Raw', @ischar); % 'Raw', 'PETH', 'Raster', 'RasterAndPETH', 'RasterAndPETHAndWaveform'
			addParameter(p, 'CopyLegend', true, @islogical);
			addParameter(p, 'CopyLabel', true, @islogical);
			parse(p, varargin{:});
			filename 	= p.Results.Filename;
			reformat 	= p.Results.Reformat;
			copyLegend 	= p.Results.CopyLegend;
			copyLabel 	= p.Results.CopyLabel;

			hButtons 	= findobj(h.Figure.Children, 'Type', 'uicontrol');
			hTexts 		= findobj(h.Figure.Children, 'Type', 'Text', 'Tag', 'HideWhenSaving');

			hWaveform 	= findobj(h.Figure, 'Tag', 'Waveform');
			hPCA 		= findobj(h.Figure, 'Tag', 'PCA');
			hRaster 	= findobj(h.Figure, 'Tag', 'Raster');
			hRaster2 	= findobj(h.Figure, 'Tag', 'Raster2');
			hPETH 		= findobj(h.Figure, 'Tag', 'PETH');

			propertiesToCopy = {'Name', 'GraphicsSmoothing', 'DefaultAxesFontSize'};

			set(hButtons, 'Visible', 'off');
			set(hTexts, 'Visible', 'off');

			switch lower(reformat)
				case lower('Raw')
					hFigurePrint = h.Figure;
				case lower('PETH')
					hFigureNew = figure('Position', [0, 0, 1776/4, 999/3]);
					for iProp = 1:length(propertiesToCopy)
						set(hFigureNew, propertiesToCopy{iProp}, get(h.Figure, propertiesToCopy{iProp}));
					end
					hPETHNew = copyobj(hPETH, hFigureNew);
					title(hPETHNew, h.Figure.Name, 'Interpreter', 'none');
					hPETHNew.OuterPosition = [0, 0, 1, 1];

					set(hPETHNew, 'XLim', get(hRaster, 'XLim'));

					if copyLegend
						legend(hPETHNew, 'Location', 'northwest');
					end

					if ~copyLabel
						xlabel(hPETHNew, '');
						ylabel(hPETHNew, '');
					end

					hFigurePrint = hFigureNew;
				case lower('Raster')
					hFigureNew = figure();
					for iProp = 1:length(propertiesToCopy)
						set(hFigureNew, propertiesToCopy{iProp}, get(h.Figure, propertiesToCopy{iProp}));
					end
					hRasterNew = copyobj(hRaster, hFigureNew);
					title(hRasterNew, h.Figure.Name, 'Interpreter', 'none');
					hRasterNew.OuterPosition = [0, 0, 1, 1];

					set(hRasterNew, 'XLim', get(hRaster, 'XLim'));

					if copyLegend
						legend(hRaster, 'Location', 'northwest');
					end

					if ~copyLabel
						xlabel(hRaster, '');
						ylabel(hRaster, '');
					end

					hFigurePrint = hFigureNew;
				case lower('RasterAndPETH')
					hFigureNew = figure('Position', [2 42 958 954]);
					for iProp = 1:length(propertiesToCopy)
						set(hFigureNew, propertiesToCopy{iProp}, get(h.Figure, propertiesToCopy{iProp}));
					end
					hRasterNew = copyobj(hRaster, hFigureNew);
					hRasterNew.OuterPosition = [0, 0.5, 1, 0.5];

					hPETHNew = copyobj(hPETH, hFigureNew);
					hPETHNew.OuterPosition = [0, 0, 1, 0.5];
					ylabel(hPETHNew, 'Spikes/s')

					set(hPETHNew, 'XLim', get(hRasterNew, 'XLim'));

					if copyLegend
						legend(hRaster, 'Location', 'northwest');
						legend(hPETHNew, 'Location', 'northwest');
					end

					% hTitle = suptitle(h.Figure.Name);
					hFigurePrint = hFigureNew;
				case lower('RasterAndPETHAndWaveform')
					hFigureNew = figure('Position', [2 42 958 954]);
					for iProp = 1:length(propertiesToCopy)
						set(hFigureNew, propertiesToCopy{iProp}, get(h.Figure, propertiesToCopy{iProp}));
					end
					hRasterNew = copyobj(hRaster, hFigureNew);
					hRasterNew.OuterPosition = [0, 0.5, 1, 0.5];

					hPETHNew = copyobj(hPETH, hFigureNew);
					hPETHNew.OuterPosition = [0, 0, 1, 0.5];
					ylabel(hPETHNew, 'Spikes/s')

					hWaveformNew = copyobj(hWaveform, hFigureNew);
					hWaveformNew.Position(1) = hRasterNew.Position(1) + 1/16*hRasterNew.Position(3);
					hWaveformNew.Position(3) = 0.25*hRasterNew.Position(3);
					hWaveformNew.Position(4) = 0.4*hRasterNew.Position(4);
					hWaveformNew.Position(2) = hRasterNew.Position(2) + hRasterNew.Position(4) - hWaveformNew.Position(3);
					hWaveformNew.Visible = 'off';

					set(hPETHNew, 'XLim', get(hRasterNew, 'XLim'));

					if copyLegend
						legend(hRasterNew, 'Location', 'north');
						legend(hPETHNew, 'Location', 'best');
					end

					% hTitle = suptitle(h.Figure.Name);
					hFigurePrint = hFigureNew;
				case lower('DualRasterAndWaveform')
					hFigureNew = figure('Position', [2 42 958 954]);
					for iProp = 1:length(propertiesToCopy)
						set(hFigureNew, propertiesToCopy{iProp}, get(h.Figure, propertiesToCopy{iProp}));
					end
					hRasterNew = copyobj(hRaster, hFigureNew);
					hRasterNew.OuterPosition = [0, hRaster2.Position(4)/(hRaster.Position(4) + hRaster2.Position(4)), 1, hRaster.Position(4)/(hRaster.Position(4) + hRaster2.Position(4))];
					xlabel(hRasterNew, '');
					title(hRasterNew, 'Lever Press')

					hRaster2New = copyobj(hRaster2, hFigureNew);
					hRaster2New.OuterPosition = [0, 0, 1, hRaster2.Position(4)/(hRaster.Position(4) + hRaster2.Position(4))];
					xlabel(hRaster2New, 'Time to movement (s)');
					title(hRaster2New, 'Lick')

					hWaveformNew = copyobj(hWaveform, hFigureNew);
					hWaveformNew.Position(1) = hRasterNew.Position(1) + 1/16*hRasterNew.Position(3);
					hWaveformNew.Position(3) = 0.25*hRasterNew.Position(3);
					hWaveformNew.Position(4) = 0.4*hRasterNew.Position(4);
					hWaveformNew.Position(2) = hRasterNew.Position(2) + hRasterNew.Position(4) - hWaveformNew.Position(3);
					hWaveformNew.Visible = 'off';

					if copyLegend
						legend(hRasterNew, 'Location', 'north');
						legend(hRaster2New, 'Location', 'north');
					end

					% hTitle = suptitle(h.Figure.Name);
					hFigurePrint = hFigureNew;

				case lower('RasterAndStimAndWaveform')
					hFigureNew = figure('Position', [2 42 958 954]);
					for iProp = 1:length(propertiesToCopy)
						set(hFigureNew, propertiesToCopy{iProp}, get(h.Figure, propertiesToCopy{iProp}));
					end
					hRasterNew = copyobj(hRaster, hFigureNew);
					hRasterNew.OuterPosition = [0, hRaster2.Position(4)/(hRaster.Position(4) + hRaster2.Position(4)), 1, hRaster.Position(4)/(hRaster.Position(4) + hRaster2.Position(4))];
					xlabel(hRasterNew, '');
					title(hRasterNew, 'Lever Press')

					hRaster2New = copyobj(hRaster2, hFigureNew);
					hRaster2New.OuterPosition = [0, 0, 1, hRaster2.Position(4)/(hRaster.Position(4) + hRaster2.Position(4))];
					xlabel(hRaster2New, 'Time to first stimOn (s)');
					title(hRaster2New, 'Stim')

					hWaveformNew = copyobj(hWaveform, hFigureNew);
					hWaveformNew.Position(1) = hRasterNew.Position(1) + 1/16*hRasterNew.Position(3);
					hWaveformNew.Position(3) = 0.25*hRasterNew.Position(3);
					hWaveformNew.Position(4) = 0.4*hRasterNew.Position(4);
					hWaveformNew.Position(2) = hRasterNew.Position(2) + hRasterNew.Position(4) - hWaveformNew.Position(3);
					hWaveformNew.Visible = 'off';

					if copyLegend
						legend(hRasterNew, 'Location', 'north');
						legend(hRaster2New, 'Location', 'north');
					end

					% hTitle = suptitle(h.Figure.Name);
					hFigurePrint = hFigureNew;

				case lower('PETHAndStim')
					hFigureNew = figure('Position', [2 42 958 954]);
					for iProp = 1:length(propertiesToCopy)
						set(hFigureNew, propertiesToCopy{iProp}, get(h.Figure, propertiesToCopy{iProp}));
					end
					hRasterNew = copyobj(hRaster2, hFigureNew);
					hRasterNew.OuterPosition = [0, 0.5, 1, 0.5];

					hPETHNew = copyobj(hPETH, hFigureNew);
					hPETHNew.OuterPosition = [0, 0, 1, 0.5];
					ylabel(hPETHNew, 'Spikes/s')

					hWaveformNew = copyobj(hWaveform, hFigureNew);
					hWaveformNew.Position(1) = hRasterNew.Position(1) + 1/16*hRasterNew.Position(3);
					hWaveformNew.Position(3) = 0.25*hRasterNew.Position(3);
					hWaveformNew.Position(4) = 0.4*hRasterNew.Position(4);
					hWaveformNew.Position(2) = hRasterNew.Position(2) + hRasterNew.Position(4) - hWaveformNew.Position(3);
					hWaveformNew.Visible = 'off';

					set(hPETHNew, 'XLim', get(hRasterNew, 'XLim'));

					if copyLegend
						legend(hRasterNew, 'Location', 'north');
						legend(hPETHNew, 'Location', 'best');
					end

					hTitle = suptitle(h.Figure.Name);
					hFigurePrint = hFigureNew;
			end

			if isempty(filename)
				print(hFigurePrint, '-clipboard', '-dbitmap')
			else
				print(hFigurePrint, filename, '-dpng')
            end
            try
    			set(hButtons, 'Visible', 'on');
        		set(hTexts, 'Visible', 'on');		
            end
			close(hFigurePrint)			
		end

		function GUIDeleteChannel(obj, hButton, evnt, iChannel, nextChn, p, h)
			obj.GUIBusy(h.Figure, true);
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
					close(h.Figure)
				else
					reference 			= p.Results.Reference;
					event 				= p.Results.Event;
					exclude 			= p.Results.Exclude;
					event2 				= p.Results.Event2;
					exclude2 			= p.Results.Exclude2;
					clusters 			= p.Results.Clusters;
					referenceCluster	= p.Results.ReferenceCluster;
					waveformWindow 		= p.Results.WaveformWindow;
					extendedWindow 		= p.Results.ExtendedWindow;
					minTrialLength 		= p.Results.MinTrialLength;
					bins 				= p.Results.Bins;
					binMethod 			= p.Results.BinMethod;
					spikeRateWindow 	= p.Results.SpikeRateWindow;
					rasterXLim 			= p.Results.RasterXLim;
					waveformYLim		= p.Results.WaveformYLim;
					fontSize 			= p.Results.FontSize;
					printMode 			= p.Results.PrintMode;
					frameRate 			= p.Results.FrameRate;
					plotStim 			= p.Results.PlotStim;

					obj.PlotChannel(nextChn, 'Reference', reference, 'Event', event, 'Exclude', exclude, 'Event2', event2, 'Exclude2', exclude2, 'Clusters', clusters, 'ReferenceCluster', referenceCluster, 'WaveformWindow', waveformWindow, 'ExtendedWindow', extendedWindow, 'MinTrialLength', minTrialLength, 'Bins', bins, 'BinMethod', binMethod, 'SpikeRateWindow', spikeRateWindow, 'RasterXLim', rasterXLim, 'WaveformYLim', waveformYLim, 'FontSize', fontSize, 'PrintMode', printMode, 'FrameRate', frameRate, 'Fig', h.Figure, 'PlotStim', plotStim);
				end
			end
			obj.GUIBusy(h.Figure, false);
		end

		function GUISelectSpikes(obj, hButton, evnt, iChannel, p, h)
			obj.GUIBusy(h.Figure, true);

			% Add onButtonDown EH for all axes
			set(h.Waveform, 'ButtonDownFcn', {@obj.GUISelectSpikesOnAxesSelected, iChannel, p, h, 'Waveform', false});
			set(h.PCA, 'ButtonDownFcn', {@obj.GUISelectSpikesOnAxesSelected, iChannel, p, h, 'Feature', false});
			set(h.Raster, 'ButtonDownFcn', {@obj.GUISelectSpikesOnAxesSelected, iChannel, p, h, 'Raster', false});
			set(h.Raster2, 'ButtonDownFcn', {@obj.GUISelectSpikesOnAxesSelected, iChannel, p, h, 'Raster', false});
			set([h.Waveform, h.PCA, h.Raster, h.Raster2], 'Box', 'on', 'LineWidth', 2, 'XColor', 'r', 'YColor', 'r');
		end

		function GUIExcludeSpikes(obj, hButton, evnt, iChannel, p, h)
			obj.GUIBusy(h.Figure, true);

			% Add onButtonDown EH for all axes
			set(h.Waveform, 'ButtonDownFcn', {@obj.GUISelectSpikesOnAxesSelected, iChannel, p, h, 'Waveform', true});
			set(h.PCA, 'ButtonDownFcn', {@obj.GUISelectSpikesOnAxesSelected, iChannel, p, h, 'Feature', true});
			set(h.Raster, 'ButtonDownFcn', {@obj.GUISelectSpikesOnAxesSelected, iChannel, p, h, 'Raster', true});
			set(h.Raster2, 'ButtonDownFcn', {@obj.GUISelectSpikesOnAxesSelected, iChannel, p, h, 'Raster', true});
			set([h.Waveform, h.PCA, h.Raster, h.Raster2], 'Box', 'on', 'LineWidth', 2, 'XColor', 'r', 'YColor', 'r');
		end

		function GUISelectAllSpikes(obj, hButton, evnt, iChannel, p, h)
			obj.GUIBusy(h.Figure, true);
			h.Figure.UserData.SelectedSampleIndex = [];
			obj.ReplotChannel(iChannel, p, h);
			obj.GUIBusy(h.Figure, false);
		end

		function GUISelectSpikesOnAxesSelected(obj, hAxes, evnt, iChannel, p, h, selectionMode, exclude)
			% Remove onButtonDown EH for all axes
			set([h.Waveform, h.PCA, h.Raster, h.Raster2], 'ButtonDownFcn', '');

			switch lower(selectionMode)
				case 'waveform'
					set([h.PCA, h.Raster], 'Box', 'off', 'LineWidth', 0.5, 'XColor', 'k', 'YColor', 'k');
					switch evnt.Button
						case 1 % LMB - hoop mode
							[x, y] = getline(hAxes);
							hoop = [x, y];
							if size(hoop, 1) == 2
								selection = obj.GetSpikesByHoop(iChannel, hoop, 'Clusters', h.Figure.UserData.SelectedClusters, 'SampleIndex', h.Figure.UserData.SelectedSampleIndex);
							else
								selection = [];
							end
						case 3 % RMB - threshold mode
							axes(hAxes);
							[~, threshold] = ginput(1);
							if ~isempty(threshold)
								selection = obj.GetSpikesByThreshold(iChannel, threshold, 'Clusters', h.Figure.UserData.SelectedClusters, 'SampleIndex', h.Figure.UserData.SelectedSampleIndex);
							else
								selection = [];
							end
					end
				case 'feature'
					set([h.Waveform, h.Raster], 'Box', 'off', 'LineWidth', 0.5, 'XColor', 'k', 'YColor', 'k');
					[x, y] = getline(hAxes, 'closed');
					polygon = [x, y];
					if size(polygon, 1) >= 3
						selection = obj.GetSpikesByFeature(iChannel, polygon, 'Clusters', h.Figure.UserData.SelectedClusters, 'SampleIndex', h.Figure.UserData.SelectedSampleIndex);
					else
						selection = [];
					end
				case 'raster'
					set([h.Waveform, h.PCA], 'Box', 'off', 'LineWidth', 0.5, 'XColor', 'k', 'YColor', 'k');
					rect = getrect(hAxes);
					if nnz(rect(3:4)) == 2
						pp = hAxes.UserData.PlotParams;
						switch lower(pp.AlignTo)
							case 'reference'
								windowReference = 'Start';
							case 'event'
								windowReference = 'End';
						end
						selection = obj.GetSpikesByTrial(iChannel, 'Reference', pp.Reference, 'Event', pp.Event, 'Exclude', pp.Exclude,...
							'Clusters', h.Figure.UserData.SelectedClusters, 'SampleIndex', h.Figure.UserData.SelectedSampleIndex, 'Window', [rect(1), rect(1) + rect(3)], 'WindowReference', windowReference);
					else
						selection = [];
					end
				otherwise
					set([h.Waveform, h.PCA, h.Raster, h.Raster2], 'Box', 'off', 'LineWidth', 0.5, 'XColor', 'k', 'YColor', 'k');
					obj.GUIBusy(h.Figure, false);
					return
			end

			if isempty(h.Figure.UserData.SelectedSampleIndex)
				originalSelection = obj.Spikes(iChannel).SampleIndex;
			else
				originalSelection = h.Figure.UserData.SelectedSampleIndex;
			end

			if ~isempty(selection)
				if exclude
					h.Figure.UserData.SelectedSampleIndex = setdiff(originalSelection, selection);
				else
					h.Figure.UserData.SelectedSampleIndex = intersect(originalSelection, selection);
				end
			end

			set([h.Waveform, h.PCA, h.Raster, h.Raster2], 'Box', 'off', 'LineWidth', 0.5, 'XColor', 'k', 'YColor', 'k');

			obj.ReplotChannel(iChannel, p, h);

			obj.GUIBusy(h.Figure, false);
		end

		function GUIDeleteSpikes(obj, hButton, evnt, iChannel, p, h)
			obj.GUIBusy(h.Figure, true);

			selectedClusters 	= h.Figure.UserData.SelectedClusters;
			selectedSampleIndex = h.Figure.UserData.SelectedSampleIndex;
			classes 			= obj.Spikes(iChannel).Cluster.Classes;
			sampleIndex 		= obj.Spikes(iChannel).SampleIndex;

			if isempty(selectedClusters)
				numSeletecdClusters = length(unique(classes));
				if isempty(selectedSampleIndex)
					numSelectedSpikes = length(obj.Spikes(iChannel).SampleIndex);
				else
					numSelectedSpikes = length(selectedSampleIndex);
				end
			else
				numSeletecdClusters = length(selectedClusters);
				if isempty(selectedSampleIndex)
					numSelectedSpikes = nnz(ismember(classes, selectedClusters));
				else
					sampleIndex = sampleIndex(ismember(classes, selectedClusters));
					numSelectedSpikes = nnz(ismember(sampleIndex, selectedSampleIndex));
				end
			end

			answer = questdlg(...
				['Permanently delete ' , num2str(numSelectedSpikes), ' spikes from ', num2str(numSeletecdClusters), ' clusters?'],...
				'Delete spikes',...
				'Delete', 'Cancel',...
				'Cancel');
			if strcmpi(answer, 'Delete')
				obj.DeleteWaveforms(iChannel, selectedSampleIndex, 'IndexType', 'SampleIndex', 'Clusters', selectedClusters);
				h.Figure.UserData.SelectedSampleIndex = [];
				h.Figure.UserData.SelectedCluster = [];
				obj.ReplotChannel(iChannel, p, h);
			end

			obj.GUIBusy(h.Figure, false);
		end

		function GUISpikesToCluster(obj, hButton, evnt, iChannel, p, h)
			obj.GUIBusy(h.Figure, true);
			liststr = {'New'};
			liststr = [liststr, cellfun(@num2str, num2cell(unique(obj.Spikes(iChannel).Cluster.Classes)), 'UniformOutput', false)];

			selectedClusters 	= h.Figure.UserData.SelectedClusters;
			selectedSampleIndex = h.Figure.UserData.SelectedSampleIndex;
			classes 			= obj.Spikes(iChannel).Cluster.Classes;
			sampleIndex 		= obj.Spikes(iChannel).SampleIndex;

			if isempty(selectedClusters)
				numSeletecdClusters = length(unique(classes));
				if isempty(selectedSampleIndex)
					selected = 1:length(obj.Spikes(iChannel).SampleIndex);
				else
					selected = find(ismember(sampleIndex, selectedSampleIndex));
				end
			else
				numSeletecdClusters = length(selectedClusters);
				if isempty(selectedSampleIndex)
					selected = find(ismember(classes, selectedClusters));
				else
					sampleIndexInCluster = sampleIndex(ismember(classes, selectedClusters));
					sampleIndexInBoth = intersect(sampleIndexInCluster, selectedSampleIndex);
					selected = find(ismember(sampleIndex, sampleIndexInBoth));
				end
			end

			numSelectedSpikes = length(selected);
			[selection, ok] = listdlg(...
				'PromptString', ['Move ' , num2str(numSelectedSpikes), ' spikes from ', num2str(numSeletecdClusters), ' clusters to selected cluster:'],...
				'SelectionMode', 'single',...
				'OKString', 'Move',...
				'ListString', liststr,...
				'InitialValue', [],...
				'ListSize', [250, 150]);

			if ok
				if strcmpi(liststr{selection}, 'New')
					newClass = 1 + max(obj.Spikes(iChannel).Cluster.Classes);
				else
					newClass = str2num(liststr{selection});
				end
				obj.Spikes(iChannel).Cluster.Classes(selected) = newClass;

				% Remove empty clusters and consolidate cluster number
				if length(unique(obj.Spikes(iChannel).Cluster.Classes)) < max(obj.Spikes(iChannel).Cluster.Classes)
					obj.ClusterMerge(iChannel, num2cell(unique(obj.Spikes(iChannel).Cluster.Classes)));
					h.Figure.UserData.SelectedClusters = [];
				else
					h.Figure.UserData.SelectedClusters = unique([h.Figure.UserData.SelectedClusters, newClass]);
				end
				h.Figure.UserData.SelectedSampleIndex = [];

				obj.ReplotChannel(iChannel, p, h);
			end
			obj.GUIBusy(h.Figure, false);
		end

		% Sort spikes and digital events into trial structure
		function varargout = GetSpikesByTrial(obj, channel, varargin)
			p = inputParser;
			addRequired(p, 'Channel', @isnumeric);
			addParameter(p, 'Reference', 'CueOn', @(x) isnumeric(x) || ischar(x)); % Timestamps for trial start events, or event name, 'CueOn'
			addParameter(p, 'Event', 'PressOn', @(x) isnumeric(x) || ischar(x)); % Timestamps for trial end events, or event name, 'PressOn'
			addParameter(p, 'Exclude', 'Lick', @ischar); % Event name, 'LickOn'
			addParameter(p, 'Clusters', [], @isnumeric);
			addParameter(p, 'SampleIndex', [], @isnumeric);
			addParameter(p, 'Window', [], @isnumeric);
			addParameter(p, 'WindowReference', 'StartAndEnd', @ischar); % Start: window + trialStart, End: window + trialEnd, StartAndEnd: window + [trialStart, trialEnd]
			parse(p, channel, varargin{:});
			channel 			= p.Results.Channel;
			reference 			= p.Results.Reference;
			event 				= p.Results.Event;
			exclude 			= p.Results.Exclude;
			clusters 			= p.Results.Clusters;
			selectedSampleIndex	= p.Results.SampleIndex;
			window				= p.Results.Window;
			windowReference 	= p.Results.WindowReference;

			if isempty(window)
				window = [0, 0];
			end

			if ischar(reference) || ischar(event)
				reference 	= obj.DigitalEvents.(reference);
				event 		= obj.DigitalEvents.(event);
				if isempty(exclude)
					exclude = [];
				else
					exclude	= obj.DigitalEvents.(exclude);
				end
				[trialStart, trialEnd] = TetrodeRecording.FindFirstInTrial(reference, event, exclude);
			else
				trialStart = reference;
				trialEnd = event;
			end

			trialStart 	= reshape(trialStart, 1, []);
			trialEnd 	= reshape(trialEnd, 1, []);

			switch lower(windowReference)
				case 'start'
					edges = [trialStart + window(1); trialStart + window(2)];
				case 'end'
					edges = [trialEnd + window(1); trialEnd + window(2)];
				case 'startandend'
					edges = [trialStart + window(1); trialEnd(1:length(trialStart)) + window(2)];
				case 'endandstart'
					edges = [trialEnd + window(1); trialStart + window(2)];
			end
			edges = edges(:);

			if isempty(edges)
				varargout = {[], [], []};
			else
				sampleIndex = obj.Spikes(channel).SampleIndex;
				timestamps 	= obj.Spikes(channel).Timestamps;
				classes 	= obj.Spikes(channel).Cluster.Classes;

				if ~isempty(clusters)
					selected 	= ismember(classes, clusters);
					sampleIndex = sampleIndex(selected);
					timestamps 	= timestamps(selected);
				end

				if ~isempty(selectedSampleIndex)
					selected 	= ismember(sampleIndex, selectedSampleIndex);
					sampleIndex = sampleIndex(selected);
					timestamps 	= timestamps(selected);
				end

				[~, ~, bins] = histcounts(timestamps, edges);
				oddBins = rem(bins, 2) ~= 0;	% Spikes in odd bins occur between reference and event, should keep these spikes

				timestamps 	= timestamps(oddBins);
				sampleIndex = sampleIndex(oddBins);
				trials 		= (bins(oddBins) + 1)/2;

				varargout = {sampleIndex, timestamps, trials};
			end
		end

		function varargout = GetSpikesByFeature(obj, channel, polygon, varargin)
			p = inputParser;
			addRequired(p, 'Channel', @isnumeric);
			addRequired(p, 'Polygon', @isnumeric);
			addParameter(p, 'SampleIndex', [], @isnumeric);
			addParameter(p, 'Clusters', [], @isnumeric);
			parse(p, channel, polygon, varargin{:});
			channel 			= p.Results.Channel;
			polygon 			= p.Results.Polygon;
			selectedSampleIndex	= p.Results.SampleIndex;
			clusters 			= p.Results.Clusters;

			sampleIndex = obj.Spikes(channel).SampleIndex;
			timestamps 	= obj.Spikes(channel).Timestamps;
			waveforms 	= obj.Spikes(channel).Waveforms;
			classes 	= obj.Spikes(channel).Cluster.Classes;
			feature 	= obj.Spikes(channel).Feature.Coeff(:, 1:2);

			if ~isempty(clusters)
				selected 	= ismember(classes, clusters);
				sampleIndex = sampleIndex(selected);
				timestamps 	= timestamps(selected);
				feature 	= feature(selected, :);
			end

			if ~isempty(selectedSampleIndex)
				selected 	= ismember(sampleIndex, selectedSampleIndex);
				sampleIndex = sampleIndex(selected);
				timestamps 	= timestamps(selected);
				feature 	= feature(selected, :);
			end

			if ~isempty(polygon)
				selected 	= inpolygon(feature(:, 1), feature(:, 2), polygon(:, 1), polygon(:, 2));
				sampleIndex = sampleIndex(selected);
				timestamps 	= timestamps(selected);
			end

			varargout = {sampleIndex, timestamps};
		end

		function varargout = GetSpikesByThreshold(obj, channel, threshold, varargin)
			p = inputParser;
			addRequired(p, 'Channel', @isnumeric);
			addRequired(p, 'Threshold', @isnumeric);
			addParameter(p, 'SampleIndex', [], @isnumeric);
			addParameter(p, 'Clusters', [], @isnumeric);
			parse(p, channel, threshold, varargin{:});
			channel 			= p.Results.Channel;
			threshold 			= p.Results.Threshold;
			selectedSampleIndex	= p.Results.SampleIndex;
			clusters 			= p.Results.Clusters;

			t 			= obj.Spikes(channel).WaveformTimestamps;
			sampleIndex = obj.Spikes(channel).SampleIndex;
			timestamps 	= obj.Spikes(channel).Timestamps;
			waveforms 	= obj.Spikes(channel).Waveforms;
			classes 	= obj.Spikes(channel).Cluster.Classes;

			% If blackrock, convert data to 'uV'
			if strcmpi(obj.System, 'Blackrock')
				waveforms = waveforms/4;
			end

			if ~isempty(clusters)
				selected 	= ismember(classes, clusters);
				sampleIndex = sampleIndex(selected);
				timestamps 	= timestamps(selected);
				waveforms 	= waveforms(selected, :);
			end

			if ~isempty(selectedSampleIndex)
				selected 	= ismember(sampleIndex, selectedSampleIndex);
				sampleIndex = sampleIndex(selected);
				timestamps 	= timestamps(selected);
				waveforms 	= waveforms(selected, :);
			end

			if ~isempty(threshold) && threshold ~= 0
				thresholded = max(sign(threshold)*waveforms, [], 2) >= abs(threshold);
				sampleIndex = sampleIndex(thresholded);
				timestamps 	= timestamps(thresholded);
			end

			varargout = {sampleIndex, timestamps};
		end

		function varargout = GetSpikesByWaveform(obj, channel, polygon, varargin)
			p = inputParser;
			addRequired(p, 'Channel', @isnumeric);
			addRequired(p, 'Polygon', @isnumeric);
			addParameter(p, 'SampleIndex', [], @isnumeric);
			addParameter(p, 'Clusters', [], @isnumeric);
			parse(p, channel, polygon, varargin{:});
			channel 			= p.Results.Channel;
			polygon 			= p.Results.Polygon;
			selectedSampleIndex	= p.Results.SampleIndex;
			clusters 			= p.Results.Clusters;

			t 			= obj.Spikes(channel).WaveformTimestamps;
			sampleIndex = obj.Spikes(channel).SampleIndex;
			timestamps 	= obj.Spikes(channel).Timestamps;
			waveforms 	= obj.Spikes(channel).Waveforms;
			classes 	= obj.Spikes(channel).Cluster.Classes;

			% If blackrock, convert data to 'uV'
			if strcmpi(obj.System, 'Blackrock')
				waveforms = waveforms/4;
			end

			if ~isempty(clusters)
				selected 	= ismember(classes, clusters);
				sampleIndex = sampleIndex(selected);
				timestamps 	= timestamps(selected);
				waveforms 	= waveforms(selected, :);
			end

			if ~isempty(selectedSampleIndex)
				selected 	= ismember(sampleIndex, selectedSampleIndex);
				sampleIndex = sampleIndex(selected);
				timestamps 	= timestamps(selected);
				waveforms 	= waveforms(selected, :);
			end

			intersectsWithPolygon = true(length(sampleIndex), 1);
			if ~isempty(polygon)
				upSampleRate = 10;
				% Upsample waveformtimestamps
				tUp = interp1(1:length(t), t, 1:(1/upSampleRate):length(t));
				for iWaveform = 1:length(sampleIndex)
					thisWaveform = waveforms(iWaveform, :);
					thisWaveformUp = interp1(1:length(t), thisWaveform, 1:(1/upSampleRate):length(t)); % Upsample waveform
					intersectsWithPolygon(iWaveform) = nnz(inpolygon(tUp, thisWaveformUp, polygon(:, 1), polygon(:, 2))) > 0;
				end
			end

			sampleIndex = sampleIndex(intersectsWithPolygon);
			timestamps 	= timestamps(intersectsWithPolygon);

			varargout = {sampleIndex, timestamps};
		end

		function varargout = GetSpikesByHoop(obj, channel, hoop, varargin)
			p = inputParser;
			addRequired(p, 'Channel', @isnumeric);
			addRequired(p, 'Hoop', @isnumeric);
			addParameter(p, 'SampleIndex', [], @isnumeric);
			addParameter(p, 'Clusters', [], @isnumeric);
			parse(p, channel, hoop, varargin{:});
			channel 			= p.Results.Channel;
			hoop 				= p.Results.Hoop;
			selectedSampleIndex	= p.Results.SampleIndex;
			clusters 			= p.Results.Clusters;

			t 			= obj.Spikes(channel).WaveformTimestamps;
			sampleIndex = obj.Spikes(channel).SampleIndex;
			timestamps 	= obj.Spikes(channel).Timestamps;
			waveforms 	= obj.Spikes(channel).Waveforms;
			classes 	= obj.Spikes(channel).Cluster.Classes;

			% If blackrock, convert data to 'uV'
			if strcmpi(obj.System, 'Blackrock')
				waveforms = waveforms/4;
			end

			if ~isempty(clusters)
				selected 	= ismember(classes, clusters);
				sampleIndex = sampleIndex(selected);
				timestamps 	= timestamps(selected);
				waveforms 	= waveforms(selected, :);
			end

			y = hoop(1:2, 2);
			if y(1)*y(2) > 0
				threshold = sign(y(1))*min(abs(y));
				selectedSampleIndex = obj.GetSpikesByThreshold(channel, threshold, 'SampleIndex', selectedSampleIndex, 'Clusters', clusters);
			end

			if ~isempty(selectedSampleIndex)
				selected 	= ismember(sampleIndex, selectedSampleIndex);
				sampleIndex = sampleIndex(selected);
				timestamps 	= timestamps(selected);
				waveforms 	= waveforms(selected, :);
			end

			if length(hoop) == 2
				intersectsWithHoop = true(length(sampleIndex), 1);
				if ~isempty(hoop)
					for iWaveform = 1:length(sampleIndex)
						intersectsWithHoop(iWaveform) = ~isempty(polyxpoly(t, waveforms(iWaveform, :), hoop(1:2, 1), hoop(1:2, 2)));
					end
				end

				sampleIndex = sampleIndex(intersectsWithHoop);
				timestamps 	= timestamps(intersectsWithHoop);
			else
				warning('Hoop can only contain two points')
			end

			varargout = {sampleIndex, timestamps};
		end
	end

	methods (Static)
		function previewObj = BatchPreview(showResults)
			previewObj = TetrodeRecording();
			dirs = uipickfiles('Prompt', 'Select (multiple) folders...');
			dirs = dirs(isfolder(dirs));
			for iDir = 1:length(dirs)
				files = dir([dirs{iDir}, '\*.rhd']);
				if ~isempty(files)
					stepSize = round(length(files)/6);
					files = {files(stepSize:stepSize:5*stepSize).name};
					sysName = 'Intan';
				else
					files = dir([dirs{iDir}, '\*.nev']);
					if ~isempty(files)
						filename = strsplit(dirs{iDir}, '\');
						filename = filename{end};
						filenameNEV = [filename, '.nev'];
						filenameNSx = [filename, '.ns5'];
						files = {filenameNEV, filenameNSx};
						sysName = 'Blackrock';
					else
						warning(['No blackrock/intan files found in folder (', dirs{iDir}, ').'])
						continue
					end
				end
				previewObj(iDir) = TetrodeRecording();
				previewObj(iDir).System = sysName;
				previewObj(iDir).Path = [dirs{iDir}, '\'];
				previewObj(iDir).Files = files;
                rig = TetrodeRecording.GetRig(previewObj(iDir).Path);
				try
					previewObj(iDir).Preview('Rig', rig, 'HideResults', true);
				catch ME
					warning(['Error when processing folder (', dirs{iDir}, ') - this one will be skipped.'])
					warning(sprintf('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message))
				end				
			end

			if nargin < 1
				showResults = true;
			end
			if showResults
				for iDir = 1:length(dirs)
					try
						previewObj(iDir).PlotAllChannels('YLim', 'auto');
						previewObj(iDir).ClearCache();
					catch ME
						warning(['Error when processing folder (', dirs{iDir}, ') - this one will be skipped.'])
						warning(sprintf('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message))
					end
				end
			end
			% TetrodeRecording.RandomWords();
		end

		function BatchProcess(previewObj, varargin)
			p = inputParser;
			addRequired(p, 'Obj', @(x) isa(x, 'TetrodeRecording'));
			addParameter(p, 'ChunkSize', 10, @isnumeric);
			addParameter(p, 'NumSigmas', 4, @isnumeric);
			addParameter(p, 'NumSigmasReturn', [], @isnumeric);
			addParameter(p, 'NumSigmasReject', [], @isnumeric);
			addParameter(p, 'WaveformWindow', [-5, 5], @isnumeric);
			addParameter(p, 'WaveformFeatureWindow', [-0.5, 0.5], @isnumeric);
			addParameter(p, 'FeatureMethod', 'WaveletTransform', @ischar);
			addParameter(p, 'ClusterMethod', 'kmeans', @ischar);
			addParameter(p, 'Dimension', 10, @isnumeric);
			addParameter(p, 'Prefix', 'tr_', @ischar);
			addParameter(p, 'SavePath', '', @ischar);
			parse(p, previewObj, varargin{:});
			previewObj 		= p.Results.Obj;
			chunkSize 		= p.Results.ChunkSize;
			numSigmas 		= p.Results.NumSigmas;
			numSigmasReturn = p.Results.NumSigmasReturn;
			numSigmasReject = p.Results.NumSigmasReject;
			waveformWindow 	= p.Results.WaveformWindow;
			waveformFeatureWindow 	= p.Results.WaveformFeatureWindow;
			featureMethod 	= p.Results.FeatureMethod;
			clusterMethod 	= p.Results.ClusterMethod;
			dimension 		= p.Results.Dimension;
			prefix 			= p.Results.Prefix;
			savePath 		= p.Results.SavePath;

			selectedChannels = {previewObj.SelectedChannels};
			allPaths = {previewObj.Path};
			for iDir = 1:length(selectedChannels)
				rig = TetrodeRecording.GetRig(allPaths{iDir});
				if rig == 1
                    if isfield(previewObj(iDir).ChannelMap, 'Rig1')
                        channelsOnRig = previewObj(iDir).ChannelMap.Rig1;
                    else
                        channelsOnRig = NaN;
                    end
				elseif rig == 2
                    if isfield(previewObj(iDir).ChannelMap, 'Rig2')
                        channelsOnRig = previewObj(iDir).ChannelMap.Rig2;
                    else
                        channelsOnRig = NaN;
                    end
				end
				channels = selectedChannels{iDir};
				if ~isempty(channels)
					try
						TetrodeRecording.TTS(['Processing folder ', num2str(iDir), '/', num2str(length(selectedChannels)), ':\n']);
						TetrodeRecording.ProcessFolder(allPaths{iDir}, savePath, chunkSize, channels, channelsOnRig, numSigmas, numSigmasReturn, numSigmasReject, waveformWindow, waveformFeatureWindow, featureMethod, clusterMethod, dimension, prefix, rig);
					catch ME
						warning(['Error when processing folder (', allPaths{iDir}, ') - this one will be skipped.'])
						warning(sprintf('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message))
					end
				end
			end
			TetrodeRecording.RandomWords();
		end

		function ProcessFolder(thisPath, savePath, chunkSize, channels, channelsOnRig, numSigmas, numSigmasReturn, numSigmasReject, waveformWindow, waveformFeatureWindow, featureMethod, clusterMethod, dimension, prefix, rig)
			tr = TetrodeRecording();
			tr.Path = thisPath;
			files = dir([tr.Path, '*.rhd']);
			if ~isempty(files)
				tr.System = 'Intan';
				tr.Files = {files.name};
			else
				files = dir([tr.Path, '*.nev']);
				if ~isempty(files)
					filename = strsplit(tr.Path, '\');
					filename = filename{end - 1};
					filenameNEV = [filename, '.nev'];
					filenameNSx = [filename, '.ns5'];
					tr.Files = {filenameNEV, filenameNSx};
					tr.System = 'Blackrock';
				else
					warning(['No blackrock/intan files found in folder (', tr.Path, ').'])
					return
				end
			end

			tr.ReadFiles(chunkSize, 'Rig', rig, 'Channels', channels, 'ChannelsOnRig', channelsOnRig, 'NumSigmas', numSigmas, 'NumSigmasReturn', numSigmasReturn, 'NumSigmasReject', numSigmasReject, 'WaveformWindow', waveformWindow);
			tr.SpikeSort([], 'FeatureMethod', featureMethod, 'ClusterMethod', clusterMethod, 'Dimension', dimension);
			TetrodeRecording.BatchSave(tr, 'Prefix', prefix, 'DiscardData', false, 'MaxChannels', 5, 'SavePath', savePath);
		end

		function varargout = BatchLoad(expNames)
			iExp = [];
			if nargin < 1
				files = uipickfiles('Prompt', 'Select .mat files containing TetrodeRecording objects to load...', 'Type', {'*.mat', 'MAT-files'});
				for iFile = 1:length(files)
					disp(['Loading file: ', files{iFile}, '...'])
					S(iFile) = load(files{iFile}, 'tr');
				end
			else
				files = {};
				for iFile = 1:length(expNames)
					if (isfile(expNames{iFile}))
						files{length(files) + 1} = expNames{iFile};
					else
						thisFile = dir(['C:\SERVER\**\SpikeSort\tr*', expNames{iFile}, '*.mat']);
						thisFile = thisFile(~[thisFile.isdir]);

						if length(thisFile) > 0
							iExp = [iExp, iFile];
							[~, iLongest] = max(cellfun(@length, {thisFile.name}));
							iLongest = find(cellfun(@(x) length(x) == length(thisFile(iLongest).name), {thisFile.name}));
							for i = transpose(iLongest(:))
								files{length(files) + 1} = [thisFile(i).folder, '\', thisFile(i).name];
							end
						end
					end
				end

				for iFile = 1:length(files)
					% Load files
					disp(['Loading file: ', files{iFile}, '...'])
					S(iFile) = load(files{iFile}, 'tr');
				end
			end
			% Merge multi-part files
			for iFile = 1:length(files)
				part = S(iFile).tr.Part;
				if nnz(part > 1) > 1
					partOne = cellfun(@(tr) strcmpi(tr.Path, S(iFile).tr.Path) && ~isempty(tr.Part) && tr.Part(1) == 1, {S.tr});
					channels = [S(iFile).tr.Spikes.Channel];
					S(partOne).tr.Spikes(channels) = S(iFile).tr.Spikes(channels);
					S(partOne).tr.Part(2) = S(partOne).tr.Part(2) - 1;
					S(iFile).tr.Part = [];
				end
			end

			partOne = cellfun(@(tr) ~isempty(tr.Part), {S.tr});

			tr = [S(partOne).tr];

			varargout = {tr, iExp};
		end

		function BatchSave(TR, varargin)
			p = inputParser;
			addParameter(p, 'Prefix', '', @ischar);
			addParameter(p, 'DiscardData', false, @islogical);
			addParameter(p, 'MaxChannels', [], @isnumeric);
			addParameter(p, 'SavePath', '', @ischar);
			parse(p, varargin{:});
			prefix = p.Results.Prefix;
			discardData = p.Results.DiscardData;
			maxChannels = p.Results.MaxChannels;
			savePath = p.Results.SavePath;

			for iTr = 1:length(TR)
				tr = TR(iTr);
				if isempty(tr.Path)
					continue
				end

				if isempty(savePath)
					thisSavePath = tr.Path;
					if ~isfolder([thisSavePath, '..\SpikeSort\'])
						mkdir([thisSavePath, '..\SpikeSort\'])
					end
					thisSavePath = [thisSavePath, '..\SpikeSort\'];
				else
					thisSavePath = savePath;
				end

				expName = tr.GetExpName();
				if discardData
					tr.Spikes = [];
					tr.DigitalEvents = [];
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
							file = [thisSavePath, prefix, expName, '(', num2str(iPart), ')', '.mat'];
							save(file, 'tr', '-v7.3');
						end
					catch ME
						tr.Spikes = spikes;
						rethrow(ME)
					end
					tr.Spikes = spikes;
				else
					file = [thisSavePath, prefix, expName, '.mat'];
					save(file, 'tr', '-v7.3');
				end
			end
		end

		%% BatchPlotSimple: function description
		function BatchPlotSimple(TR, list, varargin)
			p = inputParser;
			addParameter(p, 'RasterXLim', [-5, 0], @isnumeric);
			addParameter(p, 'RasterXLimStim', [-0.5, 0.5], @isnumeric);
			parse(p, varargin{:});

			for iPlot = 1:size(list, 1)
				thisAnimal = list{iPlot, 1};
				thisDate = list{iPlot, 2};
				thisChannel = list{iPlot, 3};
				thisCluster = list{iPlot, 4};
				for iTr = 1:length(TR)
					fprintf(1, 'Plotting: iTr = %d, %s_%d in %s\n', iTr, thisAnimal, thisDate, TR(iTr).Path)
					if ~isempty(strfind(TR(iTr).Path, thisAnimal)) && ~isempty(strfind(TR(iTr).Path, num2str(thisDate)))
						disp(iTr)
						[hFigure, ~, ~, hTitle] = TR(iTr).PlotUnitSimple(thisChannel, thisCluster,...
							'Reference', 'CueOn', 'Event', 'PressOn', 'Exclude', 'LickOn',...
							'RasterXLim', p.Results.RasterXLim, 'RasterXLimStim', p.Results.RasterXLimStim,...
							'PlotType', 'Raster', 'Position', [0, 0, 0.5, 1]);
						[hFigure2] = TR(iTr).PlotUnitSimple(thisChannel, thisCluster,...
							'Reference', 'CueOn', 'Event', 'PressOn', 'Exclude', 'LickOn',...
							'RasterXLim', p.Results.RasterXLim, 'RasterXLimStim', p.Results.RasterXLimStim,...
							'PlotType', 'PETH', 'Position', [0.5, 0, 0.5, 1]);

						disp(['Saving file:', hTitle.String, ' Raster Simple'])
						print(hFigure, [hTitle.String, ' Raster Simple'], '-dpng')
						print(hFigure2, [hTitle.String, ' PETH Simple'], '-dpng')

						% input('Type anything to continue...\n');
						try
							close(hFigure)
							close(hFigure2)
						end
						break
					end
				end
			end
		end
		

		function BatchPlot(TR, list, varargin)
			p = inputParser;
			addParameter(p, 'Reformat', 'RasterAndPETHAndWaveform', @ischar);
			addParameter(p, 'WaveformYLim', [-300, 300], @(x) isnumeric(x) || ischar(x));
			addParameter(p, 'RasterXLim', [-5, 0], @isnumeric);
			addParameter(p, 'ExtendedWindow', [-1, 0], @isnumeric);
			addParameter(p, 'CopyLegend', true, @islogical);
			addParameter(p, 'CopyLabel', true, @islogical);
			addParameter(p, 'PlotStim', false, @islogical);
			addParameter(p, 'PlotLick', false, @islogical);
			parse(p, varargin{:});
			reformat 		= p.Results.Reformat;
			waveformYLim 	= p.Results.WaveformYLim;
			rasterXLim 		= p.Results.RasterXLim;
			extendedWindow 	= p.Results.ExtendedWindow;
			copyLegend 		= p.Results.CopyLegend;
			copyLabel 		= p.Results.CopyLabel;
			plotStim 		= p.Results.PlotStim;
			plotLick 		= p.Results.PlotLick;
			for iPlot = 1:size(list, 1)
				thisAnimal = list{iPlot, 1};
				thisDate = list{iPlot, 2};
				thisChannel = list{iPlot, 3};
				thisCluster = list{iPlot, 4};

				for iTr = 1:length(TR)
					if ~isempty(strfind(TR(iTr).Path, thisAnimal)) && ~isempty(strfind(TR(iTr).Path, num2str(thisDate)))
						thisRefCluster = max(TR(iTr).Spikes(thisChannel).Cluster.Classes); % Use the last cluster is 'noise'/reference cluster
						if plotStim
							hFigure = TR(iTr).PlotChannel(thisChannel, 'PrintMode', true, 'Clusters', thisCluster, 'ReferenceCluster', thisRefCluster, 'Reference', 'CueOn', 'Event', 'PressOn', 'Exclude', 'LickOn', 'Event2', '', 'Exclude2', '', 'RasterXLim', rasterXLim, 'ExtendedWindow', extendedWindow, 'WaveformYLim', waveformYLim, 'PlotStim', true);
						elseif plotLick
							hFigure = TR(iTr).PlotChannel(thisChannel, 'PrintMode', true, 'Clusters', thisCluster, 'ReferenceCluster', thisRefCluster, 'Reference', 'CueOn', 'Event', 'PressOn', 'Exclude', 'LickOn', 'Event2', 'LickOn', 'Exclude2', 'PressOn', 'WaveformYLim', waveformYLim, 'RasterXLim', rasterXLim, 'ExtendedWindow', extendedWindow);
						else
							hFigure = TR(iTr).PlotChannel(thisChannel, 'PrintMode', true, 'Clusters', thisCluster, 'ReferenceCluster', thisRefCluster, 'Reference', 'CueOn', 'Event', 'PressOn', 'Exclude', 'LickOn', 'Event2', '', 'Exclude2', '', 'WaveformYLim', waveformYLim, 'RasterXLim', rasterXLim, 'ExtendedWindow', extendedWindow, 'PlotStim', false);
						end
						TR(iTr).GUISavePlot([], [], hFigure, 'Reformat', reformat, 'CopyLegend', copyLegend, 'CopyLabel', copyLabel, 'Filename', char(sprintf("%s Chn%d Unit%d", TR(iTr).GetExpName(), thisChannel, thisCluster)));
                        try
                            close(hFigure.Figure)
                        end
						break
					end
				end
			end
		end

		function PETH = BatchPETHistCounts(TR, list, varargin)
			p = inputParser;
			addParameter(p, 'PressOnsetCorrection', cell(size(TR)), @iscell); % Correct actual movement onset time (press only)
			addParameter(p, 'TrialLength', 6, @isnumeric);
			addParameter(p, 'AllowedTrialLength', [0, Inf], @isnumeric);
			addParameter(p, 'SpikeRateWindow', 100, @isnumeric); % in ms
			addParameter(p, 'SpikeRateWindowStim', 10, @isnumeric); % in ms
			addParameter(p, 'ExtendedWindow', 1, @isnumeric);
			addParameter(p, 'ExtendedWindowStim', [-1, 1], @isnumeric);
			addParameter(p, 'Press', true, @islogical); % Whether to process this event
			addParameter(p, 'Lick', false, @islogical); % Whether to process this event
			addParameter(p, 'Stim', false, @islogical); % Whether to process this event
			parse(p, varargin{:});
			pressOnsetCorrection	= p.Results.PressOnsetCorrection;
			trialLength 			= p.Results.TrialLength;
			allowedTrialLength 		= p.Results.AllowedTrialLength;
			spikeRateWindow 		= p.Results.SpikeRateWindow;
			spikeRateWindowStim 	= p.Results.SpikeRateWindowStim;
			extendedWindow 			= p.Results.ExtendedWindow;
			extendedWindowStim 		= p.Results.ExtendedWindowStim;
			press 					= p.Results.Press;
			lick 					= p.Results.Lick;
			stim 					= p.Results.Stim;

			PETH = [];

			for iPlot = 1:size(list, 1)
				thisAnimal 	= list{iPlot, 1};
				thisDate 	= list{iPlot, 2};
				thisChannel = list{iPlot, 3};
				thisCluster = list{iPlot, 4};
				thisExpName = [thisAnimal, '_', num2str(thisDate)];

				for iTr = 1:length(TR)
					if strcmpi(TR(iTr).GetExpName, thisExpName)
						if isempty(PETH)
							iPETH = 1;
						else
							iPETH = length(PETH) + 1;
						end

						thisRefCluster = max(TR(iTr).Spikes(thisChannel).Cluster.Classes); % Use the last cluster is 'noise'/reference cluster
						if press
							[PETH(iPETH).Press, PETH(iPETH).Time, PETH(iPETH).NumTrialsPress, PETH(iPETH).PressSingleTrial] = TR(iTr).PETHistCounts(...
								thisChannel, 'Cluster', thisCluster,...
								'Event', 'PressOn', 'Exclude', 'LickOn',...
								'TrialLength', trialLength, 'AllowedTrialLength', allowedTrialLength, 'ExtendedWindow', extendedWindow, 'SpikeRateWindow', spikeRateWindow,...
								'MoveOnsetCorrection', pressOnsetCorrection{iTr});
						end
						if lick
                            try
                                [PETH(iPETH).Lick, ~, PETH(iPETH).NumTrialsLick, PETH(iPETH).LickSingleTrial] = TR(iTr).PETHistCounts(...
                                    thisChannel, 'Cluster', thisCluster,...
                                    'Event', 'LickOn', 'Exclude', 'PressOn',...
                                    'TrialLength', trialLength, 'AllowedTrialLength', allowedTrialLength, 'ExtendedWindow', extendedWindow, 'SpikeRateWindow', spikeRateWindow);
                            catch
                            	PETH(iPETH).Lick = [];
                            	PETH(iPETH).NumTrialsLick = 0;
                            	PETH(iPETH).LickSingleTrial = [];
                            end
						end
						if stim
							PETH(iPETH).Stim = TR(iTr).PSTHistCounts(...
								thisChannel, 'Cluster', thisCluster,...
								'ExtendedWindow', extendedWindowStim, 'SpikeRateWindow', spikeRateWindowStim);
						end

						PETH(iPETH).TrialLength = trialLength;
						PETH(iPETH).SpikeRateWindow = spikeRateWindow;
						PETH(iPETH).ExtendedWindow = extendedWindow;
						PETH(iPETH).ExpName = [thisAnimal, '_', num2str(thisDate)];
						PETH(iPETH).Channel = thisChannel;
						PETH(iPETH).Cluster = thisCluster;
						PETH(iPETH).PressOnsetCorrection = pressOnsetCorrection{iTr};
						PETH(iPETH).AllowedTrialLength = allowedTrialLength;

						break
					end
				end
			end
		end

		function varargout = HeatMap(PETH, varargin)
			p = inputParser;
			addParameter(p, 'MinNumTrials', 75, @isnumeric);
			addParameter(p, 'MinSpikeRate', 15, @isnumeric);
			addParameter(p, 'Normalization', 'minmax', @ischar); % zscore, minmax, raw
			addParameter(p, 'NormalizationBaselineWindow', [-PETH(1).TrialLength, PETH(1).ExtendedWindow], @(x) isnumeric(x) && length(x) == 2); % Extend window after event
			addParameter(p, 'Sorting', 'latency', @ischar); % abs, gradient, latency
			addParameter(p, 'LatencyThreshold', 0.675, @isnumeric);
			addParameter(p, 'SortsBeforeNorms', false, @islogical); % sort before normalizing
			addParameter(p, 'Window', [-PETH(1).TrialLength, PETH(1).ExtendedWindow], @(x) isnumeric(x) && length(x) == 2); % Extend window after event
			addParameter(p, 'I', [], @isnumeric);
			addParameter(p, 'UseSameSorting', false, @islogical); % use the lever press sorting order for lick
			addParameter(p, 'Lick', false, @islogical);
			addParameter(p, 'LickMinusPress', false, @islogical);
			parse(p, varargin{:});
			minNumTrials 				= p.Results.MinNumTrials;
			minSpikeRate 				= p.Results.MinSpikeRate;
			normalization 				= p.Results.Normalization;
			normalizationBaselineWindow = p.Results.NormalizationBaselineWindow;
			sorting 					= p.Results.Sorting;
			sortsbeforeNorms 			= p.Results.SortsBeforeNorms;
			trialWindow 				= p.Results.Window;
			I 							= p.Results.I;
			useSameSorting				= p.Results.UseSameSorting;
			lick						= p.Results.Lick;
			lickMinusPress				= p.Results.LickMinusPress;

			timestamps 		= PETH(1).Time;
			selectedPress 	= [PETH.NumTrialsPress] >= minNumTrials & cellfun(@mean, {PETH.Press}) > minSpikeRate;
			if lick
				selectedLick 	= [PETH.NumTrialsLick] >= minNumTrials & cellfun(@mean, {PETH.Lick}) > minSpikeRate;
			end
			if lick && useSameSorting
				selectedPress 	= selectedPress & selectedLick;
				selectedLick 	= selectedPress;
			end
			pethPress = transpose(reshape([PETH(selectedPress).Press], length(timestamps), []));
			if lick
				pethLick = transpose(reshape([PETH(selectedLick).Lick], length(timestamps), []));
			end

			inWindow = timestamps <= trialWindow(2) & timestamps > trialWindow(1);
			pethPress = pethPress(:, inWindow);
			if lick
				pethLick  = pethLick(:, inWindow);
			end
			timestamps = timestamps(inWindow);

			baselineSamples = find(timestamps >= normalizationBaselineWindow(1) & timestamps <= normalizationBaselineWindow(2));

			% Normalize/sort PETH
			if ~sortsbeforeNorms
				pethPress = TetrodeRecording.NormalizePETH(pethPress, 'Method', normalization, 'BaselineSamples', baselineSamples);
				if lick
					pethLick  = TetrodeRecording.NormalizePETH(pethLick, 'Method', normalization, 'BaselineSamples', baselineSamples);
				end
			end

			if isempty(I)
				[pethPress, I] = TetrodeRecording.SortPETH(pethPress, 'Method', sorting, 'LatencyThreshold', p.Results.LatencyThreshold);
			else
				pethPress = pethPress(I, :);
			end	
			if lick
				if useSameSorting
					pethLick  = pethLick(I, :);
				else
					pethLick = TetrodeRecording.SortPETH(pethLick, 'Method', sorting, 'LatencyThreshold', p.Results.LatencyThreshold);
				end
			end

			IFull = find(selectedPress);
			IFull = IFull(I);
			varargout = {pethPress, I, IFull};

			if sortsbeforeNorms
				pethPress = TetrodeRecording.NormalizePETH(pethPress, 'Method', normalization, 'BaselineSamples', baselineSamples);
				if lick
					pethLick  = TetrodeRecording.NormalizePETH(pethLick, 'Method', normalization, 'BaselineSamples', baselineSamples);
				end
			end

			if (lick && lickMinusPress)
				pethLick = pethLick - pethPress;
			end

			fMargin = 0.1;
			xMargin = 0.05;
			yMargin = 0.03;
			w = 1 - 2*fMargin - 2*xMargin;
			h = 1 - 2*fMargin - 4*yMargin;
			if lick
				hPress = h*size(pethPress, 1)/(size(pethPress, 1) + size(pethLick, 1));
				hLick  = h*size(pethLick, 1)/(size(pethPress, 1) + size(pethLick, 1));
			else
				hPress = h;
			end
			hFigure = figure('DefaultAxesFontSize', 16);
			% if lick
			% 	hAxesPress = subplot('Position', [fMargin + xMargin, fMargin + 3*yMargin + hLick, w, hPress]);
			% else
			% 	hAxesPress = subplot('Position', [fMargin + xMargin, fMargin + yMargin, w, hPress]);
			% end

			if lick
				hAxesPress = subplot(1, 2, 1);
			else
				hAxesPress = axes(hFigure);
			end
			image(hAxesPress, pethPress, 'CDataMapping','scaled');

			colorbar('Peer', hAxesPress, 'Location', 'EastOutside');
			if ~strcmpi('raw', normalization)
				caxis(hAxesPress, [-6, 6]);
			end
			colormap(hAxesPress, 'jet')

			set(hAxesPress, 'XTick', find(ismember(timestamps, -100:2:100)))
			set(hAxesPress, 'XTickLabel', num2cell(timestamps(hAxesPress.XTick)))
			set(hAxesPress, 'XTickMode', 'manual')
			set(hAxesPress, 'XGrid', 'on')
			set(hAxesPress, 'GridAlpha', 1)
			set(hAxesPress, 'GridColor', 'w')
			set(hAxesPress, 'GridLineStyle', '--')
			set(hAxesPress, 'YTick', unique([1:100:sum(selectedPress), sum(selectedPress)]))
			set(hAxesPress, 'YTickMode', 'manual')
			set(hAxesPress, 'YDir', 'reverse')

			xlabel(hAxesPress, 'Time relative to lever-press (s)')
			ylabel(hAxesPress, 'Unit')
			title(hAxesPress, 'Lever Press')

			if lick
				% hAxesLick = subplot('Position', [fMargin + xMargin, fMargin + yMargin, w, hLick]);
				hAxesLick = subplot(1, 2, 2);

				image(hAxesLick, pethLick, 'CDataMapping','scaled');

				colorbar('Peer', hAxesLick, 'Location', 'EastOutside');
				if ~strcmpi('raw', normalization)
					caxis(hAxesLick, [-6, 6]);
				end
				colormap(hAxesLick, 'jet')

				set(hAxesLick, 'XTick', find(ismember(timestamps, -100:2:100)))
				set(hAxesLick, 'XTickLabel', num2cell(timestamps(hAxesPress.XTick)))
				set(hAxesLick, 'XTickMode', 'manual')
				set(hAxesLick, 'XGrid', 'on')
				set(hAxesLick, 'GridAlpha', 1)
				set(hAxesLick, 'GridColor', 'w')
				set(hAxesLick, 'GridLineStyle', '--')
				set(hAxesLick, 'YTick', unique([1:100:sum(selectedLick), sum(selectedLick)]))
				set(hAxesLick, 'YTickMode', 'manual')
				set(hAxesLick, 'YDir', 'reverse')

				xlabel(hAxesLick, 'Time relative to lick (s)')
				ylabel(hAxesLick, 'Unit')
				title(hAxesLick, 'Lick')
			end
		end

		function varargout = HeatMapStim(PETH, varargin)
			p = inputParser;
			addParameter(p, 'MinNumTrials', 50, @isnumeric);
			addParameter(p, 'MinNumTrains', 15, @isnumeric);
			addParameter(p, 'MinSpikeRate', 15, @isnumeric);
			addParameter(p, 'Sorting', 'latency', @ischar); % abs, gradient, latency, none
			addParameter(p, 'LatencyThreshold', 0.675, @isnumeric);
			addParameter(p, 'Normalization', 'zscore', @ischar); % zscore, minmax, raw
			addParameter(p, 'NormalizationBaselineWindow', [-PETH(1).TrialLength, PETH(1).ExtendedWindow], @(x) isnumeric(x) && length(x) == 2); % Extend window after event
			addParameter(p, 'NormalizationStim', 'zscore', @ischar); % zscore, minmax, raw
			addParameter(p, 'NormalizationBaselineWindowStim', [-1, 0], @(x) isnumeric(x) && length(x) == 2); % Window relative to `On
			addParameter(p, 'Window', [-PETH(1).TrialLength, PETH(1).ExtendedWindow], @(x) isnumeric(x) && length(x) == 2); % Basically XLim for Press PETH
			addParameter(p, 'WindowStim', [-1, 1], @isnumeric);
			addParameter(p, 'StimTrainTypes', [1910, 1606, 1604], @isnumeric);
			addParameter(p, 'CLimStim', [-6, 6], @(x) isnumeric(x) || ischar(x));
			addParameter(p, 'CLimPress', [-6, 6], @(x) isnumeric(x) || ischar(x));
			parse(p, varargin{:});
			minNumTrials 					= p.Results.MinNumTrials;
			minNumTrains 					= p.Results.MinNumTrains;
			minSpikeRate 					= p.Results.MinSpikeRate;
			sorting 						= p.Results.Sorting;
			latencyThreshold				= p.Results.LatencyThreshold;
			normalization 					= p.Results.Normalization;
			normalizationBaselineWindow 	= p.Results.NormalizationBaselineWindow;
			normalizationStim				= p.Results.NormalizationStim;
			normalizationBaselineWindowStim = p.Results.NormalizationBaselineWindowStim;
			windowPress						= p.Results.Window;
			windowStim 						= p.Results.WindowStim;
			stimTrainTypes 					= p.Results.StimTrainTypes;
			cLimStim 						= p.Results.CLimStim;
			cLimPress 						= p.Results.CLimPress;

			numSubplots = length(stimTrainTypes) + 1;

			fMargin = 0.05;
			spacing = 0.005;
			w = 1 - 2*fMargin - (numSubplots - 1)*spacing;
			w1 = w/numSubplots;
			h = 1 - 2*fMargin;

			hFigure = figure('DefaultAxesFontSize', 14);
			hAxesPress = subplot('Position', [fMargin, fMargin, w1, h]);
			for iTrainType = 1:length(stimTrainTypes)
				hAxesStim(iTrainType) = subplot('Position', [fMargin + (w1 + spacing) * iTrainType, fMargin, w1, h]);
				title(hAxesStim(iTrainType), num2str(stimTrainTypes(iTrainType)));
			end

			% First process PETH for press
			timestampsPress = PETH(1).Time;
			selectedPress 	= [PETH.NumTrialsPress] >= minNumTrials & cellfun(@mean, {PETH.Press}) > minSpikeRate;
			pethPress = transpose(reshape([PETH(selectedPress).Press], length(timestampsPress), []));
			inWindow = timestampsPress <= windowPress(2) & timestampsPress > windowPress(1);
			timestampsPress = timestampsPress(inWindow);
			baselineSamples = find(timestampsPress >= normalizationBaselineWindow(1) & timestampsPress <= normalizationBaselineWindow(2));

			% Normalize
			pethPress = TetrodeRecording.NormalizePETH(pethPress, 'Method', normalization, 'BaselineSamples', baselineSamples);

			% Sort
			[~, I] = TetrodeRecording.SortPETH(pethPress(:, inWindow), 'Method', sorting, 'LatencyThreshold', latencyThreshold);
			pethPress = pethPress(I, :);

			iBPL = find(selectedPress);
			iBPL = iBPL(I);

			varargout = {iBPL};

			% Plot Press PETH
			image(hAxesPress, pethPress(:, inWindow), 'CDataMapping','scaled');

			colorbar('Peer', hAxesPress, 'Location', 'EastOutside');
			if ~strcmpi('raw', normalization)
				caxis(hAxesPress, cLimPress);
			end
			colormap(hAxesPress, 'jet')
			set(hAxesPress, 'XTick', find(ismember(timestampsPress, -100:2:100)))
			set(hAxesPress, 'XTickLabel', num2cell(timestampsPress(hAxesPress.XTick)))
			title(hAxesPress, 'Lever Press')
			ylabel(hAxesPress, 'Unit')
			xlabel(hAxesPress, 'Time (Press = 0)')

			% Plot PSTH for each train type
			for iTrainType = 1:length(stimTrainTypes)
				timestampsStim{iTrainType} = [];
				numTrains{iTrainType} = [];
				trainLength{iTrainType} = round(stimTrainTypes(iTrainType)/100)/10;
				numPulses{iTrainType} = stimTrainTypes(iTrainType) - trainLength{iTrainType}*1000;
				pulseOn{iTrainType} = [];
				pulseOff{iTrainType} = [];
				psth{iTrainType} = [];
				numSamples{iTrainType} = (trainLength{iTrainType} + 2)*100;
				for iPETH = find(selectedPress)
					i = find([PETH(iPETH).Stim.TrainType] == stimTrainTypes(iTrainType));
					
					if ~isempty(i)
						Stim = PETH(iPETH).Stim(i);
					end

					hasEnoughTrains	= Stim.NumTrains >= minNumTrains;

					% Stim exists for this unit
					if hasEnoughTrains && ~isempty(i)
						% These things are filled once
						if isempty(timestampsStim{iTrainType})
							timestampsStim{iTrainType} = Stim.Timestamps;
						end
						if isempty(pulseOn{iTrainType})
							pulseOn{iTrainType} = Stim.PulseOn;
						end
						if isempty(pulseOff{iTrainType})
							pulseOff{iTrainType} = Stim.PulseOff;
						end

						% These are accumulated data
						if isempty(psth{iTrainType})
							psth{iTrainType} = Stim.SpikeRate;
						else
							psth{iTrainType} = [psth{iTrainType}; Stim.SpikeRate];
						end
						numTrains{iTrainType} = [numTrains{iTrainType}, Stim.NumTrains];
					else
						psth{iTrainType} = [psth{iTrainType}; nan(1, numSamples{iTrainType})];
						if hasEnoughTrains
							numTrains{iTrainType} = [numTrains{iTrainType}, 0];
						else
							numTrains{iTrainType} = [numTrains{iTrainType}, Stim.NumTrains];
						end
					end
				end

				% Plot this train, choo choo
				% Normalize
				inWindowStim = timestampsStim{iTrainType} <= (windowStim(2) + trainLength{iTrainType}) & timestampsStim{iTrainType} >= windowStim(1);
				baselineSamplesStim{iTrainType} = find(timestampsStim{iTrainType} >= normalizationBaselineWindowStim(1) & timestampsStim{iTrainType} <= normalizationBaselineWindowStim(2));
				psth{iTrainType} = TetrodeRecording.NormalizePSTH(psth{iTrainType}, 'Method', normalizationStim, 'BaselineSamples', baselineSamplesStim{iTrainType});
				% Plot
				image(hAxesStim(iTrainType), psth{iTrainType}(I, inWindowStim), 'CDataMapping','scaled');

				colorbar('Peer', hAxesStim(iTrainType), 'Location', 'EastOutside');
				if ~strcmpi('raw', normalizationStim)
					caxis(hAxesStim(iTrainType), cLimStim);
				end
				colormap(hAxesStim(iTrainType), 'jet')

				xlabel(hAxesStim(iTrainType), 'Time (TrainOn = 0)')
				set(hAxesStim(iTrainType), 'XTick', find(ismember(round(timestampsStim{iTrainType}(inWindowStim)*100)/100, round([pulseOn{iTrainType}, pulseOff{iTrainType}]*100)/100)))
				set(hAxesStim(iTrainType), 'XTickLabel', num2cell(timestampsStim{iTrainType}(hAxesStim(iTrainType).XTick)))
				title(hAxesStim(iTrainType), ['StimType - ', num2str(stimTrainTypes(iTrainType))])
            end

            % Common formatting
			set([hAxesPress, hAxesStim], 'XTickMode', 'manual')
			set([hAxesPress, hAxesStim], 'XGrid', 'on')
			set([hAxesPress, hAxesStim], 'GridAlpha', 1)
			set(hAxesPress, 'GridColor', 'w')
			set(hAxesStim, 'GridColor', 'g')
			set([hAxesPress, hAxesStim], 'GridLineStyle', '--')
			set([hAxesPress, hAxesStim], 'YTick', unique([1:25:sum(selectedPress), sum(selectedPress)]))
			set([hAxesPress, hAxesStim], 'YTickMode', 'manual')
			set([hAxesPress, hAxesStim], 'YDir', 'reverse')
        end

        % Ignoring stim pattern, only plot the first stim.
		function varargout = HeatMapStimSimple(PETH, varargin)
			p = inputParser;
			addParameter(p, 'MinNumTrials', 50, @isnumeric);
			addParameter(p, 'MinSpikeRate', 15, @isnumeric);
			addParameter(p, 'Sorting', 'latency', @ischar); % abs, gradient, latency, none
			addParameter(p, 'LatencyThreshold', 0.675, @isnumeric);
			addParameter(p, 'Normalization', 'zscore', @ischar); % zscore, minmax, raw
			addParameter(p, 'NormalizationBaselineWindow', [-PETH(1).TrialLength, PETH(1).ExtendedWindow], @(x) isnumeric(x) && length(x) == 2); % Extend window after event
			addParameter(p, 'NormalizationStim', 'zscore', @ischar); % zscore, minmax, raw
			addParameter(p, 'NormalizationBaselineWindowStim', [-1, 0], @(x) isnumeric(x) && length(x) == 2); % Window relative to `On
			addParameter(p, 'Window', [-PETH(1).TrialLength, PETH(1).ExtendedWindow], @(x) isnumeric(x) && length(x) == 2); % Basically XLim for Press PETH
			addParameter(p, 'WindowStim', [-1, 0.1], @isnumeric);
			addParameter(p, 'CLimStim', [-6, 6], @(x) isnumeric(x) || ischar(x));
			addParameter(p, 'CLimPress', [-6, 6], @(x) isnumeric(x) || ischar(x));
			parse(p, varargin{:});
			minNumTrials 					= p.Results.MinNumTrials;
			minSpikeRate 					= p.Results.MinSpikeRate;
			sorting 						= p.Results.Sorting;
			latencyThreshold				= p.Results.LatencyThreshold;
			normalization 					= p.Results.Normalization;
			normalizationBaselineWindow 	= p.Results.NormalizationBaselineWindow;
			normalizationStim				= p.Results.NormalizationStim;
			normalizationBaselineWindowStim = p.Results.NormalizationBaselineWindowStim;
			windowPress						= p.Results.Window;
			windowStim 						= p.Results.WindowStim;
			cLimStim 						= p.Results.CLimStim;
			cLimPress 						= p.Results.CLimPress;

			numSubplots = 2;

			fMargin = 0.05;
			spacing = 0.005;
			w = 1 - 2*fMargin - (numSubplots - 1)*spacing;
			w1 = w/numSubplots;
			h = 1 - 2*fMargin;

			hFigure = figure('DefaultAxesFontSize', 14);
			hAxesPress = subplot(1, 2, 1);
			hAxesStim = subplot(1, 2, 2);

			% First process PETH for press
			timestampsPress = PETH(1).Time;
			selectedPress 	= [PETH.NumTrialsPress] >= minNumTrials & cellfun(@mean, {PETH.Press}) > minSpikeRate;
			pethPress = transpose(reshape([PETH(selectedPress).Press], length(timestampsPress), []));
			inWindow = timestampsPress <= windowPress(2) & timestampsPress > windowPress(1);
			timestampsPress = timestampsPress(inWindow);

			% Normalize
			baselineSamples = timestampsPress >= normalizationBaselineWindow(1) & timestampsPress <= normalizationBaselineWindow(2);
			pethPress = TetrodeRecording.NormalizePETH(pethPress, 'Method', normalization, 'BaselineSamples', baselineSamples);

			% Sort
			[~, I] = TetrodeRecording.SortPETH(pethPress(:, inWindow), 'Method', sorting, 'LatencyThreshold', latencyThreshold);
			pethPress = pethPress(I, :);

			iBPL = find(selectedPress);
			iBPL = iBPL(I);

			varargout = {iBPL};

			% Plot Press PETH
			image(hAxesPress, pethPress(:, inWindow), 'CDataMapping','scaled');

			colorbar('Peer', hAxesPress, 'Location', 'EastOutside');
			if ~strcmpi('raw', normalization)
				caxis(hAxesPress, cLimPress);
			end
			colormap(hAxesPress, 'jet')
			set(hAxesPress, 'XTick', find(ismember(timestampsPress, -100:2:100)))
			set(hAxesPress, 'XTickLabel', num2cell(timestampsPress(hAxesPress.XTick)))
			title(hAxesPress, 'Lever Press')
			ylabel(hAxesPress, 'Unit')
			xlabel(hAxesPress, 'Time (Press = 0)')

			% Plot PSTH for stim
            % Calculate weighted average for all stim types.
            Stim = PETH(1).Stim(1);
            spikeRateWindowStim = (Stim.Timestamps(end) - Stim.Timestamps(1)) / (length(Stim.Timestamps) - 1);
            numSamplesPSTH = round(diff(windowStim) / spikeRateWindowStim) + 1;
            psth = zeros(length(selectedPress), numSamplesPSTH);
            psthNormalized = zeros(size(psth));
            numTrains = zeros(sum(selectedPress), 1);
            timestampsStim = zeros(size(psth));
			for iPETH = find(selectedPress)
                % Consolidate all stim types in to one PSTH.
                for iStimType = 1:length(PETH(iPETH).Stim)
                    Stim = PETH(iPETH).Stim(iStimType);
                    inStimWindow = Stim.Timestamps >= windowStim(1) & Stim.Timestamps <= windowStim(2);
                    psth(iPETH, :) = psth(iPETH, :) + Stim.NumTrains * Stim.SpikeRate(inStimWindow);
                    
                    % Normalize
                    baselineSamplesStim = Stim.Timestamps >= normalizationBaselineWindowStim(1) & Stim.Timestamps <= normalizationBaselineWindowStim(2);
                    normalizedSpikeRate = TetrodeRecording.NormalizePSTH(Stim.SpikeRate(inStimWindow), 'Method', normalizationStim, 'BaselineSamples', baselineSamplesStim);
                    psthNormalized(iPETH, :) = psthNormalized(iPETH, :) + Stim.NumTrains * normalizedSpikeRate;
                    
                    % Timestamps
                    timestampsStim(iPETH, :) = timestampsStim(iPETH, :) + Stim.NumTrains * Stim.Timestamps(inStimWindow);
                end
                numTrains(iPETH) = sum([PETH(iPETH).Stim.NumTrains]);
                psth(iPETH, :) = psth(iPETH, :) ./ numTrains(iPETH);
                psthNormalized(iPETH, :) = psthNormalized(iPETH, :) ./ numTrains(iPETH);
                timestampsStim(iPETH, :) = timestampsStim(iPETH, :) ./ numTrains(iPETH);
			end

			% Plot this train, choo choo
			% Plot
			image(hAxesStim, psthNormalized(I, :), 'CDataMapping','scaled');

			colorbar('Peer', hAxesStim, 'Location', 'EastOutside');
			if ~strcmpi('raw', normalizationStim)
				caxis(hAxesStim, cLimStim);
			end
			colormap(hAxesStim, 'jet')

			xlabel(hAxesStim, 'Time (TrainOn = 0)')
            timestampsStim = mean(timestampsStim, 1);
            [~, I] = min(abs(timestampsStim - windowStim(2)));
            xtickPos = [find(timestampsStim == 0), I];
            xtickLabels = num2cell(timestampsStim(xtickPos));
            set(hAxesStim, 'XTick', xtickPos)
            set(hAxesStim, 'XTickLabel', xtickLabels)
            
            % Common formatting
			set([hAxesPress, hAxesStim], 'XTickMode', 'manual')
			set([hAxesPress, hAxesStim], 'XGrid', 'on')
			set([hAxesPress, hAxesStim], 'GridAlpha', 1)
			set(hAxesPress, 'GridColor', 'w')
			set(hAxesStim, 'GridColor', 'g')
			set([hAxesPress, hAxesStim], 'GridLineStyle', '--')
			set([hAxesPress, hAxesStim], 'YTick', unique([1:25:sum(selectedPress), sum(selectedPress)]))
			set([hAxesPress, hAxesStim], 'YTickMode', 'manual')
			set([hAxesPress, hAxesStim], 'YDir', 'reverse')
		end

		function varargout = SortPETH(peth, varargin)
			p = inputParser;
			addParameter(p, 'Method', 'latency', @ischar); % abs, gradient, latency
			addParameter(p, 'LatencyThreshold', 0.675, @isnumeric);
			parse(p, varargin{:});
			method 				= p.Results.Method;
			latencyThreshold 	= p.Results.LatencyThreshold;

            if isempty(peth)
                varargout = {peth, [], []};
                return
            end
            
			switch lower(method)
				% Sort by max(abs(trace))
				case 'abs'
					[~, I] = sort(max(abs(peth), [], 2));
				% Sort by max(diff(trace))
				case 'gradient'
					[~, I] = sort(max(transpose(abs(diff(transpose(peth)))), [], 2));
				% Sort by how fast gradients change
				case 'latency'
					% Because the grad student is dumb
					for iCell = 1:size(peth, 1)
						thisPeth = peth(iCell, :);
						% thisSigma = nanmedian(abs(thisPeth))/0.6745;
						thisSigma = std(thisPeth);
						thisZScoredPeth = (thisPeth - mean(thisPeth))/thisSigma;
                        changeTime = find(abs(thisZScoredPeth) >= latencyThreshold*max(abs(thisZScoredPeth)), 1);
                        if (isempty(changeTime))
                            whenDidFiringRateChange(iCell) = size(peth, 2);
                        else
                            whenDidFiringRateChange(iCell) = changeTime;
                            whenDidFiringRateChange(iCell) = whenDidFiringRateChange(iCell)*sign(thisZScoredPeth(whenDidFiringRateChange(iCell)));
                        end
					end
					[whenDidFiringRateChange, I] = sort(whenDidFiringRateChange);
					I = flip(I);
                case 'raw'
                    I = 1:size(peth, 1);
                    whenDidFiringRateChange = repmat(size(peth, 2), size(peth, 1), 1);
			end

			peth = peth(I, :);

			varargout = {peth, I, whenDidFiringRateChange};
		end

		function peth = NormalizePETH(peth, varargin)
			p = inputParser;
			addParameter(p, 'Method', 'zscore', @ischar); % zscore, minmax
			addParameter(p, 'BaselineSamples', [], @(x) isnumeric(x) || islogical(x)); % zscore, minmax
			parse(p, varargin{:});
			method = p.Results.Method;
			samples = p.Results.BaselineSamples;

			if isempty(samples)
				samples = 1:length(peth, 2);
			end

			for iCell = 1:size(peth, 1)
				thisPeth = peth(iCell, :);
				switch lower(method)
					case 'zscore'
						thisSigma = std(thisPeth(samples));
						thisMean = mean(thisPeth(samples));
						peth(iCell, :) = (thisPeth - thisMean)/thisSigma;
					case 'minmax'
						peth(iCell, :) = (thisPeth - min(thisPeth(samples)))/(max(thisPeth(samples)) - min(thisPeth(samples)));
					otherwise
						peth(iCell, :) = thisPeth;
				end
			end
		end

		function psth = NormalizePSTH(psth, varargin)
			p = inputParser;
			addParameter(p, 'Method', 'zscore', @ischar); % zscore, minmax
			addParameter(p, 'BaselineSamples', [], @(x) isnumeric(x) || islogical(x));
			parse(p, varargin{:});
			method = p.Results.Method;
			samples = p.Results.BaselineSamples;

			if isempty(samples)
				samples = 1:length(psth, 2);
			end

			switch lower(method)
				case 'zscore'
					thisSigma = nanstd(psth(:, samples), 0, 2);
					thisMean = nanmean(psth(:, samples), 2);
					psth = (psth - thisMean)./thisSigma;
				case 'minmax'
					peth = (psth - min(psth(:, samples), [], 2,'omitnan'))./(max(psth(:, samples), [], 2,'omitnan') - min(psth(:, samples), [], 2,'omitnan'));
				otherwise
					psth = psth;
			end

		end

		function [thisColor, thisStyle] = GetColorAndStyle(iClass, varargin)
			p = inputParser;
			addParameter(p, 'Colors', 'rgbmyck', @ischar);
			addParameter(p, 'Styles', {'-', '--', ':', '-.'}, @iscell);
			parse(p, varargin{:});
			colors 	= p.Results.Colors;
			styles 	= p.Results.Styles;

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
			if ischar(event)
				grad = diff(['0', event]);
			else
				grad = diff([0, event]);
			end
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

		function BatchProcessStimData(TR, PTR, varargin)
			p = inputParser;
			addParameter(p, 'Window', [-20, 20], @isnumeric); % how many milliseconds before and after stim on to read.
			parse(p, varargin{:});
			readWindow = p.Results.Window * 0.001;

			for iTr = 1:length(TR)
                try
                    startTime = tic();
                    TetrodeRecording.TTS(['Processing file ', num2str(iTr), '/', num2str(length(TR)), '...']);
                    StimData = TR(iTr).ReadStimData(PTR(iTr), [TR(iTr).Spikes.Channel], 'Window', readWindow);
                    StimData = TetrodeRecording.SerializeStimData(StimData);

                    filename = strsplit(TR(iTr).Files{1}, 'nev');
                    filename = [TR(iTr).Path, '\..\SpikeSort\sd_', filename{1}, '.mat'];

                    save(filename, 'StimData', '-v7.3');

                    duration = toc(startTime);
                    TetrodeRecording.TTS([num2str(duration), ' s.\n']);
                catch
                    warning('Failed trying to preocess file');
                end
			end
		end

		% Collision test
		function StimData = SerializeStimData(StimData)
			numTrains = length(StimData.TrainOn);
			numPulses = length(StimData.StimOn);
			% maxPulseLength = ceil((max(StimData.StimOff - StimData.StimOn) + diff(StimData.Window)) * StimData.SampleRate) + 1;
			maxPulseLength = ceil(diff(StimData.Window) * StimData.SampleRate) + 1;
			numChannels = size(StimData.Data, 2);

			iPulse = 0;
			StimData.DataByPulse = zeros(numPulses, numChannels, maxPulseLength);
			StimData.TimestampByPulse = zeros(numPulses, maxPulseLength);

			for iTrain = 1:numTrains
				stimOn = StimData.StimOn;
				isPulseInTrain = (stimOn >= StimData.TrainOn(iTrain) - 0.01) & (stimOn <= StimData.TrainOff(iTrain) + 0.01);
				stimOn = stimOn(isPulseInTrain);

				for iPulseInTrain = 1:length(stimOn)
					iPulse = iPulse + 1;
					isSampleInPulse = (StimData.Timestamps(iTrain, :) >= stimOn(iPulseInTrain) + StimData.Window(1)) & (StimData.Timestamps(iTrain, :) <= stimOn(iPulseInTrain) + StimData.Window(2));
					thisPulseData = StimData.Data(iTrain, :, isSampleInPulse);
					StimData.DataByPulse(iPulse, :, 1:length(thisPulseData)) = thisPulseData;
					StimData.TimestampByPulse(iPulse, 1:length(thisPulseData)) = StimData.Timestamps(iTrain, isSampleInPulse);
				end
			end
		end

		function PlotCollisionTest(StimData, channel, varargin)
			p = inputParser;
			addParameter(p, 'StartFromTrace', 1, @isnumeric);
			addParameter(p, 'TracesPerPage', 25, @isnumeric);
			parse(p, varargin{:});
			startFromTrace	= p.Results.StartFromTrace;
			tracesPerPage 	= p.Results.TracesPerPage;

			% Normalize each pulse to this range.
			yLim = [-500, 500];
			ySpacing = 1;

			fig = figure();
			ax = axes(fig);
			grid(ax, 'on');
			hold(ax, 'on');
			xlim(ax, [-10, 10])
			
			xlabel(ax, 'Time from StimOn (ms)')
			ylabel(ax, 'Normalized voltage')

			numPlotted = 0;

			% Normalize voltage data and align time by stimOnset;
			for iPulse = startFromTrace:length(StimData.StimOn)
				isNonZero = find(StimData.TimestampByPulse(iPulse, :));

				y = (StimData.DataByPulse(iPulse, channel, isNonZero) - yLim(1)) / (yLim(2) - yLim(1)) + iPulse * ySpacing;
				t = StimData.TimestampByPulse(iPulse, isNonZero) - StimData.StimOn(iPulse);

				plot(ax, t * 1000, squeeze(y), 'k');

				numPlotted = numPlotted + 1;
				if (numPlotted >= tracesPerPage)
					w = waitforbuttonpress();
					cla(ax);
					numPlotted = 0;
				end
			end

			hold(ax, 'off');
		end

		% Find first event after reference
		function varargout = FindFirstInTrial(reference, event, eventExclude, firstOrLast)
			if nargin < 4
				firstOrLast = 'first';
			end
			if nargin < 3
				eventExclude = [];
			end
			edges = [reference(1:end - 1), max(event(end), reference(end))];
			[~, ~, bins] = histcounts(event, edges);
			event = event(bins ~= 0);
			bins = nonzeros(bins);
			[iReference, iEvent] = unique(bins, firstOrLast);
			reference = reference(iReference);
			event = event(iEvent);
            
            toRemove = [];
            
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

			varargout = {reference, event, iReference, iEvent, toRemove};
		end

		% Find last event after reference
		function varargout = FindLastInTrial(reference, event, eventExclude)
			if nargin < 3
				eventExclude = [];
			end
			[reference, event, iReference, iEvent, toRemove] = TetrodeRecording.FindFirstInTrial(reference, event, eventExclude, 'last');
			varargout = {reference, event, iReference, iEvent, toRemove};
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

		function rig = GetRig(filepath)
			if contains(filepath, {'desmond10', 'desmond11', 'desmond12', 'daisy4', 'desmond14', 'desmond16', 'desmond18', 'daisy7', 'desmond21', 'desmond22'})
				rig = 1;
			elseif contains(filepath, {'desmond13', 'daisy5', 'desmond15', 'desmond17', 'desmond19', 'desmond20', 'daisy8'})
				rig = 2;
			else
				error('This version was designed for desmond12/13 daisy4/5 only');
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
			% if speak
			% 	txt = strsplit(txt, '\');
			% 	txt = txt{1};
			% 	tts(txt);
			% end
		end

		function RandomWords()
			words = {...
				'49 times, we fought that beast. It had a chicken head with duck feet, with a woman''s face too.',...
				'Every day I worry all day about what''s waiting in the bushes of love. Something''s waiting in the bushes for us, something''s waiting in the bushes of love.',...
				'Now run, run, run, jump! I can be a backpack while you run.',...
				'Rocking, rocking and rolling. Down to the beach I''m strolling. But the seagulls poke at my head, not fun! I said "Seagulls... mmgh! Stop it now!"',...
				'Some day when you are older you could get hit by a boulder. While you''re lying there screaming come help me please, the seagulls come poke your knees.',...
				'Even though it looks like it''s the future. It''s really a long long time ago when there were knights and they got into fights using sabres of light!"',...
				'Twenty nights in the ice is a long time when there''s hostiles on the hill. It''s not about what they want, you just gotta walk your walk.',...
				'They''ve got the ultimate power in the universe. Before it gets better it''s getting worse.',...
				'It''s not a satellite but it can light up the sky when it blows. Evil comes in round shapes.',...
				'I don''t care what you said to the man with the black head. Where is the rubber band? You''re ruining my plan!'...
			};
			TetrodeRecording.TTS([words{randi(length(words))}, '\n'], true)
		end
	end
end
