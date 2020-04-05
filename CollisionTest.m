classdef CollisionTest < handle
    properties
        SampleRate = NaN
        Data = []
        Timestamps = []
        PulseOn = []
        PulseOff = []
        TrainOn = []
        TrainOff = []
        Window = [0, 0]
        Filename
        ExpName = ''
    end

    properties (Transient, Access = {})
        TR
        PTR
    end

    % Public methods
    methods
        % Constructor
        function obj = CollisionTest(varargin)
        %CollisionTest - Class for analyzing antidromic stimulation/and runnung collision tests.
        %
        % Syntax: ct = CollisionTest()
        % 
        % Optional Parameters: 
        % 'Folder'          - Skips file selection prompt if a valid path is provided.
        % 'Read'            - (Default TRUE) Whether or not to read tr/ptr/ns5 files. 
        % 'Channels'        - Which channels to look at. This is channel ID as shown in TetrodeRecording.PlotChannel([]). Leave empty to use all channels.
        % 'MaxTrains'       - Limit how many stim trains will be read.
        % 'ExtendedWindow'  - (ms) Extended window for reading data.
        % 
        % 
            p = inputParser();
            p.addParameter('Folder', '', @ischar);
            p.addParameter('Read', true, @islogical);
            p.addParameter('Channels', [], @isnumeric); % Which channels to look at. This is channel ID as shown in TetrodeRecording.PlotChannel([])
            p.addParameter('MaxTrains', [], @isnumeric); % Limit how many stim trains will be read.
            p.addParameter('ExtendedWindow', [0, 0], @(x) isnumeric(x) && length(x) == 2); % (ms) Extended window for reading data.
            p.parse(varargin{:});
            folder = p.Results.Folder;
            doRead = p.Results.Read;
            obj.Window = p.Results.ExtendedWindow;

            % File selection prompt if folder not specified during constructor call.
            if isempty(folder)
                folder = uipickfiles('Prompt', 'Select a folder containing an ns5 file.', 'NumFiles', 1);
                folder = folder{1};
            end

            % Validate folder
            [~, obj.Filename, obj.ExpName, allValid] = CollisionTest.validateFolders(folder, 'SuppressWarnings', true);
            if ~allValid
                error('Invalid folder %s', folder);
            end

            if (doRead)
                obj.read(p.Results.Channels, p.Results.MaxTrains, p.Results.ExtendedWindow ./ 1000);
            end
        end

        function plot(obj, channel, varargin)
        %plot - Plot electrode data aligned on stimulation onset events.
        % Syntax: plot(obj, channel, varargin)

            p = inputParser();
            p.addRequired('Channel', @(x) isnumeric(x));
            p.addParameter('Start', 1, @isnumeric);
            p.addParameter('TracesPerPage', 25, @isnumeric);
            p.addParameter('YLim', [-500, 500], @(x) isnumeric(x) && length(x) == 2 && diff(x) > 0);
            p.addParameter('XLim', obj.Window, @(x) isnumeric(x) && length(x) == 2 && diff(x) > 0);
            p.addParameter('YSpacing', 1, @(x) isnumeric(x) && length(x) == 1);
            p.parse(channel, varargin{:})
            channel = p.Results.Channel;
            startFromTrace = p.Results.Start;
            tracesPerPage = p.Results.TracesPerPage;
            yRange = p.Results.YLim;
            xRange = p.Results.XLim;
            ySpacing = p.Results.YSpacing;
            
            % Create figure
            fig = figure(); 
            ax = axes(fig);
            grid(ax, 'on');
            hold(ax, 'on');
            xlim(ax, xRange);
            % ylim(ax, yRange);

            xlabel(ax, 'Time from stim on (ms)')
            ylabel(ax, 'Normalized voltage (a.u.)')

            iPulseInPage = 0;
            iPulse = startFromTrace - 1;

            plotWindow = 0.001 * obj.Window;

            % Normalize voltage data and align to stimOnsetTime;
            while iPulse <= length(obj.PulseOn)
                iPulse = iPulse + 1;
                iPulseInPage = iPulseInPage + 1;

                isInPulse = obj.Timestamps > obj.PulseOn(iPulse) + plotWindow(1) & obj.Timestamps <= obj.PulseOn(iPulse) + plotWindow(2);
                pulseData = obj.Data(isInPulse, channel);
                pulseTimestamps = obj.Timestamps(isInPulse);

                % Normalize voltage to yRange.
                y = (pulseData - yRange(1)) ./ diff(yRange) + iPulse * ySpacing;
                
                % Align time to stimOn
                t = 1000 * (pulseTimestamps - obj.PulseOn(iPulse));

                % Plot trace
                plot(ax, t, y, 'k');

                % Plot stim window
                stimOnVertices(2 * iPulseInPage - 1: 2 * iPulseInPage, 1) = 0;
                stimOnVertices(2 * iPulseInPage - 1: 2 * iPulseInPage, 2) = [iPulse * ySpacing, iPulse * ySpacing + 1];
                stimOffVertices(2 * iPulseInPage - 1: 2 * iPulseInPage, 1) = 1000 * (obj.PulseOff(iPulse) - obj.PulseOn(iPulse));
                stimOffVertices(2 * iPulseInPage - 1: 2 * iPulseInPage, 2) = [iPulse * ySpacing, iPulse * ySpacing + 1];

                % Page done
                if (iPulseInPage >= tracesPerPage)
                    stimPatchVertices = vertcat(stimOnVertices, stimOffVertices(end:-1:1, :));
                    patch('XData', stimPatchVertices(:, 1), 'YData', stimPatchVertices(:, 2), 'FaceColor', '#4DBEEE', 'FaceAlpha', 0.33, 'EdgeAlpha', 0);
                    ylim([(iPulse - iPulseInPage + 1) * ySp, iPulse * ySpacing + 1])
                    title(ax, sprintf('Pulses %d - %d', iPulse - iPulseInPage + 1, iPulse))
                    [~, ~, button] = ginput(1);
                    
                    % Right click -> prev page
                    if (button == 3)
                        iPulse = iPulse - size(stimPatchVertices, 1) / 2;
                        iPulse = max(0, iPulse);
                    end
                    cla(ax);
                    stimPatchVertices = [];
                    iPulseInPage = 0;
                end
            end

            hold(ax, 'off');
        end

        function save(obj, varargin)
        %save - Save the object to a .mat file.
        %
        % Syntax: save(obj, varargin)
        % 
            p = inputParser();
            p.addOptional('Filename', '', @ischar);
            p.parse(varargin{:});
            filename = p.Results.Filename;

            if (isempty(filename))
                uisave('obj', sprintf('ct_%s.mat', obj.ExpName));
            else
                fprintf('Saving to file "%s"...', filename)
                save(filename, 'obj', '-v7.3');
                fprintf('Done.\n')
            end
        end
    end

    % Public static methods
    methods (Static)
        function batch(varargin)
        %Batch - Batch generate and save CollisionTest objects from tr/ptr/ns5 files.
        %
        % Optional Parameters: 
        % 'ExtendedWindow'  - (ms) Extended window for reading data.
        % 
        % Syntax: Batch(varargin)
        %
            p = inputParser();
            p.addParameter('ExtendedWindow', [-20, 20], @(x) isnumeric(x) && length(x) == 2 && diff(x) > 0);
            p.parse(varargin{:});

            folders = uipickfiles('Prompt', 'Select multiple folders each containing an ns5 file.');

            folders = CollisionTest.validateFolders(folders, 'SuppressWarnings', false);

            for iFolder = 1:length(folders)
                folder = folders{iFolder};
                fprintf('Processing folder %d of %d - %s:\n', iFolder, length(folders), folder);

                ct = CollisionTest('Folder', folder, 'ExtendedWindow', p.Results.ExtendedWindow);

                % Make sure save path exists
                saveFolder = sprintf('%s//..//CollisionTest', folder);
                if ~isfolder(saveFolder)
                    mkdir(saveFolder);
                end
                ct.save(sprintf('%s//ct_%s.mat', saveFolder, ct.ExpName));
            end
        end

        function load()
        %Load - Load multiple ct files.
        %
            files = uipickfiles('Prompt', 'Select (multiple) .mat files containing an CollisionTest object named "obj"', 'FilterSpec', '*.mat');
        end
    end

    % Private methodd
    methods (Access = {})
        function read(obj, channels, maxTrains, extendedWindow)
        %read - Read data from tr/ptr/ns5 files.
        % Syntax: read(obj, varargin)
            obj.PTR = TetrodeRecording.BatchLoad({obj.Filename.PTR});
            obj.TR = TetrodeRecording.BatchLoad(obj.Filename.TR);

            [obj.PulseOn, obj.PulseOff, obj.TrainOn, obj.TrainOff] = obj.readDigitalEvents();
            [obj.Data, obj.Timestamps, obj.SampleRate] = readAnalogTrains(obj, channels, maxTrains, extendedWindow, obj.TrainOn, obj.TrainOff);
        end

        function varargout = readDigitalEvents(obj)
        %readDigitalEvents - Read digital events (stimulus onsets & offsets).
            % Read digital events
			cueOn = sort(obj.TR.DigitalEvents.CueOn);
			pulseOn = sort(obj.TR.DigitalEvents.StimOn);
			pulseOff = sort(obj.TR.DigitalEvents.StimOff);

			% Get the start and end timestamps of a stim train.
			[~, trainOn] = TetrodeRecording.FindFirstInTrial(cueOn, pulseOn);
            [~, trainOff] = TetrodeRecording.FindLastInTrial(cueOn, pulseOff);
            
            varargout = {pulseOn, pulseOff, trainOn, trainOff};
        end

        function [data, timestamps, sampleRate] = readAnalogTrains(obj, channels, maxTrains, extendedWindow, trainOn, trainOff)
        %readAnalogTrains - Read analog data (only during stim trains).
            if ~isempty(maxTrains)
                trainOn = trainOn(1:maxTrains);
                trainOff = trainOff(1:maxTrains);
            end

            if isempty(channels)
                channels = [obj.TR.Spikes.Channel];
            end

            sampleRate = obj.TR.FrequencyParameters.AmplifierSampleRate;

            % Pre-allocate
            nsxChannels =  obj.mapChannels(channels);
            maxTrainLength = max(trainOff - trainOn) + diff(extendedWindow);
            maxTrainLength = ceil(sampleRate * maxTrainLength);

            dataByTrain = NaN(length(trainOn), maxTrainLength, length(channels));
            timestampsByTrain = NaN(length(trainOn), maxTrainLength);

            % Read data from each train and put them in an array.
            sysInitDelay = obj.TR.FrequencyParameters.SysInitDelay.Duration; % Get sysInitDelay to shift openNSx window to the correct value.
            if isnan(sysInitDelay)
                sysInitDelay = 0;
            end
            
            tTic = tic();
            fprintf('Reading %d trains, %d channels...', length(trainOn), length(channels));

            for iTrain = 1:length(trainOn)
                readWindow = [trainOn(iTrain) + extendedWindow(1), trainOff(iTrain) + extendedWindow(2)] + sysInitDelay;

                NSx = openNSx(obj.Filename.NSx, 'read', 'channels', nsxChannels, 'duration', readWindow, 'sec');
                numSamples = size(NSx.Data, 2);

                dataByTrain(iTrain, 1:numSamples, :) = reshape(transpose(NSx.Data), 1, numSamples, []);
                % timestampsByTrain(iTrain, 1:numSamples) = (firstSampleIndex : firstSampleIndex + numSamples - 1) / sampleRate;
                % timestampsByTrain(iTrain, 1:numSamples) = readWindow(1):1/sampleRate:readWindow(2);
                timestampsByTrain(iTrain, 1:numSamples) = readWindow(1):1/sampleRate:readWindow(1) + (numSamples - 1) * 1/sampleRate;
            end

            data = reshape(permute(dataByTrain, [2, 1, 3]), size(dataByTrain, 1) * size(dataByTrain, 2), []);
            timestamps = reshape(permute(timestampsByTrain, [2, 1]), [], 1);

            hasData = ~isnan(timestamps);
            data = data(hasData, :);
            timestamps = timestamps(hasData);

            if (~CollisionTest.validateTimestamps(timestamps))
                error('Timestamps are not monotonic increasing. Digital events data probably needs to be trimmed.')
            end

            % Train to individual pulses. Extract data by window: stimOn + extendedWindow
            % Pre-allocate.
            % maxPulseLength = ceil(diff(extendedWindow) * sampleRate) + 1;
            % data = zeros(length(pulseOn), maxPulseLength, length(channels));
            
            % for iPulse = 1:length(pulseOn)
            %     inWindow = (allTimestamps > pulseOn(iPulse) + extendedWindow(1)) & (allTimestamps <= pulseOn(iPulse) + extendedWindow(2));
            %     data(iPulse, 1:nnz(inWindow), :) = reshape(allData(inWindow, :), 1, nnz(inWindow), []);
            % end

            fprintf('Done (%.2f) s.\n', toc(tTic));
        end

        function nsxChannels = mapChannels(obj, channels)
        %mapChannels - Convert TR channel labels (1:n continuous) to NSx file channels (1:N continuous). n: number of spike sorted channels. N: number of recorded channels including both rigs.

            % Find out which rig we're using. Rig2 data is appended after rig1.
            rig = TetrodeRecording.GetRig(obj.Filename.NSx);
            if (rig == 1)
                channelMap = obj.PTR.ChannelMap.Rig1;
            else
                channelMap = obj.PTR.ChannelMap.Rig2;
            end

            % Channel label conversion:
            channelMap = channelMap(obj.PTR.SelectedChannels);
            nsxChannels = channelMap(channels);
        end
    end

    % Static private methods
    methods (Access = {}, Static)

        function isMonoIncrease = validateTimestamps(array)
            isIncrease = diff(array) > 0;
            isMonoIncrease = nnz(isIncrease) == length(array) - 1;
        end

        function varargout = validateFolders(folders, varargin)
            p = inputParser();
            p.addRequired('Folders', @(x) iscell(x) || ischar(x));
            p.addParameter('SuppressWarnings', false, @islogical);
            p.parse(folders, varargin{:});
            folders = p.Results.Folders;
            suppressWarnings = p.Results.SuppressWarnings;

            allValid = true;
            validFolders = {};

            % Convert single folder in string format to cell.
            if ischar(folders)
                folders = {folders};
            end

            % Make sure all folders are valid, filter out bad ones.
            for iFolder = 1:length(folders)
                folder = folders{iFolder};

                % Check folder for NSx files
                if (isfolder(folder))
                    nsx = dir(sprintf('%s\\*.ns5', folder));
                    if length(nsx) > 1
                        nsx = nsx(1);
                        if ~suppressWarnings
                            warning('More than one .ns5 file is detected in "%s" but only the first one is used. This may lead to unexpected results.', folder)
                        end
                    elseif length(nsx) < 1
                        if ~suppressWarnings
                            warning('No .ns5 file detected in "%s".', folder)
                        end
                        allValid = false;
                        continue
                    end
                else
                    if ~suppressWarnings
                        warning('Selected path "%s" is not a folder.', folder)
                    end
                end

                % Determine experiment name
                expName = strsplit(nsx.name, '.ns5');
                expName = expName{1};

                % Find corresponding TR/PTR files in SpikeSort folder.
                spikeSortFolder = sprintf('%s\\..\\SpikeSort', nsx.folder);
                tr = dir(sprintf('%s\\tr_%s*.mat', spikeSortFolder, expName));
                ptr = dir(sprintf('%s\\ptr_%s.mat', spikeSortFolder, expName));

                % Make sure TR/PTR files exist
                if isempty(tr) || isempty(ptr)
                    if ~suppressWarnings
                        warning('Some of the following files could not be found: \n\t"tr_%s.mat"; \n\t"ptr_%s.mat".', expName)
                    end
                    allValid = false;
                    continue
                end

                if (length(ptr) > 1) 
                    ptr = ptr(1);
                    if ~suppressWarnings
                        warning('More than one (%d) ptr file was found but only the first one is used. This may lead to unexpected results.', length(ptr));
                    end
                end

                validFolders = [validFolders, folder];

                % Store file names
                Filename(iFolder).NSx = sprintf('%s\\%s', nsx.folder, nsx.name);
                Filename(iFolder).NEV = sprintf('%s\\%s.nev', nsx.folder, expName);
                Filename(iFolder).PTR = sprintf('%s\\%s', ptr.folder, ptr.name);
                for iTr = 1:length(tr)
                    Filename(iFolder).TR{iTr} = sprintf('%s\\%s', tr(iTr).folder, tr(iTr).name);
                end
            end

            varargout = {validFolders, Filename, expName, allValid};
        end
    end
end