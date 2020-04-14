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
        Spikes
        IsSorted = false
        SysInitDelay = 0
    end

    properties (Transient, Access = {})
        TR
        PTR
    end

    %% Constructor/Read/Write
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
            p.addParameter('TRSource', 'SortedOnly', @ischar); % 'SortedOnly', 'SortedOrRaw', 'RawOnly'
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
            [~, obj.Filename, obj.ExpName, allValid, obj.IsSorted] = CollisionTest.validateFolders(folder, 'SuppressWarnings', true, 'TRSource', p.Results.TRSource);
            if ~allValid
                error('Invalid folder %s', folder);
            end

            if (doRead)
                obj.read(p.Results.Channels, p.Results.MaxTrains, p.Results.ExtendedWindow ./ 1000);
            end
        end

        function save(obj, varargin)
        %save - Save the object to a .mat file.
        %
        % Syntax: save(obj, varargin)
        % - 'Filename': Designate specific save path.
        % 
            p = inputParser();
            p.addOptional('Filename', '', @ischar);
            p.addOptional('SeparateFiles', false, @islogical);
            p.parse(varargin{:});
            filename = p.Results.Filename;
            separateFiles = p.Results.SeparateFiles;

            % Only one object, or user chose to save all objects to one file.
            if (length(obj) == 1 || ~separateFiles)
                if (isempty(filename))
                    if (length(obj) == 1)
                        defaultName = sprintf('ct_%s.mat', obj.ExpName);
                    else
                        defaultName = 'ct.mat';
                    end
                    uisave('obj', defaultName);
                else
                    tTic = tic; fprintf('Writing to file "%s"...', filename)
                    save(filename, 'obj', '-v7.3');
                    fprintf('Done (%.1f s).\n', toc(tTic))
                end
            % More than one object, user chose to save them to separate files.
            else
                selPath = uigetdir();
                allObjs = obj;
                for iObj = 1:length(obj)
                    obj = allObjs(iObj);
                    filename = sprintf('%s\\ct_%s', selPath, obj.ExpName);
                    tTic = tic; fprintf('Writing to file "%s"...', filename)
                    save(filename, 'obj', '-v7.3');
                    fprintf('Done (%.1f s).\n', toc(tTic))
                end
                obj = allObjs;
            end

        end
    end

    % Public static methods
    methods (Static)
        function varargout = batch(varargin)
        %Batch - Batch generate and save CollisionTest objects from tr/ptr/ns5 files.
        %
        % Optional Parameters: 
        % 'ExtendedWindow'  - (ms) Extended window for reading data.
        % 
        % Syntax: Batch(varargin)
        %
            p = inputParser();
            p.addParameter('ExtendedWindow', [-20, 20], @(x) isnumeric(x) && length(x) == 2 && diff(x) > 0);
            p.addParameter('OnError', 'WarningLong', @ischar); % What to do when there is an error. Can be 'WarningShort', 'WarningLong', 'Error'
            p.addParameter('TRSource', 'SortedOnly', @ischar); % 'SortedOnly', 'SortedOrRaw', 'RawOnly'
            p.parse(varargin{:});
            onError = p.Results.OnError;

            folders = uipickfiles('Prompt', 'Select multiple folders each containing an ns5 file.');

            folders = CollisionTest.validateFolders(folders, 'SuppressWarnings', false, 'TRSource', p.Results.TRSource);

            for iFolder = 1:length(folders)
                try
                    folder = folders{iFolder};
                    fprintf('Processing folder %d of %d - %s:\n', iFolder, length(folders), folder);
    
                    ct = CollisionTest('Folder', folder, 'ExtendedWindow', p.Results.ExtendedWindow, 'TRSource', p.Results.TRSource);
    
                    % Make sure save path exists
                    saveFolder = sprintf('%s//..//CollisionTest', folder);
                    if ~isfolder(saveFolder)
                        mkdir(saveFolder);
                    end
                    if ct.IsSorted
                        ct.save(sprintf('%s//ct_sorted_%s.mat', saveFolder, ct.ExpName));
                    else
                        ct.save(sprintf('%s//ct_%s.mat', saveFolder, ct.ExpName));
                    end
                catch ME
                    if (strcmpi(onError, 'Error'))
                        error('Error when processing folder "%s".', folder);
                    end
                    warning('Error when processing folder "%s". This one will be skipped.', folder);
                    if (strcmpi(onError, 'WarningLong'))
                        warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
                    end
                end
            end

            varargout = {ct}; % Only the last ct is exported.
            TetrodeRecording.RandomWords();
        end

        function ct = load()
        %Load - Load multiple ct files.
        %
            files = uipickfiles('Prompt', 'Select (multiple) .mat files containing an CollisionTest object named "obj"', 'FilterSpec', '*.mat');

            for iFile = 1:length(files)
                S(iFile) = load(files{iFile}, 'obj');
            end

            ct = [S.obj];
            ct.removeEmptySpikes();
        end

        function tr = fixMisaligned(tr, sysInitDelay)
            if tr.FrequencyParameters.SysInitDelay.DataTrimmed
                disp('Does not need trimming. Already did.')
                return
            end

            if sysInitDelay > 0
                for iChn = 1:length(tr.Spikes)
                    if ~isempty(tr.Spikes(iChn).Channel)
                        tr.Spikes(iChn).Timestamps = tr.Spikes(iChn).Timestamps - sysInitDelay;
                    end
                end
                tr.FrequencyParameters.SysInitDelay.Duration = sysInitDelay;
                tr.FrequencyParameters.SysInitDelay.DataTrimmed = true;
                disp(['Removed data in the first ', num2str(sysInitDelay), ' seconds'])
            else
                disp('Does not need fixing.');
                return
            end
        end
    end

    % Private methods
    % methods (Access = {})
    methods
        function read(obj, channels, maxTrains, extendedWindow)
        %read - Read data from tr/ptr/ns5 files.
        % Syntax: read(obj, varargin)
            obj.readTR();

            [obj.PulseOn, obj.PulseOff, obj.TrainOn, obj.TrainOff] = obj.readDigitalEvents();
            [obj.Data, obj.Timestamps, obj.SampleRate, obj.SysInitDelay] = readAnalogTrains(obj, channels, maxTrains, extendedWindow, obj.TrainOn, obj.TrainOff);



            obj.TR = CollisionTest.fixMisaligned(obj.TR, obj.SysInitDelay);
            obj.readSpikes();
            obj.removeEmptySpikes();
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

        function [data, timestamps, sampleRate, sysInitDelay] = readAnalogTrains(obj, channels, maxTrains, extendedWindow, trainOn, trainOff)
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
            nsxChannels =  obj.mapChannels(channels, 'From', 'TR', 'To', 'NSx');
            maxTrainLength = max(trainOff - trainOn) + diff(extendedWindow);
            maxTrainLength = ceil(sampleRate * maxTrainLength);

            dataByTrain = NaN(length(trainOn), maxTrainLength, length(channels));
            timestampsByTrain = NaN(length(trainOn), maxTrainLength);

            % Read data from each train and put them in an array.

            % Get sysInitDelay to shift openNSx window to the correct value.
            
            % [Deprecated] read SysInitDelay from TetrodeRecording object - we can't trust sysInitDelay from TetrodeRecording for some reason. re-read the first 20 seconds of the .NSx file to get it.
            % sysInitDelay = obj.TR.FrequencyParameters.SysInitDelay.Duration; 
            % if isnan(sysInitDelay)
            %     sysInitDelay = 0;
            % end

            NSx = openNSx(obj.Filename.NSx, 'read', 'channels', 1, 'duration', [0, 20], 'sec');
            if iscell(NSx.Data)
                if length(NSx.Data) == 2
                    sysInitDelay = size(NSx.Data{1}, 2) / sampleRate;
                else
                    error("NSx.Data is a cell array of length %d. This is not supported. It should be 2.", length(NSx.Data))
                end
            else
                sysInitDelay = 0;
            end
            
            tTic = tic();
            fprintf('Reading %d trains, %d channels...', length(trainOn), length(channels));

            for iTrain = 1:length(trainOn)
                trainWindow = [trainOn(iTrain), trainOff(iTrain)] + extendedWindow;

                % Read window for NSx file. Must shift to the right by sysInitDelay to exclude discarded data when the other rig started.
                readWindow = trainWindow + sysInitDelay;

                NSx = openNSx(obj.Filename.NSx, 'read', 'channels', nsxChannels, 'duration', readWindow, 'sec');
                numSamples = size(NSx.Data, 2);

                dataByTrain(iTrain, 1:numSamples, :) = reshape(transpose(NSx.Data), 1, numSamples, []);
                % timestampsByTrain(iTrain, 1:numSamples) = (firstSampleIndex : firstSampleIndex + numSamples - 1) / sampleRate;
                % timestampsByTrain(iTrain, 1:numSamples) = readWindow(1):1/sampleRate:readWindow(2);
                timestampsByTrain(iTrain, 1:numSamples) = trainWindow(1) : (1 / sampleRate) : (trainWindow(1) + (numSamples - 1) * 1/sampleRate);
            end

            data = reshape(permute(dataByTrain, [2, 1, 3]), size(dataByTrain, 1) * size(dataByTrain, 2), []);
            timestamps = reshape(permute(timestampsByTrain, [2, 1]), [], 1);

            hasData = ~isnan(timestamps);
            data = data(hasData, :);
            timestamps = timestamps(hasData);

            if (~CollisionTest.validateTimestamps(timestamps))
                error('Timestamps are not monotonic increasing. Digital events data probably needs to be trimmed.')
            end

            fprintf('Done (%.2f) s.\n', toc(tTic));
        end

        function readTR(obj, varargin)
        %read - Load TR/PTR files.
        % Syntax: obj.readTR('Sorted', false)
            p = inputParser();
            p.addParameter('Raw', true, @islogical);
            p.parse(varargin{:});
            readRaw = p.Results.Raw;

            obj.PTR = TetrodeRecording.BatchLoad({obj.Filename.PTR});
            obj.TR = TetrodeRecording.BatchLoad(obj.Filename.TR);
        end

        function readSpikes(obj)
        %read - Load and process sorted spikeTimes/waveforms from TR.Spikes.
        % Syntax: obj.readSpikes()
            for iChn = [obj.TR.Spikes.Channel]
                obj.Spikes(iChn).Channel = iChn;
                obj.Spikes(iChn).RawChannel = obj.mapChannels(iChn);
                obj.Spikes(iChn).WaveformWindow = obj.TR.Spikes(iChn).WaveformWindow;
                obj.Spikes(iChn).WaveformTimestamps = obj.TR.Spikes(iChn).WaveformTimestamps;
                for iUnit = unique(obj.TR.Spikes(iChn).Cluster.Classes)
                    theseIndices = obj.TR.Spikes(iChn).Cluster.Classes == iUnit;
                    theseWaveforms = obj.TR.Spikes(iChn).Waveforms(theseIndices, :);
                    obj.Spikes(iChn).Units(iUnit).Timestamps = obj.TR.Spikes(iChn).Timestamps(theseIndices);
                    obj.Spikes(iChn).Units(iUnit).Waveform.Mean = mean(theseWaveforms, 1);
                    obj.Spikes(iChn).Units(iUnit).Waveform.STD = std(theseWaveforms, 0, 1);
                    obj.Spikes(iChn).Units(iUnit).Waveform.Percentile95 = prctile(theseWaveforms, 95, 1);
                    obj.Spikes(iChn).Units(iUnit).Waveform.Percentile05 = prctile(theseWaveforms, 5, 1);
                end
            end
        end

        function removeEmptySpikes(obj)
            for iObj = 1:length(obj)
                emptyIndices = cellfun(@isempty, {obj(iObj).Spikes.Channel});
                if nnz(emptyIndices) > 0
                    obj(iObj).Spikes(emptyIndices) = [];
                    fprintf(1, 'Removed %d empty channels.\n', nnz(emptyIndices));
                end
            end
        end

        function outChannels = mapChannels(obj, channels, varargin)
        %mapChannels - Convert TR channel labels (1:n continuous) to NSx file channels (1:N continuous). n: number of spike sorted channels. N: number of recorded channels including both rigs.
            p = inputParser();
            p.addRequired('Channels', @isnumeric);
            p.addParameter('From', 'TR', @ischar);
            p.addParameter("To", "NSx", @ischar);
            p.parse(channels, varargin{:});
            channels = p.Results.Channels;
            from = p.Results.From;
            to = p.Results.To;

            switch lower(from)
                case 'tr'

                otherwise
                    warning('%s is not an acceptable value for argument "From". Only "TR" is supported.', from);
            end

            switch lower(to)
                case 'nsx'
                    % Find out which rig we're using. Rig2 data is appended after rig1.
                    rig = TetrodeRecording.GetRig(obj.Filename.NSx);
                    if (rig == 1)
                        channelMap = obj.PTR.ChannelMap.Rig1;
                    else
                        channelMap = obj.PTR.ChannelMap.Rig2;
                    end
        
                    % Channel label conversion:
                    channelMap = channelMap(obj.PTR.SelectedChannels);
                    outChannels = channelMap(channels);
                    return

                case 'data'
                    [~, outChannels] = ismember(channels, [obj.Spikes.Channel]);
                    outChannels(isnan(outChannels)) = NaN;
                    return

                otherwise
                    warning('%s is not an acceptable value for argument "To". Only "NSx" or "Data" is supported.', to);
            end
        end
    end

    % Private static methods
    methods (Access = {}, Static)

        function isMonoIncrease = validateTimestamps(array)
            isIncrease = diff(array) > 0;
            isMonoIncrease = nnz(isIncrease) == length(array) - 1;
        end

        function varargout = validateFolders(folders, varargin)
            p = inputParser();
            p.addRequired('Folders', @(x) iscell(x) || ischar(x));
            p.addParameter('SuppressWarnings', false, @islogical);
            p.addParameter('Sorted', false, @islogical);
            p.addParameter('TRSource', 'SortedOnly', @ischar); % 'SortedOnly', 'SortedOrRaw', 'RawOnly'
            p.parse(folders, varargin{:});
            folders = p.Results.Folders;
            suppressWarnings = p.Results.SuppressWarnings;
            trSource = p.Results.TRSource;

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
                switch lower(trSource)
                    case 'sortedonly'
                        tr = dir(sprintf('%s\\tr_sorted_%s*.mat', spikeSortFolder, expName));
                        isSorted = true;
                    case 'rawonly'
                        tr = dir(sprintf('%s\\tr_%s*.mat', spikeSortFolder, expName));
                        isSorted = false;
                    case 'sortedorraw'
                        tr = dir(sprintf('%s\\tr_sorted_%s*.mat', spikeSortFolder, expName));
                        isSorted = true;
                        if isempty(tr)
                            tr = dir(sprintf('%s\\tr_%s*.mat', spikeSortFolder, expName));
                            isSorted = false;
                        end
                end
                ptr = dir(sprintf('%s\\ptr_%s.mat', spikeSortFolder, expName));

                % Make sure TR/PTR files exist
                if isempty(tr) || isempty(ptr)
                    if ~suppressWarnings
                        if isSorted
                            trFilePattern = sprintf('tr_sorted_%s*.mat', expName);
                        else
                            trFilePattern = sprintf('tr_%s*.mat', expName);
                        end
                        warning('Some of the following files could not be found: \n\t%s; \n\t"ptr_%s.mat".', trFilePattern, expName)
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

            varargout = {validFolders, Filename, expName, allValid, isSorted};
        end
    end

    %% Plotting
    methods
        function plot(obj, channel, varargin)
        %plot - Plot electrode data aligned on stimulation onset events.
        % Syntax: plot(obj, channel, varargin)
        %   Press Down arrow key/Mouse wheel down to go to next page.
        %   Press Up arrow key/Mouse wheel up to go to previous page.
        %   Shift/Alt + Up/Down arrow keys to speed up.

            if isempty(channel)
                channel = obj.Spikes(1).Channel;
            end

            p = inputParser();
            p.addRequired('Channel', @(x) isnumeric(x));
            p.addParameter('Start', 1, @isnumeric);
            p.addParameter('TracesPerPage', 25, @isnumeric);
            p.addParameter('YLim', [-500, 500], @(x) isnumeric(x) && length(x) == 2 && diff(x) > 0);
            p.addParameter('XLim', obj.Window, @(x) isnumeric(x) && length(x) == 2 && diff(x) > 0);
            p.addParameter('YSpacing', 1, @(x) isnumeric(x) && length(x) == 1);
            p.parse(channel, varargin{:})
            xRange = p.Results.XLim;
            
            % Create figure
            fig = figure(); 
            ax = axes(fig);
            % ax.YDir = 'reverse';
            grid(ax, 'on');
            hold(ax, 'on');
            xlim(ax, xRange);

            xlabel(ax, 'Time from stim on (ms)')
            ylabel(ax, 'Trial number + Normalized voltage (a.u.)')

            % Plot the first page
            obj.updatePlot(ax, p, p.Results.Start);

            % Keypress listeners for updating the plot.
            fig.WindowKeyPressFcn       = {@obj.onWindowKeyPress, p};
            fig.WindowScrollWheelFcn    = {@obj.onWindowKeyPress, p};

            % Plot unit mean waveforms
            fig2 = figure();
            ax2 = axes(fig2);
            xlabel(ax2, 'Time (ms)')
            ylabel(ax2, 'Voltage (mV)')
            hold(ax2, 'on')

            channel = obj.mapChannels(channel, 'From', 'TR', 'To', 'Data');

            t = obj.Spikes(channel).WaveformTimestamps;

            colors = 'rgbcmyk';

            for iUnit = 1:length(obj.Spikes(channel).Units)
                thisColor = colors(iUnit);
                Waveform = obj.Spikes(channel).Units(iUnit).Waveform;
                line(ax2, t, Waveform.Mean, 'LineStyle', '-', 'Color', thisColor);
                patch(ax2, [t, flip(t)], [Waveform.Percentile05, flip(Waveform.Percentile95)], thisColor,...
                    'FaceAlpha', 0.15, 'EdgeColor', 'none');
                patch(ax2, [t, flip(t)], [Waveform.Mean - Waveform.STD, flip(Waveform.Mean + Waveform.STD)], thisColor,...
                    'FaceAlpha', 0.4, 'EdgeColor', 'none');
                line(ax2, t, [Waveform.Mean - Waveform.STD; Waveform.Mean + Waveform.STD], 'LineStyle', '--', 'Color', thisColor);
                line(ax2, t, [Waveform.Percentile05; Waveform.Percentile95], 'LineStyle', ':', 'Color', thisColor);
            end
        end
    end

    methods (Access = {})
        function onWindowKeyPress(obj, fig, event, p)
            if isa(event, 'matlab.ui.eventdata.ScrollWheelData')
                if (event.VerticalScrollCount < 0)
                    direction = 'up';
                else
                    direction = 'down';
                end
                turbo = 1;
            elseif isa(event, 'matlab.ui.eventdata.KeyData')
                if strcmp(event.Key, 'uparrow')
                    direction = 'up';
                elseif strcmp(event.Key, 'downarrow')
                    direction = 'down';
                else
                    return
                end
                shift = ismember('shift', event.Modifier);
                alt = ismember('alt', event.Modifier);

                if shift && alt
                    turbo = Inf;
                elseif alt
                    turbo = 10;
                elseif shift
                    turbo = 5;
                else
                    turbo = 1;
                end
            else
                return
            end

            ax = fig.findobj('Type', 'axes');

            if strcmp(direction, 'down')
                startTrace = min(ax.UserData.StartTrace + turbo * p.Results.TracesPerPage, length(obj.PulseOn));
            elseif strcmpi(direction, 'up')
                startTrace = max(ax.UserData.StartTrace - turbo * p.Results.TracesPerPage, 1);
            else
                return
            end

            obj.updatePlot(ax, p, startTrace);
        end

        function updatePlot(obj, ax, p, startTrace)
            trChannel = p.Results.Channel;
            channel = obj.mapChannels(trChannel, 'From', 'TR', 'To', 'Data');

            tracesPerPage = p.Results.TracesPerPage;
            yRange = p.Results.YLim;
            ySpacing = p.Results.YSpacing;

            % Store startTrace
            ax.UserData.StartTrace = startTrace;

            cla(ax);

            plotWindow = 0.001 * obj.Window;

            % Normalize voltage data and align to stimOnsetTime;
            for iPulse = startTrace:startTrace + tracesPerPage

                iPulseInPage = iPulse - startTrace + 1;

                isInPlotWindow = obj.Timestamps > obj.PulseOn(iPulse) + plotWindow(1) & obj.Timestamps <= obj.PulseOn(iPulse) + plotWindow(2);
                pulseData = obj.Data(isInPlotWindow, channel);
                pulseTimestamps = obj.Timestamps(isInPlotWindow);

                % Normalize voltage to yRange.
                y = (pulseData - yRange(1)) ./ diff(yRange) + iPulse * ySpacing - 0.5;
                
                % Align time to stimOn
                t = 1000 * (pulseTimestamps - obj.PulseOn(iPulse));

                % Plot trace
                plot(ax, t, y, 'k');

                % Plot stim window
                stimOnVertices(2 * iPulseInPage - 1: 2 * iPulseInPage, 1) = 0;
                stimOnVertices(2 * iPulseInPage - 1: 2 * iPulseInPage, 2) = [iPulse * ySpacing, iPulse * ySpacing + 1];
                stimOffVertices(2 * iPulseInPage - 1: 2 * iPulseInPage, 1) = 1000 * (obj.PulseOff(iPulse) - obj.PulseOn(iPulse));
                stimOffVertices(2 * iPulseInPage - 1: 2 * iPulseInPage, 2) = [iPulse * ySpacing, iPulse * ySpacing + 1];

                colors = 'rgbcmyk';

                % Plot sorted spike timestamps
                numUnits = max(1, length(obj.Spikes(channel).Units));
                for iUnit = 1:numUnits
                    unitTimestamps = obj.Spikes(channel).Units(iUnit).Timestamps;
                    isInPlotWindow = unitTimestamps > obj.PulseOn(iPulse) + plotWindow(1) & unitTimestamps <= obj.PulseOn(iPulse) + plotWindow(2);
                    t = 1000 * (unitTimestamps(isInPlotWindow) - obj.PulseOn(iPulse));
                    y = repmat(iPulse * ySpacing, [nnz(isInPlotWindow), 1]);
                    plot(t, y, sprintf('%so', colors(iUnit)));
                end

                if iPulse >= length(obj.PulseOn)
                    break
                end
            end

            % Page done
            stimPatchVertices = vertcat(stimOnVertices, stimOffVertices(end:-1:1, :));
            patch('XData', stimPatchVertices(:, 1), 'YData', stimPatchVertices(:, 2), 'FaceColor', '#4DBEEE', 'FaceAlpha', 0.33, 'EdgeAlpha', 0);
            ylim([(iPulse - iPulseInPage + 1) * ySpacing - 1, iPulse * ySpacing])
            yticks(ax, ax.YLim(1) + 1:5:ax.YLim(2)) 
            title(ax, sprintf('%s Chn%d (Pulses %d - %d)', obj.ExpName, trChannel, iPulse - iPulseInPage + 1, iPulse), 'Interpreter', 'none')
        end

        function plotUnitWaveform(obj, ax, channel)
            
        end
    end
end