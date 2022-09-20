classdef EphysUnit < handle
    properties
        ExpName = ''
        Channel = NaN
        Electrode = NaN
        Unit = NaN
        SpikeTimes = []
        Waveforms = []
        WaveformTimestamps = []
        ExtendedWindow = []
        EventTimes = struct('Cue', [], 'Press', [], 'Lick', [], 'StimOn', [], 'StimOff', [])
        Trials = struct('Press', [], 'Lick', [], 'StimTrain', [], 'Stim', [])
        SpikeCounts = uint16([])
        SpikeCountTimestamps = single([])
        SpikeRates = single([])
        SpikeRateTimestamps = single([])
        SpikeRateKernel = struct('type', [], 'params', [], 't', [], 'y', [])
        SpikeCountStats = struct('median', [], 'mad', [], 'medianITI', [], 'madITI', [], 'resolution', [])
        SpikeRateStats = struct('median', [], 'mad', [], 'medianITI', [], 'madITI', [], 'resolution', [])
    end
    
    % public methods
    methods
        % constructor
        function obj = EphysUnit(varargin)
            if nargin == 0
                return
            end

            p = inputParser();
            if isstruct(varargin{1})
                p.addRequired('PETH', @isstruct);
            else
                p.addRequired('AR', @(x) isa(x, 'AcuteRecording'));
            end
            p.addParameter('readWaveforms', false, @islogical);
            p.addParameter('cullITI', false, @islogical);
            p.addParameter('extendedWindow', [-1, 1], @(x) isnumeric(x) && length(x) == 2)
            p.parse(varargin{:});
            extendedWindow = p.Results.extendedWindow;
            
            usePETH = isfield(p.Results, 'PETH');

            % Initialize obj array
            if usePETH
                PETH = p.Results.PETH;
                obj(length(PETH)) = EphysUnit();
                uniqueExpNames = unique({PETH.ExpName});

                i = 1;
                fprintf(1, 'Reading PETH...\n')
                for iExp = 1:length(uniqueExpNames)
                    fprintf(1, 'Loading experiment %s (%i of %i)...\n', uniqueExpNames{iExp}, iExp, length(uniqueExpNames))
                    tr = TetrodeRecording.BatchLoadSimple(uniqueExpNames{iExp});
                    inExp = strcmpi({PETH.ExpName}, uniqueExpNames{iExp});
                    tTic = tic();
                    fprintf(1, 'Processing %i units...\n', nnz(inExp))
                    for e = PETH(inExp)
                         try
                            obj(i).ExpName = uniqueExpNames{iExp};
                            obj(i).Channel = e.Channel;
                            if isempty(tr.SelectedChannels)
                                obj(i).Electrode = NaN;
                            else
                                obj(i).Electrode = tr.SelectedChannels(obj(i).Channel);
                            end
                            obj(i).Unit = e.Cluster;
                            s = tr.Spikes(obj(i).Channel);
                            inCluster = s.Cluster.Classes == e.Cluster;
                            obj(i).SpikeTimes = s.Timestamps(inCluster);
                            if p.Results.readWaveforms
                                obj(i).Waveforms = int16(s.Waveforms(inCluster, :));
                                obj(i).WaveformTimestamps = s.WaveformTimestamps;
                            end
                            obj(i).EventTimes.Cue = tr.DigitalEvents.CueOn;
                            obj(i).EventTimes.Press = tr.DigitalEvents.PressOn;
                            obj(i).EventTimes.Lick = tr.DigitalEvents.LickOn;
                            obj(i).EventTimes.RewardTimes = tr.DigitalEvents.RewardOn;
                            if isfield(tr.DigitalEvents, 'StimOn')
                                obj(i).EventTimes.StimOn = tr.DigitalEvents.StimOn;
                                obj(i).EventTimes.StimOff = tr.DigitalEvents.StimOff;
                            else
                                obj(i).EventTimes.StimOn = [];
                                obj(i).EventTimes.StimOff = [];
                            end                                
    
                            if p.Results.cullITI
                                obj(i).ExtendedWindow = extendedWindow;
                            else
                                obj(i).ExtendedWindow = [];
                            end
                            
                            % Make trials
                            obj(i).Trials.Press = obj(i).makeTrials('press');
                            obj(i).Trials.Lick = obj(i).makeTrials('lick');
                            obj(i).Trials.Stim = obj(i).makeTrials('stim');
                            obj(i).Trials.StimTrain = obj(i).makeTrials('stimtrain');
    
                            % Calculate spike rates and spike counts
                            resolution_sc = 0.1;
                            resolution_sr = 1e-3;
                            
                            [sc, tsc] = obj(i).getSpikeCounts(resolution_sc);
                            [sr, tsr, kernel] = obj(i).getSpikeRates('gaussian', 0.1, resolution_sr, 'kernelWidth', 1);
                            
                            % Calculate 
                            [~, ~, scInTrial] = obj(i).cullITIData(tsc, sc, 'all', 'extendedWindow', extendedWindow);
                            [~, ~, srInTrial] = obj(i).cullITIData(tsr, sr, 'all', 'extendedWindow', extendedWindow);
                            
                            obj(i).SpikeCountStats = struct('median', median(double(sc)), 'mad', mad(double(sc), 1), 'medianITI', median(double(sc(~scInTrial))), 'madITI', mad(double(sc(~scInTrial)), 1), 'resolution', resolution_sc);
                            obj(i).SpikeRateStats = struct('median', median(sr), 'mad', mad(sr, 1), 'medianITI', median(sr(~srInTrial)), 'madITI', mad(sr(~srInTrial), 1), 'resolution', resolution_sr);
    
                            % Cull ITI spikes if told to do so.
                            if p.Results.cullITI
                                [obj(i).SpikeTimes, obj(i).Waveforms] = obj(i).cullITIData(obj(i).SpikeTimes, obj(i).Waveforms, 'all', 'extendedWindow', extendedWindow);
                                obj(i).SpikeCounts = sc(scInTrial);
                                obj(i).SpikeCountTimestamps = tsc(scInTrial);
                                obj(i).SpikeRates = sr(srInTrial);
                                obj(i).SpikeRateTimestamps = tsr(srInTrial);
                            else
                                obj(i).SpikeCounts = sc;
                                obj(i).SpikeCountTimestamps = tsc;
                                obj(i).SpikeRates = sr;
                                obj(i).SpikeRateTimestamps = tsr;
                            end
                            obj(i).SpikeRateKernel = kernel;
    
                            obj(i).save();
                            
                            i = i + 1;
                        catch ME
                            warning('Error when processing unit %i (%s) - this one will be skipped.', i, e.ExpName)
                            warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
                            i = i + 1;
                        end			
                    end
                    fprintf(1, 'Done (%.1f sec).\n', toc(tTic))
                    clear tr
                end
                fprintf(1, 'Done.\n')
            else
                ar = p.Results.AR;
                BSR = ar.bsr;
                obj(length(BSR)) = EphysUnit();
                expName = ar.expName;
                fprintf(1, 'Loading experiment %s...\n', expName)
                tr = TetrodeRecording.BatchLoadSimple(expName, true);
                tTic = tic();
                fprintf(1, 'Processing %i units...\n', length(BSR))
                i = 1;
                for bsr = BSR
                     try
                        obj(i).ExpName = expName;
                        obj(i).Channel = bsr.channel;
                        if isempty(tr.SelectedChannels)
                            obj(i).Electrode = NaN;
                        else
                            obj(i).Electrode = tr.SelectedChannels(obj(i).Channel);
                        end
                        obj(i).Unit = bsr.unit;
                        s = tr.Spikes(obj(i).Channel);
                        inCluster = s.Cluster.Classes == bsr.unit;
                        obj(i).SpikeTimes = s.Timestamps(inCluster);
                        if p.Results.readWaveforms
                            obj(i).Waveforms = int16(s.Waveforms(inCluster, :));
                            obj(i).WaveformTimestamps = s.WaveformTimestamps;
                        end
                        obj(i).EventTimes.Cue = tr.DigitalEvents.CueOn;
                        obj(i).EventTimes.Press = tr.DigitalEvents.PressOn;
                        obj(i).EventTimes.Lick = tr.DigitalEvents.LickOn;
                        obj(i).EventTimes.RewardTimes = tr.DigitalEvents.RewardOn;
                        if isfield(tr.DigitalEvents, 'StimOn')
                            obj(i).EventTimes.StimOn = tr.DigitalEvents.StimOn;
                            obj(i).EventTimes.StimOff = tr.DigitalEvents.StimOff;
                        else
                            obj(i).EventTimes.StimOn = [];
                            obj(i).EventTimes.StimOff = [];
                        end                                

                        if p.Results.cullITI
                            obj(i).ExtendedWindow = extendedWindow;
                        else
                            obj(i).ExtendedWindow = [];
                        end
                        
                        % Make trials
                        obj(i).Trials.Press = obj(i).makeTrials('press');
                        obj(i).Trials.Lick = obj(i).makeTrials('lick');
                        obj(i).Trials.Stim = obj(i).makeTrials('stim');
                        obj(i).Trials.StimTrain = obj(i).makeTrials('stimtrain');

                        % Calculate spike rates and spike counts
                        resolution_sc = 0.1;
                        resolution_sr = 1e-3;
                        
                        [sc, tsc] = obj(i).getSpikeCounts(resolution_sc);
                        [sr, tsr, kernel] = obj(i).getSpikeRates('gaussian', 0.1, resolution_sr, 'kernelWidth', 1);
                        
                        % Calculate 
                        [~, ~, scInTrial] = obj(i).cullITIData(tsc, sc, 'all', 'extendedWindow', extendedWindow);
                        [~, ~, srInTrial] = obj(i).cullITIData(tsr, sr, 'all', 'extendedWindow', extendedWindow);
                        
                        obj(i).SpikeCountStats = struct('median', median(double(sc)), 'mad', mad(double(sc), 1), 'medianITI', median(double(sc(~scInTrial))), 'madITI', mad(double(sc(~scInTrial)), 1), 'resolution', resolution_sc);
                        obj(i).SpikeRateStats = struct('median', median(sr), 'mad', mad(sr, 1), 'medianITI', median(sr(~srInTrial)), 'madITI', mad(sr(~srInTrial), 1), 'resolution', resolution_sr);

                        % Cull ITI spikes if told to do so.
                        if p.Results.cullITI
                            [obj(i).SpikeTimes, obj(i).Waveforms] = obj(i).cullITIData(obj(i).SpikeTimes, obj(i).Waveforms, 'all', 'extendedWindow', extendedWindow);
                            obj(i).SpikeCounts = sc(scInTrial);
                            obj(i).SpikeCountTimestamps = tsc(scInTrial);
                            obj(i).SpikeRates = sr(srInTrial);
                            obj(i).SpikeRateTimestamps = tsr(srInTrial);
                        else
                            obj(i).SpikeCounts = sc;
                            obj(i).SpikeCountTimestamps = tsc;
                            obj(i).SpikeRates = sr;
                            obj(i).SpikeRateTimestamps = tsr;
                        end
                        obj(i).SpikeRateKernel = kernel;

                        obj(i).save();
                        
                        i = i + 1;
                    catch ME
                        warning('Error when processing unit %i (%s) - this one will be skipped.', i, expName)
                        warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
                        i = i + 1;
                    end			
                end
                fprintf(1, 'Done (%.1f sec).\n', toc(tTic))
            end

        end
        
        function save(obj, varargin)
            p = inputParser();
            addOptional(p, 'path', 'C:\SERVER\Units\', @ischar);
            parse(p, varargin{:});
            path = p.Results.path;
            
            if path(end) ~= '\'
                path = [path, '\'];
            end
            
            for i = 1:length(obj)
                eu = obj(i);
                if ~isfolder(path)
                    mkdir(path)
                end
                file = sprintf('%s%s.mat', path, eu.getName('_'));
                
                tTic = tic();
                fprintf(1, 'Saving EpysUnit to file %s...', file)
                save(file, 'eu', '-v7.3')
                fprintf(1, 'Done! (%.2f s)\n', toc(tTic));
            end
        end
        
        function name = getAnimalName(obj)
            if length(obj) == 1
                if isempty(obj.ExpName)
                    name = '';
                    return
                else
                    name = strsplit(obj.ExpName, '_');
                    name = name{1};
                    return
                end
            else
                name = arrayfun(@getAnimalName, obj, 'UniformOutput', false);
                return
            end
        end

		function name = getName(obj, varargin)
            p = inputParser();
            addOptional(p, 'sep', ' ', @ischar);
            parse(p, varargin{:});
            sep = p.Results.sep;
            
            if length(obj) == 1
                if isnan(obj.Electrode)
                    name = sprintf('%s%sChannel%i%sUnit%i', obj.ExpName, sep, obj.Channel, sep, obj.Unit);
                else
                    name = sprintf('%s%sElectrode%i%sUnit%i', obj.ExpName, sep, obj.Electrode, sep, obj.Unit);
                end
           else
                name = arrayfun(@getName, obj, 'UniformOutput', false);
            end
        end
        
        function trials = getTrials(obj, trialType, varargin)
            p = inputParser();
            p.addRequired('trialType', @(x) all(ismember(x, {'press', 'lick', 'stim', 'stimtrain', 'stimfirstpulse'})));
            p.addOptional('sorted', true, @islogical);
            p.parse(trialType, varargin{:});
            trialType = p.Results.trialType;
            sorted = p.Results.sorted;

            if length(obj) == 1
                if ~iscell(trialType)
                    trialType = {trialType};
                end
                trials = cell(1, length(trialType));
                for itt = 1:length(trialType)
                    switch lower(trialType{itt})
                        case 'press'
                            trials{itt} = obj.Trials.Press(:);
                        case 'lick'
                            trials{itt} = obj.Trials.Lick(:);
                        case 'stim'
                            trials{itt} = obj.Trials.Stim(:);
                        case 'stimtrain'
                            trials{itt} = obj.Trials.StimTrain(:);
                        case 'stimfirstpulse'
                            trains = obj.Trials.StimTrain(:);
                            pulses = obj.Trials.Stim(:);
                            trials{itt} = Trial([trains.Start], [pulses.Stop], 'first');
                            trials{itt} = trials{itt}(:);
                    end
                end
                trials = cat(1, trials{:});
                if sorted
                    trials = trials.sortby('start', 'ascend');
                end
            else
                trials = cell(length(obj), 1);
                for i = 1:length(obj)
                    trials{i} = obj(i).getTrials(trialType, sorted);
                end
            end
        end
        
        function [X, T, N, S, B] = getBinnedTrialAverage(obj, data, varargin)
            % Group trials by length, then calculate average spike rates.
            %   X = zeros(nBins, nSamples); % nxt matrix, each row = mean spike rate across trials for that bin
            %   S = zeros(nBins, nSamples); % nxt matrix, each row = variance for that bin
            %   T = zeros(1, nSamples); % 1xt matrix = shared timestamps for all bins
            %   N = zeros(nBins, 1); % nx1 matrix, number of trials per bin
            %   B = zeros(nBins, 2); % nx2 matrix, each row representing bin edges
            p = inputParser();
            p.addRequired('data', @(x) ischar(x) && ismember(lower(x), {'rate', 'count'}))
            p.addOptional('edges', 0:2.5:10, @(x) isnumeric(x) && nnz(diff(x) <= 0) == 0)
            p.addOptional('trialType', 'press', @(x) ischar(x) && ismember(lower(x), {'press', 'lick', 'stim'}))
            p.addParameter('alignTo', 'stop', @(x) ischar(x) && ismember(lower(x), {'start', 'stop'}))
            p.addParameter('window', [], @(x) isnumeric(x) && ismember(length(x), [0, 2]))
            p.addParameter('resolution', 0.001, @(x) isnumeric(x) && x > 0)
            p.addParameter('normalize', false, @islogical)
            p.parse(data, varargin{:})
            data = lower(p.Results.data);
            edges = p.Results.edges;
            trialType = lower(p.Results.trialType);
            alignTo = lower(p.Results.alignTo);
            resolution = p.Results.resolution;
            
            if isempty(p.Results.window)
                window = [-edges(end), 1];
            else
                window = p.Results.window;
            end
            
            % Preallocate outputs
            nBins = length(edges) - 1;
            nSamples = length(window(1):resolution:window(2)) - 1;
            X = zeros(nBins, nSamples); % nxt matrix, each row = mean spike rate across trials for that bin
            S = zeros(nBins, nSamples); % nxt matrix, each row = variance for that bin
            T = zeros(1, nSamples); % 1xt matrix = shared timestamps for all bins
            N = zeros(nBins, 1); % nx1 matrix, number of trials per bin
            B = zeros(nBins, 2); % nx2 matrix, each row representing bin edges
            for iBin = 1:nBins
                B(iBin, :) = [edges(iBin), edges(iBin + 1)];
            end
            
            tTic = tic();
            for i = 1:length(obj)
                fprintf(1, 'Binning %i/%i...\n', i, length(obj))
                for iBin = 1:nBins
                    switch data
                        case 'rate'
                            stats = obj(i).SpikeRateStats;
                        case 'count'
                            stats = obj(i).SpikeCountStats;
                    end
                    
                    [xx, tt] = obj(i).getTrialAlignedData(data, window, trialType, allowedTrialDuration=[edges(iBin), edges(iBin+1)], alignTo=alignTo, resolution=resolution, includeITI=false);

                    if i == 1 && iBin == 1
                        T = tt;
                    end
                    
                    if isempty(xx)
                        continue
                    end
                    
                    % Use modified z-score if requested
                    if p.Results.normalize
                        xx = EphysUnit.normalize(xx, 'iti', stats);
                    end
                    
                    mu = mean(xx, 1, 'omitnan');
                    ss = sum((xx - mu).^2, 1, 'omitnan');
                    n = size(xx, 1);
                    
                    [mu, ss, n] = EphysUnit.combinestats(X(iBin, :), S(iBin, :), N(iBin), mu, ss, n);
                    
                    X(iBin, :) = mu;
                    S(iBin, :) = ss;
                    N(iBin) = n;
                end
            end
            
            % Convert sum of squares to (unbiased) standard deviation
            S = (S ./ (N - 1)).^0.5;
            
            fprintf(1, 'Done (%.1f sec)!\n', toc(tTic))
        end
        
        function eta = getETA(obj, data, event, varargin)
            %GETETA estimate mean spikerate around an event
            %  X - Nxt, Averaged spike rate (raw or normalized)
            %  t - tx1, common timestamps.
            %  N - Nx1, number of trials per neuron
            %  stats - Nx1 struct('mean', 'sd'), mean spike rate and sd for each neuron
            p = inputParser();
            p.addRequired('data', @(x) ischar(x) && ismember(lower(x), {'rate', 'count'}))
            p.addRequired('event', @(x) ischar(x) && ismember(lower(x), {'press', 'lick', 'stim'}))
            p.addOptional('window', [-2, 0], @(x) isnumeric(x) && length(x)>=2 && x(2) > x(1))
            p.addParameter('minTrialDuration', -Inf, @(x) isnumeric(x) && length(x)==1 && x>=0)
            p.addParameter('resolution', [], @(x) isnumeric(x) && length(x)==1 && x>=0)
            p.addParameter('normalize', 'none', @(x) isstruct(x) || isnumeric(x) || ismember(lower(x), {'none', 'iti', 'all'}))
            p.addParameter('alignTo', 'default', @(x) ismember(x, {'default', 'start', 'stop'}))
            p.addParameter('includeITI', [], @islogical)
            p.parse(data, event, varargin{:})
            data = lower(p.Results.data);
            event = p.Results.event;
            window = p.Results.window;
            minTrialDuration = p.Results.minTrialDuration;
            resolution = p.Results.resolution;
            normalize = p.Results.normalize;
            alignTo = p.Results.alignTo;
            includeITI = p.Results.includeITI;
            
            % Use default resolutions
            if isempty(resolution)
                switch data
                    case 'count'
                        resolution = obj(1).SpikeCountStats.resolution;
                    case 'rate'
                        resolution = obj(1).SpikeRateStats.resolution;
                end
            end
                
            edges = window(1):resolution:window(2);
            t = (edges(1:end-1) + edges(2:end)) / 2;
            X = zeros(length(obj), length(t));
            N = zeros(length(obj), 1);

            tTic = tic();
            fprintf(1, 'Processing ETA for %i units...', length(obj))
            for i = 1:length(obj)
                if strcmpi(alignTo, 'default') 
                    if strcmpi(event, 'stim')
                        alignTo = 'start';
                    else
                        alignTo = 'stop';
                    end
                end
                if isempty(includeITI)
                    if strcmpi(event, 'stim')
                        includeITI = true;
                    else
                        includeITI = false;
                    end
                end

                [x, ~] = obj(i).getTrialAlignedData(data, window, event, alignTo=alignTo, allowedTrialDuration=[minTrialDuration, Inf], resolution=resolution, includeITI=includeITI);
                
                if isempty(x)
                    warning('EventTriggeredAverage cannot be calculated for Unit %i (%s), likely because trial count is zero.', i, obj(i).getName())
                    X(i, :) = NaN;
                    N(i) = 0;
                    continue
                end
                
                % Normalize
                if isnumeric(normalize)
                    inNormWindow = t >= normalize(1) & t <= normalize(2);
                    stats(i).mean = mean(x(:, inNormWindow), 'all', 'omitnan');
                    stats(i).sd = std(x(:, inNormWindow), 0, 'all', 'omitnan');
                    % Average
                    n = size(x, 1);
                    x = mean(x, 1, 'omitnan');
                    x = EphysUnit.normalize(x, 'manual', stats(i));
                elseif isstruct(normalize)
                    % Average
                    n = size(x, 1);
                    x = mean(x, 1, 'omitnan');
                    x = EphysUnit.normalize(x, 'manual', normalize(i));
                elseif ~strcmpi(normalize, 'none')
                    switch data
                        case 'count'
                            stats = obj(i).SpikeCountStats;
                        case 'rate'
                            stats = obj(i).SpikeRateStats;
                    end

                    % Average
                    n = size(x, 1);
                    x = mean(x, 1, 'omitnan');
                    x = EphysUnit.normalize(x, normalize, stats);
                else
                    n = size(x, 1);
                    x = mean(x, 1, 'omitnan');
                end
                
                % Write results
                X(i, :) = x;
                N(i) = n;
            end
            fprintf(1, 'Done (%.1f sec)\n', toc(tTic))
    
            if isnumeric(normalize)
                eta = struct('X', X, 't', t, 'N', N, 'stats', stats);
            else
                eta = struct('X', X, 't', t, 'N', N);
            end
        end

        function rd = getRasterData(obj, trialType, varargin)
            p = inputParser();
            p.addRequired('trialType', @(x) all(ismember(x, {'press', 'lick', 'stim', 'stimtrain', 'stimfirstpulse'})))
            p.addOptional('window', [0, 0], @(x) isnumeric(x) && length(x) >= 2 && x(1) <= 0 && x(2) >= 0)
            p.addParameter('minTrialDuration', 0, @(x) isnumeric(x) && length(x)==1 && x>=0)
            p.addParameter('alignTo', 'default', @(x) ismember(x, {'default', 'start', 'stop'}))
            p.addParameter('sort', true, @islogical);
            p.parse(trialType, varargin{:})
            trialType = p.Results.trialType;
            window = p.Results.window;
            minTrialDuration = p.Results.minTrialDuration;
            alignTo = p.Results.alignTo;

            if strcmp(alignTo, 'default')
                switch trialType
                    case {'press', 'lick'}
                        alignTo = 'stop';
                    case {'stim', 'stimtrain', 'stimfirstpulse'}
                        alignTo = 'start';
                end
            end

            if length(obj) == 1
                trials = obj.getTrials(trialType, sorted=true);
                trials = trials(trials.duration() >= minTrialDuration);
                [~, t, I] = trials.inTrial(obj.SpikeTimes, window);
                switch alignTo
                    case 'start'
                        tRef = [trials.Start];
                        tEvent = [trials.Stop] - tRef;
                    case 'stop'
                        tRef = [trials.Stop];
                        tEvent = [trials.Start] - tRef;
                end
                t = t - tRef(I);
                IEvent = unique(I);
                tEvent = tEvent(IEvent);

                if p.Results.sort
                    [~, ISort] = sort(abs(tEvent), 'ascend');
                    tEvent = tEvent(ISort);
                    IEvent = 1:length(ISort);
                    I = changem(I, 1:length(ISort), ISort);
                end
                rd.name = obj.getName('_');
                rd.trialType = trialType;
                rd.alignTo = alignTo;
                rd.t = t;
                rd.I = I;
                rd.tEvent = tEvent;
                rd.IEvent = IEvent;
            else
                tTic = tic();
                rd(length(obj)) = struct('name', '', 'trialType', '', 'alignTo', '', 't', [], 'I', [], 'tEvent', [], 'IEvent', []);
                fprintf(1, 'Calculating raster data for %g units...\n', length(obj));
                for i = 1:length(obj)
                    try
                        rd(i) = obj(i).getRasterData(trialType, window, minTrialDuration=minTrialDuration, alignTo=alignTo);
                    catch ME
                        warning('\tError when processing unit %g', i)
                        warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
                    end
                end
                fprintf('Done and it only took %.2f sec. Easy.\n', toc(tTic));
            end
        end
    end
    
    % static methods
    methods (Static)

        function obj = load(varargin)
            p = inputParser();
            p.addOptional('path', 'C:\SERVER\Units', @isfolder);
            p.addParameter('waveforms', true, @islogical); % FALSE to cull waveforms after loading to save memory
            p.addParameter('spikecounts', true, @islogical); % FALSE to cull spikecounts after loading to save memory
            p.addParameter('spikerates', true, @islogical); % FALSE to cull spikerate after loading to save memory
            p.parse(varargin{:});
            path = p.Results.path;

            assert(isfolder(path));
            files = dir(sprintf('%s\\*.mat', path));
            S(length(files)) = struct('eu', []);
            for i = 1:length(files)
                tTic = tic();
                fprintf(1, 'Reading unit %g/%g...', i, length(files));
                S(i) = load(sprintf('%s\\%s', files(i).folder, files(i).name), 'eu');
                if ~p.Results.waveforms
                    S(i).eu.Waveforms = [];
                    S(i).eu.WaveformTimestamps = [];
                end
                if ~p.Results.spikecounts
                    S(i).eu.SpikeCounts = single([]);
                    S(i).eu.SpikeCountTimestamps = single([]);
                end
                if ~p.Results.spikerates
                    S(i).eu.SpikeRates = single([]);
                    S(i).eu.SpikeRateTimestamps = single([]);
                end
                fprintf(1, 'Done (%.2f s).\n', toc(tTic));
            end
            obj = [S.eu];
        end

        function varargout = plotETA(varargin)
            p = inputParser();
            if isgraphics(varargin{1}, 'axes')
                p.addRequired('ax', @(x) isgraphics(x, 'axes'));
            end
            p.addRequired('eta', @isstruct);
            p.addOptional('sel', [], @(x) isnumeric(x) || islogical(x))
            p.addParameter('order', [], @isnumeric);
            p.addParameter('clim', [], @isnumeric)
            p.addParameter('xlim', [], @isnumeric)
            p.addParameter('sortWindow', [-2, 0], @(x) isnumeric(x) && length(x) == 2)
            p.addParameter('signWindow', [-.5, 0], @(x) isnumeric(x) && length(x) == 2)
            p.addParameter('sortThreshold', 1, @isnumeric)
            p.addParameter('negativeSortThreshold', [], @isnumeric)
            p.addParameter('hidecolorbar', false, @islogical)
            p.parse(varargin{:})
            sortWindow = p.Results.sortWindow;
            signWindow = p.Results.signWindow;
            sortThreshold = p.Results.sortThreshold;
            negativeSortThreshold = p.Results.negativeSortThreshold;
            if isempty(negativeSortThreshold)
                negativeSortThreshold = sortThreshold;
            end
            if isfield(p.Results, 'ax')
                ax = p.Results.ax;
            else
                f = figure('Units', 'normalized', 'OuterPosition', [0, 0, 0.2, 1], 'DefaultAxesFontSize', 14);
                ax = axes(f);
            end
            X = p.Results.eta.X;
            t = p.Results.eta.t;
            N = p.Results.eta.N;
            if ~isempty(p.Results.sel)
                X = X(p.Results.sel, :);
                N = N(p.Results.sel);
            end

            if isempty(p.Results.order)
                % Sort by first significant response
                XSort = X(:, t >= sortWindow(1) & t <= sortWindow(2));
                tSort = t(t >= sortWindow(1) & t <= sortWindow(2));
                meta = mean(X(:, t >= signWindow(1) & t <= signWindow(2)), 2, 'omitnan');
                etaSign = sign(meta);
                assert(size(etaSign, 2) == 1);
                isAboveThreshold = (XSort >= sortThreshold.*etaSign & etaSign > 0) | (XSort <= negativeSortThreshold.*etaSign & etaSign < 0);
                [~, Ilate] = max(isAboveThreshold, [], 2, 'omitnan');
                [~, order] = sort(Ilate .* etaSign, 'descend');
                meta = meta(order);
                latency = tSort(Ilate);
                latency = latency(order);
                varargout = {ax, order(:), meta(:), latency(:)};
            else
                order = p.Results.order;
                if nargout < 3
                    varargout = {ax, order};
                else
                    XSort = X(:, t >= sortWindow(1) & t <= sortWindow(2));
                    tSort = t(t >= sortWindow(1) & t <= sortWindow(2));
                    meta = mean(X(:, t >= signWindow(1) & t <= signWindow(2)), 2, 'omitnan');
                    etaSign = sign(meta);
                    assert(size(etaSign, 2) == 1);
                    isAboveThreshold = (XSort >= sortThreshold.*etaSign & etaSign > 0) | (XSort <= negativeSortThreshold.*etaSign & etaSign < 0);
                    [~, Ilate] = max(isAboveThreshold, [], 2, 'omitnan');
                    meta = meta(order);
                    latency = tSort(Ilate);
                    latency = latency(order);
                    varargout = {ax, order(:), meta(:), latency(:)};
                end
            end
            
            imagesc(ax, t, 1:length(N), X(order, :))
            if ~isempty(p.Results.clim)
                clim(ax, p.Results.clim);
            end
            if ~isempty(p.Results.xlim)
                xlim(ax, p.Results.xlim);
            end
            colormap(ax, 'turbo')
            h = colorbar(ax, 'eastoutside');
            h.Label.String = 'Normalized spike rate (a.u.)';
            if p.Results.hidecolorbar
                h.Visible = 'off';
            end
            xlabel(ax, 'Time relative to movement (s)')
            ylabel(ax, 'Unit')
            title(sprintf('Event-triggered average (%g units)', length(N)))
        end

        function [ax, order, meta, latency] = plotDoubleETA(eta1, eta2, varargin)
            p = inputParser();
            p.addRequired('eta1', @isstruct)
            p.addRequired('eta2', @isstruct)
            p.addOptional('sel', [], @(x) isnumeric(x) || islogical(x))
            p.addOptional('label1', '', @ischar)
            p.addOptional('label2', '', @ischar)
            p.addParameter('clim', [], @isnumeric)
            p.addParameter('xlim', [], @(x) isnumeric(x) || (iscell(x) && all(cellfun(@isnumeric, x))))
            p.addParameter('sortWindow', [-2, 0], @(x) isnumeric(x) && length(x) == 2)
            p.addParameter('signWindow', [-.5, 0], @(x) iscell(x) || isnumeric(x) && length(x) == 2)
            p.addParameter('sortThreshold', 1, @isnumeric)
            p.addParameter('negativeSortThreshold', [], @isnumeric)
            p.parse(eta1, eta2, varargin{:});
            r = p.Results;
            eta1 = r.eta1;
            eta2 = r.eta2;
            label1 = r.label1;
            label2 = r.label2;

            f = figure(Units='normalized', OuterPosition=[0, 0, 0.5, 1], DefaultAxesFontSize=14);
            ax(1) = axes(f, Position=[0.13, 0.11, 0.3347, 0.815]);
            ax(2) = axes(f, Position=[0.5703, 0.11, 0.3347, 0.815]);
            
            if iscell(r.xlim)
                xlim1 = r.xlim{1};
                xlim2 = r.xlim{2};
            else
                xlim1 = r.xlim;
                xlim2 = r.xlim;
            end

            if iscell(r.signWindow)
                signWindow1 = r.signWindow{1};
                signWindow2 = r.signWindow{2};
            else
                signWindow1 = r.signWindow;
                signWindow2 = r.signWindow;
            end

            [~, order, meta1, latency1] = EphysUnit.plotETA(ax(1), eta1, r.sel, clim=r.clim, xlim=xlim1, sortWindow=r.sortWindow, signWindow=signWindow1, sortThreshold=r.sortThreshold, negativeSortThreshold=r.negativeSortThreshold, hidecolorbar=true);
            [~, ~, meta2, latency2] = EphysUnit.plotETA(ax(2), eta2, r.sel, order=order, clim=r.clim, xlim=xlim2, sortWindow=r.sortWindow, signWindow=signWindow2, sortThreshold=r.sortThreshold, negativeSortThreshold=r.negativeSortThreshold);

            meta = horzcat(meta1(:), meta2(:));
            latency = horzcat(latency1(:), latency2(:));

            title(ax(1), label1);
            title(ax(2), label2);
        end

        function ax = plotBinnedTrialAverage(varargin)
            p = inputParser();
            if isgraphics(varargin{1}, 'axes')
                p.addRequired('ax', @(x) isgraphics(x, 'axes'));
            end
            p.addRequired('S', @isstruct)
            p.addOptional('xlim', [0, 11], @(x) isnumeric(x) && length(x) == 2 && x(2) > x(1))
            p.parse(varargin{:})
            if isfield(p.Results, 'ax')
                ax = p.Results.ax;
            else
                f = figure();
                ax = axes(f);
            end
            X = p.Results.S.X;
            T = p.Results.S.T;
            N = p.Results.S.N;
            S = p.Results.S.S;
            B = p.Results.S.B;
            
            hold(ax, 'on')
            colors = 'rgcbmkrgcbmkrgcbmkrgcbmkrgcbmkrgcbmk';
            for iBin = 1:length(B)
                bin = B(iBin, :);
                t = T + bin(2);
                high = X(iBin, :) + 0.25*S(iBin, :); 
                low = X(iBin, :) - 0.25*S(iBin, :);
                selS = ~isnan(high);
                h(iBin) = plot(ax, t, X(iBin, :), colors(iBin), 'DisplayName', sprintf('[%.1fs, %.1fs], %i trials', bin(1), bin(2), N(iBin)));
                patch(ax, [t(selS), flip(t(selS))], [low(selS), flip(high(selS))], colors(iBin), 'FaceAlpha', 0.1, 'EdgeColor', 'none')
                yl = ax.YLim;
                plot(ax, [bin(2), bin(2)], yl, '--', 'Color', h(iBin).Color)
                ylim(ax, yl);
            end
            hold(ax, 'off')
            xlim(ax, p.Results.xlim)
            legend(ax, h)
        end
        
        function ax = plotRaster(varargin)
            p = inputParser();
            if isgraphics(varargin{1}, 'axes')
                p.addRequired('ax', @(x) isgraphics(x, 'axes'));
            end
            p.addRequired('rd', @isstruct);
            p.addParameter('xlim', [-6, 1], @(x) isnumeric(x) && length(x) == 2 && x(2) > x(1));
            p.addParameter('stim', false, @islogical)
            p.parse(varargin{:})
            rd = p.Results.rd;
            stim = p.Results.stim;

            assert(length(rd) == 1);

            if isfield(p.Results, 'ax')
                ax = p.Results.ax;
            else
                ax = axes(figure(Units='normalized', Position=[0.1, 0.1, 0.67, 0.5]));
            end

            switch rd.alignTo
                case 'start'
                    eventName = rd.trialType;
                    refName = 'Trial start';
                case 'stop'
                    eventName = 'Trial start';
                    refName = rd.trialType;
            end

            hold(ax, 'on')
            if ~stim
                h = gobjects(2, 1);
                h(1) = scatter(ax, rd.t, rd.I, 5, 'k', 'filled', DisplayName='Spikes');
                h(2) = scatter(ax, rd.tEvent, rd.IEvent, 15, 'r', 'filled', DisplayName=eventName);
            else
                h = gobjects(3, 1);
                h(1) = scatter(ax, rd.t, rd.I, 5, 'k', 'filled', DisplayName='Spikes');
                h(2) = plot(ax, [0, 0], [min(rd.IEvent), max(rd.IEvent)], 'b', DisplayName='Opto On');
                h(3) = scatter(ax, rd.tEvent, rd.IEvent, 5, 'b', 'filled', DisplayName='Opto Off');
            end

            if stim
            end
            hold(ax, 'off')
            ax.YAxis.Direction = 'reverse';
            xlim(ax, p.Results.xlim);
            ylim(ax, [min(rd.I) - 1, max(rd.I) + 1]);
            xlabel(ax, sprintf('Time relative to %s (s)', refName))
            ylabel(ax, 'Trial')
            title(ax, sprintf('Spike raster (%s)', rd.name), Interpreter="none");
            legend(ax, h, Location='northwest');
        end

    end

    % private methods
    methods %(Access = {})
        function trials = makeTrials(obj, trialType)
            if length(obj) == 1
                switch lower(trialType)
                    case 'press'
                        trials = Trial(obj.EventTimes.Cue, obj.EventTimes.Press, 'first', obj.EventTimes.Lick);
                    case 'lick'
                        trials = Trial(obj.EventTimes.Cue, obj.EventTimes.Lick, 'first', obj.EventTimes.Press);
                    case 'stim'
                        trials = Trial(obj.EventTimes.StimOn, obj.EventTimes.StimOff, 'first');
                    case 'stimtrain'
                        cueToTrainOn = Trial(obj.EventTimes.Cue, obj.EventTimes.StimOn, 'first');
                        cueToTrainOff = Trial(obj.EventTimes.Cue, obj.EventTimes.StimOff, 'last');
                        trials = Trial([cueToTrainOn.Stop], [cueToTrainOff.Stop]);
                end
            else
                trials = cell(length(obj), 1);
                for i = 1:length(obj)
                    trials{i} = obj(i).makeTrials(trialType);
                end
            end
        end

        function [t, x, inTrial] = cullITIData(obj, t, varargin)
            assert(length(obj) == 1)
            
            %CULLDISCRETEDATA is used to cull spikeTimes and spikeWaveforms
            %during ITI.
            p = inputParser();
            p.addRequired('t', @isnumeric);
            p.addOptional('x', {}, @(x) iscell(x) || isnumeric(x));
            p.addOptional('trialType', 'all', @(x) ischar(x) && ismember(lower(x), {'press', 'lick', 'stim', 'all'}));
            p.addParameter('extendedWindow', [0, 0], @(x) isnumeric(x) && length(x)>=2 && x(1)<=0 && x(2)>=0)
            p.parse(t, varargin{:})
            t = p.Results.t;
            x = p.Results.x;
            trialType = lower(p.Results.trialType);
            extendedWindow = p.Results.extendedWindow;
            
            switch trialType
                case 'press'
                    inTrial = obj.getTrials('press').inTrial(t, extendedWindow);
                case 'lick'
                    inTrial = obj.getTrials('lick').inTrial(t, extendedWindow);
                case 'stim'
                    inTrial = obj.getTrials('stimtrain').inTrial(t, extendedWindow);
                case 'all'
                    inTrialPress = obj.getTrials('press').inTrial(t, extendedWindow);
                    inTrialLick = obj.getTrials('lick').inTrial(t, extendedWindow);
                    inTrialStim = obj.getTrials('stimtrain').inTrial(t, extendedWindow);
                    inTrial = inTrialPress | inTrialLick | inTrialStim;
            end
            
            t = t(inTrial);
            if ~isempty(x)
                if isnumeric(x)
                    if numel(x) == length(x)
                        x = x(inTrial);
                    else
                        x = x(inTrial, :);
                    end
                elseif iscell(x)
                    for ic = 1:length(x)
                        xx = x{ic};
                        if numel(xx) == length(xx)
                            x{ic} = xx(inTrial);
                        else
                            x{ic} = xx(inTrial, :);
                        end
                    end
                end
            end
        end

        function [sc, t] = getSpikeCounts(obj, varargin)
            % GETSPIKECOUNTS Generates binned spike counts
            %   [sc, t] = GETSPIKECOUNTS(binWidth) uses binWidth to automatically construct bins
            %   [sc, t] = GETSPIKECOUNTS(edges) uses specific bin edges
            assert(length(obj) == 1)
            p = inputParser();
            if length(varargin{1}) == 1
                p.addRequired('binWidth', @(x) isnumeric(x) && x>0);
            else
                p.addRequired('edges', @(x) isnumeric(x) && length(x)>=2 && nnz(diff(x)<=0)==0)
            end
            p.parse(varargin{:})
            if isfield(p.Results, 'binWidth')
                binWidth = p.Results.binWidth;
                edges = [];
            else
                binWidth = [];
                edges = p.Results.edges;
            end
            
            spikes = obj.SpikeTimes;
            if isempty(edges)
                edges = spikes(1) - binWidth:binWidth:spikes(end) + binWidth;
            end
            sc = histcounts(spikes, edges);
            t = (edges(1:end-1) + edges(2:end)) / 2;

            assert(sum(sc > 2^16 - 1) == 0, 'Spike counts per bin should exceeded 2^16-1 (%i)', max(sc)) % There should not be more than 255 spikes in an 100ms bin
            sc = uint16(sc);
            t = single(t);
        end
        
        function [sr, t, kernel] = getSpikeRates(obj, varargin)
            % GETSPIKERATES Convolve discrete spikes with a Guassian (default) or Exponential kernel to get smooth spike rate estimate
            assert(length(obj) == 1)
            p = inputParser();
            defaultKernelType = 'gaussian';
            defaultResolution = 1e-3;
            p.addOptional('kernel', defaultKernelType, @(x) ischar(x) && ismember(lower(x), {'gaussian', 'exponential'}))
            if length(varargin) < 1
                kernelType = defaultKernelType;
            else
                kernelType = lower(varargin{1});
            end
            switch kernelType
                case 'gaussian'
                    p.addOptional('sigma', 0.1, @(x) isnumeric(x) && x > 0)
                case 'exponential'
                    p.addOptional('lambda1', 5, @(x) isnumeric(x) && x > 0)
                    p.addOptional('lambda2', 10, @(x) isnumeric(x) && x > 0)
            end
            p.addOptional('edgesOrResolution', defaultResolution, @isnumeric)
            p.addParameter('kernelWidth', 1.0, @(x) isnumeric(x) && x > 0)
            p.parse(varargin{:})
            if length(p.Results.edgesOrResolution) <= 1
                resolution = p.Results.edgesOrResolution;
                edges = [];
            else
                edges = p.Results.edgesOrResolution;
                assert(all(diff(edges) > 0), '''edges'' must be monotonic increasing.')
                assert(all(single(diff(edges)) == single(edges(2) - edges(1))), '''edges'' must have equal distance between neighboring elements.')
                resolution = edges(2) - edges(1);
            end
            kernelWidth = p.Results.kernelWidth;
            switch kernelType
                case 'gaussian'
                    kernelParams.sigma = p.Results.sigma;
                    kernelWindow = kernelWidth * [-0.5, 0.5];
                    tKernel = kernelWindow(1):resolution:kernelWindow(2);
                    yKernel = normpdf(tKernel, 0, kernelParams.sigma);
                case 'exponential'
                    kernelParams.lambda1 = p.Results.lambda1;
                    kernelParams.lambda2 = p.Results.lambda2;
                    kernelWindow = [0, kernelWidth];
                    tKernel = kernelWindow(1):resolution:kernelWindow(2);
                    yKernel = exp(-kernelParams.lambda1*tKernel) - exp(-kernelParams.lambda2*tKernel);
            end
            kernelParams.window = kernelWindow;
            kernelParams.resolution = resolution;
            kernelParams.width = kernelWidth;
            yKernel = yKernel / sum(yKernel) / resolution;
            kernel = struct('type', kernelType, 'params', kernelParams, 't', tKernel, 'y', yKernel);
            
            spikes = obj.SpikeTimes;
            if isempty(edges)
                edges = spikes(1) + kernelWindow(1):resolution:spikes(end) + kernelWindow(2);
                nPrepad = 0;
                nPostpad = 0;
            else
                prepad = kernelWindow(1)+edges(1):resolution:edges(1)-resolution;
                postpad = edges(end)+resolution:resolution:edges(end)+kernelWindow(2);
                nPrepad = length(prepad);
                nPostpad = length(postpad);
                edges = [prepad, edges, postpad];
            end
            spikeCounts = histcounts(spikes, edges);
            sr = conv(spikeCounts, yKernel, 'same');
            sr = sr(1+nPrepad:end-nPostpad);
            t = (edges(1:end-1) + edges(2:end)) / 2;
            t = t(1+nPrepad:end-nPostpad);
            
            sr = single(sr);
            t = single(t);
        end

        function [xAligned, tAligned] = getTrialAlignedData(obj, varargin)
            assert(length(obj) == 1)
            p = inputParser();
            if isnumeric(varargin{1}) && isnumeric(varargin{2})
                p.addRequired('x', @isnumeric)
                p.addRequired('t', @isnumeric)
                useResampleMethod = true;
            elseif ischar(varargin{1})
                p.addRequired('data', @(x) ischar(x) && ismember(lower(x), {'rate', 'count'}))
                useResampleMethod = false;
            end
            p.addOptional('window', [-4, 0], @(x) isnumeric(x) && length(x) >= 2)
            p.addOptional('trialType', 'press', @(x) ischar(x) && ismember(lower(x), {'press', 'lick', 'stim'}))
            p.addParameter('alignTo', 'stop', @(x) ischar(x) && ismember(lower(x), {'start', 'stop'}))
            p.addParameter('resolution', 0.001, @(x) isnumeric(x) && x > 0)
            p.addParameter('allowedTrialDuration', [0, Inf], @(x) isnumeric(x) && length(x) >= 2 && x(2) > x(1))
            p.addParameter('includeITI', true, @islogical) % Whether to include unaligned ITI data for averaging. When alignTo='stop', pre-trial-start data is discarded. When alignTo='start', post-trial-end data is discarded.
            p.parse(varargin{:})
            if useResampleMethod
                x = p.Results.x;
                t = p.Results.t;
            else
                data = lower(p.Results.data);
            end
            window = p.Results.window(1:2);
            trialType = lower(p.Results.trialType);
            alignTo = lower(p.Results.alignTo);
            resolution = p.Results.resolution(1);
            allowedTrialDuration = p.Results.allowedTrialDuration(1:2);
            includeITI = p.Results.includeITI;
                        
            % Filter out trials with incorrect lengths
            trials = obj.getTrials(trialType);
            durations = trials.duration();
            trials = trials(durations >= allowedTrialDuration(1) & durations <= allowedTrialDuration(2));
            if isempty(trials)
                xAligned = [];
                tAligned = [];
                return
            end

            tAligned = window(1):resolution:window(2);
            switch alignTo
                case 'stop'
                    tAlignedGlobal = tAligned + vertcat(trials.Stop);
                case 'start'
                    tAlignedGlobal = tAligned + vertcat(trials.Start);
            end

            function sel = select(tt, iTrial)
                if includeITI
                    sel = true(1, length(tt));
                elseif strcmpi(alignTo, 'stop')
                    sel = tt >= -trials(iTrial).duration();
                else % if strcmpi(alignTo, 'start')
                    sel = tt <= trials(iTrial).duration();
                end
            end
            
            % Resample by interpolating data (x, t) at new t
            if useResampleMethod
                xAligned = NaN(length(trials), length(tAligned));
                x = double(x);
                for iTrial = 1:length(trials)
                    [xx, tt] = interp1(t, x, tAlignedGlobal(iTrial, :), 'linear');
                    sel = select(tt, iTrial);
                    xAligned(iTrial, sel) = xx(sel);
                end
            % Recalculate data (count or rate) in new bins.
            else
                tAligned = (tAligned(1:end-1) + tAligned(2:end)) / 2;
                xAligned = NaN(length(trials), length(tAligned));
                switch data
                    case 'rate'
                        kernel = obj.SpikeRateKernel;
                        width = kernel.params.width;
                        if strcmpi(kernel.type, 'gaussian')
                            sigma = kernel.params.sigma;
                            for iTrial = 1:length(trials)
                                [xx, ~] = obj.getSpikeRates('gaussian', sigma, tAlignedGlobal(iTrial, :), 'kernelWidth', width);
                                sel = select(tAligned, iTrial);
                                xAligned(iTrial, sel) = xx(sel);
                            end
                        else
                            lambda1 = kernel.params.lambda1;
                            lambda2 = kernel.params.lambda2;
                            for iTrial = 1:length(trials)
                                [xx, ~] = obj.getSpikeRates('gaussian', lambda1, lambda2, tAlignedGlobal(iTrial, :), 'kernelWidth', width);
                                sel = select(tAligned, iTrial);
                                xAligned(iTrial, sel) = xx(sel);
                            end
                        end
                    case 'count'
                        for iTrial = 1:length(trials)
                            [xx, ~] = obj.getSpikeCounts(tAlignedGlobal(iTrial, :));
                            sel = select(tAligned, iTrial);
                            xAligned(iTrial, sel) = xx(sel);
                        end
                end
            end
            
        end
    end
    
    % private staic
    methods (Static, Access = {})
        function z = normalize(x, mode, stats)
            switch lower(mode)
                case 'iti'
                    m = double(stats.medianITI);
                    s = double(stats.madITI) / 0.6745;
                case 'all'
                    m = double(stats.median);
                    s = double(stats.mad) / 0.6745;
                case 'manual'
                    m = double(stats.mean);
                    s = double(stats.sd);
            end
            z = (x - m) ./ s;
        end
        
        function [mu, ss, k] = combinestats(mu_x, ss_x, m, mu_y, ss_y, n)
            k = m + n;
            mu = (m*mu_x + n*mu_y) / k;
            ss = ss_x + ss_y + m*mu_x.^2 + n*mu_y.^2 - 2*(m*mu_x + n*mu_y).*mu + (m+n)*mu.^2;
        end
    end
end
