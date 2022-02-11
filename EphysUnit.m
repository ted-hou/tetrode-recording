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
        Trials = struct('Press', [], 'Lick', [], 'Stim', [])
        SpikeCounts = uint16([])
        SpikeCountTimestamps = single([])
        SpikeRates = single([])
        SpikeRateTimestamps = single([])
        SpikeRateKernel = struct('type', [], 'window', [], 'resolution', [], 't', [], 'y', [])
        SpikeRateMedian = []
        SpikeRateMAD = []
        ITISpikeRateMedian = []
        ITISpikeRateMAD = []
    end
    
    % public methods
    methods
        % constructor
        function obj = EphysUnit(varargin)
            if nargin == 0
                return
            end

            p = inputParser();
            p.addOptional('PETH', struct(), @isstruct);
            p.addParameter('readWaveforms', false, @islogical);
            p.addParameter('cullITI', false, @islogical);
            p.addParameter('extendedWindow', [-1, 1], @(x) isnumeric(x) && length(x) == 2)
            p.parse(varargin{:});
            PETH = p.Results.PETH;
            extendedWindow = p.Results.extendedWindow;

            % Initialize obj array
            obj(length(PETH)) = EphysUnit();

            % Extract unique expNames
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

                        % Calculate spike rates and spike counts
                        obj(i).getSpikeCounts(0.1);
                        obj(i).getSpikeRates('gaussian', 0.1, 'resolution', 1e-3, 'kernelWidth', 1);

                        % Cull ITI spikes if told to do so.
                        if p.Results.cullITI
                            [obj(i).SpikeTimes, obj(i).Waveforms] = obj(i).cullITIData(obj(i).SpikeTimes, obj(i).Waveforms, 'all', 'extendedWindow', extendedWindow);
                        end

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
        
        function [X, T, N, B] = getBinnedTrialAverage(obj, data, varargin)
            p = inputParser();
            p.addRequired('data', @(x) ischar(x) && ismember(lower(x), {'spikerate', 'spikecount', 'spikerates', 'spikecounts', 'rate', 'count'}))
            p.addOptional('edges', 0:2.5:10, @(x) isnumeric(x) && nnz(diff(x) <= 0) == 0)
            p.addOptional('trialType', 'press', @(x) ischar(x) && ismember(lower(x), {'press', 'lick', 'stim'}))
            p.addParameter('alignTo', 'stop', @(x) ischar(x) && ismember(lower(x), {'start', 'stop'}))
            p.addParameter('resolution', 0.001, @(x) isnumeric(x) && x > 0)
            p.addParameter('normalize', false, @islogical)
            p.addParameter('plot', true, @islogical)
            p.parse(data, varargin{:})
            data = lower(p.Results.data);
            edges = p.Results.edges;
            trialType = lower(p.Results.trialType);
            alignTo = lower(p.Results.alignTo);
            resolution = p.Results.resolution;
            
            X = {};
            T = {};
            N = [];
            B = {};
            for i = 1:length(obj)
                for bin = 1:length(edges) - 1
                    switch data
                        case {'rate', 'spikerate', 'spikerates'}
                            x = obj(i).SpikeRates;
                            t = obj(i).SpikeRateTimestamps;
                        case {'count', 'spikecount', 'spikecounts'}
                            x = obj(i).SpikeCounts;
                            t = obj(i).SpikeCountTimestamps;
                    end
                    
                    if ~isempty(obj(i).ExtendedWindow)
                        postWindow = obj(i).ExtendedWindow(2);
                    else
                        postWindow = 1;
                    end
                    [xx, tt] = obj(i).getTrialAlignedData(x, t, [-20, postWindow], trialType, 'allowedTrialDuration', [edges(bin), edges(bin+1)], 'alignTo', alignTo, 'resolution', resolution);
                    
                    % Use modified z-score if requested
                    if p.Results.normalize
                        xx = obj(i).normalizeSpikeRate(xx);
                    end
                    
                    if i == 1
                        N(bin) = size(xx, 1);
                        X{bin} = mean(xx, 1);
                        T{bin} = tt;
                        B{bin} = [edges(bin), edges(bin+1)];
                    else
                        X{bin} = (X{bin} .* N(bin) + mean(xx, 1) .* size(xx, 1));
                        N(bin) = N(bin) + size(xx, 1);
                        X{bin} = X{bin} ./ N(bin);
                    end
                end
            end
            
            if p.Results.plot
                f = figure();
                ax = axes(f);
                hold(ax, 'on')
                for iBin = 1:length(B)
                    bin = B{iBin};
                    h(iBin) = plot(ax, T{iBin} + bin(2), X{iBin}, 'DisplayName', sprintf('[%.1fs, %.1fs], %i trials', bin(1), bin(2), N(iBin)));
                    plot([bin(2), bin(2)], [0, max(X{iBin})], '--', 'Color', h(iBin).Color)
                end
                hold(ax, 'off')
                xlim([0, 11])
                legend(ax, h)
            end
        end
        
        function [X, t, N] = getEventTriggeredAverage(obj, varargin)
            %GETEVENTTRIGGEREDAVERAGE estimate mean spikerate around an
            %event
            p = inputParser();
            p.addRequired('event', @(x) ischar(x) && ismember(lower(x), {'press', 'lick', 'stim'}))
            p.addOptional('window', [-2, 0], @(x) isnumeric(x) && length(x)>=2 && x(2) > x(1))
            p.addParameter('normalize', 'none', @(x) ischar(x) && ismember(lower(x), {'none', 'iti', 'all'}))
            p.parse(varargin{:})
            event = p.Results.event;
            window = p.Results.window;
            
            for i = 1:length(obj)
                if strcmpi(event, 'stim')
                    alignTo = 'stop';
                else
                    alignTo = 'start';
                end
                [x, t] = obj(i).getTrialAlignedData(obj(i).SpikeRates, obj(i).SpikeRateTimestamps, window, event, 'alignTo', alignTo);
                n = size(x, 1);
                x = nanmean(x, 1);
                switch lower(p.Results.normalize)
                    case 'iti'
                        x = (x - obj(i).ITISpikeRateMedian) / (obj(i).ITISpikeRateMAD / 0.6745);
                    case 'all'
                        x = (x - obj(i).SpikeRateMedian) / (obj(i).SpikeRateMAD / 0.6745);
                end
                
                % Preallocate output
                if length(obj) > 1
                    if i == 1
                        X = zeros(length(obj), size(x, 2));
                        N = zeros(length(obj), 1);
                    end
                    X(i, :) = x;
                    N(i) = n;
                else
                    X = x;
                    N = n;
                end

            end
        end
        
    end

    % private methods
    methods (Access = {})
        function trials = makeTrials(obj, trialType)
            if length(obj) == 1
                switch lower(trialType)
                    case 'press'
                        trials = Trial(obj.EventTimes.Cue, obj.EventTimes.Press, 'first', obj.EventTimes.Lick);
                    case 'lick'
                        trials = Trial(obj.EventTimes.Cue, obj.EventTimes.Lick, 'first', obj.EventTimes.Press);
                    case 'stim'
                        trials = Trial(obj.EventTimes.StimOn, obj.EventTimes.StimOff);
                end
            else
                trials = cell(length(obj), 1);
                for i = 1:length(obj)
                    trials{i} = obj(i).makeTrials(trialType);
                end
            end
        end

        function trials = getTrials(obj, trialType)
            if length(obj) == 1
                switch lower(trialType)
                    case 'press'
                        trials = obj.Trials.Press;
                    case 'lick'
                        trials = obj.Trials.Lick;
                    case 'stim'
                        trials = obj.Trials.Stim;
                    case 'all'
                        trials = [obj.Trials.Press, obj.Trials.Lick, obj.Trials.Stim];
                end
            else
                trials = cell(length(obj), 1);
                for i = 1:length(obj)
                    trials{i} = obj(i).getTrials(trialType);
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
                    inTrial = obj.getTrials('stim').inTrial(t, extendedWindow);
                case 'all'
                    inTrialPress = obj.getTrials('press').inTrial(t, extendedWindow);
                    inTrialLick = obj.getTrials('lick').inTrial(t, extendedWindow);
                    inTrialStim = obj.getTrials('stim').inTrial(t, extendedWindow);
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

        function [xAligned, tAligned] = getTrialAlignedData(obj, x, t, varargin)
            assert(length(obj) == 1)
            
            p = inputParser();
            p.addRequired('x', @isnumeric)
            p.addRequired('t', @isnumeric)
            p.addOptional('window', [-4, 0], @(x) isnumeric(x) && length(x) >= 2)
            p.addOptional('trialType', 'press', @(x) ischar(x) && ismember(lower(x), {'press', 'lick', 'stim'}))
            p.addParameter('alignTo', 'stop', @(x) ischar(x) && ismember(lower(x), {'start', 'stop'}))
            p.addParameter('resolution', 0.001, @(x) isnumeric(x) && x > 0)
            p.addParameter('allowedTrialDuration', [0, Inf], @(x) isnumeric(x) && length(x) >= 2 && x(2) > x(1))
            p.parse(x, t, varargin{:})
            x = p.Results.x;
            t = p.Results.t;
            window = p.Results.window(1:2);
            trialType = lower(p.Results.trialType);
            alignTo = lower(p.Results.alignTo);
            resolution = p.Results.resolution(1);
            allowedTrialDuration = p.Results.allowedTrialDuration(1:2);
                        
            % Filter out trials with incorrect lengths
            trials = obj.getTrials(trialType);
            durations = trials.duration();
            trials = trials(durations >= allowedTrialDuration(1) & durations < allowedTrialDuration(2));
            % durations = trials.duration();

            tAligned = window(1):resolution:window(2);
            xAligned = zeros(length(trials), length(tAligned));
            switch alignTo
                case 'stop'
                    tAlignedGlobal = tAligned + vertcat(trials.Stop);
                case 'start'
                    tAlignedGlobal = tAligned + vertcat(trials.Stop);
            end
            for iTrial = 1:length(trials)
                xAligned(iTrial, :) = interp1(t, double(x), tAlignedGlobal(iTrial, :), 'linear');
            end
        end
        
        function varargout = getSpikeCounts(obj, varargin)
            % GETSPIKECOUNTS Generates binned spike counts
            p = inputParser();
            p.addOptional('binWidth', 0.1, @(x) isnumeric(x) && x > 0 && x <= 1);
            p.parse(varargin{:})
            binWidth = p.Results.binWidth;
            
            for i = 1:length(obj)
                spikes = obj(i).SpikeTimes;
                edges = spikes(1) - binWidth:binWidth:spikes(end) + binWidth;
                spikeCounts = histcounts(spikes, edges);
                t = (edges(1:end-1) + edges(2:end)) / 2;
                
                % Cull ITI data
                if ~isempty(obj(i).ExtendedWindow)
                    [t, spikeCounts] = obj(i).cullITIData(t, spikeCounts, 'all', 'extendedWindow', obj(i).ExtendedWindow);
                end
                assert(sum(spikeCounts > 2^16 - 1) == 0, 'Spike counts per bin should exceeded 2^16-1 (%i)', max(spikeCounts)) % There should not be more than 255 spikes in an 100ms bin
                obj(i).SpikeCounts = uint16(spikeCounts);
                obj(i).SpikeCountTimestamps = single(t);
            end
            
            if length(obj) == 1
                varargout = {obj(i).SpikeCounts, obj(i).SpikeCountTimestamps};
            else
                varargout = {{obj.SpikeCounts}, {obj.SpikeCountTimestamps}};
            end
        end
        
        function varargout = getSpikeRates(obj, varargin)
            % GETSPIKERATES Convolve discrete spikes with a Guassian (default) or Exponential kernel to get smooth spike rate estimate
            p = inputParser();
            defaultKernelType = 'gaussian';
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
            p.addParameter('resolution', 1e-3, @(x) isnumeric(x) && x > 0)
            p.addParameter('kernelWidth', 1.0, @(x) isnumeric(x) && x > 0)
            p.parse(varargin{:})
            resolution = p.Results.resolution;
            kernelWidth = p.Results.kernelWidth;
            switch kernelType
                case 'gaussian'
                    sigma = p.Results.sigma;
                    kernelWindow = kernelWidth * [-0.5, 0.5];
                    tKernel = kernelWindow(1):resolution:kernelWindow(2);
                    yKernel = normpdf(tKernel, 0, sigma);
                    kernel = struct('type', kernelType, 'window', kernelWindow, 'resolution', resolution, 't', tKernel, 'y', yKernel, 'sigma', sigma);
                case 'exponential'
                    lambda1 = p.Results.lambda1;
                    lambda2 = p.Results.lambda2;
                    kernelWindow = [0, kernelWidth];
                    tKernel = kernelWindow(1):resolution:kernelWindow(2);
                    yKernel = exp(-lambda1*tKernel) - exp(-lambda2*tKernel);
                    kernel = struct('type', kernelType, 'window', kernelWindow, 'resolution', resolution, 't', tKernel, 'y', yKernel, 'lambda1', lambda1, 'lambda2', lambda2);
            end
            yKernel = yKernel / sum(yKernel) / resolution;
            
            for i = 1:length(obj)
                spikes = obj(i).SpikeTimes;
                edges = spikes(1) + kernelWindow(1):resolution:spikes(end) + kernelWindow(2);
                spikeCounts = histcounts(spikes, edges);
                sr = conv(spikeCounts, yKernel, 'same');
                
                t = (edges(1:end-1) + edges(2:end)) / 2;
                
                % Cull ITI data
                [~, ~, inTrial] = obj(i).cullITIData(t, sr, 'all', 'extendedWindow', obj(i).ExtendedWindow);
                sr_median = median(sr);
                sr_mad = mad(sr, 1);
                srITI_median = median(sr(~inTrial));
                srITI_mad = mad(sr(~inTrial), 1);
                if ~isempty(obj(i).ExtendedWindow)
                    t = t(inTrial);
                    sr = sr(inTrial);
                end
                
                obj(i).SpikeRateTimestamps = single(t);
                obj(i).SpikeRates = single(sr);
                obj(i).SpikeRateKernel = kernel;
                obj(i).SpikeRateMedian = sr_median;
                obj(i).SpikeRateMAD = sr_mad;
                obj(i).ITISpikeRateMedian = srITI_median;
                obj(i).ITISpikeRateMAD = srITI_mad;
            end
            
            if length(obj) == 1
                varargout = {obj(i).SpikeRates, obj(i).SpikeRateTimestamps, obj(i).SpikeRateKernel};
            else
                varargout = {{obj.SpikeRates}, {obj.SpikeRateTimestamps}, obj(1).SpikeRateKernel};
            end
        end

        function z = normalizeSpikeRate(obj, x, mode)
            assert(length(obj) == 1)
            if nargin < 2
                x = obj.SpikeRates;
            end
            if nargin < 3
                mode = 'ITI';
            end
            
            switch lower(mode)
                case 'iti'
                    m = obj.ITISpikeRateMedian;
                    s = obj.ITISpikeRateMAD / 0.6745;
                case 'all'
                    m = obj.SpikeRateMedial;
                    s = obj.SpikeRateMAD / 0.6745;
            end
            z = (x - m) ./ s;
        end     
    end
end
