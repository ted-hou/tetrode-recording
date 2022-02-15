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
        
        function [X, T, N, S, B] = getBinnedTrialAverage(obj, data, varargin)
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
                    
                    [xx, tt] = obj(i).getTrialAlignedData(data, window, trialType, 'allowedTrialDuration', [edges(iBin), edges(iBin+1)], 'alignTo', alignTo, 'resolution', resolution);

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
                    
                    mu = mean(xx, 1);
                    ss = sum((xx - mu).^2, 1);
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
        
        function [X, t, N] = getPETH(obj, data, event, varargin)
            %GETPETH estimate mean spikerate around an
            %event
            p = inputParser();
            p.addRequired('data', @(x) ischar(x) && ismember(lower(x), {'rate', 'count'}))
            p.addRequired('event', @(x) ischar(x) && ismember(lower(x), {'press', 'lick', 'stim'}))
            p.addOptional('window', [-2, 0], @(x) isnumeric(x) && length(x)>=2 && x(2) > x(1))
            p.addParameter('minTrialDuration', -Inf, @(x) isnumeric(x) && length(x)==1 && x>=0)
            p.addParameter('resolution', [], @(x) isnumeric(x) && length(x)==1 && x>=0)
            p.addParameter('normalize', 'none', @(x) ischar(x) && ismember(lower(x), {'none', 'iti', 'all'}))
            p.parse(data, event, varargin{:})
            data = lower(p.Results.data);
            event = p.Results.event;
            window = p.Results.window;
            minTrialDuration = p.Results.minTrialDuration;
            resolution = p.Results.resolution;
            normalize = p.Results.normalize;
            
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
                if strcmpi(event, 'stim')
                    alignTo = 'stop';
                else
                    alignTo = 'start';
                end
                
                [x, ~] = obj(i).getTrialAlignedData(data, window, event, 'alignTo', alignTo, 'allowedTrialDuration', [minTrialDuration, Inf], 'resolution', resolution);
                
                if isempty(x)
                    warning('EventTriggeredAverage cannot be calculated for Unit %i (%s), likely because trial count is zero.', i, obj(i).getName())
                    X(i, :) = NaN;
                    N(i) = 0;
                    continue
                end
                
                % Average
                n = size(x, 1);
                x = mean(x, 1);
                
                % Normalize
                if ~strcmpi(normalize, 'none')
                    switch data
                        case 'count'
                            stats = obj(1).SpikeCountStats;
                        case 'rate'
                            stats = obj(1).SpikeRateStats;
                    end
                    x = EphysUnit.normalize(x, normalize, stats);
                end
                
                % Write results
                X(i, :) = x;
                N(i) = n;
            end
            fprintf(1, 'Done (%.1f sec)\n', toc(tTic))
        end
        
    end
    
    % static methods
    methods (Static)

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
            colors = 'rgbcmkrgbcmkrgbcmkrgbcmkrgbcmkrgbcmk';
            for iBin = 1:length(B)
                bin = B(iBin, :);
                t = T + bin(2);
                high = X(iBin, :) + 0.25*S(iBin, :); 
                low = X(iBin, :) - 0.25*S(iBin, :);
                h(iBin) = plot(ax, t, X(iBin, :), colors(iBin), 'DisplayName', sprintf('[%.1fs, %.1fs], %i trials', bin(1), bin(2), N(iBin)));
                patch(ax, [t, flip(t)], [low, flip(high)], colors(iBin), 'FaceAlpha', 0.1, 'EdgeColor', 'none')
                plot(ax, [bin(2), bin(2)], [0, max(X(iBin, :))], '--', 'Color', h(iBin).Color)
            end
            hold(ax, 'off')
            xlim(ax, p.Results.xlim)
            legend(ax, h)
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
                        
            % Filter out trials with incorrect lengths
            trials = obj.getTrials(trialType);
            durations = trials.duration();
            trials = trials(durations >= allowedTrialDuration(1) & durations < allowedTrialDuration(2));
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
                    tAlignedGlobal = tAligned + vertcat(trials.Stop);
            end
            
            % Resample by interpolating data (x, t) at new t
            if useResampleMethod
                xAligned = zeros(length(trials), length(tAligned));
                x = double(x);
                for iTrial = 1:length(trials)
                    xAligned(iTrial, :) = interp1(t, x, tAlignedGlobal(iTrial, :), 'linear');
                end
            % Recalculate data (count or rate) in new bins.
            else
                tAligned = (tAligned(1:end-1) + tAligned(2:end)) / 2;
                xAligned = zeros(length(trials), length(tAligned));
                switch data
                    case 'rate'
                        kernel = obj.SpikeRateKernel;
                        width = kernel.params.width;
                        if strcmpi(kernel.type, 'gaussian')
                            sigma = kernel.params.sigma;
                            for iTrial = 1:length(trials)
                                [xAligned(iTrial, :), ~, ~] = obj.getSpikeRates('gaussian', sigma, tAlignedGlobal(iTrial, :), 'kernelWidth', width);
                            end
                        else
                            lambda1 = kernel.params.lambda1;
                            lambda2 = kernel.params.lambda2;
                            for iTrial = 1:length(trials)
                                [xAligned(iTrial, :), ~, ~] = obj.getSpikeRates('gaussian', lambda1, lambda2, tAlignedGlobal(iTrial, :), 'kernelWidth', width);
                            end
                        end
                    case 'count'
                        for iTrial = 1:length(trials)
                            [xAligned(iTrial, :), ~] = obj.getSpikeCounts(tAlignedGlobal(iTrial, :));
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
