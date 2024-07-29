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
        EventTimes = struct('Cue', [], 'Press', [], 'Lick', [], 'StimOn', [], 'StimOff', [], 'LightOn', [], 'LightOff', [])
        Trials = struct('Press', [], 'Lick', [], 'StimTrain', [], 'Stim', [], 'Light', [])
        SpikeCounts = uint16([])
        SpikeCountTimestamps = single([])
        SpikeRates = single([])
        SpikeRateTimestamps = single([])
        SpikeRateKernel = struct('type', [], 'params', [], 't', [], 'y', [])
        SpikeCountStats = struct('median', [], 'mad', [], 'medianITI', [], 'madITI', [], 'resolution', [])
        SpikeRateStats = struct('median', [], 'mad', [], 'medianITI', [], 'madITI', [], 'resolution', [])
        ProbeMap = []
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
            p.addParameter('savepath', '', @ischar)
            p.addParameter('whitelist', {}, @iscell) % If a list of unit names are provided, only these are created as EU objs.
            p.addParameter('tr', [], @(x) isa(x, 'TetrodeRecording'))
            p.addParameter('isSpontaneousPress', false, @islogical)
            p.addParameter('sessionInfo', [], @isstruct)
            p.parse(varargin{:});
            extendedWindow = p.Results.extendedWindow;
            sessionInfo = p.Results.sessionInfo;
            
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
                    if ~isempty(p.Results.tr)
                        tr = p.Results.tr;
                    else
                        tr = TetrodeRecording.BatchLoadSimple(uniqueExpNames{iExp});
                    end
                    inExp = strcmpi({PETH.ExpName}, uniqueExpNames{iExp});
                    tTic = tic();
%                     fprintf(1, 'Processing %i units...\n', nnz(inExp))
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
                            % Spontaneous
                            if isfield(tr.DigitalEvents, 'TRIAL_START')
                                obj(i).EventTimes.Cue = tr.DigitalEvents.TRIAL_START;
                            else
                                obj(i).EventTimes.Cue = tr.DigitalEvents.CueOn;
                            end
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
                            if isfield(tr.DigitalEvents, 'LampOn')
                                obj(i).EventTimes.LightOn = tr.DigitalEvents.LampOn;
                                obj(i).EventTimes.LightOff = tr.DigitalEvents.LampOff;
                            else
                                obj(i).EventTimes.LightOn = [];
                                obj(i).EventTimes.LightOff = [];
                            end
                            if isfield(tr.DigitalEvents, 'Motor1On')
                                obj(i).EventTimes.Motor1On = tr.DigitalEvents.Motor1On;
                                obj(i).EventTimes.Motor1Off = tr.DigitalEvents.Motor1Off;
                                obj(i).EventTimes.Motor2On = tr.DigitalEvents.Motor2On;
                                obj(i).EventTimes.Motor2Off = tr.DigitalEvents.Motor2Off;
                                obj(i).EventTimes.MotorBusyOn = tr.DigitalEvents.MotorBusyOn;
                                obj(i).EventTimes.MotorBusyOff = tr.DigitalEvents.MotorBusyOff;
                            else
                                obj(i).EventTimes.Motor1On = [];
                                obj(i).EventTimes.Motor1Off = [];
                                obj(i).EventTimes.Motor2On = [];
                                obj(i).EventTimes.Motor2Off = [];
                                obj(i).EventTimes.MotorBusyOn = [];
                                obj(i).EventTimes.MotorBusyOff = [];
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
                            obj(i).Trials.Light = obj(i).makeTrials('light');
    
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
        
                            if ~isempty(p.Results.savepath)
                                obj(i).save(p.Results.savepath);
                            end
                        
                            i = i + 1;
                        catch ME
                            warning('Error when processing unit %i (%s) - this one will be skipped.', i, e.ExpName)
                            warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
                            i = i + 1;
                        end			
                    end
%                     fprintf(1, 'Done (%.1f sec).\n', toc(tTic))
                    clear tr
                end
                fprintf(1, 'Done.\n')
            else
                ar = p.Results.AR;
                if ~isempty(ar.bsr)
                    BSR = ar.bsr;
                elseif ~isempty(ar.bmrPress)
                    BSR = ar.bmrPress;
                else
                    error('Both bsr and bmrPress are empty structs. What (unit list) are you trying to pull here?');
                end
                obj(length(BSR)) = EphysUnit();
                expName = ar.expName;
                if ~isempty(sessionInfo)
                    pm = ProbeMap(sessionInfo.ml, sessionInfo.ap, sessionInfo.dv, duraOffset=sessionInfo.duraOffset, facing=sessionInfo.facing, model=sessionInfo.model);
                else
                    pm = [];
                end
                fprintf(1, 'Loading experiment %s...\n', expName)
                if ~isempty(p.Results.tr)
                    tr = p.Results.tr;
                else
                    tr = TetrodeRecording.BatchLoadSimple(expName, true);
                end
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

                        if ~isempty(p.Results.whitelist) && ~ismember(obj(i).getName('_'), p.Results.whitelist)
                            fprintf(1, 'Unit %g (%s) not found in whitelist, will not be processed/saved.\n', i, obj(i).getName('_'))
                            i = i + 1;
                            continue
                        end

                        s = tr.Spikes(obj(i).Channel);
                        inCluster = s.Cluster.Classes == bsr.unit;
                        obj(i).SpikeTimes = s.Timestamps(inCluster);
                        if p.Results.readWaveforms
                            obj(i).Waveforms = int16(s.Waveforms(inCluster, :));
                            obj(i).WaveformTimestamps = s.WaveformTimestamps;
                        end
                        % Spontaneous
%                         if isfield(tr.DigitalEvents, 'TRIAL_START')
%                             obj(i).EventTimes.Cue = tr.DigitalEvents.TRIAL_START;
%                         else
%                             obj(i).EventTimes.Cue = tr.DigitalEvents.CueOn;
%                         end
                        if p.Results.isSpontaneousPress
                            obj(i).EventTimes.Cue = tr.DigitalEvents.TRIAL_START;
                            obj(i).EventTimes.PressOff = tr.DigitalEvents.PressOff;
                        elseif isfield(tr.DigitalEvents, 'CueOn')
                            obj(i).EventTimes.Cue = tr.DigitalEvents.CueOn;
                        elseif isfield(tr.DigitalEvents, 'CueLeftOn')
                            obj(i).EventTimes.CueLeftOn = tr.DigitalEvents.CueLeftOn;
                            obj(i).EventTimes.CueLeftOff = tr.DigitalEvents.CueLeftOff;
                            obj(i).EventTimes.CueRightOn = tr.DigitalEvents.CueRightOn;
                            obj(i).EventTimes.CueRightOff = tr.DigitalEvents.CueRightOff;
                        end
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
                        if isfield(tr.DigitalEvents, 'LampOn')
                            obj(i).EventTimes.LightOn = tr.DigitalEvents.LampOn;
                            obj(i).EventTimes.LightOff = tr.DigitalEvents.LampOff;
                        else
                            obj(i).EventTimes.LightOn = [];
                            obj(i).EventTimes.LightOff = [];
                        end
                        if isfield(tr.DigitalEvents, 'Motor1On')
                            obj(i).EventTimes.Motor1On = tr.DigitalEvents.Motor1On;
                            obj(i).EventTimes.Motor1Off = tr.DigitalEvents.Motor1Off;
                            obj(i).EventTimes.Motor2On = tr.DigitalEvents.Motor2On;
                            obj(i).EventTimes.Motor2Off = tr.DigitalEvents.Motor2Off;
                            obj(i).EventTimes.MotorBusyOn = tr.DigitalEvents.MotorBusyOn;
                            obj(i).EventTimes.MotorBusyOff = tr.DigitalEvents.MotorBusyOff;
                        else
                            obj(i).EventTimes.Motor1On = [];
                            obj(i).EventTimes.Motor1Off = [];
                            obj(i).EventTimes.Motor2On = [];
                            obj(i).EventTimes.Motor2Off = [];
                            obj(i).EventTimes.MotorBusyOn = [];
                            obj(i).EventTimes.MotorBusyOff = [];
                        end

                        if isfield(tr.DigitalEvents, 'Mot1HiOn')
                            obj(i).EventTimes.Mot1HiOn = tr.DigitalEvents.Mot1HiOn;
                            obj(i).EventTimes.Mot1HiOff = tr.DigitalEvents.Mot1HiOff;
                            obj(i).EventTimes.Mot1LoOn = tr.DigitalEvents.Mot1LoOn;
                            obj(i).EventTimes.Mot1LoOff = tr.DigitalEvents.Mot1LoOff;
                            obj(i).EventTimes.Mot1BusyOn = tr.DigitalEvents.Mot1BusyOn;
                            obj(i).EventTimes.Mot1BusyOff = tr.DigitalEvents.Mot1BusyOff;
                            obj(i).EventTimes.Mot2BusyOn = tr.DigitalEvents.Mot2BusyOn;
                            obj(i).EventTimes.Mot2BusyOff = tr.DigitalEvents.Mot2BusyOff;
                        else
                            obj(i).EventTimes.Motor1On = [];
                            obj(i).EventTimes.Motor1Off = [];
                            obj(i).EventTimes.Motor2On = [];
                            obj(i).EventTimes.Motor2Off = [];
                            obj(i).EventTimes.MotorBusyOn = [];
                            obj(i).EventTimes.MotorBusyOff = [];
                        end

                        if p.Results.cullITI
                            obj(i).ExtendedWindow = extendedWindow;
                        else
                            obj(i).ExtendedWindow = [];
                        end
                        
                        % Make trials
                        if p.Results.isSpontaneousPress
                            obj(i).Trials.Press = obj(i).makeTrials('press_spontaneous');
                            obj(i).Trials.PressCorrect = obj(i).makeTrials('press_spontaneous_correct');
                        else
                            obj(i).Trials.Press = obj(i).makeTrials('press');
                        end
                        obj(i).Trials.Lick = obj(i).makeTrials('lick');
                        obj(i).Trials.Stim = obj(i).makeTrials('stim');
                        obj(i).Trials.StimTrain = obj(i).makeTrials('stimtrain');
                        obj(i).Trials.Light = obj(i).makeTrials('light');

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
                        obj(i).ProbeMap = pm;

                        if ~isempty(p.Results.savepath)
                            obj(i).save(p.Results.savepath);
                        end
                        
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
            addOptional(p, 'sep', '_', @ischar);
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
            p.addRequired('trialType', @(x) all(ismember(x, {'press', 'lick', 'stim', 'stimtrain', 'stimfirstpulse', 'light', 'anylick', 'firstlick', 'circlick', 'lickbout', 'stimtwocolor', 'press_spontaneous', 'press_spontaneous2', 'press_spontaneous_medial', 'press_spontaneous_lateral'})));
            p.addOptional('sorted', true, @islogical);
            p.addParameter('minBoutCycles', 2)
            p.addParameter('maxBoutCycles', 4)
            p.addParameter('minInterval', 0.05);
            p.addParameter('maxInterval', 0.25);
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
                        case 'press_spontaneous'
                            trials{itt} = obj.Trials.PressSpontaneous(:);
                        case 'press_spontaneous_medial'
                            trials{itt} = obj.Trials.PressSpontaneousMedial(:);
                        case 'press_spontaneous_lateral'
                            trials{itt} = obj.Trials.PressSpontaneousLateral(:);
                        case 'press_spontaneous2'
                            trials{itt} = obj.makeTrials('press_spontaneous2');
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
                        case 'stimtwocolor'
                            trials{itt} = obj.makeTrials('stim');
                        case 'light'
                            trials{itt} = obj.Trials.Light(:);
                        case 'anylick'
                            trials{itt} = Trial(obj.EventTimes.Lick-0.001, obj.EventTimes.Lick);
                        case 'firstlick'
                            trials{itt} = obj.makeTrials('firstlick');
                        case 'circlick'
                            trials{itt} = obj.makeTrials('circlick');
                        case 'lickbout'
                            trials{itt} = obj.makeTrials('lickbout', ...
                                minBoutCycles=p.Results.minBoutCycles, maxBoutCycles=p.Results.maxBoutCycles, ...
                                minInterval=p.Results.minInterval, maxInterval=p.Results.maxInterval);
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
            p.addParameter('correction', [], @isnumeric)
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
                    
                    [xx, tt] = obj(i).getTrialAlignedData(data, window, trialType, allowedTrialDuration=[edges(iBin), edges(iBin+1)], alignTo=alignTo, resolution=resolution, includeInvalid=false, correction = p.Results.correction);

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
            p.addRequired('event', @(x) ischar(x) && ismember(lower(x), {'press', 'lick', 'stim', 'stimtrain', 'stimfirstpulse', 'stimtwocolor', 'anylick', ...
                'firstlick', 'circlick', 'lickbout', 'press_spontaneous', 'press_spontaneous_medial', 'press_spontaneous_lateral'}))
            p.addOptional('window', [-2, 0], @(x) isnumeric(x) && length(x)>=2 && x(2) > x(1))
            p.addParameter('minTrialDuration', 0, @(x) isnumeric(x) && length(x)==1 && x>=0)
            p.addParameter('maxTrialDuration', Inf, @(x) isnumeric(x) && length(x)==1 && x>=0)
            p.addParameter('findSingleTrialDuration', 'off', @(x) ismember(x, {'off', 'min', 'max'}))
            p.addParameter('resolution', [], @(x) isnumeric(x) && length(x)==1 && x>=0)
            p.addParameter('normalize', 'none', @(x) isstruct(x) || isnumeric(x) || ismember(lower(x), {'none', 'iti', 'all'}))
            p.addParameter('alignTo', 'default', @(x) ismember(x, {'default', 'start', 'stop'}))
            p.addParameter('includeInvalid', [], @islogical)
            p.addParameter('correction', [], @isnumeric)
            p.addParameter('trials', [], @(x) isempty(x) || iscell(x) || isa(x, 'Trial'))
            p.parse(data, event, varargin{:})
            data = lower(p.Results.data);
            event = p.Results.event;
            window = p.Results.window;
            minTrialDuration = p.Results.minTrialDuration;
            maxTrialDuration = p.Results.maxTrialDuration;
            findSingleTrialDuration = p.Results.findSingleTrialDuration;
            resolution = p.Results.resolution;
            normalize = p.Results.normalize;
            alignTo = p.Results.alignTo;
            includeInvalid = p.Results.includeInvalid;
            correction = p.Results.correction;
            trials = p.Results.trials;
            
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
            D = zeros(length(obj), 1);

            tTic = tic();
            if length(obj) > 1
                fprintf(1, 'Processing ETA for %i units...', length(obj))
            end
            for i = 1:length(obj)
                if strcmpi(alignTo, 'default') 
                    if ismember(lower(event), {'stim', 'stimtrain', 'stimfirstpulse'})
                        alignTo = 'start';
                    else
                        alignTo = 'stop';
                    end
                end
                if isempty(includeInvalid)
                    if ismember(event, {'stimtrain', 'anylick', 'firstlick'})
                        includeInvalid = true;
                    else
                        includeInvalid = false;
                    end
                end

                if iscell(trials)
                    theseTrials = trials{i};
                    if isempty(theseTrials)
                        X(i, :) = NaN;
                        N(i) = 0;
                        continue                        
                    end
                else
                    theseTrials = trials;
                end

                [x, ~, d] = obj(i).getTrialAlignedData(data, window, event, trials=theseTrials, alignTo=alignTo, allowedTrialDuration=[minTrialDuration, maxTrialDuration], ...
                    findSingleTrialDuration=findSingleTrialDuration, resolution=resolution, includeInvalid=includeInvalid, correction=correction);
                
                if isempty(x)
                    % warning('EventTriggeredAverage cannot be calculated for Unit %i (%s), likely because trial count is zero.', i, obj(i).getName())
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
                D(i) = d;
%                 fprintf(1, ' %g', i)
            end
            if length(obj) > 1
                fprintf(1, 'Done (%.1f sec)\n', toc(tTic))
            end
    
            if isnumeric(normalize)
                eta = struct('X', X, 't', t, 'N', N, 'D', D, 'stats', stats);
            else
                eta = struct('X', X, 't', t, 'N', N, 'D', D);
            end

%             if ~strcmpi(findSingleTrialDuration, 'off')
%                 error('Not implemented');
%             end
        end

        function rd = getRasterData(obj, trialType, varargin)
            p = inputParser();
            p.addRequired('trialType', @(x) all(ismember(x, {'press', 'lick', 'stim', 'stimtrain', 'stimfirstpulse', 'stimtwocolor', 'press_spontaneous', 'press_spontaneous2'})))
            p.addOptional('window', [0, 0], @(x) isnumeric(x) && length(x) >= 2 && x(1) <= 0 && x(2) >= 0)
            p.addParameter('minTrialDuration', 0, @(x) isnumeric(x) && length(x)==1 && x>=0)
            p.addParameter('maxTrialDuration', Inf, @(x) isnumeric(x) && length(x)==1 && x>=0)
            p.addParameter('durErr', 1e-3, @isnumeric)
            p.addParameter('alignTo', 'default', @(x) ismember(x, {'default', 'start', 'stop'}))
            p.addParameter('sort', true, @islogical);
            p.addParameter('trials', [], @(x) isa(x, 'Trial'))
            p.addParameter('shutterDelay', 0, @isnumeric)
            p.parse(trialType, varargin{:})
            trialType = p.Results.trialType;
            window = p.Results.window;
            minTrialDuration = p.Results.minTrialDuration;
            maxTrialDuration = p.Results.maxTrialDuration;
            durErr = p.Results.durErr;
            alignTo = p.Results.alignTo;

            if strcmp(alignTo, 'default')
                switch trialType
                    case {'press', 'lick', 'press_spontaneous'}
                        alignTo = 'stop';
                    case {'stim', 'stimtrain', 'stimfirstpulse', 'stimtwocolor'}
                        alignTo = 'start';
                end
            end

            if length(obj) == 1
                if isempty(p.Results.trials)
                    trials = obj.getTrials(trialType, sorted=true);
                else
                    trials = p.Results.trials;
                end
                dur = round(trials.duration ./ durErr) * durErr;
                trials = trials(dur >= minTrialDuration & dur <= maxTrialDuration);
                [~, t, I] = trials.inTrial(obj.SpikeTimes, window);
                switch alignTo
                    case 'start'
                        tRef = [trials.Start];
                    case 'stop'
                        tRef = [trials.Stop];
                end
                t = t - tRef(I);
                dur = round(trials.duration ./ durErr) * durErr;

                if ~strcmpi(trialType, 'stimfirstpulse')
                    iti = trials.iti();
                else
%                     For stimfirstpulse (in train), use the iti between
%                     first and second pulse in train (rather that iti
%                     between first pulse in train and first pulse in next
%                     train
%                     allTrials = obj.getTrials('stim', true);
%                     itiAll = allTrials.iti();
%                     startAll = [allTrials.Start];
%                     startFirst = [trials.Start];
% 
%                     [match, index] = ismember(round(startFirst*1000), round(startAll*1000));
% %                     assert(all(match));
%                     if ~(all(match))
%                         disp(index)
%                     end
%                     iti = itiAll(index);
                    iti = [];
                end

                if p.Results.sort
                    if strcmpi(trialType, 'stimtwocolor')
                        [tce, stimOn, stimOff, trainIndices] = obj.LoadTwoColorExperiment();
                        stimOn = stimOn + p.Results.shutterDelay;
                        stimOff = stimOff + p.Results.shutterDelay;
                        nPulses = length(stimOn);
                        trainHash = tce.getStimHash();
                        pulseHash = trainHash(trainIndices);

                        [~, trainOrder] = sort(trainHash, 'ascend');
                        [~, pulseOrder] = sort(pulseHash, 'ascend');
                        assert(length(pulseOrder) == max(pulseOrder) && min(pulseOrder) == 1)

                        trialIndices = NaN(size(obj.SpikeTimes));
                        spikesRelative = NaN(size(obj.SpikeTimes));
                        for iPulse = 1:nPulses
                            sel = obj.SpikeTimes >= stimOn(iPulse) + window(1) & obj.SpikeTimes <= stimOn(iPulse) + window(2);
                            trialIndices(sel) = iPulse;
                            spikesRelative(sel) = obj.SpikeTimes(sel) - stimOn(iPulse);
                        end
                        trialIndicesSorted = changem(trialIndices, 1:length(pulseOrder), pulseOrder);
                        
                        I = trialIndicesSorted;
                        t = spikesRelative;
                        sel = ~isnan(I);
                        I = I(sel);
                        t = t(sel);
                        dur = dur(pulseOrder);
                        iti = iti(pulseOrder);
                        rd.tce.stimLog = tce.Log;
                        rd.tce.trainOrder = trainOrder(:)';
                        rd.tce.pulseOrder = pulseOrder(:)';
                        rd.tce.trainIndices = trainIndices(:)';
                        rd.tce.params = tce.Params;
                    else
                        [dur, ISort] = sort(dur, 'ascend');
                        I = changem(I, 1:length(ISort), ISort);
                        iti = iti(ISort);
                    end
                end
                rd.name = obj.getName('_');
                rd.trialType = trialType;
                rd.alignTo = alignTo;
                rd.t = t;
                rd.I = I;
                rd.duration = dur;
                rd.iti = iti;
            else
                tTic = tic();
                if strcmpi(trialType, 'stimtwocolor')
                    rd(length(obj)) = struct('name', '', 'trialType', '', 'alignTo', '', 't', [], 'I', [], 'duration', [], 'iti', [], 'tce', []);
                else
                    rd(length(obj)) = struct('name', '', 'trialType', '', 'alignTo', '', 't', [], 'I', [], 'duration', [], 'iti', []);
                end
                fprintf(1, 'Calculating raster data for %g units...\n', length(obj));
                for i = 1:length(obj)
                    try
                        rd(i) = obj(i).getRasterData(trialType, window, ...
                            minTrialDuration=minTrialDuration, alignTo=alignTo, sort=p.Results.sort, ...
                            maxTrialDuration=maxTrialDuration, durErr=durErr, shutterDelay=p.Results.shutterDelay);
                    catch ME
                        warning('\tError when processing unit %g', i)
                        warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
                    end
                end
                fprintf('Done and it only took %.2f sec. Easy.\n', toc(tTic));
            end
        end

        function [isi, t] = getMeanPEISI(obj, trialType, trials, varargin)
            p = inputParser();
            p.addRequired('trialType', @(x) all(ismember(x, {'press', 'lick', 'stim', 'stimtrain', 'stimfirstpulse', 'stimtwocolor'})))
            p.addRequired('trials', @(x) isa(x, 'Trial'))
            p.addOptional('window', [-0.5, 0.5], @(x) isnumeric(x) && length(x) >= 2 && x(1) <= 0 && x(2) >= 0)
            p.addParameter('resolution', 1e-3, @isnumeric)
            p.addParameter('alignTo', 'default', @(x) ismember(x, {'default', 'start', 'stop'}))
            p.addParameter('shutterDelay', 0, @isnumeric)
            p.parse(trialType, trials, varargin{:})
            trialType = p.Results.trialType;
            trials = p.Results.trials;
            window = p.Results.window;
            resolution = p.Results.resolution;
            alignTo = p.Results.alignTo;
            shutterDelay = p.Results.shutterDelay;

            assert(length(obj) == 1)

            rd = obj.getRasterData(trialType, window, trials=trials, alignTo=alignTo, shutterDelay=shutterDelay, sort=false);
            isFirstSpikeInTrial = logical([1, diff(rd.I)]);
            start = [1, strfind(isFirstSpikeInTrial, [0, 1]) + 1];
            stop = [strfind(isFirstSpikeInTrial, [0, 1]), length(rd.I)];

            nTrials = length(start);
            isi = [NaN, diff(rd.t)];
            isi(isFirstSpikeInTrial) = NaN;

            t = window(1):resolution:window(2);
            ISI = NaN(nTrials, length(t));
            for iTrial = 1:nTrials
                sel = start(iTrial) + 1:stop(iTrial);
                try
                    ISI(iTrial, :) = interp1(rd.t(sel), isi(sel), t, 'linear');
                end
            end
            isi = mean(ISI, 1, 'omitnan');
        end

        function [r, lags, pairedIndices] = xcorr(obj, data, varargin)
            %% [r, lags] = obj.XCORR(_), returns full cross correlation as a function of lag, for all pairs of objects in obj array.
            p = inputParser();
            p.addOptional('data', 'rate', @(x) ismember(x, {'count', 'rate'}));
            p.addParameter('resolution', 0.001, @isnumeric);
            p.addParameter('normalize', true, @islogical)
            p.addParameter('maxlag', NaN, @isnumeric);
            p.addParameter('useCorr', false, @islogical)
            p.parse(data, varargin{:});
            data = p.Results.data;
            resolution = p.Results.resolution;
            maxlag = p.Results.maxlag;

            edges = min(arrayfun(@(obj) obj.SpikeTimes(1), obj)):resolution:max(arrayfun(@(obj) obj.SpikeTimes(end), obj));

            assert(length(obj) > 1)

            % Calculate spike rates/counts for all objects in array (saves
            % time when there are thousands of pairs)
            X = zeros(length(obj), length(edges) - 1);
            fprintf('\tCalculating binned spike counts for %i units...', length(obj))
            for i = 1:length(obj)
                switch data
                    case 'count'
                        X(i, :) = obj(i).getSpikeCounts(edges);
                    case 'rate'
                        kernel = obj(i).SpikeRateKernel;
                        width = kernel.params.width;
                        if strcmpi(kernel.type, 'gaussian')
                            sigma = kernel.params.sigma;
                            X(i, :) = obj(i).getSpikeRates('gaussian', sigma, edges, kernelWidth=width);
                        else
                            lambda1 = kernel.params.lambda1;
                            lambda2 = kernel.params.lambda2;
                            X(i, :) = obj(i).getSpikeRates('exponential', lambda1, lambda2, edges, kernelWidth=width);
                        end
                end
            end
            fprintf('Done.\n')

            % Optionally normalize (zscore) data before calculating
            % xcorr, this helps get rid of the triagle shape in r vs.
            % lag plots (due to zero padding for missing data?, so centering on zero hides this?).
            if p.Results.normalize
                fprintf('\tNormalizing spike counts...')
                X = normalize(X, 2, 'zscore', 'std');
                scaleopt = 'normalized';
                fprintf('Done.\n')
            else
                scaleopt = 'none';
            end

            % Calculate pair-wise xcorr for all units in obj array.
            pairedIndices = nchoosek(1:length(obj), 2);
            nPairs = size(pairedIndices, 1);
            r = cell(nPairs, 1);
            lags = cell(nPairs, 1);
            fprintf('\tCalculating %i pairwise cross-correlations...', nPairs)
            for iPair = 1:nPairs
                i = pairedIndices(iPair, 1);
                j = pairedIndices(iPair, 2);
    
                if ~p.Results.useCorr
                    if isinf(maxlag) || isnan(maxlag)
                        [r{iPair}, lags{iPair}] = xcorr(X(i, :), X(j, :), scaleopt);
                    else
                        [r{iPair}, lags{iPair}] = xcorr(X(i, :), X(j, :), maxlag/resolution, scaleopt);
                    end
                else
                    lags{iPair} = 0;
                    r{iPair} = corr(X(i, :)', X(j, :)');
                end
            end
            fprintf('Done.\n')

            r = cat(1, r{:});
            lags = cat(1, lags{:});
        end

        function [obj, isDuplicate] = removeDuplicates(obj, threshold)
            if nargin < 2
                threshold = 0.7;
            end

            if length(obj) <= 1
                return
            end

            [uniqueExpNames, ia, ic] = unique({obj.ExpName});
            assert(all(diff(ic) >= 0)) % So we can use ia as bin edges (where each bin is an experiment)
            isDuplicate = cell(length(uniqueExpNames), 1);
            ia(end + 1) = length(obj) + 1;
            for iExp = 1:length(uniqueExpNames)
                fprintf('Removing duplicate units in session %i of %i:\n', iExp, length(uniqueExpNames))
                iStart = ia(iExp);
                iEnd = ia(iExp + 1) - 1;
                theseObjs = obj(iStart:iEnd);
                [r, ~, pairIndices] = theseObjs.xcorr('count', resolution=0.002, normalize=false, useCorr=true);

                duplicatePairs = pairIndices(r > threshold, :);
                isDup = false(length(theseObjs), 1);
                for iPair = 1:size(duplicatePairs, 1)
                    assert(duplicatePairs(iPair, 2) > duplicatePairs(iPair, 1))
                    isDup(duplicatePairs(iPair, 2)) = true;
                end
                isDuplicate{iExp} = isDup;
                fprintf('Removed %i duplicate units out of %i.\n', nnz(isDup), length(theseObjs))
            end

            isDuplicate = cat(1, isDuplicate{:});
            obj = obj(~isDuplicate);
        end

        function [obj, isMultiUnit] = removeMultiUnits(obj, varargin)
            p = inputParser();
            p.addParameter('maxISI', 1.5e-3, @isnumeric);
            p.addParameter('maxFraction', 0.05, @isnumeric);
            p.addParameter('cullZeros', true, @islogical); % Doing this a second time should not have any additional effect
            p.parse(varargin{:});
            maxISI = p.Results.maxISI;
            maxFraction = p.Results.maxFraction;
            cullZeros = p.Results.cullZeros;

            isMultiUnit = false(size(obj));
            for i = 1:length(obj)
                if cullZeros
                    isi = [Inf, diff(obj(i).SpikeTimes)];
                    toCull = isi == 0;
                    obj(i).SpikeTimes(toCull) = [];
                    if ~isempty(obj(i).Waveforms)
                        obj(i).Waveforms(toCull) = [];
                    end
                    % Hey, remember to redo spikecounts/spikerates after this
                end
                isi = [Inf, diff(obj(i).SpikeTimes)];
                isMultiUnit(i) = (nnz(isi < maxISI) ./ length(obj(i).SpikeTimes)) > maxFraction;
            end
            obj = obj(~isMultiUnit);
        end

        function [pos, mot1, mot2, motBusy] = getMotorState(obj, t)
            name1On = 'Motor1On';
            name1Off = 'Motor1Off';
            name2On = 'Motor2On';
            name2Off  = 'Motor2Off';
            if isfield(obj.EventTimes, 'Mot1HiOn')
                name1On = 'Mot1HiOn';
                name1Off = 'Mot1HiOff';
                name2On = 'Mot1LoOn';
                name2Off  = 'Mot1LoOff';
            end
            assert(isfield(obj.EventTimes, 'Motor1On'))
            assert(length(obj.EventTimes.(name1On)) == length(obj.EventTimes.(name1Off)))
            assert(length(obj.EventTimes.(name2On)) == length(obj.EventTimes.(name2Off)))
            assert(length(obj.EventTimes.MotorBusyOn) == length(obj.EventTimes.MotorBusyOff))
            assert(all((obj.EventTimes.(name1Off) - obj.EventTimes.(name1On)) > 0))
            assert(all((obj.EventTimes.(name2Off) - obj.EventTimes.(name2On)) > 0))
            assert(all((obj.EventTimes.MotorBusyOff - obj.EventTimes.MotorBusyOn) > 0))

            if length(t) > 1
                n = length(t);
                pos = zeros(n, 1);
                mot1 = false(n, 1);
                mot2 = false(n, 1);
                motBusy = false(n, 1);
                for i = 1:n
                    [pos(i), mot1(i), mot2(i), motBusy(i)] = obj.getMotorState(t(i));
                end
                return
            end

            mot1 = false;
            mot2 = false;
            motBusy = false;
            for i = 1:length(obj.EventTimes.(name1On))
                if t >= obj.EventTimes.(name1On)(i) && t <= obj.EventTimes.(name1Off)(i)
                    mot1 = true;
                    break
                end
            end
            for i = 1:length(obj.EventTimes.(name2On))
                if t >= obj.EventTimes.(name2On)(i) && t <= obj.EventTimes.(name2Off)(i)
                    mot2 = true;
                    break
                end
            end
            for i = 1:length(obj.EventTimes.MotorBusyOn)
                if t >= obj.EventTimes.MotorBusyOn(i) && t <= obj.EventTimes.MotorBusyOff(i)
                    motBusy = true;
                    break
                end
            end

            if ~mot1 && ~mot2
                pos = 1;
            elseif ~mot1 && mot2
                pos = 2;
            elseif mot1 && ~mot2
                pos = 3;
            else
                pos = 4;
            end
        end

        function [tce, tOn, tOff, trainIndices] = LoadTwoColorExperiment(obj, varargin)
            assert(length(obj) == 1)
            p = inputParser();
            p.addParameter('pulseWidthErrorMargin', 1e-3, @isnumeric)
            p.parse(varargin{:})

            expName = obj.ExpName;
            animalName = obj.getAnimalName;
            file = dir(sprintf('Z:\\%s\\%s\\%s.mat', animalName, expName, expName));
            assert(~isempty(file))

            tce = load(sprintf('%s\\%s', file.folder, file.name));
            tce = tce.obj;
            assert(isa(tce, 'TwoColorExperiment'))
            clear file

            tOn = obj.EventTimes.StimOn;
            tOff = obj.EventTimes.StimOff;

            nPulsesPerTrain = arrayfun(@(log) log.params.nPulses, tce.Log);
            pulseWidth = arrayfun(@(log) repmat(log.params.pulseWidth, [1, log.params.nPulses]), tce.Log, UniformOutput=false);
            pulseWidth = cat(2, pulseWidth{:});
            
            % Verify congruence between Intan shutter events and TwoColorExperiment.Log
            assert(sum(nPulsesPerTrain) == nnz(tOn), 'Intan recorded %i shutterOn events while TwoColorExperiment.Log has %i.', nnz(tOn), sum(nPulsesPerTrain));
            assert(sum(nPulsesPerTrain) == nnz(tOff), 'Intan recorded %i shutterOn events while TwoColorExperiment.Log has %i.', nnz(tOn), sum(nPulsesPerTrain));
            nPulses = nnz(tOn);
            % fprintf('Intan & TwoColorExperiment.Log both recorded %i shutter pulses.\n', nnz(tOn))
            
            % Pulse durations should match between Intan and TCE as well.
            assert(all(abs(tOff - tOn - pulseWidth) < p.Results.pulseWidthErrorMargin), 'Only %i/%i pulseWidths agree (df<%gs).', ...
                nnz(abs(tOff - tOn - pulseWidth) < p.Results.pulseWidthErrorMargin), nPulses, p.Results.pulseWidthErrorMargin)
            % fprintf('All %i pulseWidths agree (df<%gs).\n', nPulses, p.Results.pulseWidthErrorMargin)
            
            % Generate pulse->train map
            iPulse = 0;
            trainIndices = NaN(1, nPulses);
            for iTrain = 1:length(tce.Log)
                trainIndices(iPulse + 1:iPulse + nPulsesPerTrain(iTrain)) = iTrain;
                iPulse = iPulse + nPulsesPerTrain(iTrain);
            end
        end

        function groups = groupTwoColorStimTrials(obj, varargin)
            p = inputParser();
            p.addOptional('groupBy', {'wavelength', 'location', 'duration', 'power'}, @(x) ~isempty(x) && all(ismember(x, {'wavelength', 'location', 'duration', 'power'})));
            p.parse(varargin{:});
            groupBy = p.Results.groupBy;

            assert(length(obj) == 1);

            [tce, ~, ~, trainIndices] = obj.LoadTwoColorExperiment();
            [groupIndices, conditions] = tce.groupStimTrains(groupBy);

            nGroups = length(conditions);
            nTrains = length(groupIndices);

            groups(nGroups).trials = [];
            groups(nGroups).pulseIndices = [];
            groups(nGroups).label = '';
            for i = 1:length(groupBy)
                switch lower(groupBy{i})
                    case 'wavelength'
                        groups(nGroups).wavelength = [];
                    case 'location'
                        groups(nGroups).location = [];
                    case 'duration'
                        groups(nGroups).duration = [];
                    case 'power'
                        groups(nGroups).power = [];
                end
            end
            trials = obj.getTrials('stim');
            trials = trials(:)';
            for iGroup = 1:nGroups
                groups(iGroup).pulseIndices = find(ismember(trainIndices, find(groupIndices == iGroup)));
                groups(iGroup).trials = trials(groups(iGroup).pulseIndices);
                for i = 1:length(groupBy)
                    switch lower(groupBy{i})
                        case 'wavelength'
                            groups(iGroup).wavelength = conditions(iGroup).wavelength;
                            groups(iGroup).label = sprintf('%s %inm', groups(iGroup).label, groups(iGroup).wavelength);
                        case 'location'
                            groups(iGroup).location = conditions(iGroup).location;
                            groups(iGroup).label = sprintf('%s %i', groups(iGroup).label, groups(iGroup).location);
                        case 'duration'
                            groups(iGroup).duration = conditions(iGroup).duration;
                            groups(iGroup).label = sprintf('%s %gms', groups(iGroup).label, 1e3*groups(iGroup).duration);
                        case 'power'
                            groups(iGroup).power = conditions(iGroup).power;
                            groups(iGroup).label = sprintf('%s %guW', groups(iGroup).label, 1e6*groups(iGroup).power);
                    end
                end
            end
        end

        
    end
    
    % static methods
    methods (Static)
        function obj = load(varargin)
            p = inputParser();
            if ischar(varargin{1})
                p.addOptional('path', 'C:\SERVER\Units', @isfolder);
            elseif iscell(varargin{1})
                p.addOptional('files', {}, @iscell);
            end
            p.addParameter('waveforms', true, @islogical); % FALSE to cull waveforms after loading to save memory
            p.addParameter('spikecounts', true, @islogical); % FALSE to cull spikecounts after loading to save memory
            p.addParameter('spikerates', true, @islogical); % FALSE to cull spikerate after loading to save memory
            p.addParameter('animalNames', {}, @iscell);
            p.parse(varargin{:});

            if isfield(p.Results, 'path')
                files = dir(sprintf('%s\\*.mat', p.Results.path));
                files = arrayfun(@(x) sprintf('%s\\%s', x.folder, x.name), files, UniformOutput=false);
            elseif isfield(p.Results, 'files')
                files = p.Results.files;
            end
            if ~isempty(p.Results.animalNames)
                files = files(contains(files, p.Results.animalNames, IgnoreCase=true));
            end

            S(length(files)) = struct('eu', []);
            for i = 1:length(files)
                tTic = tic();
                fprintf(1, 'Reading unit %g/%g...', i, length(files));
                S(i) = load(files{i}, 'eu');
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
            p.addParameter('event', 'event', @ischar);
            p.addParameter('order', [], @isnumeric);
            p.addParameter('clim', [], @isnumeric)
            p.addParameter('xlim', [], @isnumeric)
            p.addParameter('sortWindow', [-2, 0], @(x) isnumeric(x) && length(x) == 2)
            p.addParameter('signWindow', [-.5, 0], @(x) isnumeric(x) && length(x) == 2)
            p.addParameter('sortThreshold', 1, @isnumeric)
            p.addParameter('negativeSortThreshold', [], @isnumeric)
            p.addParameter('hidecolorbar', false, @islogical)
            p.addParameter('timeUnit', 's', @(x) ismember('s', 'ms'))
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
                f = figure('Units', 'normalized', 'OuterPosition', [0, 0, 0.3, 1], 'DefaultAxesFontSize', 12);
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
%                 etaSign(meta >= -negativeSortThreshold & meta <= 0) = -Inf;
%                 etaSign(meta <= sortThreshold & meta >= 0) = Inf;
                assert(size(etaSign, 2) == 1);
                isAboveThreshold = (XSort >= sortThreshold.*etaSign & etaSign > 0) | (XSort <= negativeSortThreshold.*etaSign & etaSign < 0);
                [~, iPeak] = max(flip(isAboveThreshold, 2), [], 2);
                nSignBins = nnz(t >= signWindow(1) & t <= signWindow(2) & t >= sortWindow(1) & t <= sortWindow(2));
                for i = 1:size(isAboveThreshold, 1)
                    if iPeak(i) <= nSignBins
                        isAboveThreshold(i, end-iPeak(i)+1:end) = true;
                    end
                end
                isAboveThresholdConseq = isAboveThreshold & [diff(isAboveThreshold')', zeros(size(X, 1), 1)] == 0;
                [~, Ilate] = min(flip(isAboveThresholdConseq, 2), [], 2, 'omitnan');
                Ilate = size(isAboveThresholdConseq, 2) + 2 - Ilate;
%                 [~, Ilate] = min(~isAboveThresholdConseq, [], 2, 'omitnan');
                
%                 min(Ilate)
%                 max(Ilate)

                nonSig = Ilate==size(isAboveThresholdConseq, 2) + 1;
                Ilate(nonSig) = size(isAboveThresholdConseq, 2);
%                 [~, order] = sort(Ilate .* etaSign + 0.1*abs(meta), 'descend');
                sortVal = Ilate .* etaSign;
                sortVal(nonSig) = sortVal(nonSig) + meta(nonSig);
                [~, order] = sort(sortVal, 'ascend');
                meta = meta(order);
                latency = tSort(Ilate);
                assert(all(~isnan(latency)))
%                 latency = latency(order);
                varargout = {ax, order(:), meta(:), latency(:)};
            else
                order = p.Results.order;
                if nargout < 3
                    varargout = {ax, order};
                else
                    error('Latency calc not implemented')
%                     XSort = X(:, t >= sortWindow(1) & t <= sortWindow(2));
%                     tSort = t(t >= sortWindow(1) & t <= sortWindow(2));
%                     meta = mean(X(:, t >= signWindow(1) & t <= signWindow(2)), 2, 'omitnan');
%                     etaSign = sign(meta);
%                     assert(size(etaSign, 2) == 1);
%                     isAboveThreshold = (XSort >= sortThreshold.*etaSign & etaSign > 0) | (XSort <= negativeSortThreshold.*etaSign & etaSign < 0);
%                     [~, Ilate] = max(isAboveThreshold, [], 2, 'omitnan');
%                     meta = meta(order);
%                     latency = tSort(Ilate);
%                     latency = latency(order);
%                     varargout = {ax, order(:), meta(:), latency(:)};
                end
            end
            
            switch p.Results.timeUnit
                case 's'
                    timescale = 1;
                case 'ms'
                    timescale = 1000;
            end
            imagesc(ax, t.*timescale, 1:length(N), X(order, :))
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
                delete(h)
%                 h.Visible = 'off';
            end
            xlabel(ax, sprintf('Time from %s (%s)', p.Results.event, p.Results.timeUnit))
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

            [~, order] = EphysUnit.plotETA(ax(1), eta1, r.sel, clim=r.clim, xlim=xlim1, sortWindow=r.sortWindow, signWindow=signWindow1, sortThreshold=r.sortThreshold, negativeSortThreshold=r.negativeSortThreshold, hidecolorbar=true);
            [~, ~] = EphysUnit.plotETA(ax(2), eta2, r.sel, order=order, clim=r.clim, xlim=xlim2, sortWindow=r.sortWindow, signWindow=signWindow2, sortThreshold=r.sortThreshold, negativeSortThreshold=r.negativeSortThreshold);

%             [~, order, meta1, latency1] = EphysUnit.plotETA(ax(1), eta1, r.sel, clim=r.clim, xlim=xlim1, sortWindow=r.sortWindow, signWindow=signWindow1, sortThreshold=r.sortThreshold, negativeSortThreshold=r.negativeSortThreshold, hidecolorbar=true);
%             [~, ~, meta2, latency2] = EphysUnit.plotETA(ax(2), eta2, r.sel, order=order, clim=r.clim, xlim=xlim2, sortWindow=r.sortWindow, signWindow=signWindow2, sortThreshold=r.sortThreshold, negativeSortThreshold=r.negativeSortThreshold);

%             meta = horzcat(meta1(:), meta2(:));
%             latency = horzcat(latency1(:), latency2(:));

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
            p.addParameter('nsigmas', 1, @isnumeric)
            p.addParameter('sem', false, @islogical)
            p.addParameter('showTrialNum', false, @islogical)
            p.addParameter('numFormat', '%.1f', @ischar)
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
            
            if p.Results.sem
                S = S./sqrt(N);
            end

            hold(ax, 'on')           

%             colors = 'rgcbmkrgcbmkrgcbmkrgcbmkrgcbmkrgcbmk';
            for iBin = 1:length(B)
                bin = B(iBin, :);
                t = T;% + bin(2);
                if p.Results.showTrialNum
                    h(iBin) = plot(ax, t, X(iBin, :), 'Color', hsl2rgb([0.7*(iBin-1)/(length(B)-1), 1, 0.4]), 'LineWidth', 1.5, 'DisplayName', sprintf(sprintf('%s-%ss, n=%i trials', p.Results.numFormat, p.Results.numFormat), bin(1), bin(2), N(iBin)));
                else
                    h(iBin) = plot(ax, t, X(iBin, :), 'Color', hsl2rgb([0.7*(iBin-1)/(length(B)-1), 1, 0.4]), 'LineWidth', 1.5, 'DisplayName', sprintf(sprintf('%s-%ss', p.Results.numFormat, p.Results.numFormat), bin(1), bin(2)));
                end
                if p.Results.nsigmas > 0
                    high = X(iBin, :) + p.Results.nsigmas*S(iBin, :); 
                    low = X(iBin, :) - p.Results.nsigmas*S(iBin, :);
                    selS = ~isnan(high);
                    patch(ax, [t(selS), flip(t(selS))], [low(selS), flip(high(selS))], hsl2rgb([0.7*(iBin-1)/(length(B)-1), 1, 0.4]), 'FaceAlpha', 0.1, 'EdgeColor', 'none')
                end
                %yl = ax.YLim;
                % plot(ax, [bin(2), bin(2)], yl, '--', 'Color', h(iBin).Color)
                %ylim(ax, yl);
            end
            hold(ax, 'off')
            xlim(ax, p.Results.xlim)
            if X(iBin, end) < X(iBin, 1)
                legend(ax, h, Location='southwest')
            else
                legend(ax, h, Location='northwest')
            end
                
        end
        
        function ax = plotRaster(varargin)
            p = inputParser();
            if isgraphics(varargin{1}, 'axes')
                p.addRequired('ax', @(x) isgraphics(x, 'axes'));
            end
            p.addRequired('rd', @isstruct);
            p.addParameter('xlim', [-6, 1], @(x) isnumeric(x) && length(x) == 2 && x(2) > x(1));
            p.addParameter('iti', false, @islogical);
            p.addParameter('timeUnit', 's', @(x) ismember(x, {'s', 'ms'}))
            p.addParameter('maxTrials', Inf, @isnumeric)
            p.addParameter('maxTrialsMethod', 'trim', @(x) ismember(lower(x), {'trim', 'randomsample', 'uniformsample'}))
            p.addParameter('sz', 2.5, @isnumeric)
            p.parse(varargin{:})
            rd = p.Results.rd;
            isSimpleStim = ismember(lower(rd.trialType), {'stim', 'stimtrain', 'stimfirstpulse'});
            isTwoColorStim = strcmpi(rd.trialType, 'stimtwocolor');
            timeUnit = p.Results.timeUnit;
            maxTrials = p.Results.maxTrials;

            assert(length(rd) == 1);

            if isfield(p.Results, 'ax')
                ax = p.Results.ax;
            else
%                 ax = axes(figure(Units='normalized', Position=[0.1, 0.1, 0.67, 0.5], DefaultAxesFontSize=18));
                ax = axes(figure(Units='inches', Position=[0, 0, 6.5, 4], DefaultAxesFontSize=9));
            end

            switch rd.alignTo
                case 'start'
                    eventName = rd.trialType;
                    if ~isSimpleStim && ~isTwoColorStim
                        refName = 'trial start';
                    else
                        refName = 'opto onset';
                    end
                case 'stop'
                    if ~isSimpleStim && ~isTwoColorStim
                        eventName = 'trial start';
                        refName = rd.trialType;
                    else
                        eventName = 'opto onset';
                        refName = 'opto offset';
                    end
            end

            hold(ax, 'on')
%             nTrials = length(unique(rd.I));
            nTrials = length(rd.duration);
            if strcmpi(timeUnit, 'ms')
                timescale = 1000;
            else
                timescale = 1;
            end
            sz = p.Results.sz;

            % Make sure rd is still value type, since we're about to make
            % local changes to it
            assert(isstruct(rd))
            if maxTrials < nTrials
                assert(mod(maxTrials, 1) == 0)
                switch lower(p.Results.maxTrialsMethod)
                    case 'trim'
                        sel = 1:maxTrials;
                    case 'uniformsample'
                        sel = floor(linspace(1, nTrials, maxTrials));
                    case 'randomsample'
                        sel = randsample(nTrials, maxTrials);
                    otherwise
                        error()
                end
                selSpike = ismember(rd.I, sel);
                selTrial = sel;
                rd.t = rd.t(selSpike);
                rd.I = rd.I(selSpike);
                uniqueI = unique(rd.I);
                rd.I = changem(rd.I, 1:length(uniqueI), uniqueI);
                rd.duration = rd.duration(selTrial);
                rd.iti = rd.iti(selTrial);
                nTrials = maxTrials;
            end

            switch rd.alignTo
                case 'start'
                    tEvent = rd.duration;
                case 'stop'
                    tEvent = -rd.duration;
            end

            if ~isSimpleStim && ~isTwoColorStim
%                 ax.Parent.Position(4) = ax.Parent.Position(4)*nTrials/300;
                h = gobjects(2, 1);
                h(1) = scatter(ax, rd.t .* timescale, rd.I, sz, 'k', 'filled', DisplayName='spikes');
                h(2) = scatter(ax, tEvent .* timescale, 1:length(tEvent), sz*2, 'r', 'filled', DisplayName=eventName);
            elseif isSimpleStim
                h = gobjects(2, 1);
                h(1) = scatter(ax, rd.t * timescale, rd.I, sz, 'k', 'filled', DisplayName='spikes');
                uniqueDurations = unique(rd.duration);
                if length(uniqueDurations) == 1
                    h(2) = patch(ax, [0, uniqueDurations, uniqueDurations, 0].*timescale, [0, 0, nTrials+1, nTrials+1], 'b', FaceAlpha=0.25, EdgeAlpha=0, DisplayName='opto');
                else
                    h(2) = plot(ax, [0, 0], [0, nTrials+1], 'b', DisplayName='opto on');
                    h(3) = scatter(ax, tEvent .* timescale, 1:nTrials, 5, 'r', 'filled', DisplayName='opto off');
                end
                if p.Results.iti
                    h(4) = scatter(ax, (tEvent + rd.iti) .* timescale, 1:nTrials, 5, 'g', 'filled', DisplayName='ITI end');
                end
            elseif isTwoColorStim
                h = gobjects(1, 1);
                h(1) = scatter(ax, rd.t * timescale, rd.I, sz, 'k', 'filled', DisplayName='spikes');
                stimLog = rd.tce.stimLog;
                trainOrder = rd.tce.trainOrder;
                pulseOrder = rd.tce.pulseOrder;
                trainIndices = rd.tce.trainIndices;
                nPulses = length(rd.duration);
                tceParams = rd.tce.params;
                % Draw a box to mark the stim window
                stimBoxDict = dictionary();
                for iTrain = 1:length(stimLog)
                    selPulse = trainIndices == iTrain;
                    iPulseStart = strfind(selPulse, [0, 1]) + 1;
                    iPulseEnd = strfind(selPulse, [1, 0]);
                    if isempty(iPulseStart)
                        iPulseStart = 1;
                    end
                    if isempty(iPulseEnd)
                        iPulseEnd = nPulses;
                    end
                    iPulseStart = find(pulseOrder == iPulseStart);
                    iPulseEnd = find(pulseOrder == iPulseEnd);
                    pulseWidth = stimLog(iTrain).params.pulseWidth;
                    switch stimLog(iTrain).wavelength
                        case 473
                            color = [0.2, 0.2, 0.8];
                        case 593
                            color = [0.8, 0.2, 0.2];
                    end
                    alpha = 0.2 + 0.6*((stimLog(iTrain).params.iPower - 1)./(length(tceParams.targetPowers) - 1));
                    partialHash = stimLog(iTrain).params.iPower + stimLog(iTrain).wavelength*10;
                    partialDesc = sprintf('%guW \t%inm', tceParams.targetPowers(stimLog(iTrain).params.iPower)*1e6, stimLog(iTrain).wavelength);
                    stimBoxDict(partialHash) = patch(ax, [0, pulseWidth, pulseWidth, 0], [iPulseStart, iPulseStart, iPulseEnd, iPulseEnd], color, ...
                        FaceAlpha=alpha, EdgeColor=color, EdgeAlpha=alpha, DisplayName=partialDesc);
                end
                [~, I] = sort(stimBoxDict.keys);
                handles = stimBoxDict.values;
                h = [h; handles(I)];
            end

            hold(ax, 'off')
            ax.YAxis.Direction = 'reverse';
            xlim(ax, p.Results.xlim);
            ylim(ax, [min(rd.I) - 1, max(rd.I) + 1]);
%             if nTrials <= maxTrials
%                 ylim(ax, [min(rd.I) - 1, max(rd.I) + 1]);
%             else
%                 ylim(ax, [min(rd.I) - 1, maxTrials + 1]);
%             end
            switch refName
                case 'press'
                    refName = 'touchbar-contact';
                case 'lick'
                    refName = 'spout-contact';
            end
            xlabel(ax, sprintf('Time to %s (%s)', refName, timeUnit), Interpreter='none')
            ylabel(ax, 'Trial')
            title(ax, sprintf('Spike raster (%s)', rd.name), Interpreter="none");
            if ~isTwoColorStim
                legend(ax, h, Location='northwest', FontSize=9);
            else
                legend(ax, h, Location='eastoutside', FontSize=9);
            end
        end

        function ax = plotMultiRaster(varargin)
            p = inputParser();
            if isgraphics(varargin{1}, 'axes')
                p.addRequired('ax', @(x) all(isgraphics(x, 'axes')) && length(x) > 1);
            end
            p.addRequired('rd', @isstruct);
            p.addParameter('label', string([]), @(x) isstring(x) || iscell(x)) % Ooh look at you using strings now. So fancy.
            p.addParameter('xlim', [-6, 1], @(x) isnumeric(x) && length(x) == 2 && x(2) > x(1));
            p.addParameter('iti', false, @islogical);
            p.addParameter('timeUnit', 's', @(x) ismember(x, {'s', 'ms'}))
            p.addParameter('maxTrials', Inf, @(x) isnumeric(x) || ischar(x))
            p.addParameter('sz', 2.5, @isnumeric)
            p.addParameter('figUnits', 'normalized', @ischar)
            p.addParameter('figPosition', [0, 0, 0.5, 0.75], @isnumeric)
            p.addParameter('figMargins', [0.16, 0.09], @isnumeric)
            p.parse(varargin{:})
            r = p.Results;

            rd = r.rd;
            label = r.label;
            maxTrials = r.maxTrials;

            if iscell(label)
                label = string(label);
            end
            if isempty(label)
                label = string({rd.trialType});
            end

            nTrials = arrayfun(@(rd) length(rd.duration), rd);
            if isstring(maxTrials)
                maxTrials = char(maxTrials);
            end
            if ischar(maxTrials)
                switch lower(maxTrials)
                    case 'min'
                        maxTrials = min(nTrials);
                    otherwise
                        error();
                end
            end
            nTrials = arrayfun(@(n) min(n, maxTrials), nTrials);

            if isfield(r, 'ax')
                ax = r.ax;
                fig = ax(1).Parent;
            else
                fig = figure(Units='normalized', Position=[0, 0, 0.5, 1], DefaultAxesFontSize=14);
                xMargin = r.figMargins(1);
                yMargin = r.figMargins(2);
                nS = length(rd);
                ax = gobjects(nS, 1);
                y = 1;
                for i = 1:length(rd)
                    x = xMargin;
                    w = 0.7;
                    h = nTrials(i)/sum(nTrials)*(1-yMargin*(nS + 1));
                    y = y - yMargin - h;
                    ax(i) = axes(fig, Position=[x, y, w, h]);
                end
            end
            fig.Units = r.figUnits;
            fig.Position = r.figPosition;


            assert(length(rd) == length(ax));
            assert(length(rd) == length(label));


            for i = 1:length(rd)
                EphysUnit.plotRaster(ax(i), rd(i), xlim=r.xlim, iti=r.iti, timeUnit=r.timeUnit, maxTrials=maxTrials, sz=r.sz, maxTrialsMethod='uniformsample');
                title(ax(i), label(i));
            end
        end

    end

    % private methods
    methods % (Access = {})
        function trials = makeTrials(obj, trialType, varargin)
            if length(obj) == 1
                switch lower(trialType)
                    case 'press'
                        trials = Trial(obj.EventTimes.Cue, obj.EventTimes.Press, 'first', obj.EventTimes.Lick);
                    case 'press_spontaneous'
                        trials = Trial(obj.EventTimes.PressOff, obj.EventTimes.Press, 'first');
                    case 'press_spontaneous2'
                        trials = Trial(obj.EventTimes.Press(1:end-1), obj.EventTimes.Press(2:end), advancedValidation=false);
                    case 'press_spontaneous_correct'
                        trials = Trial(obj.EventTimes.Cue, obj.EventTimes.Press, 'first');                        
                    case 'lick'
                        trials = Trial(obj.EventTimes.Cue, obj.EventTimes.Lick, 'first', obj.EventTimes.Press);
                    case 'stim'
                        trials = Trial([obj.EventTimes.StimOn, obj.EventTimes.StimOff(end)], obj.EventTimes.StimOff, 'first');
                    case 'stimtrain'
                        cueToTrainOn = Trial(obj.EventTimes.Cue, obj.EventTimes.StimOn, 'first');
                        cueToTrainOff = Trial(obj.EventTimes.Cue, obj.EventTimes.StimOff, 'last');
                        trials = Trial([cueToTrainOn.Stop], [cueToTrainOff.Stop]);
                    case 'light'
                        trials = Trial(obj.EventTimes.LightOn, obj.EventTimes.LightOff, 'first');
                    case 'firstlick'
                        trials = Trial(obj.EventTimes.Cue, obj.EventTimes.Lick, 'first');
                    case 'circlick'
                        lickTimes = unique(obj.EventTimes.Lick);
                        trials = Trial(lickTimes(1:end-1), lickTimes(2:end), advancedValidation=false);
                    case 'lickbout'
                        p = inputParser();
                        p.addParameter('minBoutCycles', 2) % shorter bouts are discarded entirely.
                        p.addParameter('maxBoutCycles', 4) % longer bouts are truncated, shorter bouts are padded with Empty Trial objects.
                        p.addParameter('maxInterval', 0.25) % interlick intervals longer than maxInterval breaks bouts
                        p.addParameter('minInterval', 0.05) % interlick intervals shorter than minInterval are considered artefacts with lick detection
                        p.parse(varargin{:})
                        minBoutCycles = p.Results.minBoutCycles;
                        maxBoutCycles = p.Results.maxBoutCycles;
                        maxInterval = p.Results.maxInterval;
                        minInterval = p.Results.minInterval;
                        
                        circTrials = obj.makeTrials('circlick');
                        isValidLick = circTrials.duration <= maxInterval & circTrials.duration >= minInterval;
                        boutStarts = strfind([false, isValidLick], [false, true]);
                        boutEnds = strfind([isValidLick, false], [true, false]);
                        boutCycles = boutEnds - boutStarts + 1;
                        isInvalidBout = boutCycles < minBoutCycles;
                        boutStarts(isInvalidBout) = [];
                        boutEnds(isInvalidBout) = [];
                        boutCycles(isInvalidBout) = [];

                        % Break down long bouts into many non-overlapping short bouts with length=maxBoutCycles
                        breakdownBoutStarts = [];
                        breakdownBoutEnds = [];
                        longBoutIndices = find(boutCycles > maxBoutCycles);
                        longBoutIndices = longBoutIndices(:)';
                        for iBout = longBoutIndices
                            breakdownBoutStarts = [breakdownBoutStarts, boutStarts(iBout):maxBoutCycles:boutEnds(iBout)-maxBoutCycles+1];
                            breakdownBoutEnds = [breakdownBoutEnds, boutStarts(iBout)+maxBoutCycles-1:maxBoutCycles:boutEnds(iBout)];
                        end
                        boutStarts(longBoutIndices) = [];
                        boutEnds(longBoutIndices) = [];
                        boutStarts = sort([boutStarts, breakdownBoutStarts]);
                        boutEnds = sort([boutEnds, breakdownBoutEnds]);
                        boutCycles = boutEnds - boutStarts + 1;
                        assert(all(boutCycles <= maxBoutCycles & boutCycles >= minBoutCycles), 'Probably messed up while splitting longer bouts into short ones.');

                        % Make trials array (2d, rows are bouts, columns are lick cycles
                        trials = repmat(Trial(), [length(boutStarts), maxBoutCycles]);
                        for iBout = 1:length(boutStarts)
                            trials(iBout, 1:boutCycles(iBout)) = circTrials(boutStarts(iBout):boutEnds(iBout));
                        end
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

            assert(sum(sc > 2^16 - 1) == 0, 'Spike counts per bin should not exceeded 2^16-1 (%i)', max(sc)) % There should not be more than 65k spikes in an 100ms bin
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
                edges = p.Results.edgesOrResolution(:)';
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

        function [xAligned, tAligned, requestedDuration] = getTrialAlignedData(obj, varargin)
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
            p.addOptional('trialType', 'press', @(x) ischar(x) && ismember(lower(x), {'press', 'lick', 'stim', 'stimtrain', 'stimfirstpulse', 'anylick', 'firstlick', 'circlick', 'lickbout', 'press_spontaneous', 'press_spontaneous2', 'press_spontaneous_medial', 'press_spontaneous_lateral'}))
            p.addParameter('alignTo', 'stop', @(x) ischar(x) && ismember(lower(x), {'start', 'stop'}))
            p.addParameter('resolution', 0.001, @(x) isnumeric(x) && x > 0)
            p.addParameter('allowedTrialDuration', [0, Inf], @(x) isnumeric(x) && length(x) >= 2 && x(2) >= x(1))
            p.addParameter('findSingleTrialDuration', 'off', @(x) ismember(x, {'off', 'min', 'max'})) % Used for opto, 'min' finds the shortest trial duration allowed by 'allowedTrialDuration', and only averages those trials.
            p.addParameter('trialDurationError', 1e-3, @isnumeric) % Used for opto, error allowed when finding identical trial durations.
            p.addParameter('includeInvalid', true, @islogical) % Whether to include unaligned ITI data for averaging. When alignTo='stop', pre-trial-start data is discarded. When alignTo='start', post-trial-end data is discarded. When 'stim', data after next opto-onset is discarded
            p.addParameter('correction', [], @isnumeric)
            p.addParameter('trials', [], @(x) isempty(x) || isa(x, 'Trial'))
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
            includeInvalid = p.Results.includeInvalid;
            err = p.Results.trialDurationError;
            correction = p.Results.correction;
                        
            % Filter out trials with incorrect lengths
            if isempty(p.Results.trials)
                trials = obj.getTrials(trialType);
            else
                trials = p.Results.trials;
            end
            if isempty(trials)
                xAligned = [];
                tAligned = [];
                requestedDuration = [];
                return
            end
            if ~isempty(correction)
                assert(length(correction) == length(trials))
                trials = Trial([trials.Start], [trials.Stop] + correction(:)');
            end
            if ~strcmpi(trialType, 'stimfirstpulse')
                iti = trials.iti();
            else
                % For stimfirstpulse (in train), use the iti between
                % first and second pulse in train (rather that iti
                % between first pulse in train and first pulse in next
                % train
                allTrials = obj.getTrials('stim', true);
                itiAll = allTrials.iti();
                startAll = [allTrials.Start];
                startFirst = [trials.Start];

                [match, index] = ismember(startFirst, startAll);

                if ~all(match)
                    warning('Not all stim pulses match stim trains. %s\n%s\n\t%s', num2str(index), num2str(match), num2str(startFirst(~match)))
                    warning('%g trials excluded', nnz(~match))
                    trials = trials(match);
                    index = index(match);
                end
                iti = itiAll(index);
            end

            if ismember(lower(trialType), {'circlick', 'lickbout'}) && strcmpi(p.Results.findSingleTrialDuration, 'off')
                requestedDuration = NaN;
            else
                durations = round(trials.duration()./err)*err;
                seltrials = durations >= allowedTrialDuration(1) & durations <= allowedTrialDuration(2);
                trials = trials(seltrials);
                durations = durations(seltrials);
                switch p.Results.findSingleTrialDuration
                    case 'min'
                        requestedDuration = allowedTrialDuration(1);
                        if ~any(durations == requestedDuration) && any(durations > requestedDuration)
                            requestedDuration = min(durations(durations > requestedDuration));
                        end
                        trials = trials(durations == requestedDuration);
                        iti = iti(durations == requestedDuration);
                    case 'max'
                        requestedDuration = allowedTrialDuration(2);
                        if ~any(durations == requestedDuration) && any(durations < requestedDuration)
                            requestedDuration = max(durations(durations < requestedDuration));
                        end
                        trials = trials(durations == requestedDuration);
                        iti = iti(durations == requestedDuration);
                    otherwise
                        requestedDuration = NaN;
                end
            end
            if isempty(trials)
                xAligned = [];
                tAligned = [];
                return
            end

            switch lower(trialType)
                case 'circlick'
                    tAligned = 0:resolution:2*pi;
                    tAlignedGlobal = zeros(length(trials), length(tAligned));
                    for iTrial = 1:length(trials)
                        tAlignedGlobal(iTrial, :) = linspace(trials(iTrial).Start, trials(iTrial).Stop, length(tAligned));
                    end
                case 'lickbout'
                    tAligned = 0:resolution:2*pi*size(trials, 2);
                    binsPerCycle = length(0:resolution:2*pi) - 1;
                    tAlignedGlobal = NaN(size(trials, 1), length(tAligned));
                    for iTrial = 1:size(trials, 1)
                        for iCycle = 1:size(trials, 2)
                            if trials(iTrial, iCycle).isEmpty
                                break;
                            end
                            tAlignedGlobal(iTrial, 1+(iCycle-1)*binsPerCycle:1+iCycle*binsPerCycle) = linspace(trials(iTrial, iCycle).Start, trials(iTrial, iCycle).Stop, binsPerCycle + 1);
                        end
                    end
                otherwise
                    tAligned = window(1):resolution:window(2);
                    switch alignTo
                        case 'stop'
                            tAlignedGlobal = tAligned + vertcat(trials.Stop);
                        case 'start'
                            tAlignedGlobal = tAligned + vertcat(trials.Start);
                    end
            end

            function sel = select(tt, iTrial)
                if includeInvalid
                    sel = true(1, length(tt));
                elseif strcmpi(alignTo, 'stop')
                    sel = tt >= -trials(iTrial).duration();
                % stim pulse: discard after next stim pulse onset
                elseif ismember(trialType, {'stim', 'stimfirstpulse'})
                    sel = tt <= trials(iTrial).duration() + iti(iTrial);
                else % discard after trial end
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
                assert(~strcmpi(trialType, 'circlick'), 'CircLick trials do not support resampling spike data.')
            % Recalculate data (count or rate) in new bins.
            else
                switch lower(trialType)
                    case {'circlick', 'lickbout'}
                        tAligned = (tAligned(1:end-1) + tAligned(2:end)) / 2;
                        xAligned = NaN(length(trials), length(tAligned));
                        switch data
                            case 'rate'
                                error('Not implemented: "rate" for circlick')
                            case 'count'
                                for iTrial = 1:length(trials)
                                    sel = ~isnan(tAlignedGlobal(iTrial, :));
                                    [xx, ~] = obj.getSpikeCounts(tAlignedGlobal(iTrial, sel));
                                    xAligned(iTrial, 1:length(xx)) = xx;
                                end
                        end
                    otherwise
                        tAligned = (tAligned(1:end-1) + tAligned(2:end)) / 2;
                        xAligned = NaN(length(trials), length(tAligned));
                        switch data
                            case 'rate'
                                kernel = obj.SpikeRateKernel;
                                width = kernel.params.width;
                                if strcmpi(kernel.type, 'gaussian')
                                    sigma = kernel.params.sigma;
        %                             sigma = 0.01;
                                    for iTrial = 1:length(trials)
                                        [xx, ~] = obj.getSpikeRates('gaussian', sigma, tAlignedGlobal(iTrial, :), 'kernelWidth', width);
                                        sel = select(tAligned, iTrial);
                                        xAligned(iTrial, sel) = xx(sel);
                                    end
                                else
                                    lambda1 = kernel.params.lambda1;
                                    lambda2 = kernel.params.lambda2;
                                    for iTrial = 1:length(trials)
                                        [xx, ~] = obj.getSpikeRates('exponential', lambda1, lambda2, tAlignedGlobal(iTrial, :), 'kernelWidth', width);
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
