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
    end
    
    % public methods
    methods
        % constructor
        function obj = EphysUnit(varargin)
            if nargin == 0
                return
            end

            p = inputParser();
            addOptional(p, 'PETH', struct(), @isstruct);
            addParameter(p, 'ReadData', true, @islogical);
            addParameter(p, 'ReadITI', false, @islogical);
            addParameter(p, 'ExtendedWindow', [-1, 1], @(x) isnumeric(x) && length(x) == 2)
            addParameter(p, 'ReadWaveforms', false, @islogical);
            parse(p, varargin{:});
            PETH = p.Results.PETH;
            extendedWindow = p.Results.ExtendedWindow;

            % Initialize obj array
            obj(length(PETH)) = EphysUnit();

            % Extract unique expNames
            uniqueExpNames = unique({PETH.ExpName});
            i = 1;
            fprintf(1, 'Reading PETH...\n')
            for iExp = 1:length(uniqueExpNames)
                fprintf(1, 'Loading experiment %s (%i of %i)...\n', uniqueExpNames{iExp}, iExp, length(uniqueExpNames))
                if p.Results.ReadData
                    tr = TetrodeRecording.BatchLoadSimple(uniqueExpNames{iExp});
                    inExp = strcmpi({PETH.ExpName}, uniqueExpNames{iExp});
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
                            if p.Results.ReadWaveforms
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

                            % Make trials
                            obj(i).Trials.Press = obj(i).makeTrials('press');
                            obj(i).Trials.Lick = obj(i).makeTrials('lick');
                            obj(i).Trials.Stim = obj(i).makeTrials('stim');
                            
                            % Cull ITI spikes if told to do so.
                            if ~p.Results.ReadITI
                                [inTrialPress, extendedWindowPress] = obj(i).findSpikesInTrial('press', extendedWindow);
                                [inTrialLick, extendedWindowLick] = obj(i).findSpikesInTrial('lick', extendedWindow);
                                [inTrialStim, extendedWindowStim] = obj(i).findSpikesInTrial('stim', extendedWindow);
                                
                                inTrial = inTrialPress | inTrialLick | inTrialStim;

                                obj(i).SpikeTimes = obj(i).SpikeTimes(inTrial);
                                if p.Results.ReadWaveforms
                                    obj(i).Waveforms = obj(i).Waveforms(inTrial, :);
                                end
                                obj(i).ExtendedWindow.Press = extendedWindowPress;
                                obj(i).ExtendedWindow.Lick = extendedWindowLick;
                                obj(i).ExtendedWindow.Stim = extendedWindowStim;
                            else
                                obj(i).ExtendedWindow.Press = [];
                                obj(i).ExtendedWindow.Lick = [];
                                obj(i).ExtendedWindow.Stim = [];
                            end
                            i = i + 1;
                        catch ME
                            warning('Error when processing unit %i (%s) - this one will be skipped.', i, e.ExpName)
                            warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
                            i = i + 1;
                        end			
                    end
                    clear tr
                else
                    inExp = strcmpi({PETH.ExpName}, uniqueExpNames{iExp});
                    for e = PETH(inExp)
                        obj(i).ExpName = uniqueExpNames{iExp};
                        obj(i).Channel = e.Channel;
                        obj(i).Unit = e.Cluster;
                        i = i + 1;
                    end
                end
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
        
        function trials = makeTrials(obj, trialType)
            switch lower(trialType)
                case 'press'
                    trials = Trial(obj.EventTimes.Cue, obj.EventTimes.Press, 'first', obj.EventTimes.Lick);
                case 'lick'
                    trials = Trial(obj.EventTimes.Cue, obj.EventTimes.Lick, 'first', obj.EventTimes.Press);
                case 'stim'
                    trials = Trial(obj.EventTimes.StimOn, obj.EventTimes.StimOff);
                case 'all'
                    for i = 1:length(obj)
                        obj(i).Trials.Press = obj(i).makeTrials('press');
                        obj(i).Trials.Lick = obj(i).makeTrials('lick');
                        obj(i).Trials.Stim = obj(i).makeTrials('stim');
                    end
            end
        end
        
        function plot(obj, varargin)
            p = inputParser();
            parse(p, varargin{:});
            
        end
    end

    % private methods
    methods (Access = {})
        function [inTrial, extendedWindow] = findSpikesInTrial(obj, trialType, extendedWindow)
            switch lower(trialType)
                case 'press'
                    [~, start, stop] = Trial.findEdges(obj.EventTimes.Cue, obj.EventTimes.Press, 'first', obj.EventTimes.Lick);
                    [edges, extendedWindow] = Trial.extendEdges(start, stop, 'window', extendedWindow);
                case 'lick'
                    [~, start, stop] = Trial.findEdges(obj.EventTimes.Cue, obj.EventTimes.Lick, 'first', obj.EventTimes.Press);
                    [edges, extendedWindow] = Trial.extendEdges(start, stop, 'window', extendedWindow);
                case 'stim'
                    % Use trainOn/trainOff as start/stop of stim trials
                    [~, ~, start] = Trial.findEdges(obj.EventTimes.Cue, obj.EventTimes.StimOn, 'first');
                    [~, ~, stop] = Trial.findEdges(obj.EventTimes.Cue, obj.EventTimes.StimOff, 'last');
                    [edges, extendedWindow] = Trial.extendEdges(start, stop, 'window', extendedWindow, 'shrinkPriority', 'postwindow');
            end
            
            if ~isempty(edges)
                % Find logical indices for spikes inside trial window.
                [~, ~, bins] = histcounts(obj.SpikeTimes, edges);
                inTrial = rem(bins, 2) ~= 0; % spikes in odd bins occur in trial window
            else
                inTrial = false(1, length(obj.SpikeTimes));
            end
        end
    end
end
