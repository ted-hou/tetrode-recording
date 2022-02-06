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
        SpikeRates = []
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

		function displayName = getDisplayName(obj)
            if length(obj) == 1
                if isnan(obj.Electrode)
                    displayName = sprintf('%s Channel%i Unit%i', obj.ExpName, obj.Channel, obj.Unit);
                else
                    displayName = sprintf('%s Electrode%i Unit%i', obj.ExpName, obj.Electrode, obj.Unit);
                end
           else
                displayName = arrayfun(@getDisplayName, obj, 'UniformOutput', false);
            end
        end
    end
    
    % static methods
    methods (Static)
    end
    
    % private methods
    methods (Access = {})
        function [inTrial, extendedWindow] = findSpikesInTrial(obj, trialType, extendedWindow)
            switch lower(trialType)
                case 'press'
                    [start, stop] = TetrodeRecording.FindFirstInTrial(obj.EventTimes.Cue, obj.EventTimes.Press, obj.EventTimes.Lick, 'first');
                case 'lick'
                    [start, stop] = TetrodeRecording.FindFirstInTrial(obj.EventTimes.Cue, obj.EventTimes.Lick, obj.EventTimes.Press, 'first');
                case 'stim'
                    % Use trainOn/trainOff as start/stop of stim trials
                    [~, start] = TetrodeRecording.FindFirstInTrial(obj.EventTimes.Cue, obj.EventTimes.StimOn);
                    [~, stop] = TetrodeRecording.FindLastInTrial(obj.EventTimes.Cue, obj.EventTimes.StimOff);
            end
            
            [edges, extendedWindow] = EphysUnit.validateTrialWindow(start, stop, extendedWindow, trialType);
            
            if ~isempty(edges)
                % Find logical indices for spikes inside trial window.
                [~, ~, bins] = histcounts(obj.SpikeTimes, edges);
                inTrial = rem(bins, 2) ~= 0; % spikes in odd bins occur in trial window
            else
                inTrial = false(1, length(obj.SpikeTimes));
            end
        end
    end
    
    % private static methods
    methods (Access = {}, Static)
        function [edges, extendedWindow] = validateTrialWindow(cue, move, extendedWindow, trialType)
            if isempty(cue) || isempty(move)
                warning('Number of valid %s trials is zero.', trialType);
                edges = [];
                return
            end
            
            edges = reshape(vertcat(cue + extendedWindow(1), move + extendedWindow(2)), [], 1);
            % Check if edges is monotonic increasing
            if nnz(diff(edges) > 0) + 1 == length(edges)
                return
            end
            
            % if trial edges not valid, retry a smaller window
            if extendedWindow(1) < 0 && extendedWindow(2) > 0
                extendedWindow = [0, extendedWindow(2)];
            elseif extendedWindow(2) > 0
                extendedWindow = [0, 0];
            else
                warning('Smallest extended window [%f, %f] would not work for %s trials.', extendedWindow(1), extendedWindow(2), trialType)
                edges = [];
                return;
            end
            
            [edges, extendedWindow] = EphysUnit.validateTrialWindow(cue, move, extendedWindow);
            
        end
    end
end
