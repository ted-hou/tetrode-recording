classdef AcuteRecording < handle
    properties
        expName = ''
        path = ''
        strain = ''
        stim = struct([])
        probeMap = struct([])
        bsr = struct([])
        bmrPress = struct([])
        bmrLick = struct([])
        conditions = struct([])
        stats = []
        statsInfo = struct([])
        statsPress = []
        statsLick = []
        statsInfoPress = struct([])
        statsInfoLick = struct([])
    end

    methods
        function obj = AcuteRecording(tr, strain)
            obj.strain = strain;

            splitPath = strsplit(tr.Path, '\');
            expName = splitPath{end-2};
            path = strjoin(splitPath(1:end-2), '\');
            if tr.Path(1) == '\'
                path = ['\', path];
            end
            obj.path = path;

            expName = strsplit(tr.GetExpName(), '_');
            obj.expName = strjoin(expName(1:2), '_');
        end

        function save(obj, varargin)
            p = inputParser();
            p.addOptional('Prefix', 'ar_', @ischar)
            p.addOptional('Path', obj.path, @ischar)
            p.parse(varargin{:})
            folder = sprintf('%s\\AcuteRecording', p.Results.Path);
            if ~isdir(folder)
                mkdir(folder)
            end
            file = sprintf('%s\\%s%s.mat', folder, p.Results.Prefix, obj.expName);
            tTic = tic();
            fprintf(1, 'Saving to file %s...', file)
            save(file, 'obj', '-v7.3');
            f = dir(file);
            fprintf(1, '(%.2fsec, %.2fMB)\n', toc(tTic), f.bytes*1e-6)
        end

        function stim = extractAllPulses(obj, tr, varargin)
            p = inputParser();
            p.addRequired('TetrodeRecording', @(x) isa(x, 'TetrodeRecording'))
            p.addOptional('FirstFiber', 'B', @(x) ismember(x, {'A', 'B'}))
            p.addOptional('FirstLight', 0.1, @(x) isnumeric(x) && length(x) == 1 && x > 0)
            p.parse(tr, varargin{:})
            tr = p.Results.TetrodeRecording;
            firstFiber = p.Results.FirstFiber;
            firstLight = p.Results.FirstLight;

            stim.calibration = obj.importCalibrationData();
            refPower = stim.calibration.(['Power_', firstFiber]);
            refLight = stim.calibration.(['Light_', firstFiber]);
            refGalvo = stim.calibration.(['Galvo_', firstFiber]);
            refPower = refPower(find(refLight == firstLight, 1));
            refGalvo = refGalvo(find(refLight == firstLight, 1));

            [stim.power, stim.galvo, ~, ~] = AcuteRecording.extractPulse(tr.DigitalEvents.LaserOn, tr.DigitalEvents.LaserOff, tr.AnalogIn.Timestamps, tr.AnalogIn.Data, 1);
            assert(stim.galvo == refGalvo);

            stim.powerCorrection = refPower - stim.power;

            nPulses = length(tr.DigitalEvents.LaserOn);
            stim.tOn = NaN(nPulses, 1);
            stim.tOff = NaN(nPulses, 1);
            stim.power = NaN(nPulses, 1);
            stim.galvo = NaN(nPulses, 1);
            stim.dvRank = NaN(nPulses, 1);
            stim.mlRank = NaN(nPulses, 1);
            stim.fiber = repmat('?', nPulses, 1);
            stim.light = NaN(nPulses, 1);
            stim.powerError = NaN(nPulses, 1);
            
            % Trim analog data where laser was off
            t = tr.AnalogIn.Timestamps;
            data = tr.AnalogIn.Data;
            sel = data(1, :) > 1000;
            t = t(sel);
            data = data(:, sel);
            
            for iPulse = 1:nPulses
                [stim.power(iPulse), stim.galvo(iPulse), stim.tOn(iPulse), stim.tOff(iPulse)] = AcuteRecording.extractPulse(tr.DigitalEvents.LaserOn, tr.DigitalEvents.LaserOff, t, data, iPulse, stim.powerCorrection);
                [stim.light(iPulse), stim.fiber(iPulse), stim.powerError(iPulse)] = AcuteRecording.findMatchingCalibration(stim.calibration, stim.power(iPulse), stim.galvo(iPulse));
            end
            
            isFiberA = stim.fiber == 'A';
            isFiberB = ~isFiberA;
            stim.mlRank(isFiberA) = 1;
            stim.mlRank(isFiberB) = 2;
            [~, stim.dvRank(isFiberA)] = ismember(abs(stim.galvo(isFiberA)), unique(abs(stim.galvo(isFiberA))));
            [~, stim.dvRank(isFiberB)] = ismember(abs(stim.galvo(isFiberB)), unique(abs(stim.galvo(isFiberB))));
            
            % Figure out iTrain and iPulseInTrain
            edges = transpose([tr.DigitalEvents.GalvoOn(:), tr.DigitalEvents.GalvoOff(:)]);
            edges = edges(:);
            [N, ~, bins] = histcounts(stim.tOn, edges); % Odd bins are in train, 1 -> 1st train, 3 -> 2nd, 5 -> 3rd, k -> (k + 1)/2
            stim.train = (bins + 1) / 2;
            N = N(1:2:end);
            assert(all(N>0), 'Not implemented: Handling trains with zero pulses.')
            stim.pulse = zeros(nPulses, 1);
            i = 0;
            for iTrain = 1:length(N)
                n = N(iTrain);
                stim.pulse(i + 1:i + n) = 1:n;
                i = i + n;
            end
            
            stim.duration = round((stim.tOff - stim.tOn) .* 1000) ./ 1000;

            obj.stim = stim;
        end

        function calibrationData = importCalibrationData(obj)
            filename = sprintf('%s\\GalvoCalibration_%s.csv', obj.path, obj.expName);

            % Set up the Import Options and import the data
            opts = delimitedTextImportOptions("NumVariables", 7);

            % Specify range and delimiter
            opts.DataLines = [2, Inf];
            opts.Delimiter = ",";

            % Specify column names and types
            opts.VariableNames = ["Light_A", "Power_A", "Galvo_A", "Light_B", "Power_B", "Galvo_B", "Var7"];
            opts.SelectedVariableNames = ["Light_A", "Power_A", "Galvo_A", "Light_B", "Power_B", "Galvo_B"];
            opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "string"];

            % Specify file level properties
            opts.ExtraColumnsRule = "ignore";
            opts.EmptyLineRule = "read";

            % Specify variable properties
            opts = setvaropts(opts, "Var7", "WhitespaceRule", "preserve");
            opts = setvaropts(opts, "Var7", "EmptyFieldRule", "auto");

            % Import the data
            calibrationData = readtable(filename, opts);
        end      

        function plotPowerError(obj)
            fig = figure();
            ax = axes(fig);
            hold(ax, 'on');
            plot(obj.stim.power, 'k', 'DisplayName', 'Recorded');
            plot(obj.stim.power - obj.stim.powerError, 'r', 'DisplayName', 'Corrected');
            hold(ax, 'off');
            legend(ax);
        end

        function probeMap = importProbeMap(obj, frontOrBack, varargin)
            p = inputParser();
            p.addRequired('FrontOrBack', @(x) ischar(x) && ismember(x, {'front', 'back'})); % Whether probe is anterior (front) or posterior (back) to the headstage. Usually this is back, unless recording in anterior SNr, and probe needed to be reversed (front) to make space for fiber optics.
            p.addOptional('ML', 1300, @isnumeric); % Center of probe. Negative for left hemisphere.
            p.addOptional('DV', -4600, @isnumeric); % From surface of brain, not bregma. Bregma is usually 200um above surface here.
            p.addOptional('AP', -3280, @isnumeric);
            p.addParameter('DVOffset', -200, @isnumeric); % Distance from bregma to surface of brain at SNr.
            p.parse(frontOrBack, varargin{:});
            front = strcmpi(p.Results.FrontOrBack, 'front');
            back = ~front;
            ml = p.Results.ML;
            dv = p.Results.DV + p.Results.DVOffset;
            ap = p.Results.AP;
            left = ml < 0;
            right = ~left;
            
            load('128DN_bottom.mat', 's');
            % Determine whether to flip the probe x positions/shank number
            % For the probe, when viewed from the front (black epoxy visible), shank=1, x=0 is left-most shank, z=0 is deepest in brain.
            % For our convention, we want to convert to brain coordinates (ml/dv) (shank=1 is medial)
            flipDir = 1;
            if right
                flipDir = flipDir * -1;
            end
            if back
                flipDir = flipDir * -1;
            end
            
            if flipDir == -1
                s.shaft = 5 - s.shaft;
                s.x = 450 - s.x;
            end
            
            map = zeros(32, 4);
            for iShaft = 1:4
                inShaft = s.shaft == iShaft;
                [~, I] = sort(s.z(inShaft), 'descend');
                channels = find(inShaft);
                channels = channels(I);
                map(:, iShaft) = channels;
            end
            
            probeMap.channel = s.channels + 1;
            probeMap.ml = ml - 225 + s.x;
            probeMap.dv = dv + s.tipelectrode + s.z;
            probeMap.ap = repmat(ap, size(s.x));
            probeMap.shaft = s.shaft;
            probeMap.map = map;
            
            obj.probeMap = probeMap;
        end

        function varargout = binMoveResponse(obj, tr, moveType, varargin)
            p = inputParser();
            p.addRequired('TetrodeRecording', @(x) isa(x, 'TetrodeRecording'))
            p.addRequired('MoveType', @(x) ismember(lower(x), {'press', 'lick'}))
            p.addOptional('Channels', [], @isnumeric);
            p.addOptional('Units', [], @isnumeric);
            p.addParameter('BinWidth', 0.1, @isnumeric); % Bin width in seconds
            p.addParameter('Window', [-6, 1]); % Additional seconds before and after tMove
            p.addParameter('BaselineWindow', [-6, -2]); % Window used for normalization, relative to move time
            p.addParameter('Store', false, @islogical);
            p.parse(tr, moveType, varargin{:})
            tr = p.Results.TetrodeRecording;
            moveType = p.Results.MoveType;
            channels = p.Results.Channels;
            units = p.Results.Units;
            binWidth = p.Results.BinWidth;
            window = p.Results.Window;
            baselineWindow = p.Results.BaselineWindow;
            store = p.Results.Store;

            % Do all channels
            if isempty(channels)
                channels = [tr.Spikes.Channel];
            end
            
            nChannels = length(channels);
            nUnitsInChannel = zeros(nChannels, 1);
            for i = 1:nChannels
                nUnitsInChannel(i) = max(tr.Spikes(channels(i)).Cluster.Classes) - 1;
            end
            assert(isempty(units) || max(nUnitsInChannel) >= max(units), 'Requested units (%s) exceeds total number of units (%i) in data.', num2str(units), max(units));
            
            if isempty(units)
                nUnits = sum(nUnitsInChannel);
            else
                nUnits = length(units) * nChannels;
            end
            
            bmr(nUnits) = struct('expName', '', 'channel', [], 'unit', [], 't', [], 'tRight', [], 'spikeRates', [], 'normalizedSpikeRates', [], 'trialLength', []);
            m = zeros(nUnits, 1);
            s = zeros(nUnits, 1);

            tRef = tr.DigitalEvents.CueOn;
            switch lower(moveType)
                case 'press'
                    tMove = tr.DigitalEvents.PressOn;
                    tExclude = sort([tr.DigitalEvents.LickOn, tr.DigitalEvents.RewardOn]);
                case 'lick'
                    tMove = tr.DigitalEvents.LickOn;
                    tExclude = sort([tr.DigitalEvents.PressOn, tr.DigitalEvents.RewardOn]);
                otherwise
                    error('Unrecognized move type: %s', lower(moveType))
            end

            [tRef, tMove] = tr.FindFirstInTrial(tRef, tMove, tExclude);

            % Bin spikes
            nTrials = length(tMove);
            t = window(1):binWidth:window(2);
            tBaseline = baselineWindow(1):binWidth:baselineWindow(2);
            nBins = length(t) - 1;
            nBaselineBins = length(tBaseline) - 1;

            iUnit = 0;
            for iChn = 1:nChannels
                channel = channels(iChn);
                
                if isempty(units)
                    selUnits = 1:nUnitsInChannel(iChn);
                else
                    selUnits = units;
                end
                
                for iUnitInChannel = 1:length(selUnits)
                    iUnit = iUnit + 1;
                    unit = selUnits(iUnitInChannel);

                    % Extract unit spike times
                    sel = tr.Spikes(channel).Cluster.Classes == unit;
                    spikeTimes = tr.Spikes(channel).Timestamps(sel);

                    % Stats for selecting units (session-mean firing rate 
                    spikeRates = histcounts(spikeTimes, spikeTimes(1):binWidth:spikeTimes(end)) ./ binWidth;
                    m(iUnit) = mean(spikeRates);
                    s(iUnit) = mad(spikeRates, 0) / 0.6745;
                    bmr(iUnit).expName = obj.expName;
                    bmr(iUnit).channel = channel;
                    bmr(iUnit).unit = unit;
                    bmr(iUnit).spikeRates = zeros(nTrials, nBins);
                    bmr(iUnit).normalizedSpikeRates = zeros(nTrials, nBins);
                    bmr(iUnit).trialLength = tMove - tRef;
                    for iTrial = 1:nTrials
                        edges = tMove(iTrial) + window(1):binWidth:tMove(iTrial) + window(2);
                        baselineEdges = tMove(iTrial) + baselineWindow(1):binWidth:tMove(iTrial) + baselineWindow(2);
                        bmr(iUnit).tRight = edges(2:end) - tMove(iTrial);
                        bmr(iUnit).t = (edges(1:end-1) + edges(2:end))/2 - tMove(iTrial);
                        bmr(iUnit).spikeRates(iTrial, :) = histcounts(spikeTimes, edges) ./ binWidth;
                        baselineSpikeRates = histcounts(spikeTimes, baselineEdges) ./ binWidth;
                        % Normalize by substracting mean norm window and dividing by global MAD
                        bmr(iUnit).normalizedSpikeRates(iTrial, :) = (bmr(iUnit).spikeRates(iTrial, :) - mean(baselineSpikeRates)) ./ s(iUnit);
                    end
                end
            end

            if store
                switch lower(moveType)
                    case 'press'
                        obj.bmrPress = bmr;
                    case 'lick'
                        obj.bmrLick = bmr;
                end
            end

            varargout = {bmr, m, s};
        end

        function stats = summarizeMoveResponse(obj, bmrOrType, varargin)
            % stats = ar.summarizeMoveResponse('press', 'peak', [-1, 0], 'AllowedTrialLength', [2, Inf])
            p = inputParser();
            if ischar(bmrOrType)
                p.addRequired('MoveType', @(x) ismember(lower(x), {'press', 'lick'}))
            else
                p.addRequired('BinnedMoveResponse', @isstruct)
            end
            p.addOptional('Method', 'peak', @(x) ismember(x, {'peak', 'mean'}));
            p.addOptional('Window', [-1, 0], @(x) isnumeric(x) && length(x) == 2);
            p.addParameter('AllowedTrialLength', [0, Inf], @(x) isnumeric(x) && length(x) == 2);
            p.addParameter('Store', false, @islogical);
            p.parse(bmrOrType, varargin{:});

            if ischar(bmrOrType)
                switch lower(p.Results.MoveType)
                    case 'press'
                        bmr = obj.bmrPress;
                    case 'lick'
                        bmr = obj.bmrLick;
                    otherwise
                        error()
                end
            else
                bmr = p.Results.BinnedMoveResponse;
            end

            nUnits = length(bmr);
            stats = NaN(nUnits, 1);
            for i = 1:length(bmr)
                sel = bmr(i).trialLength >= p.Results.AllowedTrialLength(1) & bmr(i).trialLength <= p.Results.AllowedTrialLength(2);
                meanNSR = mean(bmr(i).normalizedSpikeRates(sel, :), 1);
                stats(i, :) = AcuteRecording.summarizeMeanSpikeRates(meanNSR, bmr(i).t, p.Results.Method, p.Results.Window);
            end

            if p.Results.Store && isfield(p.Results, 'MoveType')
                switch lower(p.Results.MoveType)
                    case 'press'
                        obj.statsPress = stats;
                        obj.statsInfoPress = struct('Method', p.Results.Method, 'Window', p.Results.Window, 'AllowedTrialLength', p.Results.AllowedTrialLength);
                    case 'lick'
                        obj.statsLick = stats;
                        obj.statsInfoLick = struct('Method', p.Results.Method, 'Window', p.Results.Window, 'AllowedTrialLength', p.Results.AllowedTrialLength);
                end
            end
        end

        function varargout = binStimResponse(obj, tr, varargin)
            p = inputParser();
            p.addRequired('TetrodeRecording', @(x) isa(x, 'TetrodeRecording'))
            p.addOptional('Channels', @isnumeric);
            p.addOptional('Units', [], @isnumeric);
            p.addParameter('BinWidth', 0.01, @isnumeric); % Bin width in seconds
            p.addParameter('Window', [-0.2, 0.5]); % Additional seconds before and after tOn
            p.addParameter('Store', false, @islogical);
            p.parse(tr, varargin{:});
            tr = p.Results.TetrodeRecording;
            channels = p.Results.Channels;
            units = p.Results.Units;
            binWidth = p.Results.BinWidth;
            window = p.Results.Window;
            store = p.Results.Store;
            
            % Do all channels
            if isempty(channels)
                channels = [tr.Spikes.Channel];
            end
            
            nChannels = length(channels);
            nUnitsInChannel = zeros(nChannels, 1);
            for i = 1:nChannels
                nUnitsInChannel(i) = max(tr.Spikes(channels(i)).Cluster.Classes) - 1;
            end
            assert(isempty(units) || max(nUnitsInChannel) >= max(units), 'Requested units (%s) exceeds total number of units (%i) in data.', num2str(units), max(units));
            
            if isempty(units)
                nUnits = sum(nUnitsInChannel);
            else
                nUnits = length(units) * nChannels;
            end
            
            bsr(nUnits) = struct('expName', '', 'channel', [], 'unit', [], 't', [], 'tRight', [], 'spikeRates', [], 'normalizedSpikeRates', []);
            m = zeros(nUnits, 1);
            s = zeros(nUnits, 1);

            % Bin spikes
            nPulses = length(obj.stim.tOn);
            t = window(1):binWidth:window(2);
            nBins = length(t) - 1;

            iUnit = 0;
            for iChn = 1:nChannels
                channel = channels(iChn);
                
                if isempty(units)
                    selUnits = 1:nUnitsInChannel(iChn);
                else
                    selUnits = units;
                end
                
                for iUnitInChannel = 1:length(selUnits)
                    iUnit = iUnit + 1;
                    unit = selUnits(iUnitInChannel);

                    % Extract unit spike times
                    sel = tr.Spikes(channel).Cluster.Classes == unit;
                    spikeTimes = tr.Spikes(channel).Timestamps(sel);

                    % Stats for normalizing spike rates
                    spikeRates = histcounts(spikeTimes, spikeTimes(1):binWidth:spikeTimes(end)) ./ binWidth;
                    m(iUnit) = mean(spikeRates);
                    s(iUnit) = mad(spikeRates, 0) / 0.6745;
                    bsr(iUnit).expName = obj.expName;
                    bsr(iUnit).channel = channel;
                    bsr(iUnit).unit = unit;
                    bsr(iUnit).spikeRates = zeros(nPulses, nBins);
                    bsr(iUnit).normalizedSpikeRates = zeros(nPulses, nBins);
                    isPreWindow = t < 0;
                    for iPulse = 1:nPulses
                        edges = obj.stim.tOn(iPulse) + window(1):binWidth:obj.stim.tOn(iPulse) + window(2);
                        % bsr(iUnit).t = (edges(1:end - 1) + edges(2:end))/2 - obj.stim.tOn(iPulse);
                        bsr(iUnit).tRight = edges(2:end) - obj.stim.tOn(iPulse);
                        bsr(iUnit).t = (edges(1:end-1) + edges(2:end))/2 - obj.stim.tOn(iPulse);
                        bsr(iUnit).spikeRates(iPulse, :) = histcounts(spikeTimes, edges) ./ binWidth;
                        % Normalize by substracting pre stim window and dividing by global MAD
                        bsr(iUnit).normalizedSpikeRates(iPulse, :) = (bsr(iUnit).spikeRates(iPulse, :) - mean(bsr(iUnit).spikeRates(iPulse, isPreWindow))) ./ s(iUnit);
                    end
                end
            end

            if store
                obj.bsr = bsr;
            end
            
            varargout = {bsr, m, s};
        end

        function varargout = selectStimResponse(obj, varargin)
            p = inputParser();
            p.addOptional('BinnedStimResponse', obj.bsr, @isstruct);
            p.addParameter('Light', [], @isnumeric);
            p.addParameter('Duration', [], @isnumeric);
            p.addParameter('MLRank', [], @isnumeric); % 1: most medial, 4: most lateral
            p.addParameter('DVRank', [], @isnumeric); % 1: most ventral, 4: most dorsal
            p.addParameter('Fiber', {}, @(x) ischar(x) || iscell(x));
            p.addParameter('Galvo', [], @isnumeric);
            p.parse(varargin{:});
            bsr = p.Results.BinnedStimResponse;
            crit.light = p.Results.Light;
            crit.duration = p.Results.Duration;
            crit.mlRank = p.Results.MLRank;
            crit.dvRank = p.Results.DVRank;
            crit.fiber = p.Results.Fiber;
            crit.galvo = p.Results.Galvo;
            
            sel = true(size(obj.stim.tOn));
            if ~isempty(crit.light)
                sel = sel & ismember(obj.stim.light, crit.light);
            end
            if ~isempty(crit.duration)
                sel = sel & ismember(obj.stim.duration, crit.duration);
            end
            if ~isempty(crit.mlRank)
                sel = sel & ismember(obj.stim.mlRank, crit.mlRank);
            end
            if ~isempty(crit.dvRank)
                sel = sel & ismember(obj.stim.dvRank, crit.dvRank);
            end
            if ~isempty(crit.fiber)
                sel = sel & ismember(obj.stim.fiber, crit.fiber);
            end
            if ~isempty(crit.galvo)
                sel = sel & ismember(obj.stim.galvo, crit.galvo);
            end
            
            selBSR = bsr;
            selPulse = find(sel);
            for i = 1:length(bsr)
                selBSR(i).spikeRates = selBSR(i).spikeRates(sel, :);
                selBSR(i).normalizedSpikeRates = selBSR(i).normalizedSpikeRates(sel, :);
                selBSR(i).selPulse = selPulse;
            end
            
            fprintf(1, '%i trials selected with provided criteria.\n', length(selPulse))
            
            varargout = {selBSR, selPulse};
        end

        function [stats, conditions] = summarize(obj, bsr, varargin)
            p = inputParser();
            p.addRequired('BinnedStimResponse', @isstruct)
            p.addOptional('Method', 'peak', @(x) ismember(x, {'peak', 'mean', 'firstPeak'}))
            p.addOptional('Window', [0, 0.05], @(x) isnumeric(x) && length(x) == 2)
            p.addOptional('FirstPeakThreshold', 0, @isnumeric)
            p.addParameter('Store', false, @islogical)
            p.parse(bsr, varargin{:})
            bsr = p.Results.BinnedStimResponse;

            [conditions, ~, ~, ~] = obj.groupByConditions(bsr(1));
            nConditions = length(conditions);
            nUnits = length(bsr);
            stats = NaN(nUnits, nConditions);
            for i = 1:length(bsr)
                [~, ~, condNSR, ~] = obj.groupByConditions(bsr(i));
                stats(i, :) = AcuteRecording.summarizeMeanSpikeRates(condNSR, bsr(i).t, p.Results.Method, p.Results.Window, p.Results.FirstPeakThreshold);
            end

            if p.Results.Store
                obj.conditions = conditions;
                obj.stats = stats;
                obj.statsInfo = struct('Method', p.Results.Method, 'Window', p.Results.Window, 'FirstPeakThreshold', p.Results.FirstPeakThreshold);
            end
        end

        function varargout = groupByConditions(obj, bsr)
            assert(length(bsr) == 1)

            spikeRates = bsr.spikeRates;
            normalizedSpikeRates = bsr.normalizedSpikeRates;
            if isfield(bsr, 'selPulse')
                I = bsr.selPulse;
            else
                I = 1:length(obj.stim.tOn);
            end

            % Find total number of unique conditions
            light = obj.stim.light(I);
            duration = obj.stim.duration(I);
            [~, lightRank] = ismember(light, unique(light));
            [~, durationRank] = ismember(duration, unique(duration));
            mlRank = obj.stim.mlRank(I);
            dvRank = obj.stim.dvRank(I);
            fiber = obj.stim.fiber(I);
            galvo = obj.stim.galvo(I);
            
            base = max([max(lightRank), max(durationRank), max(mlRank), max(dvRank)]);
            condId = lightRank*(base^3) + durationRank*(base^2) + mlRank*(base^1) + dvRank*(base^0);
            [uniqueConditions, ia] = unique(condId);
            nConditions = length(uniqueConditions);
            nColorsA = length(unique(dvRank(mlRank==1)));
            nColorsB = length(unique(dvRank(mlRank==2)));
            conditions(nConditions) = struct('id', [], 'label', '', 'numTrials', [], 'light', [], 'duration', [], 'fiber', '', 'galvo', [], 'mlRank', [], 'dvRank', [] , 'linewidth', '');
            condSR = NaN(nConditions, size(spikeRates, 2));
            condNSR = NaN(nConditions, size(spikeRates, 2));
            for iCond = 1:nConditions
                i = ia(iCond);
                conditions(iCond).light = light(i);
                conditions(iCond).duration = duration(i);
                conditions(iCond).mlRank = mlRank(i);
                conditions(iCond).dvRank = dvRank(i);
                conditions(iCond).fiber = fiber(i);
                conditions(iCond).galvo = galvo(i);
                conditions(iCond).id = condId(i);
                switch mlRank(i)
                    case 1
                        mlText = 'mStr';
                    case 2
                        mlText = 'lStr';
                end
                switch dvRank(i)
                    case 1
                        dvText = '-4.15';
                    case 2
                        dvText = '-3.48';
                    case 3
                        dvText = '-2.81';
                    case 4
                        dvText = '-2.15';
                end
                conditions(iCond).label = sprintf('%.1fmW, %.0fms (%s %s)', light(i), duration(i)*1000, mlText, dvText);
                conditions(iCond).linewidth = AcuteRecording.lerp(1, 2, (lightRank(i) - 1)/(max(lightRank) - 1));
                if mlRank(i) == 1
                    conditions(iCond).linecolor = [0.9, AcuteRecording.lerp(0.1, 0.9, (dvRank(i)-1)/(nColorsA-1)), 0.1];
                elseif mlRank(i) == 2
                    conditions(iCond).linecolor = [AcuteRecording.lerp(0.1, 0.9, (dvRank(i)-1)/(nColorsB-1)), 0.9, 0.1];
                else
                    error('ml rank must be 1 or 2, got %i instead', mlrank(i));
                end
                condSel = condId == condId(ia(iCond));
                conditions(iCond).numTrials = nnz(condSel);
                condSR(iCond, :) = mean(spikeRates(condSel, :), 1);
                condNSR(iCond, :) = mean(normalizedSpikeRates(condSel, :), 1);
            end
            
            varargout = {conditions, condSR, condNSR, condId};
        end

        function plotStimResponse(obj, bsr, varargin)
            p = inputParser();
            p.addRequired('BinnedStimResponse', @isstruct);
            p.addParameter('Mode', 'both', @(x) ischar(x) && ismember(x, {'heatmap', 'line', 'both'}));
            p.addParameter('CLim', [], @isnumeric);
            p.parse(bsr, varargin{:});
            bsr = p.Results.BinnedStimResponse;

            if length(bsr) > 1
                for i = 1:length(bsr)
                    obj.plotStimResponse(bsr(i), 'Mode', p.Results.Mode, 'CLim', p.Results.CLim);
                end
                return
            end
            
            [conditions, condSR, condNSR, ~] = groupByConditions(obj, bsr);
            nConditions = length(conditions);

            fig = figure('Units', 'Normalized', 'Position', [0.1, 0, 0.8, AcuteRecording.lerp(0.125, 0.3, nConditions/10)]);
            
            switch p.Results.Mode
                case 'line'
                    ax(1) = axes(fig);
                case 'heatmap'
                    ax(2) = axes(fig);
                case 'both'
                    ax(1) = subplot(2, 1, 1);
                    ax(2) = subplot(2, 1, 2);
                    fig.Position(4) = fig.Position(4)*2;
            end
            
            if ismember(p.Results.Mode, {'line', 'both'})
                hold(ax(1), 'on')
                for iCond = 1:nConditions
                    plot(ax(1), bsr.tRight, condSR(iCond, :), 'Color', conditions(iCond).linecolor, 'LineWidth', conditions(iCond).linewidth, 'DisplayName', conditions(iCond).label);
                end
                hold(ax(1), 'off')
                xlabel(ax(1), 'Time (s)')
                ylabel(ax(1), 'Spike Rate (sp/s)')
                l = legend(ax(1), 'Location', 'west');
                l.Position(1) = l.Position(1)-l.Position(3);
            end
            if ismember(p.Results.Mode, {'heatmap', 'both'})
                imagesc(ax(2), bsr.t, 1:nConditions, condNSR, [-max(abs(condNSR(:))), max(abs(condNSR(:)))])
                colormap(ax(2), 'jet');
                cb = colorbar(ax(2));
                cb.Label.String = 'Normalized Spike Rate (a.u.)';
                xlabel(ax(2), 'Time (s)')
                yticks(ax(2), 1:nConditions)
                yticklabels(ax(2), {conditions.label})
                ylabel(ax(2), 'Stim Condition')
                title(ax(2), sprintf('%s Chn %i Unit %i', bsr.expName, bsr.channel, bsr.unit), 'Interpreter', 'none')

                if ~isempty(p.Results.CLim)
                    ax(2).CLim = p.Results.CLim;
                end
            end
        end

        function titleText = plotStimResponseMap(obj, bsr, varargin)
            p = inputParser();
            p.addRequired('BinnedStimResponse', @(x) isstruct(x) || iscell(x));
            p.addOptional('SRange', [0.25, 3], @(x) isnumeric(x) && length(x)==2);
            p.addOptional('Threshold', 0.25, @isnumeric);
            p.addOptional('Method', 'peak', @(x) ismember(x, {'peak', 'firstPeak', 'mean'}));
            p.addOptional('Window', [0, 0.05], @(x) isnumeric(x) && length(x)==2);
            p.addOptional('FirstPeakThreshold', [], @isnumeric);
            p.addParameter('UseSignedML', false, @islogical);
            p.addParameter('HideFlatUnits', false, @islogical);
            p.addParameter('ConditionBase', 4, @isnumeric); % Max number of stim conditions per category (e.g. if there are 2 durations, 3 light levels, 2 fibers and 4 galvo voltages, then use base 4 = max([2, 3, 2, 4])). Suggest 4, or [] to autocalculate.            
            p.parse(bsr, varargin{:});
            srange = p.Results.SRange;
            threshold = p.Results.Threshold;
            method = p.Results.Method;
            window = p.Results.Window;
            firstPeakThreshold = p.Results.FirstPeakThreshold;
            if isempty(firstPeakThreshold)
                firstPeakThreshold = threshold;
            end

            % Single experiment
            if length(obj) == 1
                coords = obj.getProbeCoords([bsr.channel]);
                [stats, conditions] = obj.summarize(bsr, method, window, firstPeakThreshold);

                nConditions = length(conditions);
                nCols = max(1, floor(sqrt(nConditions)));
                nRows = max(4, ceil(nConditions / nCols));
                fig = figure('Units', 'normalized', 'Position', [0, 0, 0.4, 1]);
                ax = gobjects(nConditions, 1);
                methodLabel = method;
                methodLabel(1) = upper(method(1));
                for iCond = 1:nConditions
                    [i, j] = ind2sub([nRows, nCols], iCond);
                    iSubplot = sub2ind([nCols, nRows], j, nRows + 1 - i);
                    ax(iCond) = subplot(nRows, nCols, iSubplot);
                    h = AcuteRecording.plotMap(ax(iCond), coords, stats(:, iCond), srange, threshold, [bsr.channel], methodLabel, 'UseSignedML', p.Results.UseSignedML);
                    title(ax(iCond), conditions(iCond).label)
                    axis(ax(iCond), 'equal')
                end
                figure(fig);
                titleText = sprintf('%s (%s)', obj.expName, obj.strain);
                suptitle(titleText);
            elseif length(obj) > 1
                assert(iscell(bsr) && length(obj) == length(bsr))
                
                % Pool experiments
                for i = 1:length(obj)
                    coords{i} = obj(i).getProbeCoords([bsr{i}.channel]);
                    [stats{i}, conditions{i}] = obj(i).summarize(bsr{i}, method, window, firstPeakThreshold);
                end
                pooledConditions = AcuteRecording.poolConditions(conditions, p.Results.ConditionBase);
                nConditions = length(pooledConditions);
                pooledCoords = cell(1, nConditions);
                pooledStats = cell(1, nConditions);
                pooledChannels = cell(1, nConditions);
                for iCond = 1:nConditions
                    id = pooledConditions(iCond).id;
                    for iExp = 1:length(obj)
                        iCondInExp = find([conditions{iExp}.id] == id);
                        if ~isempty(iCondInExp)
                            pooledCoords{iCond} = vertcat(pooledCoords{iCond}, coords{iExp});
                            pooledStats{iCond} = vertcat(pooledStats{iCond}, stats{iExp}(:, iCondInExp));
                            pooledChannels{iCond} = vertcat(pooledChannels{iCond}, [bsr{iExp}.channel]');
                        end
                    end
                end
                
                % Plot
                nCols = max(1, floor(sqrt(nConditions)));
                nRows = max(4, ceil(nConditions / nCols));
                fig = figure('Units', 'normalized', 'Position', [0, 0, 0.4, 1]);
                ax = gobjects(nConditions, 1);
                methodLabel = method;
                methodLabel(1) = upper(method(1));
                for iCond = 1:nConditions
                    [i, j] = ind2sub([nRows, nCols], iCond);
                    iSubplot = sub2ind([nCols, nRows], j, nRows + 1 - i);
                    ax(iCond) = subplot(nRows, nCols, iSubplot);
                    h = AcuteRecording.plotMap(ax(iCond), pooledCoords{iCond}, pooledStats{iCond}, srange, threshold, pooledChannels{iCond}, methodLabel, 'UseSignedML', p.Results.UseSignedML, 'HideFlatUnits', p.Results.HideFlatUnits);
                    title(ax(iCond), pooledConditions(iCond).label)
                    axis(ax(iCond), 'equal')
%                     xlim(ax(iCond), [0.9, 1.7])
                end
                
                strain = unique({obj.strain});
                if length(strain) == 1
                    strain = strain{1};
                else
                    strain = 'Multiple Strains';
                end
                nUnits = sum(cellfun(@length, bsr));
                nSessions = length(obj);
                animalNames = cellfun(@toAnimalName, {obj.expName}, 'UniformOutput', false);
                nAnimals = length(unique(animalNames));
                figure(fig);
                titleText = sprintf('%s (%i animals, %i sessions)', strain, nAnimals, nSessions);
                suptitle(titleText);
            end

            function animalName = toAnimalName(expName)
                animalName = strsplit(expName, '_');
                animalName = animalName{1};
            end
        end

        function varargout = plotStimVsMoveResponse(obj, bsr, moveType, varargin)
            function m = mean2(x)
                m = mean(x, 2);
            end

            function m = max2(x)
                [~, I] = max(abs(x), [], 2);
                m = diag(x(:, I));
            end
            
            p = inputParser();
            p.addRequired('BinnedStimResponse', @(x) isstruct(x) || iscell(x));
            p.addRequired('MoveType', @(x) ismember(lower(x), {'press', 'lick'}));
            p.addParameter('StimThreshold', 0.5, @isnumeric);
            p.addParameter('MoveThreshold', 1, @isnumeric);
            p.addParameter('Highlight', 'intersect', @(x) ismember(lower(x), {'move', 'stim', 'union', 'intersect'}))
            p.addParameter('PoolConditions', 'off', @(x) ismember(lower(x), {'off', 'mean', 'max'})) % Should multiple stim conditions be merged into one?
            p.addParameter('ConditionBase', 4, @isnumeric); % Max number of stim conditions per category (e.g. if there are 2 durations, 3 light levels, 2 fibers and 4 galvo voltages, then use base 4 = max([2, 3, 2, 4])). Suggest 4, or [] to autocalculate.
            p.parse(bsr, moveType, varargin{:});
            bsr = p.Results.BinnedStimResponse;
            moveType = p.Results.MoveType;

            if length(obj) == 1
                % Get single-numeric move/stim stats
                switch lower(moveType)
                    case 'press'
                        bmr = obj.bmrPress;
                    case 'lick'
                        bmr = obj.bmrLick;
                    otherwise
                        error()
                end
                moveStats = obj.summarizeMoveResponse(bmr, 'peak', [-1, 0], 'AllowedTrialLength', [1, Inf]);
                [stimStats, stimConditions] = obj.summarize(bsr, 'peak', [0, 0.05]);

                % Figure layout, one subplot per stim condition
                if strcmpi(p.Results.PoolConditions, 'off')
                    nConditions = length(stimConditions);
                    nCols = max(1, floor(sqrt(nConditions)));
                    nRows = max(4, ceil(nConditions / nCols));
                    fig = figure('Units', 'normalized', 'Position', [0, 0, 0.4, 1]);
                    ax = gobjects(nConditions, 1);
                    for iCond = 1:nConditions
                        [i, j] = ind2sub([nRows, nCols], iCond);
                        iSubplot = sub2ind([nCols, nRows], j, nRows + 1 - i);
                        ax(iCond) = subplot(nRows, nCols, iSubplot);
                        AcuteRecording.plotScatter(ax(iCond), moveStats, stimStats(:, iCond), p.Results.MoveThreshold, p.Results.StimThreshold, ...
                            'MoveType', moveType, 'Highlight', p.Results.Highlight, 'Title', stimConditions(iCond).label);
                    end
                    figure(fig);
                    titleText = sprintf('%s (%s)', obj.expName, obj.strain);
                    suptitle(titleText);
                else
                    fig = figure('Units', 'Normalized', 'Position', [0, 0, 0.4, 0.4]);
                    ax = axes(fig);
                    if strcmpi(p.Results.PoolConditions, 'mean')
                        stimStats = mean2(stimStats);
                        titleText = sprintf('Mean stim effect vs. %s', moveType);
                    elseif strcmpi(p.Results.PoolConditions, 'max')
                        stimStats = max2(stimStats);
                        titleText = sprintf('Max stim effect vs. %s', moveType);
                    end
                    AcuteRecording.plotScatter(ax, moveStats, stimStats, p.Results.MoveThreshold, p.Results.StimThreshold, ...
                        'MoveType', moveType, 'Highlight', p.Results.Highlight, 'Title', titleText);
                end
            elseif length(obj) > 1
                assert(iscell(bsr) && length(obj) == length(bsr))
                
                % Get POOLED single-numeric move/stim stats
                switch lower(moveType)
                    case 'press'
                        bmr = {obj.bmrPress};
                    case 'lick'
                        bmr = {obj.bmrLick};
                    otherwise
                        error()
                end

                for i = 1:length(obj)
                    moveStats{i} = obj(i).summarizeMoveResponse(bmr{i}, 'peak', [-1, 0], 'AllowedTrialLength', [1, Inf]);
                    [stimStats{i}, stimConditions{i}] = obj(i).summarize(bsr{i}, 'peak', [0, 0.05]);
                end

                if strcmpi(p.Results.PoolConditions, 'off')
                    pooledConditions = AcuteRecording.poolConditions(stimConditions, p.Results.ConditionBase);
                    nConditions = length(pooledConditions);
                    pooledStimStats = cell(1, nConditions);
                    pooledMoveStats = cell(1, nConditions);
                    for iCond = 1:nConditions
                        id = pooledConditions(iCond).id;
                        for iExp = 1:length(obj)
                            iCondInExp = find([stimConditions{iExp}.id] == id);
                            if ~isempty(iCondInExp)
                                pooledStimStats{iCond} = vertcat(pooledStimStats{iCond}, stimStats{iExp}(:, iCondInExp));
                                pooledMoveStats{iCond} = vertcat(pooledMoveStats{iCond}, moveStats{iExp});
                            end
                        end
                    end
    
                    % Plot
                    nCols = max(1, floor(sqrt(nConditions)));
                    nRows = max(4, ceil(nConditions / nCols));
                    fig = figure('Units', 'normalized', 'Position', [0, 0, 0.4, 1]);
                    ax = gobjects(nConditions, 1);
                    methodLabel = 'Peak';
                    for iCond = 1:nConditions
                        stimStats = pooledStimStats{iCond};
                        moveStats = pooledMoveStats{iCond};
        
                        [i, j] = ind2sub([nRows, nCols], iCond);
                        iSubplot = sub2ind([nCols, nRows], j, nRows + 1 - i);
                        ax(iCond) = subplot(nRows, nCols, iSubplot, 'Tag', 'scatter');
                        AcuteRecording.plotScatter(ax(iCond), moveStats, stimStats, p.Results.MoveThreshold, p.Results.StimThreshold, ...
                            'MoveType', moveType, 'Highlight', p.Results.Highlight, 'Title', pooledConditions(iCond).label);
                    end

                    varargout = {fig, pooledStimStats, pooledMoveStats, pooledConditions};
                else
                    switch lower(p.Results.PoolConditions)
                        case 'mean'
                            stimStats = cellfun(@mean2, stimStats, 'UniformOutput', false);
                            titleText = sprintf('Mean stim effect vs. %s', moveType);
                        case 'max'
                            stimStats = cellfun(@max2, stimStats, 'UniformOutput', false);
                            titleText = sprintf('Max stim effect vs. %s', moveType);
                    end
                    stimStats = cat(1, stimStats{:});
                    moveStats = cat(1, moveStats{:});

                    fig = figure('Units', 'Normalized', 'Position', [0, 0, 0.4, 0.4]);
                    ax = axes(fig, 'Tag', 'scatter');
                    AcuteRecording.plotScatter(ax, moveStats, stimStats, p.Results.MoveThreshold, p.Results.StimThreshold, ...
                        'MoveType', moveType, 'Highlight', p.Results.Highlight, 'Title', titleText);

                    varargout = {fig, stimStats, moveStats};
                end
                
                strain = unique({obj.strain});
                if length(strain) == 1
                    strain = strain{1};
                else
                    strain = 'Multiple Strains';
                end
                nUnits = sum(cellfun(@length, bsr));
                nSessions = length(obj);
                animalNames = cellfun(@toAnimalName, {obj.expName}, 'UniformOutput', false);
                nAnimals = length(unique(animalNames));
                figure(fig);
                titleText = sprintf('%s (%i animals, %i sessions)', strain, nAnimals, nSessions);
                suptitle(titleText);
            end

            function animalName = toAnimalName(expName)
                animalName = strsplit(expName, '_');
                animalName = animalName{1};
            end
        end

        function coords = getProbeCoords(obj, channels)
            map = obj.probeMap;
            coords = zeros(length(channels), 3);
            for i = 1:length(channels)
                I = find(map.channel == channels(i), 1);
                assert(~isempty(I))
                coords(i, 1) = map.ml(I);
                coords(i, 2) = map.dv(I);
                coords(i, 3) = map.ap(I);
            end
        end
    end

    methods (Static)
        function obj = load(varargin)
            p = inputParser();
            p.addOptional('FilesOrDirs', '', @(x) ischar(x) || iscell(x))
            p.parse(varargin{:})
            filesOrDirs = p.Results.FilesOrDirs;

            if isempty(filesOrDirs)
                filesOrDirs = uipickfiles('FilterSpec', '*.mat');
            end
            
            % Read all files in folder (string, dir path)
            files = {};
            if ~iscell(filesOrDirs)
                assert(ischar(filesOrDirs));
                filesOrDirs = {filesOrDirs};
            end


            for i = 1:length(filesOrDirs)
                thisPath = filesOrDirs{i};
                if isfile(thisPath)
                    files = [files, thisPath];
                elseif isdir(thisPath)
                    if thisPath(end) == '\'
                        thisPath = thisPath(1:end-1);
                    end
                    filesInPath = dir(sprintf('%s\\ar_*.mat', thisPath));
                    filesInPath = cellfun(@(x) sprintf('%s\\%s', thisPath, x), {filesInPath.name}, 'UniformOutput', false);
                    files = [files, filesInPath];
                else
                    error('Path %i of %i (''%s'') cannot be read because it is neither a file or an directory.', i, length(filesOrDirs), thisPath)
                end
            end

            tTic = tic();
            S(length(files)) = struct('obj', []);
            for i = 1:length(files)
                fprintf(1, 'Loading file %i of %i (%s)...\n', i, length(files), files{i});
                S(i) = load(files{i}, 'obj');
            end
            obj = [S.obj];
            fprintf(1, 'Loaded %i files in %.2f seconds.\n', length(files), toc(tTic));
        end

        % TODO: Return iPulseInTrain
        function [power, galvo, tOn, tOff] = extractPulse(laserOn, laserOff, analogTimestamps, analogData, iPulse, offset, threshold)
            if nargin < 6
                offset = 0;
            end
            if nargin < 7
                threshold = 1000;
            end

            sel = analogTimestamps >= laserOn(iPulse) & analogTimestamps <= laserOff(iPulse) + 0.005;
            t = analogTimestamps(sel);
            power = analogData(1, sel);
            galvo = round(analogData(2, sel) / 100) * 100;

            iOn = find(power > threshold, 1);
            iOff = find(power > threshold, 1, 'last');
            power = mean(power(iOn:iOff)) + offset;
            galvo = mean(galvo);
            tOn = t(iOn);
            tOff = t(iOff);
        end

        function [light, fiber, powerOffset] = findMatchingCalibration(data, power, galvo)
            % A
            sel_A = data.Galvo_A == galvo;
            sel_B = data.Galvo_B == galvo;

            if any(sel_A) && any(sel_B)
                power_A = data.Power_A(sel_A);
                light_A = data.Light_A(sel_A);
                [df_A, i] = min(abs(power_A - power));
                light_A = light_A(i);
                powerOffset_A = power - power_A(i);

                power_B = data.Power_B(sel_B);
                light_B = data.Light_B(sel_B);
                [df_B, i] = min(abs(power_B - power));
                light_B = light_B(i);
                powerOffset_B = power - power_B(i);

                if df_A < df_B
                    fiber = 'A';
                    light = light_A;
                    powerOffset = powerOffset_A;
                else
                    fiber = 'B';
                    light = light_B;
                    powerOffset = powerOffset_B;
                end
            elseif any(sel_A)
                fiber = 'A';
                power_A = data.Power_A(sel_A);
                light_A = data.Light_A(sel_A);
                [~, i] = min(abs(power_A - power));
                light = light_A(i);
                powerOffset = power - power_A(i);
            elseif any(sel_B)
                fiber = 'B';
                power_B = data.Power_B(sel_B);
                light_B = data.Light_B(sel_B);
                [~, i] = min(abs(power_B - power));
                light = light_B(i);
                powerOffset = power - power_B(i);
            else
                error('No galvo voltage matches %f in calibration data.', galvo);
            end
        end

        function [stats, info] = summarizeMeanSpikeRates(msr, t, varargin)
            p = inputParser();
            p.addRequired('MeanSpikeRates', @isnumeric)
            p.addRequired('Timestamps', @isnumeric)
            p.addOptional('Method', 'peak', @(x) ismember(x, {'peak', 'mean', 'firstPeak'}))
            p.addOptional('Window', [0, 0.05], @(x) isnumeric(x) && length(x) == 2)
            p.addOptional('FirstPeakThreshold', 0, @isnumeric)
            p.parse(msr, t, varargin{:});
            msr = p.Results.MeanSpikeRates;
            t = p.Results.Timestamps;
            window = p.Results.Window;

            N = size(msr, 1);
            sel = t <= window(2) & t >= window(1);
            msr = msr(:, sel); % Truncate by window
            switch p.Results.Method
                case 'mean'
                    stats = mean(msr, 2);
                case 'peak'
                    [~, I] = max(abs(msr), [], 2);
                    stats = diag(msr(:, I));
                case 'firstPeak'
                    stats = NaN(N, 1);
                    for i = 1:N
                        [peaks, ~] = AcuteRecording.findpeaks(msr(i, :), p.Results.FirstPeakThreshold);
                        stats(i) = peaks(1);
                    end
                otherwise
                    error('Not implemented method %s', p.Results.Method)
            end
            info.name = p.Results.Method;
            info.window = window;
        end

        function h = plotMap(varargin)
            p = inputParser();
            if isgraphics(varargin{1})
                p.addRequired('ax', @isgraphics)
                ax = varargin{1};
            else
                ax = gca();
            end
            p.addRequired('coords', @isnumeric)
            p.addRequired('stats', @isnumeric)
            p.addOptional('srange', [0.25, 3], @(x) isnumeric(x) && length(x) == 2)
            p.addOptional('threshold', 0.25, @isnumeric)
            p.addOptional('channels', [], @isnumeric)
            p.addOptional('method', 'Stat', @ischar)
            p.addParameter('SLim', [9, 72], @isnumeric)
            p.addParameter('ARange', [0.25, 1], @isnumeric)
            p.addParameter('HideFlatUnits', false, @islogical)
            p.addParameter('UseSignedML', false, @islogical);
            p.parse(varargin{:})
            
            coords = p.Results.coords;
            stats = p.Results.stats;
            srange = p.Results.srange;
            threshold = p.Results.threshold;
            channels = p.Results.channels;
            method = p.Results.method;

            t = AcuteRecording.inverseLerp(srange(1), srange(2), abs(stats));
            C = zeros(length(stats), 3);
            isUp = stats>=threshold;
            isDown = stats<=-threshold;
            isFlat = ~(isUp | isDown);
            C(isUp, 1) = 1;
            C(isDown, 3) = 1;
            C(isFlat, :) = 0.5;
            S = AcuteRecording.lerp(p.Results.SLim(1), p.Results.SLim(2), t);
            A = ones(size(S)) * 40;
            if p.Results.HideFlatUnits
                A(isFlat) = 1;
            end
            ml = coords(:, 1) / 1000;
            dv = coords(:, 2) / 1000;
            ap = coords(:, 3) / 1000;
            if ~p.Results.UseSignedML
                ml = abs(ml);
            end
            h = scatter3(ax, ml, dv, ap, S, C, 'filled', 'MarkerFaceAlpha', 'flat', 'AlphaData', A, 'AlphaDataMapping', 'direct');
            view(ax, 0, 90)
            if p.Results.UseSignedML
                mlLabel = 'ML';
            else
                mlLabel = 'ML (abs)';
            end
            xlabel(ax, mlLabel)
            ylabel(ax, 'DV')
            zlabel(ax, 'AP')
            
            h.DataTipTemplate.DataTipRows(1).Label = mlLabel;
            h.DataTipTemplate.DataTipRows(2).Label = 'DV';
            h.DataTipTemplate.DataTipRows(3).Label = 'AP';
            h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Channel', channels);
            h.DataTipTemplate.DataTipRows(end+2) = dataTipTextRow(method, stats);
            h.DataTipTemplate.DataTipRows(end+3) = dataTipTextRow('Index', 1:length(stats));
        end

        function varargout = plotScatter(varargin)
            p = inputParser();
            if isgraphics(varargin{1})
                p.addRequired('ax', @isgraphics)
                ax = varargin{1};
            else
                ax = gca();
            end
            p.addRequired('moveStats', @isnumeric)
            p.addRequired('stimStats', @isnumeric)
            p.addOptional('moveThreshold', 1, @isnumeric)
            p.addOptional('stimThreshold', 0.5, @isnumeric)
            p.addParameter('Highlight', 'intersect', @(x) ismember(lower(x), {'stim', 'move', 'union', 'intersect'}))
            p.addParameter('Title', '', @ischar)
            p.addParameter('MoveType', 'Move', @(x) ismember(lower(x), {'press', 'lick'}))
            p.parse(varargin{:})
            
            moveStats = p.Results.moveStats;
            stimStats = p.Results.stimStats;
            moveThreshold = p.Results.moveThreshold;
            stimThreshold = p.Results.stimThreshold;

            hold(ax, 'on')
    
            switch lower(p.Results.Highlight)
                case 'stim'
                    sel = abs(stimStats) >= stimThreshold; 
                    labelSig = sprintf('N=%i (stim resp)', nnz(sel));
                    labelInsig = sprintf('N=%i', nnz(~sel));
                case 'move'
                    sel = abs(moveStats) >= moveThreshold;    
                    labelSig = sprintf('N=%i (move resp)', nnz(sel));
                    labelInsig = sprintf('N=%i', nnz(~sel));
                case 'union'
                    sel = abs(moveStats) >= moveThreshold | abs(stimStats) >= stimThreshold;
                    labelSig = sprintf('N=%i (move&stim resp)', nnz(sel));
                    labelInsig = sprintf('N=%i', nnz(~sel));
                case 'intersect'
                    sel = abs(moveStats) >= moveThreshold & abs(stimStats) >= stimThreshold;
                    labelSig = sprintf('N=%i (move/stim resp)', nnz(sel));
                    labelInsig = sprintf('N=%i', nnz(~sel));
            end

            hSig = scatter(ax, moveStats(sel), stimStats(sel), 35, 'blue', 'filled', 'DisplayName', labelSig);
            hInsig = scatter(ax, moveStats(~sel), stimStats(~sel), 25, 'black', 'DisplayName', labelInsig);
            axis(ax, 'equal')
            hY = plot(ax, ax.XLim, [0, 0], 'k--');
            hX = plot(ax, [0, 0], ax.YLim, 'k--');
            title(ax, p.Results.Title)
            xlabel(ax, p.Results.MoveType)
            ylabel(ax, 'Stim')
            hold(ax, 'off')
            hLegend = legend(ax, [hSig, hInsig], 'Location', 'northoutside', 'Orientation', 'horizontal');
            varargout = {ax, hSig, hInsig, hX, hY, hLegend};
        end

        function id = calculateConditionID(conditions, base)
            light = [conditions.light];
            duration = [conditions.duration];
            [~, lightRank] = ismember(light, unique(light));
            [~, durationRank] = ismember(duration, unique(duration));
            mlRank = [conditions.mlRank];
            dvRank = [conditions.dvRank];
            
            if nargin < 2
                base = max([max(lightRank), max(durationRank), max(mlRank), max(dvRank)]);
            end
            id = lightRank*(base^3) + durationRank*(base^2) + mlRank*(base^1) + dvRank*(base^0);
        end

        function pc = poolConditions(conditions, base)
            if nargin < 2
                base = []; % Suggest use base 4, or leave empty to autocalculate.
            end

            assert(iscell(conditions) && length(conditions) > 1)
            for i = 1:length(conditions)
                id = AcuteRecording.calculateConditionID(conditions{i}, base);
                for j = 1:length(conditions{i})
                    conditions{i}(j).id = id(j);
                end
            end

            allConditions = conditions{1};
            for i = 2:length(conditions)
                allConditions = [allConditions, conditions{i}];
            end
            uniqueIDs = unique([allConditions.id]);
            pc(length(uniqueIDs)) = struct('id', [], 'label', [], 'numTrials', [], 'light', [], 'duration', [], 'fiber', [], 'galvo', [], 'mlRank', [], 'dvRank', [], 'linewidth', [], 'linecolor', []);
            for i = 1:length(uniqueIDs)
                id = uniqueIDs(i);
                sel = [allConditions.id] == id;
                pc(i).id = id;
                pc(i).numTrials = range([allConditions(sel).numTrials]);
                pc(i).light = range([allConditions(sel).light]);
                pc(i).duration = range([allConditions(sel).duration]);
                pc(i).fiber = unique([allConditions(sel).fiber]);
                pc(i).galvo = unique([allConditions(sel).galvo]);
                pc(i).mlRank = unique([allConditions(sel).mlRank]);
                pc(i).dvRank = unique([allConditions(sel).dvRank]);
                pc(i).linewidth = mean(unique([allConditions(sel).linewidth]));
                pc(i).linecolor = mean(vertcat(allConditions(sel).linecolor), 1);
                
                assert(length(pc(i).fiber) == 1)
                assert(length(pc(i).mlRank) == 1)
                assert(length(pc(i).dvRank) == 1)
                
                switch pc(i).mlRank
                    case 1
                        mlText = 'mStr';
                    case 2
                        mlText = 'lStr';
                end
                switch pc(i).dvRank
                    case 1
                        dvText = '-4.15';
                    case 2
                        dvText = '-3.48';
                    case 3
                        dvText = '-2.81';
                    case 4
                        dvText = '-2.15';
                end
                lightText = range2text(pc(i).light, '%.1f');
                durationText = range2text(pc(i).duration*1000, '%.0f');
                numTrialsText = range2text(pc(i).numTrials, '%.0f');
                
                pc(i).label = sprintf('%smW, %sms (%s %s) (%s trials)', lightText, durationText, mlText, dvText, numTrialsText);
            end
            
            function r = range(x)
                r = [min(x), max(x)];
                if r(1) == r(2)
                    r = r(1);
                end
            end
            
            function s = range2text(r, formatSpec)
                if length(r) == 1
                    s = sprintf(formatSpec, r);
                else
                    s = sprintf([formatSpec, '-', formatSpec], r(1), r(2));
                end
            end
        end
    end

    methods (Static, Access = {})
        function x = lerp(a, b, t)
            if isnan(t)
                t = 0;
            end
            t = max(0, min(1, t));
            x = a + (b - a).*t;
        end

        function t = inverseLerp(a, b, x)
            t = (x - a) / (b - a);
            t = max(0, min(1, t));
        end

        function [peaks, I] = findpeaks(x, varargin)
            % Find local maxima and minima (maginitue must be greater than
            % threshold (default 0). If no peak found, use max(abs).
            p = inputParser();
            p.addRequired('X', @isnumeric)
            p.addOptional('Threshold', 0, @isnumeric)
            p.parse(x, varargin{:})
            x = p.Results.X;
            threshold = p.Results.Threshold;
            
            x = [0, x(:)'];
            
            df = [0, diff(x)];
            df1 = df > 0;
            df2 = df <= 0;
            df3 = df < 0;
            df4 = df >= 0;
            I = find(abs(x) > threshold & ((df1 & circshift(df2, -1)) | (df3 & circshift(df4, -1))));
            if isempty(I)
                [~, I] = max(abs(x));
            end
            peaks = x(I);
            I = I - 1;
        end
    end
end