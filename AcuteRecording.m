classdef AcuteRecording < handle
    properties
        tr
        stim = struct([])
        probeMap = struct([])
        strain = ''
    end
    
    methods
        function obj = AcuteRecording(tr, strain)
            obj.tr = tr;
            obj.strain = strain;
        end
        
        function stim = extractAllPulses(obj, firstFiber, firstLight)
            if nargin < 2
                firstFiber = 'B';
            end
            if nargin < 3
                firstLight = 0.1;
            end

            stim.calibration = obj.importCalibrationData();
            refPower = stim.calibration.(['Power_', firstFiber]);
            refLight = stim.calibration.(['Light_', firstFiber]);
            refGalvo = stim.calibration.(['Galvo_', firstFiber]);
            refPower = refPower(find(refLight == firstLight, 1));
            refGalvo = refGalvo(find(refLight == firstLight, 1));

            [stim.power, stim.galvo, ~, ~] = AcuteRecording.extractPulse(obj.tr.DigitalEvents.LaserOn, obj.tr.DigitalEvents.LaserOff, obj.tr.AnalogIn.Timestamps, obj.tr.AnalogIn.Data, 1);
            assert(stim.galvo == refGalvo);

            stim.powerCorrection = refPower - stim.power;

            nPulses = length(obj.tr.DigitalEvents.LaserOn);
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
            t = obj.tr.AnalogIn.Timestamps;
            data = obj.tr.AnalogIn.Data;
            sel = data(1, :) > 1000;
            t = t(sel);
            data = data(:, sel);
            
            for iPulse = 1:nPulses
                [stim.power(iPulse), stim.galvo(iPulse), stim.tOn(iPulse), stim.tOff(iPulse)] = AcuteRecording.extractPulse(obj.tr.DigitalEvents.LaserOn, obj.tr.DigitalEvents.LaserOff, t, data, iPulse, stim.powerCorrection);
                [stim.light(iPulse), stim.fiber(iPulse), stim.powerError(iPulse)] = AcuteRecording.findMatchingCalibration(stim.calibration, stim.power(iPulse), stim.galvo(iPulse));
            end
            
            isFiberA = stim.fiber == 'A';
            isFiberB = ~isFiberA;
            stim.mlRank(isFiberA) = 1;
            stim.mlRank(isFiberB) = 2;
            [~, stim.dvRank(isFiberA)] = ismember(abs(stim.galvo(isFiberA)), unique(abs(stim.galvo(isFiberA))));
            [~, stim.dvRank(isFiberB)] = ismember(abs(stim.galvo(isFiberB)), unique(abs(stim.galvo(isFiberB))));
            
            % Figure out iTrain and iPulseInTrain
            edges = transpose([obj.tr.DigitalEvents.GalvoOn(:), obj.tr.DigitalEvents.GalvoOff(:)]);
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
            splitPath = strsplit(obj.tr.Path, '\');
            expName = splitPath{end-2};
            path = strjoin(splitPath(1:end-2), '\');
            if obj.tr.Path(1) == '\'
                path = ['\', path];
            end
            filename = sprintf('%s\\GalvoCalibration_%s.csv', path, expName);

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
            p.addRequired('FrontOrBack', @(x) ischar(x) && ismember(x, {'front', 'back'})); % Whether probe is anterior or posterior to the headstage
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
            if right
                probeMap.ml = ml - 225 + s.x;
            else
                probeMap.ml = ml + 225 + s.x;
            end
            probeMap.dv = dv + s.tipelectrode + s.z;
            probeMap.ap = repmat(ap, size(s.x));
            probeMap.shaft = s.shaft;
            probeMap.map = map;
            
            obj.probeMap = probeMap;
        end
        
        function [bsr, m, s] = binStimResponse(obj, channels, varargin)
            p = inputParser();
            p.addRequired('Channels', @isnumeric);
            p.addOptional('Units', [], @isnumeric);
            p.addParameter('BinWidth', 0.01, @isnumeric); % Bin width in seconds
            p.addParameter('Window', [-0.2, 0.5]); % Additional seconds before and after tOn
            p.parse(channels, varargin{:});
            channels = p.Results.Channels;
            units = p.Results.Units;
            binWidth = p.Results.BinWidth;
            window = p.Results.Window;
            
            % Do all channels
            if isempty(channels)
                channels = [obj.tr.Spikes.Channel];
            end
            
            nChannels = length(channels);
            nUnitsInChannel = zeros(nChannels, 1);
            for i = 1:nChannels
                nUnitsInChannel(i) = max(obj.tr.Spikes(channels(i)).Cluster.Classes) - 1;
            end
            assert(isempty(units) || max(nUnitsInChannel) >= max(units), 'Requested units (%s) exceeds total number of units (%i) in data.', num2str(units), max(units));
            
            if isempty(units)
                nUnits = sum(nUnitsInChannel);
            else
                nUnits = length(units) * nChannels;
            end
            
            bsr(nUnits) = struct('expName', '', 'channel', [], 'unit', [], 't', [], 'spikeRates', [], 'normalizedSpikeRates', []);
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
                    sel = obj.tr.Spikes(channel).Cluster.Classes == unit;
                    spikeTimes = obj.tr.Spikes(channel).Timestamps(sel);

                    % Stats for normalizing spike rates
                    spikeRates = histcounts(spikeTimes, spikeTimes(1):binWidth:spikeTimes(end)) ./ binWidth;
                    m(iUnit) = mean(spikeRates);
                    s(iUnit) = mad(spikeRates, 0) / 0.6745;
                    expName = strsplit(obj.tr.GetExpName(), '_');
                    expName = strjoin(expName(1:2), '_');
                    bsr(iUnit).expName = expName;
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
        end
        
        function [selBSR, selPulse] = selectStimResponse(obj, bsr, varargin)
            p = inputParser();
            p.addRequired('BinnedStimResponse', @isstruct);
            p.addParameter('Light', [], @isnumeric);
            p.addParameter('Duration', [], @isnumeric);
            p.addParameter('MLRank', [], @isnumeric); % 1: most medial, 4: most lateral
            p.addParameter('DVRank', [], @isnumeric); % 1: most ventral, 4: most dorsal
            p.addParameter('Fiber', {}, @(x) ischar(x) || iscell(x));
            p.addParameter('Galvo', [], @isnumeric);
            p.parse(bsr, varargin{:});
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
        end

        function plotPSTHByStimCondition(obj, bsr, varargin)
            p = inputParser();
            p.addRequired('BinnedStimResponse', @isstruct);
            p.addParameter('Mode', 'both', @(x) ischar(x) && ismember(x, {'heatmap', 'line', 'both'}));
            p.addParameter('CLim', [], @isnumeric);
            p.parse(bsr, varargin{:});
            bsr = p.Results.BinnedStimResponse;
            
            if length(bsr) > 1
                for i = 1:length(bsr)
                    obj.plotPSTHByStimCondition(bsr(i), 'Mode', p.Results.Mode, 'CLim', p.Results.CLim);
                end
                return
            end
            
            % [bsr, I] = obj.selectStimResponse(bsr, 'Light', p.Results.Light, 'Duration', p.Results.Duration, 'MLRank', p.Results.MLRank, 'DVRank', p.Results.DVRank, 'Fiber', p.Results.Fiber, 'Galvo', p.Results.Galvo);
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
            
            base = max([max(lightRank), max(durationRank), max(mlRank), max(dvRank)]);
            condInt = lightRank*(base^3) + durationRank*(base^2) + mlRank*(base^1) + dvRank*(base^0);
            [uniqueConditions, ia] = unique(condInt);
            nConditions = length(uniqueConditions);
            nColorsA = length(unique(dvRank(mlRank==1)));
            nColorsB = length(unique(dvRank(mlRank==2)));
            condLabel = cell(nConditions, 1);
            condColor = cell(nConditions, 1);
            condWidth = ones(nConditions, 1);
            condSR = NaN(nConditions, size(spikeRates, 2));
            condNSR = NaN(nConditions, size(spikeRates, 2));
            for iCond = 1:nConditions
                i = ia(iCond);
                condLabel{iCond} = sprintf('%.1fmW, %.0fms (ML %i, DV %i))', light(i), duration(i)*1000, mlRank(i), dvRank(i));
                condWidth(iCond) = AcuteRecording.lerp(1, 2, (lightRank(i) - 1)/(max(lightRank) - 1));
                if mlRank(i) == 1
                    condColor{iCond} = [0.9, AcuteRecording.lerp(0.1, 0.9, (dvRank(i)-1)/(nColorsA-1)), 0.1];
                elseif mlRank(i) == 2
                    condColor{iCond} = [AcuteRecording.lerp(0.1, 0.9, (dvRank(i)-1)/(nColorsB-1)), 0.9, 0.1];
                else
                    error('ml rank must be 1 or 2, got %i instead', mlrank(i));
                end
                condSel = condInt == condInt(ia(iCond));
                condSR(iCond, :) = mean(spikeRates(condSel, :), 1);
                condNSR(iCond, :) = mean(normalizedSpikeRates(condSel, :), 1);
            end
            
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
                    plot(ax(1), bsr.tRight, condSR(iCond, :), 'Color', condColor{iCond}, 'LineWidth', condWidth(iCond), 'DisplayName', condLabel{iCond});
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
                yticklabels(ax(2), condLabel)
                ylabel(ax(2), 'Stim Condition')
                title(ax(2), sprintf('%s Chn %i Unit %i', bsr.expName, bsr.channel, bsr.unit), 'Interpreter', 'none')

                if ~isempty(p.Results.CLim)
                    ax(2).CLim = p.Results.CLim;
                end
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

        function [peaks, I] = findpeaks(x)
            df = [0, diff(x)];
            df1 = df > 0;
            df2 = df <= 0;
            df3 = df < 0;
            df4 = df >= 0;
            I = find((df1 & circshift(df2, -1)) | (df3 & circshift(df4, -1)));
            peaks = x(I);
        end

        function [stats, bsr] = summarizeStimResponse(bsr, varargin)
            p = inputParser();
            p.addRequired('BinnedStimResponse', @isstruct)
            p.addParameter('Window', [0, 0.05], @(x) isnumeric(x) && length(x) == 2)
            p.addParameter('Method', 'max', @(x) ismember(x, {'peak', 'mean', 'firstPeak'}))
            p.addParameter('Normalized', true, @islogical)
            p.parse(bsr, varargin{:});
            bsr = p.Results.BinnedStimResponse;
            window = p.Results.Window;

            stats = NaN(length(bsr), 1);
            for i = 1:length(bsr)
                sel = bsr(i).t <= window(2) & bsr(i).t >= window(1);
                if p.Results.Normalized
                    sr = bsr(i).normalizedSpikeRates(:, sel);
                else
                    sr = bsr(i).spikeRates(:, sel);
                end
                msr = mean(sr, 1);  % Mean (normalized/raw) spike rate across trials
                switch p.Results.Method
                    case 'mean'
                        stats(i) = mean(msr);
                    case 'peak'
                        [M, I] = max(abs(msr), [], 2);
                        stats(i) = M * sign(msr(I));
                    case 'firstPeak'
                        [peaks, I] = AcuteRecording.findpeaks(msr);
                        stats(i) = peaks(1);
                    otherwise
                        error('Not implemented method %s', p.Results.Method)
                end
                bsr(i).stat.value = stats(i);
                bsr(i).stat.name = p.Results.Method;
                bsr(i).stat.normalized = p.Results.Normalized;
                bsr(i).stat.window = window;
            end
        end

        function plotMap(coords, stats)
            fig = figure();
            ax = axes(fig);
            h = scatter3(ax, coords(:, 1), coords(:, 2), coords(:, 3), stats);
            xlabel('ML')
            ylabel('DV')
            zlabel('AP')
        end
    end
    
    methods (Static, Access={})
        function x = lerp(a, b, t)
            if isnan(t)
                t = 0;
            end
            t = max(0, min(1, t));
            x = a + (b - a)*t;
        end
    end
end

