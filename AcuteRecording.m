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
            
            stim.duration = (stim.tOff - stim.tOn);

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
        
        function [psth] = getPSTH(obj, channel, unit)
        end
        
        function plotPSTHByLocation(obj)
        end
        
        function plotPSTHByStim(obj)
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
    end
end

