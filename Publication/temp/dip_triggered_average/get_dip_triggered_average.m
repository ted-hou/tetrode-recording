folders = [ ...
    "C:\SERVER\daisy26\daisy26_20250424", ...
    "C:\SERVER\daisy26\daisy26_20250425", ...
    "C:\SERVER\desmond38\desmond38_20250401", ...
    "C:\SERVER\desmond38\desmond38_20250402", ...
    "C:\SERVER\desmond38\desmond38_20250407", ...
    "C:\SERVER\desmond38\desmond38_20250417", ...
    "C:\SERVER\desmond39\desmond39_20250416", ...
    "C:\SERVER\desmond39\desmond39_20250423", ...
    "C:\SERVER\daisy27\daisy27_20250624", ...
    "C:\SERVER\daisy27\daisy27_20250626", ...
    "C:\SERVER\daisy27\daisy27_20250707", ...
    "C:\SERVER\daisy27\daisy27_20250715", ...
    "C:\SERVER\daisy27\daisy27_20250717", ...
    "C:\SERVER\daisy27\daisy27_20250721", ...
    "C:\SERVER\daisy27\daisy27_20250724", ...
    "C:\SERVER\daisy28\daisy28_20250701", ...
    "C:\SERVER\daisy28\daisy28_20250702", ...
    "C:\SERVER\daisy28\daisy28_20250716", ...
    "C:\SERVER\daisy28\daisy28_20250718", ...
    "C:\SERVER\daisy28\daisy28_20250723", ...
    "C:\SERVER\daisy28\daisy28_20250725", ...
    "C:\SERVER\daisy28\daisy28_20250729", ...
    "C:\SERVER\daisy28\daisy28_20250714", ...
    "C:\SERVER\daisy29\daisy29_20251023", ...
    "C:\SERVER\daisy29\daisy29_20251024", ...
    "C:\SERVER\daisy29\daisy29_20251025", ...
    "C:\SERVER\daisy29\daisy29_20251027", ...
    "C:\SERVER\daisy29\daisy29_20251028", ...
    "C:\SERVER\daisy29\daisy29_20251031", ...
    "C:\SERVER\daisy29\daisy29_20251118", ...
    "C:\SERVER\daisy29\daisy29_20251120", ...
    "C:\SERVER\daisy29\daisy29_20251103", ...
    "C:\SERVER\daisy30\daisy30_20251029", ...
    "C:\SERVER\daisy30\daisy30_20251030", ...
    "C:\SERVER\daisy30\daisy30_20251119", ...
    "C:\SERVER\daisy30\daisy30_20251121", ...
    "C:\SERVER\desmond41\desmond41_20251028", ...
    "C:\SERVER\desmond41\desmond41_20251029", ...
    "C:\SERVER\desmond42\desmond42_20251030", ...
    "C:\SERVER\desmond42\desmond42_20251031", ...
    "C:\SERVER\desmond41\desmond41_20251104", ...
    "C:\SERVER\desmond41\desmond41_20251117", ...
    "C:\SERVER\desmond42\desmond42_20251121", ...
    "C:\SERVER\desmond41\desmond41_20251124", ...
    ];

euFolders = [ ...
    "C:\SERVER\Units\TwoColor_SNr_SCRetro\SingleUnit_NonDuplicate_NonDrift_SNr", ... daisy26, desmond38, desmond39
    "C:\SERVER\Units\TwoColor_SNr_SCRetro\ReverseInjection\SingleUnit_NonDuplicate_NonDrift_SNr_withTrials", ... daisy27, 28
    "C:\SERVER\Units\TwoColor_Striatonigral\SingleUnit_NonDuplicate_NonDrift_SNr", ... daisy29, 30, desmond41, 42
];

dlcResultsPath = 'C:\SERVER\DeepLabCut\Results\FourPawsTongueJawSpine';

%% Load ephysunits
eu = cell(length(euFolders), 1);
for iDir = 1:length(euFolders)
    eu{iDir} = EphysUnit.load(euFolders(iDir));
end

eu = cat(2, eu{:});


%% Make CompleteExperiment3

exp = CompleteExperiment3(eu, cameras='lr', deeplabcutPath=dlcResultsPath);

exp.alignTimestamps(refEventNameArduino={'OPTO1_ON', 'OPTO2_ON', 'REWARD_ON'}, refEventNameEphys={'StimOn', 'RewardOn'}, trialDurationTolerance=2);

clear results
results(length(exp)) = struct(name=[], varsL=[], hasTimestampsL=[], varsR=[], hasTimestampsR=[], isValid=[]);
for iExp = 1:length(exp)
    results(iExp).name = exp(iExp).name;
    results(iExp).varsR = string(exp(iExp).vtdR.Properties.VariableNames)';
    results(iExp).hasTimestampsR = ismember('Timestamp', exp(iExp).vtdR.Properties.VariableNames);

    if ~isempty(exp(iExp).vtdL)
        results(iExp).varsL = string(exp(iExp).vtdL.Properties.VariableNames)';
        results(iExp).hasTimestampsL = ismember('Timestamp', exp(iExp).vtdL.Properties.VariableNames);
    else
        results(iExp).hasTimestampsL = false;
    end

    results(iExp).isValid = length(results(iExp).varsL) == 23 && length(results(iExp).varsR) == 23;
end
clear iExp

% Cull bad experiments (some are still pending DLC)
assert(isequal({exp.name}, {results.name}));
exp = exp([results.isValid]);
results = results([results.isValid]);
eu = [exp.eu];

% Rename variables
from = ["HandIpsiCam", "HandContraCam", "LegIpsiCam", "LegContraCam", "Jaw", "Tongue", "Spine"];
toL = ["HandL", "HandR", "FootL", "FootR", "Jaw", "Tongue", "Spine"];
toR = ["HandR", "HandL", "FootR", "FootL", "Jaw", "Tongue", "Spine"];
for iExp = 1:length(exp)
    if results(iExp).hasTimestampsL
        vtd = exp(iExp).vtdL;
        for i = 1:length(toL)
            if ismember(sprintf("%s_X", from(i)), results(iExp).varsL)
                vtd = renamevars(vtd, ...
                    [sprintf("%s_X", from(i)), sprintf("%s_Y", from(i)), sprintf("%s_Likelihood", from(i))], ...
                    [sprintf("%s_X", toL(i)), sprintf("%s_Y", toL(i)), sprintf("%s_Likelihood", toL(i))]);
            end
        end
        exp(iExp).vtdL = vtd;
    end
    if results(iExp).hasTimestampsR
        vtd = exp(iExp).vtdR;
        for i = 1:length(toR)
            if ismember(sprintf("%s_X", from(i)), results(iExp).varsR)
                vtd = renamevars(vtd, ...
                    [sprintf("%s_X", from(i)), sprintf("%s_Y", from(i)), sprintf("%s_Likelihood", from(i))], ...
                    [sprintf("%s_X", toR(i)), sprintf("%s_Y", toR(i)), sprintf("%s_Likelihood", toR(i))]);
            end
        end
        exp(iExp).vtdR = vtd;
    end
end

%% Reset results
clearvars -except eu exp results

expIndices = cellfun(@(name) find(strcmpi(name, {exp.name}), 1, 'first'), {eu.ExpName});

%% Detect dips in firing
p.features = ["HandR", "HandL", "FootR", "FootL", "Tongue", "Jaw", "Spine"];
p.featureStats = ["xPos", "xPos", "xPos", "xPos", "likelihood", "yPos", "yPos"]; % xPos, yPos, xVel, yVel, likelihood, displacement, speed
p.vtdNames = ["vtdR", "vtdL", "vtdR", "vtdL", "both", "both", "both"];
% p.smoothWindow = [10, 10, 10, 10, 5, 10, 10];
p.smoothWindow = [20, 20, 20, 20, 5, 20, 20];
p.minL = [0.5, 0.5, 0.5, 0.5, 0.2, 0.5, 0.5];
p.spikeRes = 0.1;
p.xta.res = 1/30;
p.xta.window = [-1, 1];
p.xta.meanWindow = [-0.3, 0.3];

p.xta.dip.samples = [2, 8];
p.xta.dip.thresholdQuantile = 0.25;
p.xta.dip.thresholdSubQuantile = 0.25;
p.xta.dip.pattern = arrayfun(@(n) [0, 0, 0, ones(1, n), 0, 0, 0] , p.xta.dip.samples(1):p.xta.dip.samples(2), UniformOutput=false); % 100-300ms dips
p.xta.dip.patternOnset = cellfun(@(pat) find(pat, 1, 'first') - 1, p.xta.dip.pattern); % finds the onset


p.xta.rise.samples = [2, 8];
p.xta.rise.thresholdQuantile = 1 - p.xta.dip.thresholdQuantile;
p.xta.rise.thresholdSubQuantile = 1 - p.xta.dip.thresholdSubQuantile;
p.xta.rise.pattern = arrayfun(@(n) [0, 0, 0, ones(1, n), 0, 0, 0] , p.xta.rise.samples(1):p.xta.rise.samples(2), UniformOutput=false);
p.xta.rise.patternOnset = cellfun(@(pat) find(pat, 1, 'first') - 1, p.xta.rise.pattern);

p.blank(1).event = "StimOn";
p.blank(1).window = [-1, 1];

p.nBoot = 100;
if p.nBoot < 1000
    warning("Running bootstrap with nBoot=%i<1000 is only recommended for testing purposes. Run a real bootstrap pls you lazy bum.", p.nBoot)
end
p.bootAlpha = 0.05;

selUnits = 1:length(eu);

% Process vtdFeatues
clear kinematics
kinematics(length(exp)) = struct(HandR=[], HandL=[], FootR=[], FootL=[], Tongue=[]);
for iExp = 1:length(exp)
    for iFeature = 1:length(p.features)
        fn = p.features(iFeature);
        vn = p.vtdNames(iFeature);
        sn = p.featureStats(iFeature);
        if ismember(vn, ["vtdL", "vtdR"])
            vtd = exp(iExp).(vn);
            if isempty(vtd)
                continue
            end
        end
        % Tongue uses bilateral likelihood
        if fn == "Tongue"
            assert(sn == "likelihood");
            iSide = 0;
            t = 0:p.xta.res:max(exp(iExp).vtdL.Timestamp(end), exp(iExp).vtdR.Timestamp(end));
            L = NaN(length(t), 1);
            for side = ["vtdL", "vtdR"]
                iSide = iSide + 1;
                l = exp(iExp).(side).(sprintf("%s_Likelihood", fn));
                l(l<p.minL(iFeature)) = 0;
                l(l>p.minL(iFeature)) = 1;
                % l = smoothdata(l, 'gaussian', 7);
                L(:, iSide) = interp1(exp(iExp).(side).Timestamp, l, t, 'linear');
            end
            L = sum(L, 2);
            L = single(L > 0);
            clear iSide side l
            L = smoothdata(L, 'gaussian', p.smoothWindow(iFeature));
            kinematics(iExp).(fn) = struct(X=L, t=t);
            clear t
        % Licks from digital events
        elseif fn == "Lick"
            assert(sn == "likelihood");
            L = exp(iExp).eu(1).EventTimes.LickOn;
            t = 0:p.xta.res:exp(iExp).vtdR.Timestamp(end);
            edges = [t - p.xta.res/2, t(end) + p.xta.res/2];
            L = histcounts(L, edges);
            L = smoothdata(L, 'gaussian', p.smoothWindow(iFeature));
            kinematics(iExp).(fn) = struct(X=L, t=t);
            clear edges t
        % Average displacement from both sides
        elseif vn == "both"
            iSide = 0;
            t = 0:p.xta.res:max(exp(iExp).vtdL.Timestamp(end), exp(iExp).vtdR.Timestamp(end));
            S = NaN(length(t), 1);
            for side = ["vtdL", "vtdR"]
                iSide = iSide + 1;
                x = exp(iExp).(side).(sprintf("%s_X", fn));
                y = exp(iExp).(side).(sprintf("%s_Y", fn));
                l = exp(iExp).(side).(sprintf("%s_Likelihood", fn));
                x(l<p.minL(iFeature)) = NaN;
                y(l<p.minL(iFeature)) = NaN;
                x = (x - mean(x, 'all', 'omitnan')) ./ std(x, 0, 'all', 'omitnan');
                y = (y - mean(y, 'all', 'omitnan')) ./ std(y, 0, 'all', 'omitnan');
                tt = exp(iExp).(side).Timestamp;

                switch sn
                    case "xPos"
                        if vn == "vtdR"
                            s = x;
                        else
                            s = -x;
                        end
                    case "yPos"
                        s = y;
                    case "xVel"
                        if vn == "vtdR"
                            s = [NaN; diff(smoothdata(x, 'gaussian', p.smoothWindow(iFeature)))]./[NaN; diff(tt)];
                        else
                            s = -[NaN; diff(smoothdata(x, 'gaussian', p.smoothWindow(iFeature)))]./[NaN; diff(tt)];
                        end
                    case "yVel"
                        s = [NaN; diff(smoothdata(y, 'gaussian', p.smoothWindow(iFeature)))]./[NaN; diff(tt)];
                    case "displacement"
                        s = sqrt(x.^2 + y.^2);
                    case "speed"
                        x = [NaN; diff(smoothdata(x, 'gaussian', p.smoothWindow(iFeature)))]./[NaN; diff(tt)];
                        y = [NaN; diff(smoothdata(y, 'gaussian', p.smoothWindow(iFeature)))]./[NaN; diff(tt)];
                        s = sqrt(x.^2 + y.^2);
                    otherwise
                        error("Unsupported stat '%s' for feature '%s'", sn, fn)
                end

                % l = smoothdata(l, 'gaussian', 7);
                S(:, iSide) = interp1(tt, s, t, 'linear');
            end
            clear iSide side x y l d
            S = mean(S, 2, 'omitnan');
            if ~ismember(sn, ["xVel", "yVel", "speed"])
                S = smoothdata(S, 'gaussian', p.smoothWindow(iFeature));
            end
            kinematics(iExp).(fn) = struct(X=S, t=t);
            clear t
        % Other tracking points use position or speed
        elseif ismember(sprintf("%s_X", fn), vtd.Properties.VariableNames)
            x = vtd.(sprintf("%s_X", fn));
            y = vtd.(sprintf("%s_Y", fn));
            l = vtd.(sprintf("%s_Likelihood", fn));
            t = vtd.Timestamp;
            x(l<p.minL(iFeature)) = NaN;
            y(l<p.minL(iFeature)) = NaN;
            x = (x - mean(x, 'all', 'omitnan')) ./ std(x, 0, 'all', 'omitnan');
            y = (y - mean(y, 'all', 'omitnan')) ./ std(y, 0, 'all', 'omitnan');
            switch sn
                case "xPos"
                    if vn == "vtdR"
                        s = x;
                    else
                        s = -x;
                    end
                case "yPos"
                    s = y;
                case "xVel"
                    if vn == "vtdR"
                        s = [NaN; diff(smoothdata(x, 'gaussian', p.smoothWindow(iFeature)))]./[NaN; diff(t)];
                    else
                        s = -[NaN; diff(smoothdata(x, 'gaussian', p.smoothWindow(iFeature)))]./[NaN; diff(t)];
                    end
                case "yVel"
                    s = [NaN; diff(smoothdata(y, 'gaussian', p.smoothWindow(iFeature)))]./[NaN; diff(t)];
                case "displacement"
                    s = sqrt(x.^2 + y.^2);
                case "speed"
                    x = [NaN; diff(smoothdata(x, 'gaussian', p.smoothWindow(iFeature)))]./[NaN; diff(t)];
                    y = [NaN; diff(smoothdata(y, 'gaussian', p.smoothWindow(iFeature)))]./[NaN; diff(t)];
                    s = sqrt(x.^2 + y.^2);
                otherwise
                    error("Unsupported stat '%s' for feature '%s'", sn, fn)
            end
            % selnan = isnan(D);
            % D(selnan) = interp1(vtd.Timestamp(~selnan), D(~selnan), vtd.Timestamp(selnan), 'linear');
            if ~ismember(sn, ["xVel", "yVel", "speed"])
                s = smoothdata(s, 'gaussian', p.smoothWindow(iFeature));
            end
            kinematics(iExp).(fn) = struct(X=s, t=vtd.Timestamp);
        end
    end
end
clear iExp iFeature vn fn vtd L X Y S selnan

% Process dip/rise-triggered kinematics
tLocal = p.xta.window(1):p.xta.res:p.xta.window(2);
clear xta
xta.dip(length(eu)) = struct(iExp=[], params=[], t0=[], spikerate=[]);
xta.rise(length(eu)) = struct(iExp=[], params=[], t0=[], spikerate=[]);
lineLength = 0;
tTicTotal = tic();
rng(42)
if isempty(gcp('nocreate'))
    parpool('Processes');
end
for iEu = selUnits
    iExp = expIndices(iEu);

    % Get z-scored whole-session spike rates
    [x, t] = eu(iEu).getSpikeCounts(p.spikeRes);
    x = double(x)./p.spikeRes;
    mu = mean(x);
    sd = std(x, 0);
    x = (x-mu)/sd;

    t0 = struct(dip=[], rise=[]);
    nTotal = struct(dip=[], rise=[]);
    threshold = struct(dip=[], rise=[]);
    % Find onset of dips and rises
    for dir = ["dip", "rise"]
        switch dir
            case "dip"
                threshold.(dir) = quantile(x(x<0), p.xta.(dir).thresholdQuantile);
                i0 = arrayfun(@(i) strfind(x<=threshold.(dir), p.xta.(dir).pattern{i}) + p.xta.(dir).patternOnset(i), 1:length(p.xta.(dir).pattern), UniformOutput=false);
            case "rise"
                threshold.(dir) = quantile(x(x>0), p.xta.(dir).thresholdQuantile);
                i0 = arrayfun(@(i) strfind(x>=threshold.(dir), p.xta.(dir).pattern{i}) + p.xta.(dir).patternOnset(i), 1:length(p.xta.(dir).pattern), UniformOutput=false);
        end
        i0 = cat(2, i0{:});
        t0.(dir) = t(i0);
        nTotal.(dir) = length(t0.(dir));

        % Blank out dips around certain behavioral/stim events (opto, etc.)
        for iEvent = 1:length(p.blank)
            tEvent = eu(iEu).EventTimes.(p.blank(iEvent).event);
            windows = tEvent(:) + p.blank(iEvent).window;
            for i = 1:length(tEvent)
                t0.(dir)(t0.(dir)>=windows(i, 1) & t0.(dir)<=windows(i, 2)) = [];
            end
        end
        clear iEvent tEvent windows i i0
    end

    fprintf(repmat('\b', [1, lineLength]))
    lineLength = fprintf('Unit %i/%i; %i(-%i) dips (x<%.2f), %i(-%i) rises (x>%.2f)... %.1fs elapsed...', iEu, length(eu), length(t0.dip), nTotal.dip-length(t0.dip), threshold.dip, length(t0.rise), nTotal.rise-length(t0.rise), threshold.rise, toc(tTicTotal));

    for dir = ["dip", "rise"]
        xta.(dir)(iEu).iExp = iExp;
        xta.(dir)(iEu).params = p;
        xta.(dir)(iEu).t0 = t0.(dir);

        t0Temp = t0.(dir);
        if ~isempty(t0Temp)
            T = tLocal + t0Temp';
            % Spike rates
            xta.(dir)(iEu).spikerate = struct(t=tLocal, X=NaN(length(t0Temp), length(tLocal)));
            for iDip = 1:length(t0Temp)
                xta.(dir)(iEu).spikerate.X(iDip, :) = interp1(t, x, T(iDip, :), 'linear');
            end
            % Subselect dips/rises by quantile
            if ~isnan(p.xta.(dir).thresholdSubQuantile)
                magnitude = mean(xta.(dir)(iEu).spikerate.X(:, xta.(dir)(iEu).spikerate.t > p.xta.meanWindow(1) & xta.(dir)(iEu).spikerate.t < p.xta.meanWindow(2)), 2, 'omitnan');
                switch dir
                    case "dip"
                        sel = magnitude < quantile(magnitude, p.xta.(dir).thresholdSubQuantile);
                    case "rise"
                        sel = magnitude > quantile(magnitude, p.xta.(dir).thresholdSubQuantile);
                end
                t0Temp = t0Temp(sel);
                T = T(sel, :);
                xta.(dir)(iEu).spikerate.X = xta.(dir)(iEu).spikerate.X(sel, :);
                xta.(dir)(iEu).t0 = t0Temp;
                lineLength = lineLength + fprintf('\tSubselecting %i/%i %ss...', nnz(sel), length(sel), dir);
            end

            % dip/rise-triggered kinematics
            for fn = p.features
                xta.(dir)(iEu).(fn) = struct(t=tLocal, X=NaN(length(t0Temp), length(tLocal)));
            end
            for fn = p.features
                if isempty(kinematics(iExp).(fn))
                    xta.(dir)(iEu).(fn) = [];
                    continue
                end
                for iDip = 1:length(t0Temp)
                    xta.(dir)(iEu).(fn).X(iDip, :) = interp1(kinematics(iExp).(fn).t, kinematics(iExp).(fn).X, T(iDip, :), 'linear'); % Consider doing shifts instead of interp1 for faster bootstraping
                end
            end

            % Bootstrap dip-triggered kinematics
            if p.nBoot > 0
                XBoot = NaN(p.nBoot, length(tLocal), length(p.features));
                IT = 1:length(tLocal);
                IFEATURE = 1:length(p.features);
                maxT = kinematics(iExp).HandR.t(end);
                kinematicsTemp = kinematics(iExp);
                tTic = tic();
                parfor iBoot = 1:p.nBoot   
                    t0Boot = rand([length(t0Temp), 1]) * maxT;
                    T = tLocal + t0Boot;
                    for iFeature = IFEATURE
                        fn = p.features(iFeature);
                        if isempty(kinematicsTemp.(fn))
                            continue
                        end
                        X = NaN(length(t0Boot), length(tLocal));
                        for iDip = 1:length(t0Boot)
                            X(iDip, :) = interp1(kinematicsTemp.(fn).t, kinematicsTemp.(fn).X, T(iDip, :), 'linear');
                        end
                        XBoot(iBoot, IT, iFeature) = mean(X, 1, 'omitnan');
                    end
                end
                for iFeature = 1:length(p.features)
                    fn = p.features(iFeature);
                    if isempty(kinematicsTemp.(fn))
                        continue
                    end
                    xta.(dir)(iEu).(fn).XBoot = XBoot(:, :, iFeature);
                end
                lineLength = lineLength + fprintf('%.1fs;', toc(tTic));
                clear iBoot t0Boot T X XBoot IT IFEATURE maxT kinematicsTemp tTic
            end
        end        
    end
    lineLength = lineLength + fprintf('\n');
end

clear iEu iExp lineLength tLocal tTicTotal x t mu sd threshold iRise t0Temp riseThreshold iRise tRise T magnitude sel riseMagnitude selRises
clear lineLength2 I2 iDip iFeature fn nDipsTotal nRisesTotal


% Post-hoc do a bootstrap for xta.dip/xta.rise spikerate
lineLength = 0;
tTicTotal = tic();
tLocal = p.xta.window(1):p.xta.res:p.xta.window(2);
if p.nBoot > 0
    for iEu = selUnits
        fprintf(repmat('\b', [1, lineLength]))
        lineLength = fprintf('Unit %i/%i... %.1fs;', iEu, length(eu), toc(tTicTotal));

        % Whole session spike rates
        [x, t] = eu(iEu).getSpikeCounts(p.spikeRes);
        x = double(x)./p.spikeRes;
        mu = mean(x);
        sd = std(x, 0);
        x = (x-mu)/sd;

        maxT = t(end);
        IT = 1:length(tLocal);
        XBoot = NaN(p.nBoot, length(tLocal));

        for dir = ["dip", "rise"]
            if ~isempty(xta.(dir)(iEu).t0)
                t0Temp = xta.(dir)(iEu).t0;
                lineLength = lineLength + fprintf(' %i %ss...', length(t0Temp), dir);
                tTic = tic();
                parfor iBoot = 1:p.nBoot   
                    t0Boot = rand([length(t0Temp), 1]) * maxT;
                    T = tLocal + t0Boot;
                    X = NaN(length(t0Boot), length(tLocal));
                    for i0 = 1:length(t0Boot)
                        X(i0, :) = interp1(t, x, T(i0, :), 'linear');
                    end
                    XBoot(iBoot, IT) = mean(X, 1, 'omitnan');
                end
                xta.(dir)(iEu).spikerate.XBoot = XBoot(:, :);
                lineLength = lineLength + fprintf('%.1fs;', toc(tTic));
            end
        end
    end
end
fprintf('\n')
clear tLocal iEu x t mu sd maxT IT XBoot t0Temp tRise tTic iBoot t0Boot T X i0 iRise lineLength tTic tTicTotal


% Post-hoc calculate the bootstrpped 95%CI of the STD of the displacement traces
if isfield(p, 'rms')
    p = rmfield(p, 'rms');
end
p.std.features = ["HandR", "HandL", "FootR", "FootL", "Tongue", "Jaw", "Spine"];
p.std.window = [-0.3, 0.6];

for iUnit = selUnits
    for dir = ["dip", "rise"]
        for fn = p.std.features
            if isempty(xta.(dir)(iUnit).(fn))
                continue
            end
            selT = xta.(dir)(iUnit).(fn).t >= p.std.window(1) & xta.(dir)(iUnit).(fn).t <= p.std.window(2);
            xta.(dir)(iUnit).(fn).stats.std = std(mean(xta.(dir)(iUnit).(fn).X(:, selT), 1, 'omitnan'), 0, 2, 'omitnan');
            if p.nBoot > 0
                xta.(dir)(iUnit).(fn).stats.stdBoot = std(xta.(dir)(iUnit).(fn).XBoot(:, selT), 0, 2, 'omitnan');
            else
                xta.(dir)(iUnit).(fn).stats.stdBoot = [];
            end
        end
    end
end

clear iUnit fn selT

% Save results
exportPath = fullfile("C:\SERVER\LickVsReach_DTA_RTA_boot\NewData", sprintf("LickVsReach_DLC_dta_rta_%i_%i_%ito%ims_units%ito%i_%iboots.mat", 100*p.xta.dip.thresholdQuantile, 100*p.xta.dip.thresholdSubQuantile, 100*p.xta.dip.samples(1), 100*p.xta.rise.samples(2), selUnits(1), selUnits(end), p.nBoot));
save(exportPath, 'xta', 'kinematics', 'p', '-v7.3')
fprintf("Saved to %s\n", exportPath);
