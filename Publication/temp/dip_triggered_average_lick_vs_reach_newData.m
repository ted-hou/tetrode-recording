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
save(exportPath, 'xta', 'xta', 'kinematics', 'p', '-v7.3')

%% Load data
load('\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\LickVsReach_DTA_RTA_boot\NewData\LickVsReach_DLC_dta_rta_25_25_200to800ms_units1to1138_100boots');

%% Make a metaDTA/metaRTA
for fn = ["spikerate", "HandR", "HandL", "FootR", "FootL", "Spine", "Jaw", "Tongue"]
    X = arrayfun(@(xta) xta.(fn).X, xta.dip(selUnits), UniformOutput=false);
    X = cat(1, X{:});
    xta.dip(length(eu) + 1).(fn) = struct(X=mean(X, 1, 'omitnan'), t=xta.dip(1).(fn).t);
    X = arrayfun(@(xta) xta.(fn).X, xta.rise(selUnits), UniformOutput=false);
    X = cat(1, X{:});
    xta.rise(length(eu) + 1).(fn) = struct(X=mean(X, 1, 'omitnan'), t=xta.rise(1).(fn).t);
    clear X
end
xta.dip(length(eu) + 1).iExp = 0;
xta.rise(length(eu) + 1).iExp = 0;
xta.dip(length(eu) + 1).t0 = [xta.dip(selUnits).t0];
xta.rise(length(eu) + 1).t0 = [xta.rise(selUnits).t0];


%% Plot dip-triggered average kinematics (STD Version
close all
exportPath = fullfile("C:\SERVER\LickVsReach_DTA_RTA_boot\NewData\Figures", sprintf("LickVsReach_DLC_dta_rta_%i_%i_%ito%ims_std", 100*p.xta.dip.thresholdQuantile, 100*p.xta.dip.thresholdSubQuantile, 100*p.xta.dip.samples(1), 100*p.xta.rise.samples(2)));
if ~exist(exportPath, 'dir')
    mkdir(exportPath)
end
features = ["spikerate", "HandR", "HandL", "FootR", "FootL", "Spine", "Jaw", "Tongue"];
featureUnits = ["spike rate (a.u.)", "AP pos (a.u.)", "AP pos (a.u.)", "AP pos (a.u.)", "AP pos (a.u.)", "DV pos (a.u.)", "DV pos (a.u.)", "prob"];
statFeatures = ["HandR", "HandL", "FootR", "FootL", "Spine", "Jaw", "Tongue"];
dirs = ["dip", "rise"];
dimensionReduction = "tsne"; % pca, tsne, umap
% yl = {[-2.5, 5], [-0.1, 5.1], [-0.1, 5.1], [-0.1, 5.1], [-0.1, 5.1], [-0.1, 5.1], [-0.1, 5.1], [-0.1, 1.1]};
yl = {[-2.5, 5], [-2, 2], [-2, 2], [-2, 2], [-2, 2], [-2, 2], [-2, 2], [0, 1]};
fig = figure(Units='inches', InnerPosition=[2, 2, 1.5*(2+length(features)), 5]);
tlp = tiledlayout(fig, 2, 1, TileSpacing='compact', Padding='compact');
tl = gobjects(2, 1);
tl(1) = tiledlayout(tlp, 1, length(features) + 2, TileSpacing='compact', Padding='compact');
tl(2) = tiledlayout(tlp, 1, length(features) + 2, TileSpacing='compact', Padding='compact');
tl(1).Layout.Tile = 1;
tl(2).Layout.Tile = 2;

ax = gobjects(2, length(features) + 2);
for iDir = 1:2
    for iAx = 1:length(features) + 2
        ax(iDir, iAx) = nexttile(tl(iDir));
    end
end
for iUnit = 680:length(eu)
% for iUnit = [length(eu)+1, selUnits]
% for iUnit = length(eu)+1
    for iDir = 1:2
        for iAx = 1:length(features) + 2
            cla(ax(iDir, iAx))
        end
    end
    for iDir = 1:2
        dir = dirs(iDir);
        % Check existence
        if isempty(xta.(dir)(iUnit).t0) && xta.(dir)(iUnit).iExp > 0
            for iAx = 1:length(features)
                cla(ax(iDir, iAx))
                ax(iDir, iAx).Visible = false;
            end
            continue
        end

        % Movement diversity scatter plot
        if p.nBoot > 0 && xta.(dir)(iUnit).iExp > 0
            iAx = length(features)+1;
            if length(xta.(dir)(iUnit).t0) <= 2
                idx = ones(length(xta.(dir)(iUnit).t0), 1);
                nClusters = 1;
            else
                hold(ax(iDir, iAx), 'on')
                mr = NaN(length(xta.(dir)(iUnit).t0), length(p.std.features));% movement range: nTrials x nFeatures
                for i = 1:length(p.std.features)
                    fn = p.std.features(statFeatureOrder(i));
                    if isempty(xta.(dir)(iUnit).(fn))
                        continue
                    end
                    t = xta.(dir)(iUnit).(fn).t;
                    selT = t >= p.std.window(1) & t <= p.std.window(2);
                    xx = xta.(dir)(iUnit).(fn).X(:, selT);
                    mr(:, i) = std(xx, 0, 2, 'omitnan');
                end
                mr(isnan(mr)) = 0;
                [~, pcScore, ~, ~, explained] = pca(mr);
                switch dimensionReduction   
                    case "pca"
                        score = pcScore;
                    case "umap"
                        score = umap(mr, NumDimensions=2);
                    case "tsne"
                        score = tsne(mr);
                end
                eva = evalclusters(pcScore(:, 1:3), 'kmeans', 'CalinskiHarabasz', KList=1:3);
                nClusters = eva.OptimalK;
                clear eva
                idx = kmeans(pcScore(:, 1:3), nClusters);
                for k = 1:nClusters
                    sel = idx==k;
                    scatter(ax(iDir, iAx), score(sel, 1), score(sel, 2), 10, getColor(k, 2))
                end
                switch dimensionReduction   
                    case "pca"
                        xlabel(ax(iDir, iAx), sprintf("PC%i (%.1f%%)", 1, explained(1)))
                        ylabel(ax(iDir, iAx), sprintf("PC%i (%.1f%%)", 2, explained(2)))
                        title(ax(iDir, iAx), "PCA")
                    case "tsne"
                        xlabel(ax(iDir, iAx), sprintf("PC%i", 1))
                        ylabel(ax(iDir, iAx), sprintf("PC%i", 2))
                        title(ax(iDir, iAx), "t-SNE")
                    case "umap"
                        xlabel(ax(iDir, iAx), sprintf("PC%i", 1))
                        ylabel(ax(iDir, iAx), sprintf("PC%i", 2))
                        title(ax(iDir, iAx), "UMAP")
                end
                xticks(ax(iDir, iAx), [])
                yticks(ax(iDir, iAx), [])
            end
        end

        % Movement diversity matrix
        if p.nBoot > 0 && xta.(dir)(iUnit).iExp > 0
            iAx = length(features)+2;
            mdm = NaN(length(xta.(dir)(iUnit).t0), length(p.std.features));
            for i = 1:length(p.std.features)
                fn = p.std.features(statFeatureOrder(i));
                if isempty(xta.(dir)(iUnit).(fn))
                    continue
                end
                t = xta.(dir)(iUnit).(fn).t;
                selT = t >= p.std.window(1) & t <= p.std.window(2);
                stdObs = std(xta.(dir)(iUnit).(fn).X(:, selT), 0, 2, 'omitnan');
                mdm(:, i) = arrayfun(@(data) nnz(xta.(dir)(iUnit).(fn).stats.stdBoot < data) ./ length(xta.(dir)(iUnit).(fn).stats.stdBoot), stdObs, UniformOutput=true);
            end
            mdm(isnan(mdm)) = 0;
            hash = sum((mdm > 0.95) .* 2.^(size(mdm, 2)-1:-1:0), 2);
            hash = hash + (idx-1) .* 2.^(size(mdm, 2));
            [~, I] = sort(hash, 'ascend');
            idxSorted = idx(I);
            sepHash = arrayfun(@(idx) find(idxSorted==idx, 1, 'last'), 1:max(idx)-1);
            imagesc(ax(iDir, iAx), mdm(I, :))
            if ~isempty(sepHash)
                yline(ax(iDir, iAx), 0.5+sepHash, 'k--')
                yticks(ax(iDir, iAx), 0.5+unique([1, sepHash, length(idx)]))
                yticklabels(ax(iDir, iAx), string(unique([1, sepHash, length(idx)])))
            end
            % colormap(ax(iType, iAx), [1, 1, 1; 0, 0, 0])
            % applyCustomColormap(ax(iType, iAx), [-1, 1], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);
            % applyCustomColormap(ax(iDir, iAx), [0, 1], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.025, h0=0.33);
            applyCustomColormap(ax(iDir, iAx), [-1, 1], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.05, h0=0.33);
            % ax(iType, iAx).ColorScale = 'log';
            xticks(ax(iDir, iAx), 1:length(p.std.features))
            xticklabels(ax(iDir, iAx), p.std.features(statFeatureOrder))
            colorbar(ax(iDir, iAx), Orientation='horizontal', Location='southoutside')
            ax(iDir, iAx).XAxisLocation = 'top';
            ylabel(ax(iDir, iAx), 'Trial')
            clear mdm i fn t selT stdObs pObs hash I
        end
        
        for iAx = 1:length(features)
            fn = features(iAx);
            if isempty(xta.(dir)(iUnit).(fn))
                ax(iDir, iAx).Visible = false;
                continue
            end

            if p.nBoot > 0 && ismember(fn, p.std.features) && isfield(xta.(dir)(iUnit).(fn), 'stats')
                prcSTD = quantile(xta.(dir)(iUnit).(fn).stats.stdBoot, [0.95, 0.99, 0.999]);
                nStarsSTD = sum(xta.(dir)(iUnit).(fn).stats.std > prcSTD);
            else
                nStarsSTD = 0;
            end

            ax(iDir, iAx).Visible = true;

            hold(ax(iDir, iAx), 'on')
            hold(ax(2+1-iDir, iAx), 'on')
            t = 1e3*xta.(dir)(iUnit).(fn).t;
            X = mean(xta.(dir)(iUnit).(fn).X, 1, 'omitnan');
            plot(ax(iDir, iAx), t, X, Color=[0.15, 0.15, 0.15, 1], LineWidth=1.5, LineStyle=':');
            plot(ax(2+1-iDir, iAx), t, X, Color=[0.15, 0.15, 0.15, 0.1], LineWidth=1.5, LineStyle=':');
            for k = 1:nClusters
                % c = getColor(iAx, length(features), 0.7);
                c = getColor(k, nClusters, 0.7);
                % c = [0.15, 0.15, 0.15];
    
                X = xta.(dir)(iUnit).(fn).X(idx==k, :);
                mu = mean(X, 1, 'omitnan');
                err = std(X, 0, 1, 'omitnan')./sqrt(size(xta.(dir)(iUnit).(fn).X, 1));
    
                plot(ax(iDir, iAx), t, mu, Color=c, LineWidth=1.5);
            end
            if p.nBoot > 0 && isfield(xta.(dir)(iUnit).(fn), 'XBoot')
                prc = quantile(xta.(dir)(iUnit).(fn).XBoot, [p.bootAlpha/2, 1-p.bootAlpha/2], 1);
                patch(ax(iDir, iAx), [t, flip(t)], [prc(1, :), flip(prc(2, :))], [0.15, 0.15, 0.15], FaceAlpha=0.05, EdgeColor=[0.15, 0.15, 0.15], EdgeAlpha=0.5);
            end
            % xline(ax(iRow, iAx), 1e3*p.xta.meanWindow, 'k--', Alpha=0.1)
            xline(ax(iDir, iAx), 1e3*p.std.window, 'k--', Alpha=0.1)
            xline(ax(iDir, iAx), 0, 'k-', Alpha=0.1)
            % xticks(ax(iRow, iAx), 1e3*p.xta.meanWindow)
            xticks(ax(iDir, iAx), [-300, 0, 600])
            xtickangle(ax(iDir, iAx), 0)

            xlabel(ax(iDir, iAx), 'time (ms)')
            ylabel(ax(iDir, iAx), featureUnits(iAx))

            fnDisp = sprintf("%s %s", fn, repmat('*', [1, nStarsSTD]));
            title(ax(iDir, iAx), fnDisp, Interpreter='none')
            ylim(ax(iDir, iAx), yl{iAx})
            hold(ax(iDir, iAx), 'off')
            hold(ax(2+1-iDir, iAx), 'off')
        end

        % Correlegram
        % [lia, statFeatureOrder] = ismember(statFeatures, p.std.features);
        % assert(all(lia), 'Some members of statFeatureOrder are not found.')
        % 
        % iAx = iAx + 1;
        % r = NaN(length(p.std.features));
        % for i = 1:length(p.std.features)
        %     fni = p.std.features(statFeatureOrder(i));
        %     if isempty(xta.(dir)(iUnit).(fni))
        %         continue
        %     end
        %     for j = 1:length(p.std.features)
        %         fnj = p.std.features(statFeatureOrder(j));
        %         if isempty(xta.(dir)(iUnit).(fnj))
        %             continue
        %         end
        %         selT = xta.(dir)(iUnit).(fni).t >= p.std.window(1) & xta.(dir)(iUnit).(fni).t <= p.std.window(2);
        %         r(i, j) = corr(std(xta.(dir)(iUnit).(fni).X(:, selT), 0, 2, 'omitnan'), std(xta.(dir)(iUnit).(fnj).X(:, selT), 0, 2, 'omitnan'), Rows='complete');
        %     end
        % end
        % r(isnan(r)) = 0;
        % imagesc(ax(iDir, iAx), r);
        % xticks(ax(iDir, iAx), 1:length(p.std.features))
        % yticks(ax(iDir, iAx), 1:length(p.std.features))
        % xticklabels(ax(iDir, iAx), p.std.features(statFeatureOrder))
        % yticklabels(ax(iDir, iAx), p.std.features(statFeatureOrder))
        % xtickangle(ax(iDir, iAx), 90)
        % ax(iDir, iAx).XAxisLocation = 'top';
        % axis(ax(iDir, iAx), 'image')
        % ax(iDir, iAx).XAxis.Direction = 'normal';
        % clim(ax(iDir, iAx), [0, 1])
        % colormap(ax(iDir, iAx), 'gray')
        % colorbar(ax(iDir, iAx), 'eastoutside')
        % applyCustomColormap(ax(iDir, iAx), [-1, 1], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);
    end
    for iDir = 1:2
        dir = dirs(iDir);
        if ~isempty(xta.(dir)(iUnit).HandR)
            if xta.(dir)(iUnit).iExp > 0
                title(tl(iDir), sprintf("Unit %i (n=%i %ss)", iUnit, length(xta.(dir)(iUnit).t0), dir), FontWeight='bold')
            else
                title(tl(iDir), sprintf("%i units (n=%i %ss)", length(eu), length(xta.(dir)(iUnit).t0), dir), FontWeight='bold')
            end
        else
            xlabel(tl(iDir), '')
            title(tl(iDir), '')
        end
    end
    xlim(ax(:, 1:end-2), 1e3*[-0.5, 1])
    fontsize(fig, 9, 'points')
    print(fig, fullfile(exportPath, sprintf("dta_rta_unit_%03i", iUnit)), '-dpng', '-r0')
end
clear features featureUnits fig tlp tl ax iAx iUnit fn c t mu err iDir iFeature tlp yl fnDisp
clear prcSTD nStarsSTD i j fni fnj r selT



%% Do scatter plots of dip-triggered/rise-triggered movement magnitudes
statFeatures = ["spikerate", "HandR", "HandL", "FootR", "FootL", "Spine", "Jaw", "Tongue"];
dipWindow = [0, 0.6];
yRange = [-2, 2];
% Normalize movement magnitude by boot-strapped range
selUnits = 1:length(eu);
clear dX YRise YDip
dX(max(selUnits)) = struct();
YRise(max(selUnits)) = struct();
YDip(max(selUnits)) = struct();
for iUnit = selUnits
    for fn = statFeatures
        try
            t = xta.dip(iUnit).(fn).t;
            selT = t>=dipWindow(1) & t<=dipWindow(2);
            X = mean(xta.dip(iUnit).(fn).X(:, selT), 1, 'omitnan');
            XBoot = xta.dip(iUnit).(fn).XBoot(:, selT);
            [XAbs, I] = max(abs(X), [], 'all');
            XSign = sign(X(I));
            Xdip = XSign .* XAbs ./ range(XBoot, 'all');
            YDip(iUnit).(fn) = Xdip;
        catch
            YDip(iUnit).(fn) = NaN;
        end

        try
            t = xta.rise(iUnit).(fn).t;
            selT = t>=dipWindow(1) & t<=dipWindow(2);
            X = mean(xta.rise(iUnit).(fn).X(:, selT), 1, 'omitnan');
            XBoot = xta.rise(iUnit).(fn).XBoot(:, selT);
            [XAbs, I] = max(abs(X), [], 'all');
            XSign = sign(X(I));
            Xrise = XSign .* XAbs ./ range(XBoot, 'all');
            YRise(iUnit).(fn) = Xrise;
        catch
            YRise(iUnit).(fn) = NaN;
        end

        % dX(iUnit).(fn) = Xrise - Xdip; % dX positive -> rises gives bigger forward/upwards movement
    end
end

close all
fig = figure(Units='inches', Position=[0.5 0.5 10 10]);
tl = tiledlayout(fig, length(statFeatures), length(statFeatures), TileSpacing='tight', Padding='tight');
for i = 1:length(statFeatures)
    fni = statFeatures(i);
    for j = 1:length(statFeatures)
        fnj = statFeatures(j);
        ax = nexttile(tl);
        hold(ax, 'on')
        x = [YRise.(fnj)];
        y = [YRise.(fni)];
        sel = abs(x)>0.5 & abs(y)>0.5;
        scatter(ax, x(sel), y(sel), 2, 'black')
        xlim(ax, yRange)
        ylim(ax, yRange)
        plot(ax, yRange, yRange, 'k--')
        xline(ax, 0, 'k--')
        yline(ax, 0, 'k--')
        hold(ax, 'off')
        if i == length(statFeatures)
            xlabel(ax, fnj)
        end
        if j == 1
           ylabel(ax, fni)
        end
        axis(ax, 'square')
    end
end
title(tl, 'Rise')

fig = figure(Units='inches', Position=[10.5 0.5 10 10]);
tl = tiledlayout(fig, length(statFeatures), length(statFeatures), TileSpacing='tight', Padding='tight');
for i = 1:length(statFeatures)
    fni = statFeatures(i);
    for j = 1:length(statFeatures)
        fnj = statFeatures(j);
        ax = nexttile(tl);
        hold(ax, 'on')
        x = [YDip.(fnj)];
        y = [YDip.(fni)];
        sel = abs(x)>0.5 & abs(y)>0.5;
        scatter(ax, x(sel), y(sel), 2, 'black')
        xlim(ax, yRange)
        ylim(ax, yRange)
        plot(ax, yRange, yRange, 'k--')
        xline(ax, 0, 'k--')
        yline(ax, 0, 'k--')
        hold(ax, 'off')
        if i == length(statFeatures)
            xlabel(ax, fnj)
        end
        if j == 1
           ylabel(ax, fni)
        end
        axis(ax, 'square')
    end
end
title(tl, 'Dip')