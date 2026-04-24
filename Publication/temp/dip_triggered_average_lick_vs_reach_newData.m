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
p.dtaRes = 1/30;

p.dipSamples = [2, 8];
p.dipThresholdQuantile = 0.25;
p.dipThresholdSubQuantile = 0.25;
p.dipPattern = arrayfun(@(n) [0, 0, 0, ones(1, n), 0, 0, 0] , p.dipSamples(1):p.dipSamples(2), UniformOutput=false); % 100-300ms dips
p.dipPatternOnset = cellfun(@(pat) find(pat, 1, 'first') - 1, p.dipPattern); % finds the onset
p.dtaWindow = [-1, 1];
p.mdtaWindow = [-0.3, 0.3];

p.riseSamples = [2, 8];
p.riseThresholdQuantile = 1 - p.dipThresholdQuantile;
p.riseThresholdSubQuantile = 1 - p.dipThresholdSubQuantile;
p.risePattern = arrayfun(@(n) [0, 0, 0, ones(1, n), 0, 0, 0] , p.riseSamples(1):p.riseSamples(2), UniformOutput=false);
p.risePatternOnset = cellfun(@(pat) find(pat, 1, 'first') - 1, p.risePattern);

p.blank(1).event = "StimOn";
p.blank(1).window = [-1, 1];

p.nBoot = 0;
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
            t = 0:p.dtaRes:max(exp(iExp).vtdL.Timestamp(end), exp(iExp).vtdR.Timestamp(end));
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
            t = 0:p.dtaRes:exp(iExp).vtdR.Timestamp(end);
            edges = [t - p.dtaRes/2, t(end) + p.dtaRes/2];
            L = histcounts(L, edges);
            L = smoothdata(L, 'gaussian', p.smoothWindow(iFeature));
            kinematics(iExp).(fn) = struct(X=L, t=t);
            clear edges t
        % Average displacement from both sides
        elseif vn == "both"
            iSide = 0;
            t = 0:p.dtaRes:max(exp(iExp).vtdL.Timestamp(end), exp(iExp).vtdR.Timestamp(end));
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
tLocal = p.dtaWindow(1):p.dtaRes:p.dtaWindow(2);
clear dta rta
dta(length(eu)) = struct(iExp=[], params=[], tDip=[], spikerate=[]);
rta(length(eu)) = struct(iExp=[], params=[], tRise=[], spikerate=[]);
lineLength = 0;
tTicTotal = tic();
rng(42)
if isempty(gcp('nocreate'))
    parpool('Processes');
end
for iEu = selUnits
    iExp = expIndices(iEu);

    % Get z-scored whole-session spike rates
    [x, t] = eu(iEu).getSpikeCounts(0.1);
    x = double(x)./0.1;
    mu = mean(x);
    sd = std(x, 0);
    x = (x-mu)/sd;

    % Find dips
    dipThreshold = quantile(x(x<0), p.dipThresholdQuantile);
    iDip = arrayfun(@(i) strfind(x<=dipThreshold, p.dipPattern{i}) + p.dipPatternOnset(i), 1:length(p.dipPattern), UniformOutput=false);
    iDip = cat(2, iDip{:});
    tDip = t(iDip);
    nDipsTotal = length(tDip);

    for iEvent = 1:length(p.blank)
        tEvent = eu(iEu).EventTimes.(p.blank(iEvent).event);
        windows = tEvent(:) + p.blank(iEvent).window;
        for i = 1:length(tEvent)
            tDip(tDip>=windows(i, 1) & tDip<=windows(i, 2)) = [];
        end
    end
    clear iEvent tEvent windows i iDip

    % Find rises
    riseThreshold = quantile(x(x>0), p.riseThresholdQuantile);
    iRise = arrayfun(@(i) strfind(x>=riseThreshold, p.risePattern{i}) + p.risePatternOnset(i), 1:length(p.risePattern), UniformOutput=false);
    iRise = cat(2, iRise{:});
    tRise = t(iRise);
    nRisesTotal = length(tRise);

    for iEvent = 1:length(p.blank)
        tEvent = eu(iEu).EventTimes.(p.blank(iEvent).event);
        windows = tEvent(:) + p.blank(iEvent).window;
        for i = 1:length(tEvent)
            tRise(tRise>=windows(i, 1) & tRise<=windows(i, 2)) = [];
        end
    end
    clear iEvent tEvent windows i iRise


    fprintf(repmat('\b', [1, lineLength]))
    lineLength = fprintf('Unit %i/%i; %i(-%i) dips (x<%.2f), %i(-%i) rises (x>%.2f)... %.1fs elapsed...', iEu, length(eu), length(tDip), nDipsTotal-length(tDip), dipThreshold, length(tRise), nRisesTotal-length(tRise), riseThreshold, toc(tTicTotal));

    dta(iEu).iExp = iExp;
    rta(iEu).iExp = iExp;
    dta(iEu).params = p;
    rta(iEu).params = p;
    dta(iEu).tDip = tDip;
    rta(iEu).tRise = tRise;

    if ~isempty(tDip)
        T = tLocal + tDip';
        % Spike rates
        dta(iEu).spikerate = struct(t=tLocal, X=NaN(length(tDip), length(tLocal)));
        for iDip = 1:length(tDip)
            dta(iEu).spikerate.X(iDip, :) = interp1(t, x, T(iDip, :), 'linear');
        end
        % Subselect dips by quantile
        if ~isnan(p.dipThresholdSubQuantile)
            dipMagnitude = mean(dta(iEu).spikerate.X(:, dta(iEu).spikerate.t > p.mdtaWindow(1) & dta(iEu).spikerate.t < p.mdtaWindow(2)), 2, 'omitnan');
            selDips = dipMagnitude < quantile(dipMagnitude, p.dipThresholdSubQuantile);
            tDip = tDip(selDips);
            T = T(selDips, :);
            dta(iEu).spikerate.X = dta(iEu).spikerate.X(selDips, :);
            dta(iEu).tDip = tDip;
            lineLength = lineLength + fprintf('\tSubselecting %i/%i dips...', nnz(selDips), length(selDips));
        end

        % dip-triggered kinematics
        for fn = p.features
            dta(iEu).(fn) = struct(t=tLocal, X=NaN(length(tDip), length(tLocal)));
        end
        for fn = p.features
            if isempty(kinematics(iExp).(fn))
                dta(iEu).(fn) = [];
                continue
            end
            for iDip = 1:length(tDip)
                dta(iEu).(fn).X(iDip, :) = interp1(kinematics(iExp).(fn).t, kinematics(iExp).(fn).X, T(iDip, :), 'linear'); % Consider doing shifts instead of interp1 for faster bootstraping
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
                tDipBoot = rand([length(tDip), 1]) * maxT;
                T = tLocal + tDipBoot;
                for iFeature = IFEATURE
                    fn = p.features(iFeature);
                    if isempty(kinematicsTemp.(fn))
                        continue
                    end
                    X = NaN(length(tDipBoot), length(tLocal));
                    for iDip = 1:length(tDipBoot)
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
                dta(iEu).(fn).XBoot = XBoot(:, :, iFeature);
            end
            lineLength = lineLength + fprintf('%.1fs;', toc(tTic));
            clear iBoot tDipBoot T X XBoot IT IFEATURE maxT kinematicsTemp tTic
        end
    end
   
    if ~isempty(tRise)
        T = tLocal + tRise';
        % Spike rates
        rta(iEu).spikerate = struct(t=tLocal, X=NaN(length(tRise), length(tLocal)));
        for iRise = 1:length(tRise)
            rta(iEu).spikerate.X(iRise, :) = interp1(t, x, T(iRise, :), 'linear');
        end
        % Subselect rises by quantile
        if ~isnan(p.riseThresholdSubQuantile)
            riseMagnitude = mean(rta(iEu).spikerate.X(:, rta(iEu).spikerate.t >= p.mdtaWindow(1) & rta(iEu).spikerate.t <= p.mdtaWindow(2)), 2, 'omitnan');
            selRises = riseMagnitude > quantile(riseMagnitude, p.riseThresholdSubQuantile);
            tRise = tRise(selRises);
            rta(iEu).spikerate.X = rta(iEu).spikerate.X(selRises, :);
            rta(iEu).tRise = tRise;
            lineLength = lineLength + fprintf('\tSubselecting %i/%i rises...', nnz(selRises), length(selRises));
        end

        % rise-triggered kinematics
        for fn = p.features
            rta(iEu).(fn) = struct(t=tLocal, X=NaN(length(tRise), length(tLocal)));
        end
        for fn = p.features
            if isempty(kinematics(iExp).(fn))
                rta(iEu).(fn) = [];
                continue
            end
            for iRise = 1:length(tRise)
                rta(iEu).(fn).X(iRise, :) = interp1(kinematics(iExp).(fn).t, kinematics(iExp).(fn).X, T(iRise, :), 'linear');
            end
        end

        % Bootstrap rise-triggered kinematics
        if p.nBoot > 0
            XBoot = NaN(p.nBoot, length(tLocal), length(p.features));
            IT = 1:length(tLocal);
            IFEATURE = 1:length(p.features);
            maxT = kinematics(iExp).HandR.t(end);
            kinematicsTemp = kinematics(iExp);
            tTic = tic();
            parfor iBoot = 1:p.nBoot   
                tRiseBoot = rand([length(tRise), 1]) * maxT;
                T = tLocal + tRiseBoot;
                for iFeature = IFEATURE
                    fn = p.features(iFeature);
                    if isempty(kinematicsTemp.(fn))
                        continue
                    end
                    X = NaN(length(tRiseBoot), length(tLocal));
                    for iRise = 1:length(tRiseBoot)
                        X(iRise, :) = interp1(kinematicsTemp.(fn).t, kinematicsTemp.(fn).X, T(iRise, :), 'linear');
                    end
                    XBoot(iBoot, IT, iFeature) = mean(X, 1, 'omitnan');
                end
            end
            for iFeature = 1:length(p.features)
                fn = p.features(iFeature);
                if isempty(kinematicsTemp.(fn))
                    continue
                end
                rta(iEu).(fn).XBoot = XBoot(:, :, iFeature);
            end
            lineLength = lineLength + fprintf('%.1fs;', toc(tTic));
            clear iBoot tRiseBoot T X XBoot IT IFEATURE maxT kinematicsTemp tTic
        end
    end
    lineLength = lineLength + fprintf('\n');
end

clear iEu iExp lineLength tLocal tTicTotal x t mu sd dipThreshold iRise tDip riseThreshold iRise tRise T dipMagnitude selDips riseMagnitude selRises
clear lineLength2 I2 iDip iFeature fn nDipsTotal nRisesTotal


% Post-hoc do a bootstrap for dta/rta spikerate
lineLength = 0;
tTicTotal = tic();
tLocal = p.dtaWindow(1):p.dtaRes:p.dtaWindow(2);
if p.nBoot > 0
    for iEu = selUnits
        fprintf(repmat('\b', [1, lineLength]))
        lineLength = fprintf('Unit %i/%i... %.1fs;', iEu, length(eu), toc(tTicTotal));

        % Whole session spike rates
        [x, t] = eu(iEu).getSpikeCounts(0.1);
        x = double(x)./0.1;
        mu = mean(x);
        sd = std(x, 0);
        x = (x-mu)/sd;

        maxT = t(end);
        IT = 1:length(tLocal);
        XBoot = NaN(p.nBoot, length(tLocal));

        if ~isempty(dta(iEu).tDip)
            tDip = dta(iEu).tDip;
            lineLength = lineLength + fprintf(' %i dips...', length(tDip));
            tTic = tic();
            parfor iBoot = 1:p.nBoot   
                tDipBoot = rand([length(tDip), 1]) * maxT;
                T = tLocal + tDipBoot;
                X = NaN(length(tDipBoot), length(tLocal));
                for iDip = 1:length(tDipBoot)
                    X(iDip, :) = interp1(t, x, T(iDip, :), 'linear');
                end
                XBoot(iBoot, IT) = mean(X, 1, 'omitnan');
            end
            dta(iEu).spikerate.XBoot = XBoot(:, :);
            lineLength = lineLength + fprintf('%.1fs;', toc(tTic));
        end

        if ~isempty(rta(iEu).tRise)
            tRise = rta(iEu).tRise;
            lineLength = lineLength + fprintf(' %i rises...', length(tRise));
            tTic = tic();
            parfor iBoot = 1:p.nBoot   
                tRiseBoot = rand([length(tRise), 1]) * maxT;
                T = tLocal + tRiseBoot;
                X = NaN(length(tRiseBoot), length(tLocal));
                for iRise = 1:length(tRiseBoot)
                    X(iRise, :) = interp1(t, x, T(iRise, :), 'linear');
                end
                XBoot(iBoot, IT) = mean(X, 1, 'omitnan');
            end
            rta(iEu).spikerate.XBoot = XBoot(:, :);
            lineLength = lineLength + fprintf('%.1fs;', toc(tTic));
        end
    end
end
fprintf('\n')
clear tLocal iEu x t mu sd maxT IT XBoot tDip tRise tTic iBoot tDipBoot T X iDip iRise lineLength tTic tTicTotal


% Post-hoc calculate the bootstrpped 95%CI of the STD of the displacement traces
if isfield(p, 'rms')
    p = rmfield(p, 'rms');
end
p.std.features = ["HandR", "HandL", "FootR", "FootL", "Tongue", "Jaw", "Spine"];
p.std.window = [-0.3, 0.6];

for iUnit = selUnits
    for fn = p.std.features
        if isempty(dta(iUnit).(fn))
            continue
        end
        selT = dta(iUnit).(fn).t >= p.std.window(1) & dta(iUnit).(fn).t <= p.std.window(2);
        dta(iUnit).(fn).stats.std = std(mean(dta(iUnit).(fn).X(:, selT), 1, 'omitnan'), 0, 2, 'omitnan');
        if p.nBoot > 0
            dta(iUnit).(fn).stats.stdBoot = std(dta(iUnit).(fn).XBoot(:, selT), 0, 2, 'omitnan');
        else
            dta(iUnit).(fn).stats.stdBoot = [];
        end
    end
    for fn = p.std.features
        if isempty(rta(iUnit).(fn))
            continue
        end
        selT = rta(iUnit).(fn).t >= p.std.window(1) & rta(iUnit).(fn).t <= p.std.window(2);
        rta(iUnit).(fn).stats.std = std(mean(rta(iUnit).(fn).X(:, selT), 1, 'omitnan'), 0, 2, 'omitnan');
        if p.nBoot > 0
            rta(iUnit).(fn).stats.stdBoot = std(rta(iUnit).(fn).XBoot(:, selT), 0, 2, 'omitnan');
        else
            rta(iUnit).(fn).stats.stdBoot = [];
        end
    end
end

clear iUnit fn selT

%% Save results
exportPath = fullfile("C:\SERVER\LickVsReach_DTA_RTA_boot\NewData", sprintf("LickVsReach_DLC_dta_rta_%i_%i_%ito%ims_units%ito%i_%iboots.mat", 100*p.dipThresholdQuantile, 100*p.dipThresholdSubQuantile, 100*p.dipSamples(1), 100*p.riseSamples(2), selUnits(1), selUnits(end), p.nBoot));
save(exportPath, 'dta', 'rta', 'kinematics', 'p', '-v7.3')

%% Make a metaDTA/metaRTA
for fn = ["spikerate", "HandR", "HandL", "FootR", "FootL", "Spine", "Jaw", "Tongue"]
    X = arrayfun(@(xta) xta.(fn).X, dta(selUnits), UniformOutput=false);
    X = cat(1, X{:});
    dta(length(eu) + 1).(fn) = struct(X=mean(X, 1, 'omitnan'), t=dta(1).(fn).t);
    X = arrayfun(@(xta) xta.(fn).X, rta(selUnits), UniformOutput=false);
    X = cat(1, X{:});
    rta(length(eu) + 1).(fn) = struct(X=mean(X, 1, 'omitnan'), t=rta(1).(fn).t);
    clear X
end
dta(length(eu) + 1).iExp = 0;
rta(length(eu) + 1).iExp = 0;
dta(length(eu) + 1).tDip = [dta(selUnits).tDip];
rta(length(eu) + 1).tRise = [rta(selUnits).tRise];

%% Plot dip-triggered average kinematics (STD Version
close all
exportPath = fullfile("C:\SERVER\LickVsReach_DTA_RTA_boot\NewData\Figures", sprintf("LickVsReach_DLC_dta_rta_%i_%i_%ito%ims_std", 100*p.dipThresholdQuantile, 100*p.dipThresholdSubQuantile, 100*p.dipSamples(1), 100*p.riseSamples(2)));
if ~exist(exportPath, 'dir')
    mkdir(exportPath)
end
features = ["spikerate", "HandR", "HandL", "FootR", "FootL", "Spine", "Jaw", "Tongue"];
featureUnits = ["a.u.", "AP pos (a.u.)", "AP pos (a.u.)", "AP pos (a.u.)", "AP pos (a.u.)", "DV pos (a.u.)", "DV pos (a.u.)", "prob"];
statFeatures = ["HandR", "HandL", "FootR", "FootL", "Spine", "Jaw", "Tongue"];
% yl = {[-2.5, 5], [-0.1, 5.1], [-0.1, 5.1], [-0.1, 5.1], [-0.1, 5.1], [-0.1, 5.1], [-0.1, 5.1], [-0.1, 1.1]};
yl = {[-2.5, 5], [-2, 2], [-2, 2], [-2, 2], [-2, 2], [-2, 2], [-2, 2], [0, 1]};
fig = figure(Units='inches', InnerPosition=[2, 2, 2*(2+length(features)), 5]);
tlp = tiledlayout(fig, 2, 1, TileSpacing='compact', Padding='compact');
tl = gobjects(2, 1);
tl(1) = tiledlayout(tlp, 1, length(features) + 2, TileSpacing='compact', Padding='compact');
tl(2) = tiledlayout(tlp, 1, length(features) + 2, TileSpacing='compact', Padding='compact');
tl(1).Layout.Tile = 1;
tl(2).Layout.Tile = 2;

ax = gobjects(2, length(features) + 2);
for iType = 1:2
    for iAx = 1:length(features) + 2
        ax(iType, iAx) = nexttile(tl(iType));
    end
end
for iUnit = [length(eu)+1, selUnits]
% for iUnit = length(eu)+1
    for iType = 1:2
        for iAx = 1:length(features)
            cla(ax(iType, iAx))
        end
    end
    for iType = 1:2
        % Check existence
        if iType == 1
            if isempty(dta(iUnit).tDip) && dta(iUnit).iExp > 0
                for iAx = 1:length(features)
                    cla(ax(iType, iAx))
                    ax(iType, iAx).Visible = false;
                end
                continue
            end
        else
            if isempty(rta(iUnit).tRise) && dta(iUnit).iExp > 0
                for iAx = 1:length(features)
                    cla(ax(iType, iAx))
                    ax(iType, iAx).Visible = false;
                end
                continue
            end
        end
        for iAx = 1:length(features)
            fn = features(iAx);
            if iType == 1 && isempty(dta(iUnit).(fn))
                ax(iType, iAx).Visible = false;
                continue
            elseif iType == 2 && isempty(rta(iUnit).(fn))
                ax(iType, iAx).Visible = false;
                continue
            end

            ax(iType, iAx).Visible = true;

            hold(ax(iType, iAx), 'on')
            c = getColor(iAx, length(features), 0.7);
            if iType == 1
                t = 1e3*dta(iUnit).(fn).t;
                X = dta(iUnit).(fn).X;
                mu = mean(dta(iUnit).(fn).X, 1, 'omitnan');
                err = std(dta(iUnit).(fn).X, 0, 1, 'omitnan')./sqrt(size(dta(iUnit).(fn).X, 1));
                if p.nBoot > 0
                    prc = quantile(dta(iUnit).(fn).XBoot, [p.bootAlpha/2, 1-p.bootAlpha/2], 1);
                end

                if p.nBoot > 0 && ismember(fn, p.std.features)
                    prcSTD = quantile(dta(iUnit).(fn).stats.stdBoot, [0.95, 0.99, 0.999]);
                    nStarsSTD = sum(dta(iUnit).(fn).stats.std > prcSTD);
                else
                    nStarsSTD = 0;
                end
            else
                t = 1e3*rta(iUnit).(fn).t;
                X = rta(iUnit).(fn).X;
                mu = mean(rta(iUnit).(fn).X, 1, 'omitnan');
                err = std(rta(iUnit).(fn).X, 0, 1, 'omitnan')./sqrt(size(rta(iUnit).(fn).X, 1));
                if p.nBoot > 0
                    prc = quantile(rta(iUnit).(fn).XBoot, [p.bootAlpha/2, 1-p.bootAlpha/2], 1);
                end

                if p.nBoot > 0 && ismember(fn, p.std.features)
                    prcSTD = quantile(rta(iUnit).(fn).stats.stdBoot, [0.95, 0.99, 0.999]);
                    nStarsSTD = sum(rta(iUnit).(fn).stats.std > prcSTD);
                else
                    nStarsSTD = 0;
                end
            end
            % plot(ax(iType, iAx), t, X, Color=[0.15, 0.15, 0.15, 0.1]);
            plot(ax(2+1-iType, iAx), t, mu, Color=[0.15, 0.15, 0.15, 0.1], LineWidth=1.5);
            plot(ax(iType, iAx), t, mu, Color=c, LineWidth=1.5);
            if p.nBoot > 0
                patch(ax(iType, iAx), [t, flip(t)], [prc(1, :), flip(prc(2, :))], c, FaceAlpha=0.05, EdgeColor=c, EdgeAlpha=0.5);
            else
                patch(ax(iType, iAx), [t, flip(t)], [mu-err, flip(mu+err)], c, FaceAlpha=0.05, EdgeColor=c, EdgeAlpha=0.5);
            end
            % xline(ax(iRow, iAx), 1e3*p.mdtaWindow, 'k--', Alpha=0.1)
            xline(ax(iType, iAx), 1e3*p.std.window, 'k--', Alpha=0.1)
            xline(ax(iType, iAx), 0, 'k-', Alpha=0.1)
            % xticks(ax(iRow, iAx), 1e3*p.mdtaWindow)
            xticks(ax(iType, iAx), [-300, 0, 600])
            xtickangle(ax(iType, iAx), 0)
    
            ylabel(ax(iType, iAx), featureUnits(iAx))
            
            fnDisp = sprintf("%s %s", fn, repmat('*', [1, nStarsSTD]));
            title(ax(iType, iAx), fnDisp, Interpreter='none')
            ylim(ax(iType, iAx), yl{iAx})
        end

        % Correlegram
        [lia, statFeatureOrder] = ismember(statFeatures, p.std.features);
        assert(all(lia), 'Some members of statFeatureOrder are not found.')

        iAx = iAx + 1;
        r = NaN(length(p.std.features));
        for i = 1:length(p.std.features)
            fni = p.std.features(statFeatureOrder(i));
            if isempty(dta(iUnit).(fni))
                continue
            end
            for j = 1:length(p.std.features)
                fnj = p.std.features(statFeatureOrder(j));
                if isempty(dta(iUnit).(fnj))
                    continue
                end
                if iType == 1
                    selT = dta(iUnit).(fni).t >= p.std.window(1) & dta(iUnit).(fni).t <= p.std.window(2);
                    r(i, j) = corr(std(dta(iUnit).(fni).X(:, selT), 0, 2, 'omitnan'), std(dta(iUnit).(fnj).X(:, selT), 0, 2, 'omitnan'), Rows='complete');
                else
                    selT = rta(iUnit).(fni).t >= p.std.window(1) & rta(iUnit).(fni).t <= p.std.window(2);
                    r(i, j) = corr(std(rta(iUnit).(fni).X(:, selT), 0, 2, 'omitnan'), std(rta(iUnit).(fnj).X(:, selT), 0, 2, 'omitnan'), Rows='complete');
                end
            end
        end
        r(isnan(r)) = 0;
        imagesc(ax(iType, iAx), r);
        xticks(ax(iType, iAx), 1:length(p.std.features))
        yticks(ax(iType, iAx), 1:length(p.std.features))
        xticklabels(ax(iType, iAx), p.std.features(statFeatureOrder))
        yticklabels(ax(iType, iAx), p.std.features(statFeatureOrder))
        xtickangle(ax(iType, iAx), 90)
        ax(iType, iAx).XAxisLocation = 'top';
        axis(ax(iType, iAx), 'image')
        ax(iType, iAx).XAxis.Direction = 'normal';
        clim(ax(iType, iAx), [0, 1])
        colormap(ax(iType, iAx), 'gray')
        colorbar(ax(iType, iAx), 'eastoutside')
        applyCustomColormap(ax(iType, iAx), [-1, 1], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);

        % Movement diversity matrix
        if p.nBoot > 0
            iAx = iAx + 1;
            if iType == 1
                xta = dta(iUnit);
                mdm = NaN(length(xta.tDip), length(p.std.features));
            else
                xta = rta(iUnit);
                mdm = NaN(length(xta.tRise), length(p.std.features));
            end
            for i = 1:length(p.std.features)
                fn = p.std.features(statFeatureOrder(i));
                if isempty(xta.(fn))
                    continue
                end
                t = xta.(fn).t;
                selT = t >= p.std.window(1) & t <= p.std.window(2);
                stdObs = std(xta.(fn).X(:, selT), 0, 2, 'omitnan');
                mdm(:, i) = arrayfun(@(data) nnz(xta.(fn).stats.stdBoot < data) ./ length(xta.(fn).stats.stdBoot), stdObs, UniformOutput=true);
            end
            mdm(isnan(mdm)) = 0;
            hash = sum((mdm > 0.95) .* 2.^(size(mdm, 2)-1:-1:0), 2);
            [~, I] = sort(hash, 'descend');
            imagesc(ax(iType, iAx), mdm(I, :))
            % colormap(ax(iType, iAx), [1, 1, 1; 0, 0, 0])
            % applyCustomColormap(ax(iType, iAx), [-1, 1], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);
            applyCustomColormap(ax(iType, iAx), [0, 1], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.025, h0=0.33);
            % ax(iType, iAx).ColorScale = 'log';
            xticks(ax(iType, iAx), 1:length(p.std.features))
            xticklabels(ax(iType, iAx), p.std.features(statFeatureOrder))
            colorbar(ax(iType, iAx), Orientation='horizontal', Location='southoutside')
            ax(iType, iAx).XAxisLocation = 'top';
            ylabel(ax(iType, iAx), 'Trial')
            clear iType xta mdm i fn t selT stdObs pObs hash I
        end        
    end
    if ~isempty(dta(iUnit).HandR)
        xlabel(tl(1), 'time from dip onset (ms)')
        if dta(iUnit).iExp > 0
            title(tl(1), sprintf("Unit %i (n=%i dips)", iUnit, length(dta(iUnit).tDip)), FontWeight='bold')
        else
            title(tl(1), sprintf("%i units (n=%i dips)", length(eu), length(dta(iUnit).tDip)), FontWeight='bold')
        end
    else
        xlabel(tl(1), '')
        title(tl(1), '')
    end
    if ~isempty(rta(iUnit).HandR)
        xlabel(tl(2), 'time from rise onset (ms)')
        if rta(iUnit).iExp > 0
            title(tl(2), sprintf("Unit %i (n=%i rises)", iUnit, length(rta(iUnit).tRise)), FontWeight='bold')
        else
            title(tl(2), sprintf("%i units (n=%i rises)", length(eu), length(rta(iUnit).tRise)), FontWeight='bold')
        end
    else
        xlabel(tl(2), '')
        title(tl(2), '')
    end
    % xlim(ax(:, 1:end-1), 1e3*p.std.window)
    xlim(ax(:, 1:end-2), 1e3*[-0.5, 1])
    fontsize(fig, 9, 'points')
    print(fig, fullfile(exportPath, sprintf("dta_rta_unit_%03i", iUnit)), '-dpng', '-r0')
end
clear features featureUnits fig tlp tl ax iAx iUnit fn c t mu err iType iFeature tlp yl fnDisp
clear prcSTD nStarsSTD i j fni fnj r selT
