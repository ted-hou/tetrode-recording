if ~exist('euArtiFree', 'var')
    if exist('E:\Data\Units\PressVsLick_ArtifactsRemoved_Full\FixedEventsAndTrials', 'dir')
        euArtiFree = EphysUnit.load('E:\Data\Units\PressVsLick_ArtifactsRemoved_Full\FixedEventsAndTrials');
        % load('E:\Data\Units\meta_PressVsLick_ArtifactsRemoved_Full_20260107.mat');
        % etaArtiFree = metaArtiFree.eta;
    else
        euArtiFree = EphysUnit.load('C:\SERVER\Units\PressVsLick_ArtifactsRemoved_Full\FixedEventsAndTrials');
        % load('C:\SERVER\Units\meta_PressVsLick_ArtifactsRemoved_Full_20260107.mat');
        % etaArtiFree = metaArtiFree.eta;
    end
end
% Daisy 2, 3, 8, 9, 10, 13, 14, 15, desmond10, 11, 22, 23, 24, 25, 26, 27
% Load DLC from 
% paths = [ ...
%     "\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\DeepLabCut\Results", ...
%     "\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\DeepLabCut\Projects\July14", ...
%     "\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\DeepLabCut\Projects\July3-DT-2025-07-03" ...
%     ];

paths = [ ...
    "\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\DeepLabCut\Results" ...
    % "\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\DeepLabCut\Projects\July14", ...
    % "\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\DeepLabCut\Projects\July3-DT-2025-07-03" ...
    ];

clear sessions
sessionNames = unique(string({euArtiFree.ExpName}));
nFound = 0;
for iSession = 1:length(sessionNames)
    sessions(iSession).name = sessionNames(iSession);
    for iPath = 1:length(paths)
        files = dir(sprintf("%s\\%s*.csv", paths(iPath), sessionNames(iSession)));
        if ~isempty(files)
            sessions(iSession).files = files;
            sessions(iSession).path = paths(iPath);
            sessions(iSession).eu = euArtiFree(ismember(string({euArtiFree.ExpName}), sessionNames(iSession)));
            break
        end
    end
end
clear iSession iPath files

% Make CompleteExperiment3 objects
clc
exp = CompleteExperiment3([sessions.eu], cameras='lr', deeplabcutPath='\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\DeepLabCut\Results');
exp.alignTimestamps(refEventNameArduino={'CUE_ON'}, refEventNameEphys={'Cue'}, trialDurationTolerance=2);

clear results
results(length(exp)) = struct(name=[], varsL=[], hasTimestampsL=[], varsR=[], hasTimestampsR=[]);
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
end
clear iExp

% Rename variables
from = ["HandCameraSide", "handIpsi", "handCont", "footIpsi", "footCont", "tongue"];
toL = ["HandL", "HandL", "HandR", "FootL", "FootR", "Tongue"];
toR = ["HandR", "HandR", "HandL", "FootR", "FootL", "Tongue"];
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

clear results
results(length(exp)) = struct(name=[], varsL=[], hasTimestampsL=[], varsR=[], hasTimestampsR=[]);
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
end
clear iExp toL toR vtd i from 

% Remove bad sessions
sel = [results.hasTimestampsL] | [results.hasTimestampsR];
exp = exp(sel);
results = results(sel);
eu = [exp.eu];

clear sel

% Reset
clearvars -except euArtiFree paths sessionNames sessions exp results

% Get the onset of arm reach/offset of arm retraction via video
p.useContinuousSequence = true;
% Get ReachStart
for iExp = 1:length(exp)
    if p.useContinuousSequence
        theseTrials = exp(iExp).eu(1).Trials.PressCorrect;
    else
        theseTrials = exp(iExp).eu(1).Trials.PressValid;
    end
    [X, Y, L, t] = exp(iExp).getTrajectoryByTrial('r', 'HandR', trialType='press', trials=theseTrials, alignTo='stop', window=[-4, 0], includeInvalid=true, likelihoodThreshold=0.5);
    selBaseline = t >= -4 & t <= -2;
    X = (X - mean(X(:, selBaseline), 2, 'omitnan')) ./ std(X(:, selBaseline), 0, 2, 'omitnan');
    Y = (Y - mean(Y(:, selBaseline), 2, 'omitnan')) ./ std(Y(:, selBaseline), 0, 2, 'omitnan');
    theta = 3;
    B = abs(X) >= theta | abs(Y) >= theta;
    tOnset = NaN(1, size(B, 1));
    for iTrial = 1:size(B, 1)
        iOnset = strfind(B(iTrial, :), [0, 1, 1]) + 1;
        if isempty(iOnset)
            tOnset(1, iTrial) = NaN;
        else
            iOnset = iOnset(end);
            tOnset(1, iTrial) = t(iOnset);
        end
    end
    fprintf('Found reach onset for %.0f%% (%i/%i) of trials, NaNs are replaced with the session median of %.0f ms.\n', 100*nnz(~isnan(tOnset))./length(tOnset), nnz(~isnan(tOnset)), length(tOnset), 1e3*median(tOnset, 'all', 'omitnan'));
    tOnset(isnan(tOnset)) = median(tOnset, 'all', 'omitnan');
    cueTrials = Trial([theseTrials.Start], [theseTrials.Stop] + tOnset, advancedValidation=false);
    theseTrials = Trial([theseTrials.Stop] + tOnset, [theseTrials.Stop], advancedValidation=false);
    for iEu = 1:length(exp(iExp).eu)
        if p.useContinuousSequence
            exp(iExp).eu(iEu).Trials.CorrectReachStartToReachEnd = theseTrials;
            exp(iExp).eu(iEu).Trials.CorrectCueToReachStart = cueTrials;
        else
            exp(iExp).eu(iEu).Trials.ValidReachStartToReachEnd = theseTrials;
            exp(iExp).eu(iEu).Trials.ValidCueToReachStart = cueTrials;
        end
    end
end
clear iExp X Y L t selBaseline theta B tOnset iTrial iOnset theseTrials

% Get RetractEnd
for iExp = 1:length(exp)
    if p.useContinuousSequence
        theseTrials = exp(iExp).eu(1).Trials.CueToLeverReleaseCorrect;
    else
        theseTrials = exp(iExp).eu(1).Trials.CueToLeverRelease;
    end
    [X, Y, L, t] = exp(iExp).getTrajectoryByTrial('r', 'HandR', trialType='press', trials=theseTrials, alignTo='stop', window=[0, 10], includeInvalid=true, likelihoodThreshold=0.5);
    selBaseline = t >= 4 & t <= 10;
    X = (X - mean(X(:, selBaseline), 2, 'omitnan')) ./ std(X(:, selBaseline), 0, 2, 'omitnan');
    Y = (Y - mean(Y(:, selBaseline), 2, 'omitnan')) ./ std(Y(:, selBaseline), 0, 2, 'omitnan');
    theta = 2;
    B = abs(X) <= theta & abs(Y) <= theta;
    tOffset = NaN(1, size(B, 1));
    for iTrial = 1:size(B, 1)
        iOffset = strfind(B(iTrial, :), [0, 1, 1]) + 1;
        if isempty(iOffset)
            tOffset(1, iTrial) = NaN;
        else
            iOffset = iOffset(1);
            tOffset(1, iTrial) = t(iOffset);
        end
    end
    fprintf('Found retract offset for %.0f%% (%i/%i) of trials, NaNs are replaced with the session median of %.0f ms.\n', 100*nnz(~isnan(tOffset))./length(tOffset), nnz(~isnan(tOffset)), length(tOffset), 1e3*median(tOffset, 'all', 'omitnan'));
    tOffset(isnan(tOffset)) = median(tOffset, 'all', 'omitnan');
    cueTrials = Trial([theseTrials.Start], [theseTrials.Stop] + tOffset, advancedValidation=false);
    theseTrials = Trial([theseTrials.Stop], [theseTrials.Stop] + tOffset, advancedValidation=false);
    for iEu = 1:length(exp(iExp).eu)
        if p.useContinuousSequence
            exp(iExp).eu(iEu).Trials.CorrectRetractStartToRetractEnd = theseTrials;
            exp(iExp).eu(iEu).Trials.CorrectCueToRetractEnd = cueTrials;
        else
            exp(iExp).eu(iEu).Trials.ValidRetractStartToRetractEnd = theseTrials;
            exp(iExp).eu(iEu).Trials.ValidCueToRetractEnd = cueTrials;
        end
    end
end
clear iExp X Y L t selBaseline theta B tOffset iTrial iOffset theseTrials cueTrials

eu = [exp.eu];
expIndices = cellfun(@(name) find(strcmpi(name, {exp.name}), 1, 'first'), {eu.ExpName});


%% Detect dips in firing
p.features = ["HandR", "HandL", "FootR", "FootL", "Tongue", "Lick"];
p.vtdNames = ["vtdR", "vtdL", "vtdR", "vtdL", "vtdR", "vtdR"];
p.smoothWindow = [10, 10, 10, 10, 5, 5];
p.minL = [0.5, 0.5, 0.5, 0.5, 0.2, NaN];
p.dtaRes = 1/30;

p.dipSamples = [2, 8];
p.dipThresholdQuantile = 0.05;
p.dipThresholdSubQuantile = 0.05;
p.dipPattern = arrayfun(@(n) [0, ones(1, n), 0] , p.dipSamples(1):p.dipSamples(2), UniformOutput=false); % 100-300ms dips
p.dipPatternOnset = cellfun(@(pat) find(pat, 1, 'first') - 1, p.dipPattern); % finds the onset
p.dtaWindow = [-1, 1];
p.mdtaWindow = [-0.3, 0.3];

p.riseSamples = [2, 8];
p.riseThresholdQuantile = 1 - p.dipThresholdQuantile;
p.riseThresholdSubQuantile = 1 - p.dipThresholdSubQuantile;
p.risePattern = arrayfun(@(n) [0, ones(1, n), 0] , p.riseSamples(1):p.riseSamples(2), UniformOutput=false);
p.risePatternOnset = cellfun(@(pat) find(pat, 1, 'first') - 1, p.risePattern);

p.nBoot = 1000;
p.bootAlpha = 0.05;

% Process vtdFeatues
clear kinematics
kinematics(length(exp)) = struct(HandR=[], HandL=[], FootR=[], FootL=[], Tongue=[]);
for iExp = 1:length(exp)
    for iFeature = 1:length(p.features)
        vn = p.vtdNames(iFeature);
        fn = p.features(iFeature);
        vtd = exp(iExp).(vn);
        if isempty(vtd)
            continue
        end
        % Tongue uses bilateral likelihood
        if fn == "Tongue"
            if isempty(exp(iExp).vtdL)
                L = exp(iExp).vtdR.(sprintf("%s_Likelihood", fn));
                L(L<p.minL(iFeature)) = 0;
                L(L>p.minL(iFeature)) = 1;
                L = single(L > 0);
                t = vtd.Timestamp;
            else
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
            end
            L = smoothdata(L, 'gaussian', p.smoothWindow(iFeature));
            kinematics(iExp).(fn) = struct(X=L, t=t);
            clear t
        % Licks from digital events
        elseif fn == "Lick"
            L = exp(iExp).eu(1).EventTimes.LickOn;
            t = 0:p.dtaRes:exp(iExp).vtdR.Timestamp(end);
            edges = [t - p.dtaRes/2, t(end) + p.dtaRes/2];
            L = histcounts(L, edges);
            L = smoothdata(L, 'gaussian', p.smoothWindow(iFeature));
            kinematics(iExp).(fn) = struct(X=L, t=t);
            clear edges t
        % Other tracking points use position or speed
        elseif ismember(sprintf("%s_X", fn), vtd.Properties.VariableNames)
            X = vtd.(sprintf("%s_X", fn));
            Y = vtd.(sprintf("%s_Y", fn));
            L = vtd.(sprintf("%s_Likelihood", fn));
            t = vtd.Timestamp;
            X(L<p.minL(iFeature)) = NaN;
            Y(L<p.minL(iFeature)) = NaN;
            X = (X - mean(X, 'all', 'omitnan')) ./ std(X, 0, 'all', 'omitnan');
            Y = (Y - mean(Y, 'all', 'omitnan')) ./ std(Y, 0, 'all', 'omitnan');
            % X = [NaN; diff(X)]./[NaN; diff(t)];
            % Y = [NaN; diff(Y)]./[NaN; diff(t)];
            D = sqrt(X.^2 + Y.^2);
            % selnan = isnan(D);
            % D(selnan) = interp1(vtd.Timestamp(~selnan), D(~selnan), vtd.Timestamp(selnan), 'linear');
            D = smoothdata(D, 'gaussian', p.smoothWindow(iFeature));
            kinematics(iExp).(fn) = struct(X=D, t=vtd.Timestamp);
        end
    end
end
clear iExp iFeature vn fn vtd L X Y D selnan

%%
tLocal = p.dtaWindow(1):p.dtaRes:p.dtaWindow(2);
clear dta rta
lineLength = 0;
lineLength2 = 0;
tTicTotal = tic();
rng(42)
for iEu = 1:length(eu)
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

    % Find rises
    riseThreshold = quantile(x(x>0), p.riseThresholdQuantile);
    iRise = arrayfun(@(i) strfind(x>=riseThreshold, p.risePattern{i}) + p.risePatternOnset(i), 1:length(p.risePattern), UniformOutput=false);
    iRise = cat(2, iRise{:});
    tRise = t(iRise);

    fprintf(repmat('\b', [1, lineLength]))
    lineLength = fprintf('Unit %i/%i;; %i dips (x<%.2f), %i rises (x>%.2f)... %.1fs elapsed...', iEu, length(eu), length(tDip), dipThreshold, length(tRise), riseThreshold, toc(tTicTotal));

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
clear lineLength2 I2 iDip
exportPath = fullfile("C:\SERVER\LickVsReach_DTA_RTA_boot", sprintf("LickVsReach_DLC_dta_rta_%i_%i_%ito%ims.mat", 100*p.dipThresholdQuantile, 100*p.dipThresholdSubQuantile, 100*p.dipSamples(1), 100*p.riseSamples(2)));
save(exportPath, 'dta', 'rta', 'kinematics', 'p', '-v7.3')


%% Post-hoc do a bootstrap for dta/rta spikerate
lineLength = 0;
tTicTotal = tic();
tLocal = p.dtaWindow(1):p.dtaRes:p.dtaWindow(2);
if p.nBoot > 0
    for iEu = 1:length(dta)
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

%% Post-hoc calculate the bootstrpped 95%CI of the RMS of the displacement traces
p.rms.features = ["HandR", "HandL", "FootR", "FootL", "Tongue", "Lick"];
p.rms.window = [-0.3, 0.6];

for iUnit = 1:length(dta)
    for fn = p.rms.features
        if isempty(dta(iUnit).(fn))
            continue
        end
        selT = dta(iUnit).(fn).t >= p.rms.window(1) & dta(iUnit).(fn).t <= p.rms.window(2);
        dta(iUnit).(fn).stats.rms = rms(mean(dta(iUnit).(fn).X(:, selT), 1, 'omitnan'), 2, 'omitnan');
        dta(iUnit).(fn).stats.rmsBoot = rms(dta(iUnit).(fn).XBoot(:, selT), 2, 'omitnan');
    end
    for fn = p.rms.features
        if isempty(rta(iUnit).(fn))
            continue
        end
        selT = rta(iUnit).(fn).t >= p.rms.window(1) & rta(iUnit).(fn).t <= p.rms.window(2);
        rta(iUnit).(fn).stats.rms = rms(mean(rta(iUnit).(fn).X(:, selT), 1, 'omitnan'), 2, 'omitnan');
        rta(iUnit).(fn).stats.rmsBoot = rms(rta(iUnit).(fn).XBoot(:, selT), 2, 'omitnan');
    end
end

clear iUnit fn selT

%% Post-hoc calculate the bootstrpped 95%CI of the STD of the displacement traces
p = rmfield(p, 'rms');

p.std.features = ["HandR", "HandL", "FootR", "FootL", "Tongue", "Lick"];
p.std.window = [-0.3, 0.6];

for iUnit = 1:length(dta)
    for fn = p.std.features
        if isempty(dta(iUnit).(fn))
            continue
        end
        selT = dta(iUnit).(fn).t >= p.std.window(1) & dta(iUnit).(fn).t <= p.std.window(2);
        dta(iUnit).(fn).stats.std = std(mean(dta(iUnit).(fn).X(:, selT), 1, 'omitnan'), 0, 2, 'omitnan');
        dta(iUnit).(fn).stats.stdBoot = std(dta(iUnit).(fn).XBoot(:, selT), 0, 2, 'omitnan');
    end
    for fn = p.std.features
        if isempty(rta(iUnit).(fn))
            continue
        end
        selT = rta(iUnit).(fn).t >= p.std.window(1) & rta(iUnit).(fn).t <= p.std.window(2);
        rta(iUnit).(fn).stats.std = std(mean(rta(iUnit).(fn).X(:, selT), 1, 'omitnan'), 0, 2, 'omitnan');
        rta(iUnit).(fn).stats.stdBoot = std(rta(iUnit).(fn).XBoot(:, selT), 0, 2, 'omitnan');
    end
end

clear iUnit fn selT


%% Plot dip-triggered average kinematics (RMS VERSION)
% close all
% exportPath = fullfile("C:\SERVER\LickVsReach_DTA_RTA_boot\Figures", sprintf("LickVsReach_DLC_dta_rta_%i_%i_%ito%ims", 100*p.dipThresholdQuantile, 100*p.dipThresholdSubQuantile, 100*p.dipSamples(1), 100*p.riseSamples(2)));
% if ~exist(exportPath, 'dir')
%     mkdir(exportPath)
% end
% features = ["spikerate", "HandR", "HandL", "FootR", "FootL", "Tongue", "Lick"];
% featureUnits = ["a.u.", "a.u.", "a.u.", "a.u.", "a.u.", "prob", "prob"];
% yl = {[-2.5, 5], [-0.1, 5.1], [-0.1, 5.1], [-0.1, 5.1], [-0.1, 5.1], [-0.1, 1.1], [-0.1, 1.1]};
% fig = figure(Units='inches', InnerPosition=[2, 2, 2*(2+length(features)), 5]);
% tlp = tiledlayout(fig, 2, 1, TileSpacing='compact', Padding='compact');
% tl = gobjects(2, 1);
% tl(1) = tiledlayout(tlp, 1, length(features) + 2, TileSpacing='compact', Padding='compact');
% tl(2) = tiledlayout(tlp, 1, length(features) + 2, TileSpacing='compact', Padding='compact');
% tl(1).Layout.Tile = 1;
% tl(2).Layout.Tile = 2;
% 
% ax = gobjects(2, length(features) + 2);
% for iType = 1:2
%     for iAx = 1:length(features) + 2
%         ax(iType, iAx) = nexttile(tl(iType));
%     end
% end
% for iUnit = 1:length(dta)
%     for iType = 1:2
%         % Check existence
%         if iType == 1
%             if isempty(dta(iUnit).tDip)
%                 for iAx = 1:length(features)
%                     cla(ax(iType, iAx))
%                     ax(iType, iAx).Visible = false;
%                 end
%                 continue
%             end
%         else
%             if isempty(rta(iUnit).tRise)
%                 for iAx = 1:length(features)
%                     cla(ax(iType, iAx))
%                     ax(iType, iAx).Visible = false;
%                 end
%                 continue
%             end
%         end
%         for iAx = 1:length(features)
%             cla(ax(iType, iAx))
%             fn = features(iAx);
%             if iType == 1 && isempty(dta(iUnit).(fn))
%                 ax(iType, iAx).Visible = false;
%                 continue
%             elseif iType == 2 && isempty(rta(iUnit).(fn))
%                 ax(iType, iAx).Visible = false;
%                 continue
%             end
% 
%             ax(iType, iAx).Visible = true;
% 
%             hold(ax(iType, iAx), 'on')
%             c = getColor(iAx, length(features), 0.7);
%             if iType == 1
%                 t = 1e3*dta(iUnit).(fn).t;
%                 X = dta(iUnit).(fn).X;
%                 mu = mean(dta(iUnit).(fn).X, 1, 'omitnan');
%                 err = std(dta(iUnit).(fn).X, 0, 1, 'omitnan')./sqrt(size(dta(iUnit).(fn).X, 1));
%                 prc = quantile(dta(iUnit).(fn).XBoot, [p.bootAlpha/2, 1-p.bootAlpha/2], 1);
% 
%                 if ismember(fn, p.rms.features)
%                     prcRMS = quantile(dta(iUnit).(fn).stats.rmsBoot, [0.95, 0.99, 0.999]);
%                     nStarsRMS = sum(dta(iUnit).(fn).stats.rms > prcRMS);
%                 else
%                     nStarsRMS = 0;
%                 end
%             else
%                 t = 1e3*rta(iUnit).(fn).t;
%                 X = rta(iUnit).(fn).X;
%                 mu = mean(rta(iUnit).(fn).X, 1, 'omitnan');
%                 err = std(rta(iUnit).(fn).X, 0, 1, 'omitnan')./sqrt(size(rta(iUnit).(fn).X, 1));
%                 prc = quantile(rta(iUnit).(fn).XBoot, [p.bootAlpha/2, 1-p.bootAlpha/2], 1);
% 
%                 if ismember(fn, p.rms.features)
%                     prcRMS = quantile(rta(iUnit).(fn).stats.rmsBoot, [0.95, 0.99, 0.999]);
%                     nStarsRMS = sum(rta(iUnit).(fn).stats.rms > prcRMS);
%                 else
%                     nStarsRMS = 0;
%                 end
%             end
%             plot(ax(iType, iAx), t, X, Color=[c, 0.1]);
%             plot(ax(iType, iAx), t, mu, Color=c, LineWidth=1.5);
%             % patch(ax(iRow, iAx), [t, flip(t)], [mu-err, flip(mu+err)], c, FaceAlpha=0.15, EdgeColor=c);
%             patch(ax(iType, iAx), [t, flip(t)], [prc(1, :), flip(prc(2, :))], c, FaceAlpha=0.05, EdgeColor=c, EdgeAlpha=0.5);
%             % xline(ax(iRow, iAx), 1e3*p.mdtaWindow, 'k--', Alpha=0.1)
%             xline(ax(iType, iAx), 1e3*p.rms.window, 'k--', Alpha=0.1)
%             xline(ax(iType, iAx), 0, 'k-', Alpha=0.1)
%             % xticks(ax(iRow, iAx), 1e3*p.mdtaWindow)
%             xticks(ax(iType, iAx), [-300, 0, 600])
%             xtickangle(ax(iType, iAx), 0)
% 
%             ylabel(ax(iType, iAx), featureUnits(iAx))
% 
%             fnDisp = sprintf("%s %s", fn, repmat('*', [1, nStarsRMS]));
%             title(ax(iType, iAx), fnDisp, Interpreter='none')
%             ylim(ax(iType, iAx), yl{iAx})
%         end
% 
%         % Correlegram
%         iAx = iAx + 1;
%         r = NaN(length(p.rms.features));
%         for i = 1:length(p.rms.features)
%             fni = p.rms.features(i);
%             if isempty(dta(iUnit).(fni))
%                 continue
%             end
%             for j = 1:length(p.rms.features)
%                 fnj = p.rms.features(j);
%                 if isempty(dta(iUnit).(fnj))
%                     continue
%                 end
%                 if iType == 1
%                     selT = dta(iUnit).(fni).t >= p.rms.window(1) & dta(iUnit).(fni).t <= p.rms.window(2);
%                     r(i, j) = corr(rms(dta(iUnit).(fni).X(:, selT), 2, 'omitnan'), rms(dta(iUnit).(fnj).X(:, selT), 2, 'omitnan'), Rows='complete');
%                 else
%                     selT = rta(iUnit).(fni).t >= p.rms.window(1) & rta(iUnit).(fni).t <= p.rms.window(2);
%                     r(i, j) = corr(rms(rta(iUnit).(fni).X(:, selT), 2, 'omitnan'), rms(rta(iUnit).(fnj).X(:, selT), 2, 'omitnan'), Rows='complete');
%                 end
%             end
%         end
%         r(isnan(r)) = 0;
%         imagesc(ax(iType, iAx), r);
%         xticks(ax(iType, iAx), 1:length(p.rms.features))
%         yticks(ax(iType, iAx), 1:length(p.rms.features))
%         xticklabels(ax(iType, iAx), p.rms.features)
%         yticklabels(ax(iType, iAx), p.rms.features)
%         xtickangle(ax(iType, iAx), 90)
%         ax(iType, iAx).XAxisLocation = 'top';
%         axis(ax(iType, iAx), 'image')
%         ax(iType, iAx).XAxis.Direction = 'normal';
%         clim(ax(iType, iAx), [0, 1])
%         colormap(ax(iType, iAx), 'gray')
%         colorbar(ax(iType, iAx), 'eastoutside')
%         applyCustomColormap(ax(iType, iAx), [-1, 1], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);
% 
%         % Movement diversity matrix
%         iAx = iAx + 1;
%         if iType == 1
%             xta = dta(iUnit);
%             mdm = NaN(length(xta.tDip), length(p.rms.features));
%         else
%             xta = rta(iUnit);
%             mdm = NaN(length(xta.tRise), length(p.rms.features));
%         end
%         for i = 1:length(p.rms.features)
%             fn = p.rms.features(i);
%             if isempty(xta.(fn))
%                 continue
%             end
%             t = xta.(fn).t;
%             selT = t >= p.rms.window(1) & t <= p.rms.window(2);
%             rmsObs = rms(xta.(fn).X(:, selT), 2, 'omitnan');
%             mdm(:, i) = arrayfun(@(data) nnz(xta.(fn).stats.rmsBoot < data) ./ length(xta.(fn).stats.rmsBoot), rmsObs, UniformOutput=true);
%         end
%         mdm(isnan(mdm)) = 0;
%         hash = sum((mdm > 0.95) .* 2.^(size(mdm, 2)-1:-1:0), 2);
%         [~, I] = sort(hash, 'descend');
%         imagesc(ax(iType, iAx), mdm(I, :))
%         % colormap(ax(iType, iAx), [1, 1, 1; 0, 0, 0])
%         % applyCustomColormap(ax(iType, iAx), [-1, 1], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);
%         applyCustomColormap(ax(iType, iAx), [0, 1], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.025, h0=0.33);
%         % ax(iType, iAx).ColorScale = 'log';
%         xticks(ax(iType, iAx), 1:length(p.rms.features))
%         xticklabels(ax(iType, iAx), p.rms.features)
%         colorbar(ax(iType, iAx), Orientation='horizontal', Location='southoutside')
%         ax(iType, iAx).XAxisLocation = 'top';
%         ylabel(ax(iType, iAx), 'Trial')
%         clear iType xta mdm i fn t selT rmsObs pObs hash I
% 
%     end
%     if ~isempty(dta(iUnit).HandR)
%         xlabel(tl(1), 'time from dip onset (ms)')
%         title(tl(1), sprintf("Unit %i (n=%i dips)", iUnit, length(dta(iUnit).tDip)), FontWeight='bold')
%     else
%         xlabel(tl(1), '')
%         title(tl(1), '')
%     end
%     if ~isempty(rta(iUnit).HandR)
%         xlabel(tl(2), 'time from rise onset (ms)')
%         title(tl(2), sprintf("Unit %i (n=%i rises)", iUnit, length(rta(iUnit).tRise)), FontWeight='bold')
%     else
%         xlabel(tl(2), '')
%         title(tl(2), '')
%     end
%     % xlim(ax(:, 1:end-1), 1e3*p.rms.window)
%     xlim(ax(:, 1:end-2), 1e3*[-0.5, 1])
%     fontsize(fig, 9, 'points')
%     print(fig, fullfile(exportPath, sprintf("dta_rta_unit_%03i", iUnit)), '-dpng', '-r0')
% end
% clear features featureUnits fig tlp tl ax iAx iUnit fn c t mu err iType iFeature tlp yl fnDisp
% clear prcRMS nStarsRMS i j fni fnj r selT

%% Plot dip-triggered average kinematics (STD Version
close all
exportPath = fullfile("C:\SERVER\LickVsReach_DTA_RTA_boot\Figures", sprintf("LickVsReach_DLC_dta_rta_%i_%i_%ito%ims_std", 100*p.dipThresholdQuantile, 100*p.dipThresholdSubQuantile, 100*p.dipSamples(1), 100*p.riseSamples(2)));
if ~exist(exportPath, 'dir')
    mkdir(exportPath)
end
features = ["spikerate", "HandR", "HandL", "FootR", "FootL", "Tongue", "Lick"];
featureUnits = ["a.u.", "a.u.", "a.u.", "a.u.", "a.u.", "prob", "prob"];
yl = {[-2.5, 5], [-0.1, 5.1], [-0.1, 5.1], [-0.1, 5.1], [-0.1, 5.1], [-0.1, 1.1], [-0.1, 1.1]};
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
for iUnit = 1:length(dta)
    for iType = 1:2
        % Check existence
        if iType == 1
            if isempty(dta(iUnit).tDip)
                for iAx = 1:length(features)
                    cla(ax(iType, iAx))
                    ax(iType, iAx).Visible = false;
                end
                continue
            end
        else
            if isempty(rta(iUnit).tRise)
                for iAx = 1:length(features)
                    cla(ax(iType, iAx))
                    ax(iType, iAx).Visible = false;
                end
                continue
            end
        end
        for iAx = 1:length(features)
            cla(ax(iType, iAx))
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
                prc = quantile(dta(iUnit).(fn).XBoot, [p.bootAlpha/2, 1-p.bootAlpha/2], 1);

                if ismember(fn, p.std.features)
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
                prc = quantile(rta(iUnit).(fn).XBoot, [p.bootAlpha/2, 1-p.bootAlpha/2], 1);

                if ismember(fn, p.std.features)
                    prcSTD = quantile(rta(iUnit).(fn).stats.stdBoot, [0.95, 0.99, 0.999]);
                    nStarsSTD = sum(rta(iUnit).(fn).stats.std > prcSTD);
                else
                    nStarsSTD = 0;
                end
            end
            plot(ax(iType, iAx), t, X, Color=[0.15, 0.15, 0.15, 0.1]);
            plot(ax(iType, iAx), t, mu, Color=c, LineWidth=1.5);
            % patch(ax(iRow, iAx), [t, flip(t)], [mu-err, flip(mu+err)], c, FaceAlpha=0.15, EdgeColor=c);
            patch(ax(iType, iAx), [t, flip(t)], [prc(1, :), flip(prc(2, :))], c, FaceAlpha=0.05, EdgeColor=c, EdgeAlpha=0.5);
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
        iAx = iAx + 1;
        r = NaN(length(p.std.features));
        for i = 1:length(p.std.features)
            fni = p.std.features(i);
            if isempty(dta(iUnit).(fni))
                continue
            end
            for j = 1:length(p.std.features)
                fnj = p.std.features(j);
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
        xticklabels(ax(iType, iAx), p.std.features)
        yticklabels(ax(iType, iAx), p.std.features)
        xtickangle(ax(iType, iAx), 90)
        ax(iType, iAx).XAxisLocation = 'top';
        axis(ax(iType, iAx), 'image')
        ax(iType, iAx).XAxis.Direction = 'normal';
        clim(ax(iType, iAx), [0, 1])
        colormap(ax(iType, iAx), 'gray')
        colorbar(ax(iType, iAx), 'eastoutside')
        applyCustomColormap(ax(iType, iAx), [-1, 1], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);

        % Movement diversity matrix
        iAx = iAx + 1;
        if iType == 1
            xta = dta(iUnit);
            mdm = NaN(length(xta.tDip), length(p.std.features));
        else
            xta = rta(iUnit);
            mdm = NaN(length(xta.tRise), length(p.std.features));
        end
        for i = 1:length(p.std.features)
            fn = p.std.features(i);
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
        xticklabels(ax(iType, iAx), p.std.features)
        colorbar(ax(iType, iAx), Orientation='horizontal', Location='southoutside')
        ax(iType, iAx).XAxisLocation = 'top';
        ylabel(ax(iType, iAx), 'Trial')
        clear iType xta mdm i fn t selT stdObs pObs hash I
        
    end
    if ~isempty(dta(iUnit).HandR)
        xlabel(tl(1), 'time from dip onset (ms)')
        title(tl(1), sprintf("Unit %i (n=%i dips)", iUnit, length(dta(iUnit).tDip)), FontWeight='bold')
    else
        xlabel(tl(1), '')
        title(tl(1), '')
    end
    if ~isempty(rta(iUnit).HandR)
        xlabel(tl(2), 'time from rise onset (ms)')
        title(tl(2), sprintf("Unit %i (n=%i rises)", iUnit, length(rta(iUnit).tRise)), FontWeight='bold')
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

%% Plot scatter of body movements
close all
exportPath = fullfile("E:\DATA\LickVsReach_DTA_RTA_boot\Figures", sprintf("LickVsReach_DLC_dta_rta_%i_%i_%ito%ims", 100*p.dipThresholdQuantile, 100*p.dipThresholdSubQuantile, 100*p.dipSamples(1), 100*p.riseSamples(2)));
if ~exist(exportPath, 'dir')
    mkdir(exportPath)
end
features = ["HandR", "HandL", "FootR", "FootL", "Tongue", "Lick"];
featureUnits = ["a.u.", "a.u.", "a.u.", "a.u.", "prob", "prob"];
yl = {[0, 2], [0, 2], [0, 2], [0, 2], [0, 1], [0, 1]};
fig = figure(Units='inches', InnerPosition=[0.1, 0.1, 2*length(features), 1*length(features)]);
tlp = tiledlayout(fig, 1, 2, TileSpacing='compact', Padding='compact');
tl = gobjects(2, 1);
tl(1) = tiledlayout(tlp, length(features), length(features), TileSpacing='tight', Padding='tight');
tl(2) = tiledlayout(tlp, length(features), length(features), TileSpacing='tight', Padding='tight');
tl(1).Layout.Tile = 1;
tl(2).Layout.Tile = 2;

ax = gobjects(2, length(features), length(features));
for iType = 1:2
    for i = 1:length(features)
        for j = 1:length(features)
            ax(iType, i, j) = nexttile(tl(iType));
            ax(iType, i, j).XAxis.Direction = 'reverse';
            % ax(iType, i, j).XAxisLocation = 'top';
            ax(iType, i, j).YAxisLocation = 'right';
            axis(ax(iType, i, j), 'square')
        end
    end
end
set(ax, Visible=false)
set(ax, Box='on')

for iUnit = 1:length(dta)
    set(ax, Visible=false)
    cla(ax)
    for iType = 1:2
        switch iType
            case 1
                xta = dta(iUnit);
            case 2
                xta = rta(iUnit);
        end
        for j = 1:length(features)
            fnj = features(j);
            if isempty(xta.(fnj))
                continue
            end
            for i = 1:j
                fni = features(i);
                if isempty(xta.(fni))
                    continue
                end
                axTemp = ax(iType, i, j);
                set(axTemp, Visible=true)
                hold(axTemp, 'on')
                t = xta.(fni).t;
                selT = t >= p.rms.window(1) & t <= p.rms.window(2);
                rmsI = rms(xta.(fni).X(:, selT), 2, 'omitnan');
                rmsJ = rms(xta.(fnj).X(:, selT), 2, 'omitnan');
                prcI = arrayfun(@(data) nnz(xta.(fni).stats.rmsBoot <= data) ./ length(xta.(fni).stats.rmsBoot), rmsI, UniformOutput=true);
                prcJ = arrayfun(@(data) nnz(xta.(fnj).stats.rmsBoot <= data) ./ length(xta.(fnj).stats.rmsBoot), rmsJ, UniformOutput=true);
                if i == j
                    col = [0.75, 0.15, 0.15, 1];
                else
                    col = [0.15, 0.15, 0.15, 1];
                end
                plot(axTemp, prcI, prcJ, 'o', Color=col, MarkerSize=2)
                xlim(axTemp, [-0.05, 1.05])
                ylim(axTemp, [-0.05, 1.05])
                xlabel(axTemp, features(i))
                ylabel(axTemp, features(j))
                if i == j
                    col = [0.75, 0.15, 0.15, 0.33];
                else
                    col = [0.15, 0.15, 0.15, 0.33];
                end
                plot(axTemp, [0, 1], [0, 1], '--', Color=col)
                hold(axTemp, 'off')
                xticks(axTemp, [0, 1])
                yticks(axTemp, [0, 1])
            end
        end
        if iType == 1
            title(tl(iType), sprintf('Dips (n=%i)', length(xta.tDip)))
        else
            title(tl(iType), sprintf('Rises (n=%i)', length(xta.tRise)))
        end
    end

    fontsize(fig, 9, 'points')
    print(fig, fullfile(exportPath, sprintf("dta_rta_unit_%03i_movement_scatter", iUnit)), '-dpng', '-r0')
end
clear features featureUnits fig tlp tl ax iType i j iUnit fni fnj t selT rmsI rmsJ prcI prcJ
