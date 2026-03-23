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
save(fullfile("E:\Data", sprintf("LickVsReach_DLC_dta_rta_%i_%i_%ito%ims.mat", 100*p.dipThresholdQuantile, 100*p.dipThresholdSubQuantile, 100*p.dipSamples, 100*p.riseSamples)), 'dta', 'rta', '-v7.3')


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

%% Plot dip-triggered average kinematics
close all
exportPath = fullfile("C:\SERVER\LickVsReach_DTA_RTA_boot", sprintf("LickVsReach_DLC_dta_rta_%i_%i_%ito%ims.mat", 100*p.dipThresholdQuantile, 100*p.dipThresholdSubQuantile, 100*p.dipSamples(1), 100*p.riseSamples(2)));
if ~exist(exportPath, 'dir')
    mkdir(exportPath)
end
features = ["spikerate", "HandR", "HandL", "FootR", "FootL", "Tongue", "Lick"];
featureUnits = ["a.u.", "a.u.", "a.u.", "a.u.", "a.u.", "prob", "prob"];
yl = {[-2.5, 5], [-0.1, 5.1], [-0.1, 5.1], [-0.1, 5.1], [-0.1, 5.1], [-0.1, 1.1], [-0.1, 1.1]};
fig = figure(Units='inches', InnerPosition=[2, 2, 2*length(features), 5]);
tlp = tiledlayout(fig, 2, 1, TileSpacing='compact', Padding='compact');
tl = gobjects(2, 1);
tl(1) = tiledlayout(tlp, 1, length(features), TileSpacing='compact', Padding='compact');
tl(2) = tiledlayout(tlp, 1, length(features), TileSpacing='compact', Padding='compact');
tl(1).Layout.Tile = 1;
tl(2).Layout.Tile = 2;

ax = gobjects(2, length(features));
for iRow = 1:2
    for iAx = 1:length(features)
        ax(iRow, iAx) = nexttile(tl(iRow));
    end
end
for iUnit = 1:length(dta)
    for iRow = 1:2
        % Check existence
        if iRow == 1
            if isempty(dta(iUnit).tDip)
                for iAx = 1:length(features)
                    cla(ax(iRow, iAx))
                    ax(iRow, iAx).Visible = false;
                end
                continue
            end
        else
            if isempty(rta(iUnit).tRise)
                for iAx = 1:length(features)
                    cla(ax(iRow, iAx))
                    ax(iRow, iAx).Visible = false;
                end
                continue
            end
        end
        for iAx = 1:length(features)
            cla(ax(iRow, iAx))
            fn = features(iAx);
            if iRow == 1 && isempty(dta(iUnit).(fn))
                ax(iRow, iAx).Visible = false;
                continue
            elseif iRow == 2 && isempty(rta(iUnit).(fn))
                ax(iRow, iAx).Visible = false;
                continue
            end

            ax(iRow, iAx).Visible = true;

            hold(ax(iRow, iAx), 'on')
            c = getColor(iAx, length(features), 0.7);
            if iRow == 1
                t = 1e3*dta(iUnit).(fn).t;
                mu = mean(dta(iUnit).(fn).X, 1, 'omitnan');
                err = std(dta(iUnit).(fn).X, 0, 1, 'omitnan')./sqrt(size(dta(iUnit).(fn).X, 1));
                prc = quantile(dta(iUnit).(fn).XBoot, [p.bootAlpha/2, 1-p.bootAlpha/2], 1);
            else
                t = 1e3*rta(iUnit).(fn).t;
                mu = mean(rta(iUnit).(fn).X, 1, 'omitnan');
                err = std(rta(iUnit).(fn).X, 0, 1, 'omitnan')./sqrt(size(rta(iUnit).(fn).X, 1));
                prc = quantile(rta(iUnit).(fn).XBoot, [p.bootAlpha/2, 1-p.bootAlpha/2], 1);
            end
            plot(ax(iRow, iAx), t, mu, Color=c, LineWidth=1.5);
            % patch(ax(iRow, iAx), [t, flip(t)], [mu-err, flip(mu+err)], c, FaceAlpha=0.15, EdgeColor=c);
            patch(ax(iRow, iAx), [t, flip(t)], [prc(1, :), flip(prc(2, :))], c, FaceAlpha=0.05, EdgeColor=c, EdgeAlpha=0.5);
            xline(ax(iRow, iAx), 1e3*p.mdtaWindow, 'k--', Alpha=0.1)
            xline(ax(iRow, iAx), 0, 'k-', Alpha=0.1)
            xticks(ax(iRow, iAx), 1e3*p.mdtaWindow)
    
            ylabel(ax(iRow, iAx), featureUnits(iAx))
            title(ax(iRow, iAx), fn, Interpreter='none')
            ylim(ax(iRow, iAx), yl{iAx})
        end
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
    fontsize(fig, 9, 'points')
    print(fig, fullfile(exportPath, sprintf("dta_rta_unit_%03i", iUnit)), '-dpng', '-r0')
end
clear features featureUnits fig tl ax iAx iUnit fn c t mu err iRow iFeature tlp yl
