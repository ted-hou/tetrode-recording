
%% Load EU
euReachDir = EphysUnit.load('C:\SERVER\Units\acute_3cam_reach_direction');

%% Make CompleteExperiment3 from EU
expReachDir = CompleteExperiment3(euReachDir);
expReachDir.alignTimestamps();
posOrder = zeros(length(expReachDir), 4);
posNames = {'contra-out', 'contra-front', 'contra-in', 'ipsi-front'};
expIndex = zeros(length(euReachDir), 1);

i = 1;
% Correct feature names
for iExp = 1:length(expReachDir)
    expIndex(i:i+length(expReachDir(iExp).eu)-1) = iExp;
    i = i + length(expReachDir(iExp).eu);
    switch expReachDir(iExp).animalName
        case {'desmond28', 'desmond30'}
            expReachDir(iExp).vtdL = renamevars(expReachDir(iExp).vtdL, ...
                {'handIpsiCam_X', 'handIpsiCam_Y', 'handIpsiCam_Likelihood', 'footIpsiCam_X', 'footIpsiCam_Y', 'footIpsiCam_Likelihood'}, ...
                {'handContra_X', 'handContra_Y', 'handContra_Likelihood', 'footContra_X', 'footContra_Y', 'footContra_Likelihood'});
            expReachDir(iExp).vtdR = renamevars(expReachDir(iExp).vtdR, ...
                {'handIpsiCam_X', 'handIpsiCam_Y', 'handIpsiCam_Likelihood', 'footIpsiCam_X', 'footIpsiCam_Y', 'footIpsiCam_Likelihood'}, ...
                {'handIpsi_X', 'handIpsi_Y', 'handIpsi_Likelihood', 'footIpsi_X', 'footIpsi_Y', 'footIpsi_Likelihood'});
            posOrder(iExp, :) = [4, 3, 2, 1];

        case 'desmond29'
            expReachDir(iExp).vtdF = renamevars(expReachDir(iExp).vtdF, ...
                {'handIpsi_X', 'handIpsi_Y', 'handIpsi_Likelihood', 'footIpsi_X', 'footIpsi_Y', 'footIpsi_Likelihood', ...
                'handContra_X', 'handContra_Y', 'handContra_Likelihood', 'footContra_X', 'footContra_Y', 'footContra_Likelihood'}, ...
                {'handContra_X', 'handContra_Y', 'handContra_Likelihood', 'footContra_X', 'footContra_Y', 'footContra_Likelihood', ...
                'handIpsi_X', 'handIpsi_Y', 'handIpsi_Likelihood', 'footIpsi_X', 'footIpsi_Y', 'footIpsi_Likelihood'});
            expReachDir(iExp).vtdL = renamevars(expReachDir(iExp).vtdL, ...
                {'handIpsiCam_X', 'handIpsiCam_Y', 'handIpsiCam_Likelihood', 'footIpsiCam_X', 'footIpsiCam_Y', 'footIpsiCam_Likelihood'}, ...
                {'handIpsi_X', 'handIpsi_Y', 'handIpsi_Likelihood', 'footIpsi_X', 'footIpsi_Y', 'footIpsi_Likelihood'});
            expReachDir(iExp).vtdR = renamevars(expReachDir(iExp).vtdR, ...
                {'handIpsiCam_X', 'handIpsiCam_Y', 'handIpsiCam_Likelihood', 'footIpsiCam_X', 'footIpsiCam_Y', 'footIpsiCam_Likelihood'}, ...
                {'handContra_X', 'handContra_Y', 'handContra_Likelihood', 'footContra_X', 'footContra_Y', 'footContra_Likelihood'});
            posOrder(iExp, :) = [1, 2, 3, 4];            
    end
end
clear iExp i

%% Gather 2d trajctories from front cam
pReachDir.resolution = 1/30;
pReachDir.window = [-10, 0.5];
pReachDir.features = {'handIpsi', 'handContra', 'jaw', 'lever'};
pReachDir.llThreshold = 0.90;
pReachDir.nTrajectorySamples = 100;

clear trajectories X Y L XX YY LL motPos;
motPos = cell(length(expReachDir), 1);

for iExp = 1:length(expReachDir)
    trials = expReachDir(iExp).eu(1).getTrials('press');
    motPos{iExp} = expReachDir(iExp).eu(1).getMotorState([trials.Start]);
    nTrials = histcounts(motPos{iExp}, [0.5, 1.5, 2.5, 3.5, 4.5]);
    t = pReachDir.window(1):pReachDir.resolution:pReachDir.window(2);
    for iFeat = 1:length(pReachDir.features)
        X = zeros(4, length(t));
        Y = X;
        L = X;
        
        for iPos = 1:4
            [XX, YY, LL, ~] = expReachDir(iExp).getTrajectoryByTrial('f', pReachDir.features{iFeat}, trials=trials(motPos{iExp}==iPos), window=pReachDir.window, likelihoodThreshold=pReachDir.llThreshold, includeInvalid=false);
            X(iPos, :) = mean(XX, 1, 'omitnan');
            Y(iPos, :) = mean(YY, 1, 'omitnan');
            L(iPos, :) = mean(LL, 1, 'omitnan');
        end
        [XAll, YAll, LLAll] = expReachDir(iExp).getTrajectoryByTrial('f', pReachDir.features{iFeat}, trials=trials, window=pReachDir.window, likelihoodThreshold=pReachDir.llThreshold, includeInvalid=false);


        trajectories(iExp).(pReachDir.features{iFeat}).X = X;
        trajectories(iExp).(pReachDir.features{iFeat}).Y = Y;
        trajectories(iExp).(pReachDir.features{iFeat}).L = L;
        trajectories(iExp).(pReachDir.features{iFeat}).XAll= XAll;
        trajectories(iExp).(pReachDir.features{iFeat}).YAll= YAll;
        trajectories(iExp).(pReachDir.features{iFeat}).LAll= LLAll;
        trajectories(iExp).(pReachDir.features{iFeat}).t = t;
        trajectories(iExp).(pReachDir.features{iFeat}).n = nTrials;
    end
end

% Get TrueStartTime
pReachDir.tstThreshold = 2.5;
tstReachDir = cell(length(expReachDir), 1);
for iExp = 1:length(expReachDir)
    trials = expReachDir(iExp).eu(1).getTrials('press');
    trueStartTime = NaN(length(trials), 1);
    t = trajectories(iExp).handContra.t;
    for iTrial = 1:length(trials)
        data = [ ...
            normalize(trajectories(iExp).handContra.XAll(iTrial, :), 'zscore', 'robust'); ...
            normalize(trajectories(iExp).handContra.YAll(iTrial, :), 'zscore', 'robust'); ...
            normalize(trajectories(iExp).handIpsi.XAll(iTrial, :), 'zscore', 'robust'); ...
            normalize(trajectories(iExp).handIpsi.YAll(iTrial, :), 'zscore', 'robust'); ...
            ];
        trueStartIndex = find(all(abs(data) < pReachDir.tstThreshold, 1), 1, 'last');
        if ~isempty(trueStartIndex)
            trueStartTime(iTrial) = t(trueStartIndex) - t(end);
        end
    end
    tstReachDir{iExp} = trueStartTime;
end

%% Get resampled trajectories from reach onset to touch

for iExp = 1:length(expReachDir)
    trials = expReachDir(iExp).eu(1).getTrials('press');
    % Only include trials with negative tst (no zeros or NaNs)
    sel = tstReachDir{iExp} < 0;
    stopTime = reshape([trials(sel).Stop], [], 1);
    trueStartTime = stopTime + reshape(tstReachDir{iExp}(sel), [], 1);
    trajectoryTrials = Trial(trueStartTime, stopTime);
    motPos = expReachDir(iExp).eu(1).getMotorState([trajectoryTrials.Start]);
    nTrials = histcounts(motPos, [0.5, 1.5, 2.5, 3.5, 4.5]);

    for iFeat = 1:length(pReachDir.features)
        XAll = zeros(length(trajectoryTrials), pReachDir.nTrajectorySamples);
        YAll = XAll;
        LAll = XAll;
        TAll = XAll;

        X = zeros(4, pReachDir.nTrajectorySamples);
        Y = X;
        L = X;
        
        for iTrial = 1:length(trajectoryTrials)
            t = linspace(trajectoryTrials(iTrial).Start, trajectoryTrials(iTrial).Stop, pReachDir.nTrajectorySamples);
            [XAll(iTrial, :), YAll(iTrial, :), LAll(iTrial, :)] = expReachDir(iExp).getTrajectory(t, 'f', pReachDir.features{iFeat}, likelihoodThreshold=pReachDir.llThreshold);
            TAll(iTrial, :) = t;
        end
        trajectories(iExp).(pReachDir.features{iFeat}).Resampled.XAll = XAll;
        trajectories(iExp).(pReachDir.features{iFeat}).Resampled.YAll = YAll;
        trajectories(iExp).(pReachDir.features{iFeat}).Resampled.LAll = LAll;
        trajectories(iExp).(pReachDir.features{iFeat}).Resampled.TAll = TAll;

        for iPos = 1:4
            X(iPos, :) = mean(XAll(motPos==iPos, :), 1, 'omitnan');
            Y(iPos, :) = mean(YAll(motPos==iPos, :), 1, 'omitnan');
            L(iPos, :) = mean(LAll(motPos==iPos, :), 1, 'omitnan');
        end
        trajectories(iExp).(pReachDir.features{iFeat}).Resampled.X = X;
        trajectories(iExp).(pReachDir.features{iFeat}).Resampled.Y = Y;
        trajectories(iExp).(pReachDir.features{iFeat}).Resampled.L = L;
        trajectories(iExp).(pReachDir.features{iFeat}).Resampled.n = nTrials;
    end
end

%% Calculate ETA, grouped by lever pos, correct with true start time
euReachDir = [expReachDir.eu];
clear etaReachDir

etaReachDir(length(euReachDir), 4) = struct('X', [], 't', [], 'N', [], 'D', [], 'stats', []);
nTrials = zeros(length(euReachDir), 4);

for iEu = 1:length(euReachDir)
    trials = euReachDir(iEu).getTrials('press');
    tTouch = [trials.Stop];
    motPos = euReachDir(iEu).getMotorState(tTouch);
    nTrials(iEu, :) = histcounts(motPos, [0.5, 1.5, 2.5, 3.5, 4.5]);
    iExp = expIndex(iEu);
    for iPos = 1:4
        rawPos = posOrder(iExp, iPos);
        etaReachDir(iEu, iPos) = euReachDir(iEu).getETA('count', 'press', window=[-4, 2], resolution=0.1, normalize=[-4,-2], ...
            alignTo='stop', includeInvalid=true, trials=trials(motPos==rawPos), correction=tstReachDir{iExp}(motPos==rawPos));
    end
    
    clear trials tTouch motPos iPos rawPos iExp
end
clear iEu

% Merge ETAs
clear etaReachDirMerged
etaReachDirMerged(1, 4) = struct('X', [], 't', [], 'N', [], 'D', [], 'stats', []);
for iPos = 1:4
    etaReachDirMerged(1, iPos).X = vertcat(etaReachDir(:, iPos).X);
    etaReachDirMerged(1, iPos).t = etaReachDir(1, iPos).t;
    etaReachDirMerged(1, iPos).N = NaN(length(euReachDir), 1);
    etaReachDirMerged(1, iPos).D = NaN(length(euReachDir), 1);
end

% Diff ETAs
clear etaReachDirDiff
etaReachDirDiff(3, 4) = struct('X', [], 't', [], 'N', [], 'D', [], 'stats', []);
for i = 1:3
    for j = i+1:4
        etaReachDirDiff(i, j).X = etaReachDirMerged(j).X - etaReachDirMerged(i).X;
        etaReachDirDiff(i, j).t = etaReachDirMerged(i).t;
        etaReachDirDiff(i, j).N = etaReachDirMerged(i).N;
        etaReachDirDiff(i, j).D = etaReachDirMerged(i).D;
    end
end

% Calculate META
clear metaReachDir
metaReachDir(1, 4) = struct('pressDirectional', []);
for iPos = 1:4
    t = etaReachDirMerged(1, iPos).t;
    sel = t <= 0.5 & t >= -0.5;
    metaReachDir(iPos).pressDirectional = mean(etaReachDirMerged(iPos).X(:, sel), 2, 'omitnan');
end