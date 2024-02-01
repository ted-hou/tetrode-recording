
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

clear trajectories X Y L XX YY LL motPosRaw;
motPosRaw = cell(length(expReachDir), 1);

for iExp = 1:length(expReachDir)
    trials = expReachDir(iExp).eu(1).getTrials('press');
    motPosRaw{iExp} = expReachDir(iExp).eu(1).getMotorState([trials.Start]);
    nTrials = histcounts(motPosRaw{iExp}, [0.5, 1.5, 2.5, 3.5, 4.5]);
    t = pReachDir.window(1):pReachDir.resolution:pReachDir.window(2);
    for iFeat = 1:length(pReachDir.features)
        X = zeros(4, length(t));
        Y = X;
        L = X;
        
        for iPos = 1:4
            [XX, YY, LL, ~] = expReachDir(iExp).getTrajectoryByTrial('f', pReachDir.features{iFeat}, trials=trials(motPosRaw{iExp}==iPos), window=pReachDir.window, likelihoodThreshold=pReachDir.llThreshold, includeInvalid=false);
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

% Get resampled trajectories from reach onset to touch

for iExp = 1:length(expReachDir)
    trials = expReachDir(iExp).eu(1).getTrials('press');
    % Only include trials with negative tst (no zeros or NaNs)
    sel = tstReachDir{iExp} < 0;
    stopTime = reshape([trials(sel).Stop], [], 1);
    trueStartTime = stopTime + reshape(tstReachDir{iExp}(sel), [], 1);
    trajectoryTrials = Trial(trueStartTime, stopTime);
    motPosRaw = expReachDir(iExp).eu(1).getMotorState([trajectoryTrials.Start]);
    nTrials = histcounts(motPosRaw, [0.5, 1.5, 2.5, 3.5, 4.5]);

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
            X(iPos, :) = mean(XAll(motPosRaw==iPos, :), 1, 'omitnan');
            Y(iPos, :) = mean(YAll(motPosRaw==iPos, :), 1, 'omitnan');
            L(iPos, :) = mean(LAll(motPosRaw==iPos, :), 1, 'omitnan');
        end
        trajectories(iExp).(pReachDir.features{iFeat}).Resampled.X = X;
        trajectories(iExp).(pReachDir.features{iFeat}).Resampled.Y = Y;
        trajectories(iExp).(pReachDir.features{iFeat}).Resampled.L = L;
        trajectories(iExp).(pReachDir.features{iFeat}).Resampled.n = nTrials;
    end
end

% Calculate ETA, grouped by lever pos, correct with true start time
euReachDir = [expReachDir.eu];
clear etaReachDir

etaReachDir(length(euReachDir), 4) = struct('X', [], 't', [], 'N', [], 'D', [], 'stats', []);
nTrials = zeros(length(euReachDir), 4);

% Calculate ETA per motPos
for iEu = 1:length(euReachDir)
    trials = euReachDir(iEu).getTrials('press');
    tTouch = [trials.Stop];
    motPosRaw = euReachDir(iEu).getMotorState(tTouch);
    nTrials(iEu, :) = histcounts(motPosRaw, [0.5, 1.5, 2.5, 3.5, 4.5]);
    iExp = expIndex(iEu);
    for iPos = 1:4
        thisMotPosRaw = posOrder(iExp, iPos);
        etaReachDir(iEu, iPos) = euReachDir(iEu).getETA('count', 'press', window=[-4, 2], resolution=0.1, normalize=[-4,-2], ...
            alignTo='stop', includeInvalid=true, trials=trials(motPosRaw==thisMotPosRaw), correction=tstReachDir{iExp}(motPosRaw==thisMotPosRaw));
    end
    
    clear trials tTouch motPosRaw iPos thisMotPosRaw iExp
end
clear iEu

% Merge ETAs (across sessions for grand average)
clear etaReachDirMerged
etaReachDirMerged(1, 4) = struct('X', [], 't', [], 'N', [], 'D', [], 'stats', []);
for iPos = 1:4
    etaReachDirMerged(1, iPos).X = vertcat(etaReachDir(:, iPos).X);
    etaReachDirMerged(1, iPos).t = etaReachDir(1, iPos).t;
    etaReachDirMerged(1, iPos).N = NaN(length(euReachDir), 1);
    etaReachDirMerged(1, iPos).D = NaN(length(euReachDir), 1);
end

% Diff ETAs (pairwise diff matrix)
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

%% Bootstrap directional selectivity
% p.bootAlphaReachDir = 0.05;
% bootReachDirPre_005 = bootstrapReachDirETA(10000, posOrder, expIndex, euReachDir, distWindow=[-2, 0], alpha=p.bootAlphaReachDir, correction=tstReachDir);
% bootReachDirPeri_005 = bootstrapReachDirETA(10000, posOrder, expIndex, euReachDir, distWindow=[-2, 0.5], alpha=p.bootAlphaReachDir, correction=tstReachDir);
% 
% fprintf('%i/%i (%.2f%%) units had significant pre-reach direction selectivity (p < %g)\n', nnz(bootReachDirPre_005.distH), length(euReachDir), 100*nnz(bootReachDirPre_005.distH)/length(euReachDir), p.bootAlphaReachDir)
% fprintf('%i/%i (%.2f%%) units had significant peri-reach direction selectivity (p < %g)\n', nnz(bootReachDirPeri_005.distH), length(euReachDir), 100*nnz(bootReachDirPeri_005.distH)/length(euReachDir), p.bootAlphaReachDir)

p.bootAlphaReachDir = 0.01;
bootReachDirPre_001 = bootstrapReachDirETA(10000, posOrder, expIndex, euReachDir, distWindow=[-2, 0], alpha=p.bootAlphaReachDir, correction=tstReachDir);
bootReachDirPeri_001 = bootstrapReachDirETA(10000, posOrder, expIndex, euReachDir, distWindow=[-2, 0.5], alpha=p.bootAlphaReachDir, correction=tstReachDir);

fprintf('%i/%i (%.2f%%) units had significant pre-reach direction selectivity (p < %g)\n', nnz(bootReachDirPre_001.distH), length(euReachDir), 100*nnz(bootReachDirPre_001.distH)/length(euReachDir), p.bootAlphaReachDir)
fprintf('%i/%i (%.2f%%) units had significant peri-reach direction selectivity (p < %g)\n', nnz(bootReachDirPeri_001.distH), length(euReachDir), 100*nnz(bootReachDirPeri_001.distH)/length(euReachDir), p.bootAlphaReachDir)

% etaReachDir struct array sized: (nUnits, nMotPos)
function boot = bootstrapReachDirETA(nboot, posOrder, expIndex, eu, varargin)
    p = inputParser();
    p.addRequired('nboot', @isnumeric)
    p.addRequired('posOrder', @isnumeric)
    p.addRequired('expIndex', @isnumeric)
    p.addRequired('eu', @(x) isa(x, 'EphysUnit'))
    p.addOptional('sel', [], @(x) islogical(x) || isnumeric(x))
    p.addParameter('method', 'perm', @(x) ismember(lower(x), {'perm', 'boot'}))
    p.addParameter('distWindow', [-2, 0], @(x) isnumeric(x) && length(x) == 2)
    p.addParameter('alpha', 0.05, @isnumeric)
    p.addParameter('correction', {}, @iscell)

    p.parse(nboot, posOrder, expIndex, eu, varargin{:})
    r = p.Results;
    nboot = r.nboot;
    posOrder = r.posOrder;
    correction = r.correction;
    expIndex = r.expIndex;
    eu = r.eu;
    nMotPos = size(posOrder, 2); % num of motor positions
    assert(nMotPos == 4)

    assert(strcmpi(r.method, 'perm'), 'Only permutation test is implemented. If you want to do bootstrap use the builtin matlab function.')
    
    if isempty(r.sel)
        euIndices = 1:length(eu);
    elseif islogical(r.sel)
        euIndices = reshape(find(r.sel), 1, []);
    else
        euIndices = reshape(r.sel, 1, []);
    end

    ii = 0;
    boot.distH = NaN(length(eu), 1);
    boot.distCI = NaN(length(eu), 2);
    boot.distObs = NaN(length(eu), 1);

    for iEu = euIndices
        ii = ii + 1;
        fprintf(1, '%d/%d ', ii, length(euIndices))
        if mod(ii, 10) == 0
            fprintf(1, '\n')
        end
        iExp = expIndex(iEu);
        trials = eu(iEu).getTrials('press');
        if ~isempty(correction)
            trials = Trial([trials.Start], [trials.Stop] + reshape(correction{iExp}, 1, []));
        end
        [xx, tt] = eu(iEu).getTrialAlignedData('count', r.distWindow, 'press', allowedTrialDuration=[0, Inf], alignTo='stop', resolution=0.1, includeInvalid=false, trials=trials);
        dd = trials.duration;
        assert(length(dd) == size(xx, 1))
    
        motPosRaw = eu(iEu).getMotorState([trials.Stop]);
        motPosRectified = posOrder(iExp, motPosRaw);
        N = histcounts(motPosRectified, [0.5, 1.5, 2.5, 3.5, 4.5]);

        xxMean = NaN(nMotPos, length(tt));
        for iMotPos = 1:nMotPos
            xxMean(iMotPos, :) = mean(xx(motPosRectified == iMotPos, :), 1, 'omitnan');
        end

        
        pairs = nchoosek(1:nMotPos, 2);
        nPairs = size(pairs, 1);
        dist = NaN(nPairs, 1);
        for iPair = 1:nPairs
            x1 = xxMean(pairs(iPair, 1), :);
            x2 = xxMean(pairs(iPair, 2), :);
            sqrs = (x1 - x2).^2;
            dist(iPair) = mean(sqrs, 'omitnan');
        end
        dist = mean(dist);
        
        bsample = zeros(sum(N), nboot);
        for iboot = 1:nboot
            bsample(:, iboot) = randperm(sum(N));
        end
        
        bootDivEdges = [0, cumsum(N)];
        xxMeanBoot = NaN(nMotPos, length(tt), nboot);
        for iMotPos = 1:nMotPos
            for iboot = 1:nboot
                selTrialsBoot = bsample(bootDivEdges(iMotPos)+1:bootDivEdges(iMotPos+1), iboot);
                xxMeanBoot(iMotPos, :, iboot) = mean(xx(selTrialsBoot, :), 1, 'omitnan');
            end
        end
        distBoot = NaN(nPairs, nboot);
        for iPair = 1:nPairs
            x1 = squeeze(xxMeanBoot(pairs(iPair, 1), :, :));
            x2 = squeeze(xxMeanBoot(pairs(iPair, 2), :, :));
            sqrs = (x1 - x2).^2;
            distBoot(iPair, :) = mean(sqrs, 1, 'omitnan');
        end
        distBoot = mean(distBoot, 1, 'omitnan');
        distCI = prctile(distBoot, [50*r.alpha, 100-50*r.alpha]);
        distH = dist >= distCI(2) || dist <= distCI(1);
        boot.distH(iEu) = distH;
        boot.distCI(iEu, 1:2) = distCI;
        boot.distObs(iEu) = dist;
    end

end
