%% Load EphysUnits
euReachDir2Tgt = EphysUnit.load('C:\SERVER\Units\acute_3cam_reach_direction_2tgts\SingleUnits_NonDuplicate');
for iEu = 1:length(euReachDir2Tgt)
    trials = euReachDir2Tgt(iEu).makeTrials('press_spontaneous2');
    goodTrials = trials(trials.duration() > 4);
    euReachDir2Tgt(iEu).Trials.PressSpontaneous = goodTrials;
end

pa2tgt = Pawnalyzer2(euReachDir2Tgt, refEvent='press');


% Make trajectories
pa2tgt.getClips(noImage=true, nFramesBefore=15, nFramesAfter=0, keepData=false, trials='PressSpontaneous');
pa2tgt.load('auto')

%% Analysis

%% Select good sessions (remove when sofia finishes labelling other sessions)
goodExpNames = ["daisy23_20240627", "daisy24_20240626", "daisy25_20240701"];


%% Calculate trajectories for each session (traj), then combine across sessions (tarjCombined)
close all
nExp = pa2tgt.getLength('exp');
clear traj2tgt
traj2tgt(nExp) = struct(contra=[], ipsi=[], target=[], t=[], pca=[], kmeans=[]);
nt = 10;
for iExp = 1:nExp
    if ~ismember(pa2tgt.exp(iExp).name, goodExpNames)
        continue
    end
    nFrames = 16;
    [traj2tgt(iExp).contra, traj2tgt(iExp).ipsi, traj2tgt(iExp).target, traj2tgt(iExp).t] = pa2tgt.getTrajectories(iExp, zero=nFrames - nt + 1);
    if ~ismember(pa2tgt.exp(iExp).animalName, {'desmond29', 'daisy25'})
        traj2tgt(iExp).contra.x = -traj2tgt(iExp).contra.x;
        traj2tgt(iExp).ipsi.x = -traj2tgt(iExp).ipsi.x;
    end
end

%% Combine trajectories across sessions
selExp = find(ismember({pa2tgt.exp.name}, goodExpNames));

contraX = arrayfun(@(t) t.contra.x, traj2tgt(selExp), 'UniformOutput', false);
contraY = arrayfun(@(t) t.contra.y, traj2tgt(selExp), 'UniformOutput', false);
contraZ = arrayfun(@(t) t.contra.z, traj2tgt(selExp), 'UniformOutput', false);

contra.x = cat(1, contraX{:});
contra.y = cat(1, contraY{:});
contra.z = cat(1, contraZ{:});

ipsiX = arrayfun(@(t) t.ipsi.x, traj2tgt(selExp), 'UniformOutput', false);
ipsiY = arrayfun(@(t) t.ipsi.y, traj2tgt(selExp), 'UniformOutput', false);
ipsiZ = arrayfun(@(t) t.ipsi.z, traj2tgt(selExp), 'UniformOutput', false);

ipsi.x = cat(1, ipsiX{:});
ipsi.y = cat(1, ipsiY{:});
ipsi.z = cat(1, ipsiZ{:});

trajCombined2tgt = struct(contra=contra, ipsi=ipsi, target=[traj2tgt(selExp).target], t=traj2tgt(selExp(1)).t);

%% Calculate ETA using movement intiation correction
% Calculate trueStartTime
targetNames = {'contra-out', 'contra-in'};
trueStartTime2tgts = cell(nExp, 1);
for iExp = 1:nExp
    if ~ismember(pa2tgt.exp(iExp).name, goodExpNames)
        continue
    end
    nTrials = length(traj2tgt(iExp).target);
    trueStartTime2tgts{iExp} = NaN(nTrials, 1);
    t = traj2tgt(iExp).t;
    for iTrial = 1:nTrials
        data = [ ...
            traj2tgt(iExp).contra.x(iTrial, :); ...
            traj2tgt(iExp).contra.y(iTrial, :); ...
            traj2tgt(iExp).contra.z(iTrial, :); ...
            traj2tgt(iExp).ipsi.x(iTrial, :); ...
            traj2tgt(iExp).ipsi.y(iTrial, :); ...
            traj2tgt(iExp).ipsi.z(iTrial, :); ...
            ];
        dist = sqrt(sum(data.^2, 1));
        normDist = normalize(dist, 'zscore', 'robust')*0.6745;
        if any(isnan(normDist))
            continue
        end
        lastIndexAboveThreshold = find(normDist > 2.5, 1, 'last'); % We don't do abs since dist is already positive, although z-scoring will generate negatives, true movement should be positive.
        if isempty(lastIndexAboveThreshold)
            continue
        end
        trueStartIndex = find((normDist < 2.5) & ((1:nFrames) <= lastIndexAboveThreshold), 1, 'last');
        if ~isempty(trueStartIndex)
            trueStartTime2tgts{iExp}(iTrial) = t(trueStartIndex) - t(end);
        end
    end
    [~, targetIndex] = ismember(traj2tgt(iExp).target, targetNames);
    targetIndex = targetIndex';
    for iTarget = 1:2
        sel = targetIndex==iTarget;
        traj2tgt(iExp).eta(iTarget) = pa2tgt.exp(iExp).eu.getETA('count', 'press_spontaneous', window=[-4, 0.5], resolution=0.1, normalize=[-4, -0.5], alignTo='stop', includeInvalid=true, ...
            trials=pa2tgt.exp(iExp).eu(1).Trials.Press(sel), correction=trueStartTime2tgts{iExp}(sel));
    end
end

%% Combine ETAs across sessions
for iTarget = 1:2
    X = arrayfun(@(traj) traj.eta(iTarget).X, traj2tgt(selExp), 'UniformOutput', false);
    N = arrayfun(@(traj) traj.eta(iTarget).N, traj2tgt(selExp), 'UniformOutput', false);
    stats = arrayfun(@(traj) traj.eta(iTarget).stats, traj2tgt(selExp), 'UniformOutput', false);
    trajCombined2tgt.eta(iTarget).X = cat(1, X{:});
    trajCombined2tgt.eta(iTarget).t = traj2tgt(selExp(1)).eta(iTarget).t;
    trajCombined2tgt.eta(iTarget).N = cat(1, N{:});
    trajCombined2tgt.eta(iTarget).stats = cat(2, stats{:});
    trajCombined2tgt.eta(iTarget).target = targetNames{iTarget};
end

