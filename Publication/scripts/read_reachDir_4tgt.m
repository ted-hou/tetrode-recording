%% Load EphysUnits
euReachDir4Tgt = EphysUnit.load('C:\SERVER\Units\acute_3cam_reach_direction\SingleUnits_NonDuplicate');
% pa = Pawnalyzer2(eu, refEvent='press');


% Make trajectories
paReachDir4Tgt = Pawnalyzer2(euReachDir4Tgt, refEvent='cue');
paReachDir4Tgt.getClips(noImage=true, nFramesBefore=15, nFramesAfter=0, keepData=false, trials='Press');
paReachDir4Tgt.load('C:\SERVER\Units\acute_3cam_reach_direction\Pawnalyzer2\pa.mat')

%% Analysis
%% Calculate trajectories for each session (traj), then combine across sessions (tarjCombined)
close all
nExp = paReachDir4Tgt.getLength('exp');
clear traj
traj(nExp) = struct(contra=[], ipsi=[], target=[], t=[], pca=[], kmeans=[]);
nt = 10;
for iExp = 1:nExp
    nFrames = paReachDir4Tgt.getLength('frame', exp=iExp, trial=1);
    [traj(iExp).contra, traj(iExp).ipsi, traj(iExp).target, traj(iExp).t] = paReachDir4Tgt.getTrajectories(iExp, zero=nFrames - nt + 1);
    if ~ismember(paReachDir4Tgt.exp(iExp).animalName, {'desmond29', 'daisy25'})
        traj(iExp).contra.x = -traj(iExp).contra.x;
        traj(iExp).ipsi.x = -traj(iExp).ipsi.x;
    end
end

% Combine trajectories across sessions
contraX = arrayfun(@(t) t.contra.x, traj, 'UniformOutput', false);
contraY = arrayfun(@(t) t.contra.y, traj, 'UniformOutput', false);
contraZ = arrayfun(@(t) t.contra.z, traj, 'UniformOutput', false);

contra.x = cat(1, contraX{:});
contra.y = cat(1, contraY{:});
contra.z = cat(1, contraZ{:});

ipsiX = arrayfun(@(t) t.ipsi.x, traj, 'UniformOutput', false);
ipsiY = arrayfun(@(t) t.ipsi.y, traj, 'UniformOutput', false);
ipsiZ = arrayfun(@(t) t.ipsi.z, traj, 'UniformOutput', false);

ipsi.x = cat(1, ipsiX{:});
ipsi.y = cat(1, ipsiY{:});
ipsi.z = cat(1, ipsiZ{:});

trajCombined = struct(contra=contra, ipsi=ipsi, target=[traj.target], t=traj(1).t);

%% Calculate ETA using movement intiation correction
% Calculate trueStartTime
targetNames = {'contra-out', 'contra-front', 'contra-in', 'ipsi-front'};
trueStartTime4tgts = cell(nExp, 1);
for iExp = 1:nExp
    nTrials = length(traj(iExp).target);
    trueStartTime4tgts{iExp} = NaN(nTrials, 1);
    t = traj(iExp).t;
    for iTrial = 1:nTrials
        data = [ ...
            traj(iExp).contra.x(iTrial, :); ...
            traj(iExp).contra.y(iTrial, :); ...
            traj(iExp).contra.z(iTrial, :); ...
            traj(iExp).ipsi.x(iTrial, :); ...
            traj(iExp).ipsi.y(iTrial, :); ...
            traj(iExp).ipsi.z(iTrial, :); ...
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
            trueStartTime4tgts{iExp}(iTrial) = t(trueStartIndex) - t(end);
        end
    end
    [~, targetIndex] = ismember(traj(iExp).target, targetNames);
    targetIndex = targetIndex';
    for iTarget = 1:4
        sel = targetIndex==iTarget;
        traj(iExp).eta(iTarget) = paReachDir4Tgt.exp(iExp).eu.getETA('count', 'press', window=[-4, 0.5], resolution=0.1, normalize=[-4, -0.5], alignTo='stop', includeInvalid=true, ...
            trials=paReachDir4Tgt.exp(iExp).eu(1).Trials.Press(sel), correction=trueStartTime4tgts{iExp}(sel));
    end
end

%% Combine ETAs across sessions
for iTarget = 1:4
    X = arrayfun(@(traj) traj.eta(iTarget).X, traj, 'UniformOutput', false);
    N = arrayfun(@(traj) traj.eta(iTarget).N, traj, 'UniformOutput', false);
    stats = arrayfun(@(traj) traj.eta(iTarget).stats, traj, 'UniformOutput', false);
    trajCombined.eta(iTarget).X = cat(1, X{:});
    trajCombined.eta(iTarget).t = traj(1).eta(iTarget).t;
    trajCombined.eta(iTarget).N = cat(1, N{:});
    trajCombined.eta(iTarget).stats = cat(2, stats{:});
    trajCombined.eta(iTarget).target = targetNames{iTarget};
end

