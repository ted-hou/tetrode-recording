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

%% Calculate trajectories for each session (traj), then combine across sessions (tarjCombined)
close all
nExp = pa2tgt.getLength('exp');
clear traj2tgt
traj2tgt(nExp) = struct(contra=[], ipsi=[], target=[], t=[], pca=[], kmeans=[]);
nt = 16;
for iExp = 1:nExp
    nFrames = 16;
    [traj2tgt(iExp).contra, traj2tgt(iExp).ipsi, traj2tgt(iExp).target, traj2tgt(iExp).t] = pa2tgt.getTrajectories(iExp, zero=nFrames - nt + 1);
    if ~ismember(pa2tgt.exp(iExp).animalName, {'desmond29', 'daisy25'})
        traj2tgt(iExp).contra.x = -traj2tgt(iExp).contra.x;
        traj2tgt(iExp).ipsi.x = -traj2tgt(iExp).ipsi.x;
    end
end

%% Combine trajectories across sessions
contraX = arrayfun(@(t) t.contra.x, traj2tgt, 'UniformOutput', false);
contraY = arrayfun(@(t) t.contra.y, traj2tgt, 'UniformOutput', false);
contraZ = arrayfun(@(t) t.contra.z, traj2tgt, 'UniformOutput', false);

contra.x = cat(1, contraX{:});
contra.y = cat(1, contraY{:});
contra.z = cat(1, contraZ{:});

ipsiX = arrayfun(@(t) t.ipsi.x, traj2tgt, 'UniformOutput', false);
ipsiY = arrayfun(@(t) t.ipsi.y, traj2tgt, 'UniformOutput', false);
ipsiZ = arrayfun(@(t) t.ipsi.z, traj2tgt, 'UniformOutput', false);

ipsi.x = cat(1, ipsiX{:});
ipsi.y = cat(1, ipsiY{:});
ipsi.z = cat(1, ipsiZ{:});

trajCombined2tgt = struct(contra=contra, ipsi=ipsi, target=[traj2tgt.target], t=traj2tgt(1).t);

%% Determine which paw was used
close all
clear r
r.contra = [...
    range(trajCombined2tgt.contra.x, 2), ...
    range(trajCombined2tgt.contra.y, 2), ...
    range(trajCombined2tgt.contra.z, 2), ...
    ];
r.ipsi = [...
    range(trajCombined2tgt.ipsi.x, 2), ...
    range(trajCombined2tgt.ipsi.y, 2), ...
    range(trajCombined2tgt.ipsi.z, 2), ...
    ];
r.contra = sqrt(sum(r.contra.^2, 2));
r.ipsi = sqrt(sum(r.ipsi.^2, 2));


ax = axes(figure);
hold(ax, 'on')
h = gobjects(2, 1);
h(1) = scatter(ax, r.contra(trajCombined2tgt.target=="contra-out"), r.ipsi(trajCombined2tgt.target=="contra-out"), 10, getColor(1, 4, 0.8), DisplayName="contra-out");
h(2) = scatter(ax, r.contra(trajCombined2tgt.target=="contra-in"), r.ipsi(trajCombined2tgt.target=="contra-in"), 10, getColor(3, 4, 0.8), DisplayName="contra-in");
xlabel(ax, 'contra range')
ylabel(ax, 'ipsi range')
axis(ax, 'equal')
set(ax, XLim=[0, 150], YLim=[0, 150])

mdl = fitlm(r.contra, r.ipsi);
plot(ax, [0, 300], mdl.predict([0; 300]), 'k--', DisplayName='LM')

legend(ax)

trajCombined2tgt.usedIpsiPaw = r.ipsi' > mdl.predict(r.contra)';
fprintf('Out of %i trials:\n', length(trajCombined2tgt.target));
selPaw = trajCombined2tgt.usedIpsiPaw;
selTarget = trajCombined2tgt.target == "contra-out";
fprintf('%i trials when target was at "contra-out", we estimate %i (%.1f%%) used contra paw, %i (%.1f%%) used ipsi paw.\n', nnz(selTarget), nnz(selTarget & ~selPaw), nnz(selTarget & ~selPaw)./nnz(selTarget).*100, nnz(selTarget & selPaw), nnz(selTarget & selPaw)./nnz(selTarget).*100);
selTarget = trajCombined2tgt.target == "contra-in";
fprintf('%i trials when target was at "contra-in", we estimate %i (%.1f%%) used contra paw, %i (%.1f%%) used ipsi paw.\n', nnz(selTarget), nnz(selTarget & ~selPaw), nnz(selTarget & ~selPaw)./nnz(selTarget).*100, nnz(selTarget & selPaw), nnz(selTarget & selPaw)./nnz(selTarget).*100);

clear b ax h

%% Calculate ETA using movement intiation correction
% Calculate trueStartTime
nExp = pa2tgt.getLength('exp');
onsetThreshold = 0.25;
targetNames = {'contra-out', 'contra-in'};
trueStartTime2tgts = cell(nExp, 1);
sessionEdges = cumsum([0, arrayfun(@(traj) length(traj.target), traj2tgt)]); %Right inclusive trial indices for each session
for iExp = 1:nExp
    nTrials = length(traj2tgt(iExp).target);
    trueStartTime2tgts{iExp} = NaN(nTrials, 1);
    t = traj2tgt(iExp).t;
    dist = zeros(nTrials, length(t));
    for iTrial = 1:nTrials
        pos = [ ...
            traj2tgt(iExp).contra.x(iTrial, :); ...
            traj2tgt(iExp).contra.y(iTrial, :); ...
            traj2tgt(iExp).contra.z(iTrial, :); ...
            ];
        dist(iTrial, :) = sqrt(sum(pos.^2, 1));
    end
    distNorm = (dist - median(dist(:), 'omitnan'))./(mad(dist(:), 1)./0.6745); % Normalize using session median/mad
    for iTrial = 1:nTrials
        if any(isnan(distNorm(iTrial, :)))
            continue
        end
        
        isAbove = distNorm(iTrial, :) >= onsetThreshold;
        iLastAbove = find(isAbove, 1, 'last'); % We don't do abs since dist is already positive, although z-scoring will generate negatives, true movement should be positive.
        if isempty(iLastAbove)
            continue
        end
        iOnset = strfind(isAbove(1:iLastAbove), [0 1 1]) + 1;
        if isempty(iOnset)
            continue
        end
        iOnset = iOnset(end);
        trueStartTime2tgts{iExp}(iTrial) = t(iOnset) - t(end);
    end

    % Calculate ETA for this session
    usedIpsiPawInSession = trajCombined2tgt.usedIpsiPaw(sessionEdges(iExp) + 1:sessionEdges(iExp + 1));

    [~, targetIndex] = ismember(traj2tgt(iExp).target, targetNames);
    targetIndex = targetIndex';
    for iTarget = 1:2
        selTrials = targetIndex==iTarget & ~usedIpsiPawInSession';
        traj2tgt(iExp).eta(iTarget) = pa2tgt.exp(iExp).eu.getETA('count', 'press_spontaneous', window=[-4, 0.5], resolution=0.1, normalize=[-4, -2], alignTo='stop', includeInvalid=true, ...
            trials=pa2tgt.exp(iExp).eu(1).Trials.Press(selTrials), correction=trueStartTime2tgts{iExp}(selTrials));
        traj2tgt(iExp).trials{iTarget} = pa2tgt.exp(iExp).eu(1).Trials.Press(selTrials);
        traj2tgt(iExp).correction{iTarget} = trueStartTime2tgts{iExp}(selTrials);

        selTrials = targetIndex==iTarget & usedIpsiPawInSession';
        traj2tgt(iExp).etaIpsiPaw(iTarget) = pa2tgt.exp(iExp).eu.getETA('count', 'press_spontaneous', window=[-4, 0.5], resolution=0.1, normalize=[-4, -2], alignTo='stop', includeInvalid=true, ...
            trials=pa2tgt.exp(iExp).eu(1).Trials.Press(selTrials), correction=trueStartTime2tgts{iExp}(selTrials)); 
        traj2tgt(iExp).trialsIpsiPaw{iTarget} = pa2tgt.exp(iExp).eu(1).Trials.Press(selTrials);
        traj2tgt(iExp).correctionIpsiPaw{iTarget} = trueStartTime2tgts{iExp}(selTrials);       
    end
end
clear usedIpsiPawInSession targetIndex iTarget selTrials iExp sessionEdges dist distNorm iTrial isAbove iLastAbove iOnset

% Combine ETAs across sessions
[~, euExpIndices] = ismember({euReachDir2Tgt.ExpName}, {pa2tgt.exp.name});

for iTarget = 1:2
    X = arrayfun(@(traj) traj.eta(iTarget).X, traj2tgt, 'UniformOutput', false);
    N = arrayfun(@(traj) traj.eta(iTarget).N, traj2tgt, 'UniformOutput', false);
    stats = arrayfun(@(traj) traj.eta(iTarget).stats, traj2tgt, 'UniformOutput', false);
    trajCombined2tgt.eta(iTarget).X = cat(1, X{:});
    trajCombined2tgt.eta(iTarget).t = traj2tgt(1).eta(iTarget).t;
    trajCombined2tgt.eta(iTarget).N = cat(1, N{:});
    trajCombined2tgt.eta(iTarget).stats = cat(2, stats{:});
    trajCombined2tgt.eta(iTarget).target = targetNames{iTarget};
    trialsByExp = arrayfun(@(traj) traj.trials{iTarget}, traj2tgt, UniformOutput=false);
    correctionByExp = arrayfun(@(traj) traj.correction{iTarget}, traj2tgt, UniformOutput=false);
    trajCombined2tgt.trials{iTarget} = trialsByExp(euExpIndices);
    trajCombined2tgt.correction{iTarget} = correctionByExp(euExpIndices);
end

for iTarget = 1:2
    X = arrayfun(@(traj) traj.etaIpsiPaw(iTarget).X, traj2tgt, 'UniformOutput', false);
    N = arrayfun(@(traj) traj.etaIpsiPaw(iTarget).N, traj2tgt, 'UniformOutput', false);
    stats = arrayfun(@(traj) traj.etaIpsiPaw(iTarget).stats, traj2tgt, 'UniformOutput', false);
    trajCombined2tgt.etaIpsiPaw(iTarget).X = cat(1, X{:});
    trajCombined2tgt.etaIpsiPaw(iTarget).t = traj2tgt(1).etaIpsiPaw(iTarget).t;
    trajCombined2tgt.etaIpsiPaw(iTarget).N = cat(1, N{:});
    trajCombined2tgt.etaIpsiPaw(iTarget).stats = cat(2, stats{:});
    trajCombined2tgt.etaIpsiPaw(iTarget).target = targetNames{iTarget};
end

