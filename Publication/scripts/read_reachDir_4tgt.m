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
nt = 16;
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

trajCombined = struct(contra=contra, ipsi=ipsi, target=[traj.target]', t=traj(1).t);

%% Determine which paw was used
close all
clear r
r.contra = [...
    range(trajCombined.contra.x, 2), ...
    range(trajCombined.contra.y, 2), ...
    range(trajCombined.contra.z, 2), ...
    ];
r.ipsi = [...
    range(trajCombined.ipsi.x, 2), ...
    range(trajCombined.ipsi.y, 2), ...
    range(trajCombined.ipsi.z, 2), ...
    ];
r.contra = sqrt(sum(r.contra.^2, 2));
r.ipsi = sqrt(sum(r.ipsi.^2, 2));


% Draw quadrants to separate trials into (neither, ipsi-only, contra-only,
% both)
% Visualize the thresholds
target = trajCombined.target;
contra0 = quantile(r.contra(target=="ipsi-front"), 0.5);
ipsi0 = quantile(r.ipsi(ismember(target, ["contra-front", "contra-in", "contra-out"])), 0.5);

figure
ax = scatterhist(r.contra, r.ipsi, Group=target, Color=getColor(1:4, 4, 0.8), ...
    NBins=50, Marker='o', MarkerSize=4, LineStyle='-');
ax(1).Legend.AutoUpdate = false;
hold(ax, 'on')
plot(ax(1), [contra0, contra0], [0, 150], 'k--')
plot(ax(1), [0, 150], [ipsi0, ipsi0], 'k--')

plot(ax(2), [contra0, contra0], [0, 1], 'k--')
plot(ax(3), [ipsi0, ipsi0], [0, 1], 'k--')

xlim(ax(1:3), [0, 150])
ylim(ax(1), [0, 150])
xlabel(ax, 'contra paw movement')
ylabel(ax, 'ipsi paw movement')

% Categorize trials
nTrials = length(target);
trajCombined.paw = repmat("", nTrials, 1);
trajCombined.paw(r.contra >= contra0 & r.ipsi < ipsi0) = "contra";
trajCombined.paw(r.contra >= contra0 & r.ipsi >= ipsi0) = "both";
trajCombined.paw(r.contra < contra0 & r.ipsi >= ipsi0) = "ipsi";
trajCombined.paw(r.contra < contra0 & r.ipsi < ipsi0) = "neither";

paw = trajCombined.paw;

PAWS = ["contra", "both", "ipsi"];
TARGETS = ["contra-out", "contra-front", "contra-in", "ipsi-front"];
T = nan(length(TARGETS), length(PAWS));
for iTarget = 1:length(TARGETS)
    for iPaw = 1:length(PAWS)
        T(iTarget, iPaw) = nnz(target == TARGETS(iTarget) & paw == PAWS(iPaw));
    end
end
T = array2table(T, VariableNames=PAWS, RowNames=TARGETS);
disp(T)
clear contra0 ipsi0 ax PAWS TARGETS T
clear target paw

%% Calculate movement initiation time for each paw
onsetThreshold = 0.25;
onset = struct(contra=NaN(nTrials, 1), ipsi=NaN(nTrials, 1));
t = trajCombined.t;
for iTrial = 1:nTrials
    posContra = [ ...
        trajCombined.contra.x(iTrial, :); ...
        trajCombined.contra.y(iTrial, :); ...
        trajCombined.contra.z(iTrial, :); ...
        ];
    posIpsi = [ ...
        trajCombined.ipsi.x(iTrial, :); ...
        trajCombined.ipsi.y(iTrial, :); ...
        trajCombined.ipsi.z(iTrial, :); ...
        ];
    distContra(iTrial, :) = sqrt(sum(posContra.^2, 1));
    distIpsi(iTrial, :) = sqrt(sum(posIpsi.^2, 1));
    distContraNorm = (distContra - median(distContra(:), 'omitnan'))./(mad(distContra(:), 1)./0.6745);
    distIpsiNorm = (distIpsi - median(distIpsi(:), 'omitnan'))./(mad(distIpsi(:), 1)./0.6745);
end

for iTrial = 1:nTrials
    if any(isnan(distContraNorm(iTrial, :)))
        continue
    end
    
    isAbove = distContraNorm(iTrial, :) >= onsetThreshold;
    iLastAbove = find(isAbove, 1, 'last'); % We don't do abs since dist is already positive, although z-scoring will generate negatives, true movement should be positive.
    if isempty(iLastAbove)
        continue
    end
    iOnset = strfind(isAbove(1:iLastAbove), [0 1 1]) + 1;
    if isempty(iOnset)
        continue
    end
    iOnset = iOnset(end);
    onset.contra(iTrial) = t(iOnset) - t(end);
end

for iTrial = 1:nTrials
    if any(isnan(distIpsiNorm(iTrial, :)))
        continue
    end
    
    isAbove = distIpsiNorm(iTrial, :) >= onsetThreshold;
    iLastAbove = find(isAbove, 1, 'last'); % We don't do abs since dist is already positive, although z-scoring will generate negatives, true movement should be positive.
    if isempty(iLastAbove)
        continue
    end
    iOnset = strfind(isAbove(1:iLastAbove), [0 1 1]) + 1;
    if isempty(iOnset)
        continue
    end
    iOnset = iOnset(end);
    onset.ipsi(iTrial) = t(iOnset) - t(end);
end
onset.both = min(onset.contra, onset.ipsi, 'omitnan');


clear t iTrial posContra posIpsi distContra distIpsi distContraNorm distIpsiNorm isAbove iLastAbove iOnset

%% Calculate ETA using movement intiation correction
minNumTrials = 4;
targetNames = ["contra-out", "contra-front", "contra-in", "ipsi-front"];
pawNames = ["contra", "both", "ipsi"];
sessionEdges = cumsum([0, cellfun(@length, {traj.target})]);
for iExp = 1:nExp
    iTrialsInExp = sessionEdges(iExp)+1:sessionEdges(iExp+1);
    targetExp = trajCombined.target(iTrialsInExp);
    pawExp = trajCombined.paw(iTrialsInExp);
    onsetExp = struct( ...
        contra=onset.contra(iTrialsInExp), ...
        ipsi=onset.ipsi(iTrialsInExp), ...
        both=onset.both(iTrialsInExp) ...
        );
    for iTarget = 1:4
        for iPaw = 1:3
            selTrials = targetExp==targetNames(iTarget) & pawExp==pawNames(iPaw);
            if nnz(selTrials) >= minNumTrials
                traj(iExp).eta(iTarget, iPaw) = paReachDir4Tgt.exp(iExp).eu.getETA('count', 'press', window=[-4, 0.5], resolution=0.1, normalize=[-4, -2], alignTo='stop', includeInvalid=true, ...
                    trials=paReachDir4Tgt.exp(iExp).eu(1).Trials.Press(selTrials), correction=onsetExp.(pawNames(iPaw))(selTrials));
            else
                t = (-4:0.1:0.4)+0.05;
                traj(iExp).eta(iTarget, iPaw) = struct(X=NaN(length(paReachDir4Tgt.exp(iExp).eu), length(t)), t=t, N=repmat(nnz(selTrials), [length(paReachDir4Tgt.exp(iExp).eu), 1]), D=[], stats=[]);
            end
        end
    end
end

% Combine ETAs across sessions
for iTarget = 1:4
    for iPaw = 1:3
        X = arrayfun(@(traj) traj.eta(iTarget, iPaw).X, traj, 'UniformOutput', false);
        N = arrayfun(@(traj) traj.eta(iTarget, iPaw).N, traj, 'UniformOutput', false);
        stats = arrayfun(@(traj) traj.eta(iTarget, iPaw).stats, traj, 'UniformOutput', false);
        trajCombined.eta(iTarget, iPaw).X = cat(1, X{:});
        trajCombined.eta(iTarget, iPaw).t = traj(1).eta(iTarget, iPaw).t;
        trajCombined.eta(iTarget, iPaw).N = cat(1, N{:});
        trajCombined.eta(iTarget, iPaw).stats = cat(2, stats{:});
        trajCombined.eta(iTarget, iPaw).target = targetNames{iTarget};
        trajCombined.eta(iTarget, iPaw).paw = pawNames{iPaw};
    end
end

%% Plot ETAs
fig = figure(Units='inches', Position=[0 0 7.5 10]);
for iTarget = 1:4
    for iPaw = 1:3
        ax = subplot(4, 3, 3*(iTarget-1)+iPaw);
        if iTarget == 1 && iPaw == 1
            [~, order] = EphysUnit.plotETA(ax, trajCombined.eta(iTarget, iPaw), event='reach onset', clim=[-1, 1], xlim=[-4, 0.5], ...
                sortWindow=[-2, 0.5], signWindow=[-0.2, 0.2], sortThreshold=0.25, negativeSortThreshold=0.25);
        else
            EphysUnit.plotETA(ax, trajCombined.eta(iTarget, iPaw), event='reach onset', clim=[-1, 1], xlim=[-4, 0.5], ...
                order=order);
        end
        ylabel(ax, targetNames(iTarget))
        title(ax, pawNames(iPaw))
    end
end

%% Calc METAs
trajCombined.meta = cell(4, 3);
for iTarget = 1:4
    for iPaw = 1:3
        eta = trajCombined.eta(iTarget, iPaw);
        X = eta.X;
        t = eta.t;
        trajCombined.meta{iTarget, iPaw} = mean(X(:, t>=-0.2 & t<= 0.2), 2, 'omitnan');
    end
end

%% Plot METAs
close all
skippedPairs = [13, 23, 41, 12, 22, 32, 42];
fig = figure(Units='inches', Position=[0   0   10.6771   10.4896]);
N = 12 - length(skippedPairs);
tl = tiledlayout(fig, N, N, TileSpacing='none', Padding=['none' ...
    '']);
AX = gobjects(1, N.^2);

iax = 0;
for ity = 1:4
    for ipy = 1:3
        for itx = 1:4
            for ipx = 1:3
                if any(ismember([itx*10+ipx, ity*10+ipy], skippedPairs))
                    continue
                end
                iax = iax + 1;
                ax = nexttile(tl);
                AX(iax) = ax;
                hold(ax, 'on')
                scatter(ax, trajCombined.meta{itx, ipx}, trajCombined.meta{ity, ipy}, 5, 'k');
                plot(ax, [-1, 3], [-1, 3], 'k:')
                plot(ax, [0, 0], [-1, 3], 'k:')
                plot(ax, [-1, 3], [0, 0], 'k:')
                hold(ax, 'off')
                if mod(iax, N) == 1
                    ylabel(ax, sprintf('%s\n(%s)', targetNames(ity), pawNames(ipy)))
                end
                if floor(iax./N) == N-1 || iax == N^2
                    xlabel(ax, sprintf('%s\n(%s)', targetNames(itx), pawNames(ipx)))
                end
            end
        end
    end
end
yticks(AX, [])
xticks(AX, [])
yticks(AX(1:N:N^2-N+1), [0, 2])
xticks(AX(N^2-N+1:end), [0, 2])

xlim(AX, [-1, 3])
ylim(AX, [-1, 3])
axis(AX, 'equal')