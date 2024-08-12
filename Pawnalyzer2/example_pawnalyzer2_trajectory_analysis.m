
eu = EphysUnit.load('C:\SERVER\Units\acute_3cam_reach_direction');

% pa = Pawnalyzer2(eu, refEvent='press');

%% Select good sessions
expNames = ["daisy23_20240626", "daisy23_20240627", "daisy24_20240626", "daisy25_20240701"];
selEu = ismember({eu.ExpName}, expNames);
eu = eu(selEu);
clearvars -except eu expNames

%% Make trajectories
pa = Pawnalyzer2(eu, refEvent='cue');
pa.getClips(nFramesBefore=15, nFramesAfter=0, keepData=false, trials='Press');
pa.load('\\research.files.med.harvard.edu\Neurobio\NEUROBIOLOGY SHARED\Assad Lab\Lingfeng\Data\Units\acute_3cam_reach_direction\Pawnalyzer2\pa.mat')

%% Analysis
%% PLOT TRAJECTORIES (ACROSS SESSIONS)
p.fontSize = 9;
% p.view = [-20, 25]; %%
p.view = [0, 90]; %%
close all
nExp = pa.getLength('exp');
clear traj
traj(nExp) = struct(contra=[], ipsi=[], target=[], t=[], pca=[], kmeans=[]);
nt = 10;
for iExp = 1:nExp
    nFrames = pa.getLength('frame', exp=iExp, trial=1);
    [traj(iExp).contra, traj(iExp).ipsi, traj(iExp).target, traj(iExp).t] = pa.getTrajectories(iExp, zero=nFrames - nt + 1);
    if ~ismember(pa.exp(iExp).animalName, {'desmond29', 'daisy25'})
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

%% Draw trajectories (aggregate across sessions)
ax = axes(figure(DefaultAxesFontSize=p.fontSize, Units='inches', Position=[0, 0, 3, 3]));
axis(ax, 'image');
hold(ax, 'on')
nTargets = 4;
nFrames = length(trajCombined.t);
targetNames = ["contra-out", "contra-front", "contra-in", "ipsi-front"];
h = gobjects(nTargets, 1);
for iTarget = 1:nTargets
    selFrames = nFrames - nt + 1:nFrames;
    sel = trajCombined.target == targetNames(iTarget);
    x = mean(trajCombined.contra.x(sel, selFrames), 1, 'omitnan');
    y = mean(trajCombined.contra.y(sel, selFrames), 1, 'omitnan');
    z = mean(trajCombined.contra.z(sel, selFrames), 1, 'omitnan');
    h(iTarget) = plot3(ax, x, y, z, LineWidth=1.5, Color=getColor(iTarget, 4, 0.8), DisplayName=targetNames(iTarget));
    scatter3(ax, x, y, z, 2*(selFrames-selFrames(1)+1).^1.6, getColor(iTarget, 4, 0.8), 'filled', Marker='o', DisplayName=targetNames(iTarget));
end
ax.ZAxis.Direction = 'reverse';
ax.YAxis.Direction = 'reverse';
ax.View = p.view;

xl = [-18, 5];
yl = [-5, 30];
zl = [-27, 2];
set(ax, XLim=xl, YLim=yl, ZLim=zl)
xrange = diff(ax.XLim);
yrange = diff(ax.YLim);
zrange = diff(ax.ZLim);

xticks(ax, ax.XLim + xrange*[0.125, -0.125])
yticks(ax, ax.YLim + yrange*[0.125, -0.125])
zticks(ax, ax.ZLim + zrange*[0.125, -0.125])
xticklabels(ax, ["med", "lat"])
yticklabels(ax, ["back", "front"])
zticklabels(ax, ["up", "down"])
set(ax, XMinorGrid='on', YMinorGrid='on', ZMinorGrid='on', Box='off')
ax.XAxis.MinorTickValues=ax.XLim(1) + xrange*[0.375, 0.625];
ax.YAxis.MinorTickValues=ax.YLim(1) + yrange*[0.375, 0.625];
ax.ZAxis.MinorTickValues=ax.ZLim(1) + zrange*[0.375, 0.625];

% Draw trajectories (ipsi)
ax = axes(figure(DefaultAxesFontSize=p.fontSize, Units='inches', Position=[0, 0, 3, 3]));
axis(ax, 'image');
hold(ax, 'on')
nTargets = 4;
nFrames = length(trajCombined.t);
targetNames = ["contra-out", "contra-front", "contra-in", "ipsi-front"];
h = gobjects(nTargets, 1);
for iTarget = 1:nTargets
    selFrames = nFrames - nt + 1:nFrames;
    sel = trajCombined.target == targetNames(iTarget);
    x = mean(trajCombined.ipsi.x(sel, selFrames), 1, 'omitnan');
    y = mean(trajCombined.ipsi.y(sel, selFrames), 1, 'omitnan');
    z = mean(trajCombined.ipsi.z(sel, selFrames), 1, 'omitnan');
    h(iTarget) = plot3(ax, x, y, z, LineWidth=1.5, Color=getColor(iTarget, 4, 0.8), DisplayName=targetNames(iTarget));
    scatter3(ax, x, y, z, 2*(selFrames-selFrames(1)+1).^1.6, getColor(iTarget, 4, 0.8), 'filled', Marker='o', DisplayName=targetNames(iTarget));
end
ax.ZAxis.Direction = 'reverse';
ax.YAxis.Direction = 'reverse';
ax.View = p.view;

set(ax, XLim=flip(-xl), YLim=yl, ZLim=zl)

xticks(ax, ax.XLim + xrange*[0.125, -0.125])
yticks(ax, ax.YLim + yrange*[0.125, -0.125])
zticks(ax, ax.ZLim + zrange*[0.125, -0.125])
ax.XAxis.MinorTickValues=ax.XLim(1) + xrange*[0.375, 0.625];
ax.YAxis.MinorTickValues=ax.YLim(1) + yrange*[0.375, 0.625];
ax.ZAxis.MinorTickValues=ax.ZLim(1) + zrange*[0.375, 0.625];
xticklabels(ax, ["med", "lat"])
yticklabels(ax, ["back", "front"])
zticklabels(ax, ["up", "down"])
set(ax, XMinorGrid='on', YMinorGrid='on', ZMinorGrid='on', Box='off')

%% Draw trajectories (each session)
fig = figure(DefaultAxesFontSize=p.fontSize, Units='inches', Position=[0, 0, 6, 3*nExp]);

for iExp = 1:nExp
    ax = subplot(nExp, 2, 2*(iExp - 1) + 1);
    axis(ax, 'image');
    hold(ax, 'on')
    nTargets = 4;
    nFrames = length(traj(iExp).t);
    targetNames = ["contra-out", "contra-front", "contra-in", "ipsi-front"];
    h = gobjects(nTargets, 1);
    for iTarget = 1:nTargets
        selFrames = nFrames - nt + 1:nFrames;
        sel = traj(iExp).target == targetNames(iTarget);
        x = mean(traj(iExp).contra.x(sel, selFrames), 1, 'omitnan');
        y = mean(traj(iExp).contra.y(sel, selFrames), 1, 'omitnan');
        z = mean(traj(iExp).contra.z(sel, selFrames), 1, 'omitnan');
        h(iTarget) = plot3(ax, x, y, z, LineWidth=1.5, Color=getColor(iTarget, 4, 0.8), DisplayName=targetNames(iTarget));
        scatter3(ax, x, y, z, 2*(selFrames-selFrames(1)+1).^1.6, getColor(iTarget, 4, 0.8), 'filled', Marker='o', DisplayName=targetNames(iTarget));
    end
    ax.ZAxis.Direction = 'reverse';
    ax.YAxis.Direction = 'reverse';
    ax.View = p.view;
    
    xlim(ax, 'auto')
    ylim(ax, 'auto')
    zlim(ax, 'auto')
    xl = ax.XLim; % [-18, 5];
    yl = ax.YLim; % [-5, 30];
    zl = ax.ZLim; % [-27, 2];
    set(ax, XLim=xl, YLim=yl, ZLim=zl)
    xrange = diff(ax.XLim);
    yrange = diff(ax.YLim);
    zrange = diff(ax.ZLim);
    
    xticks(ax, ax.XLim + xrange*[0.125, -0.125])
    yticks(ax, ax.YLim + yrange*[0.125, -0.125])
    zticks(ax, ax.ZLim + zrange*[0.125, -0.125])
    xticklabels(ax, ["med", "lat"])
    yticklabels(ax, ["back", "front"])
    zticklabels(ax, ["up", "down"])
    set(ax, XMinorGrid='on', YMinorGrid='on', ZMinorGrid='on', Box='off')
    ax.XAxis.MinorTickValues=ax.XLim(1) + xrange*[0.375, 0.625];
    ax.YAxis.MinorTickValues=ax.YLim(1) + yrange*[0.375, 0.625];
    ax.ZAxis.MinorTickValues=ax.ZLim(1) + zrange*[0.375, 0.625];
    
    % Draw trajectories (ipsi)
    ax = subplot(nExp, 2, 2*(iExp - 1) + 2);
    axis(ax, 'image');
    hold(ax, 'on')
    nTargets = 4;
    nFrames = length(traj(iExp).t);
    targetNames = ["contra-out", "contra-front", "contra-in", "ipsi-front"];
    h = gobjects(nTargets, 1);
    for iTarget = 1:nTargets
        selFrames = nFrames - nt + 1:nFrames;
        sel = traj(iExp).target == targetNames(iTarget);
        x = mean(traj(iExp).ipsi.x(sel, selFrames), 1, 'omitnan');
        y = mean(traj(iExp).ipsi.y(sel, selFrames), 1, 'omitnan');
        z = mean(traj(iExp).ipsi.z(sel, selFrames), 1, 'omitnan');
        h(iTarget) = plot3(ax, x, y, z, LineWidth=1.5, Color=getColor(iTarget, 4, 0.8), DisplayName=targetNames(iTarget));
        scatter3(ax, x, y, z, 2*(selFrames-selFrames(1)+1).^1.6, getColor(iTarget, 4, 0.8), 'filled', Marker='o', DisplayName=targetNames(iTarget));
    end
    ax.ZAxis.Direction = 'reverse';
    ax.YAxis.Direction = 'reverse';
    ax.View = p.view;

    set(ax, XLim=flip(-xl), YLim=yl, ZLim=zl)
    
    xticks(ax, ax.XLim + xrange*[0.125, -0.125])
    yticks(ax, ax.YLim + yrange*[0.125, -0.125])
    zticks(ax, ax.ZLim + zrange*[0.125, -0.125])
    ax.XAxis.MinorTickValues=ax.XLim(1) + xrange*[0.375, 0.625];
    ax.YAxis.MinorTickValues=ax.YLim(1) + yrange*[0.375, 0.625];
    ax.ZAxis.MinorTickValues=ax.ZLim(1) + zrange*[0.375, 0.625];
    xticklabels(ax, ["med", "lat"])
    yticklabels(ax, ["back", "front"])
    zticklabels(ax, ["up", "down"])
    set(ax, XMinorGrid='on', YMinorGrid='on', ZMinorGrid='on', Box='off')
end

%% Plot PETH for each sub-clustered trajectories
% Calculate ETA
targetNames = {'contra-out', 'contra-front', 'contra-in', 'ipsi-front'};
for iExp = 1:nExp
    nTrials = length(traj(iExp).target);
    trueStartTime = NaN(nTrials, 1);
    t = traj(iExp).t;
    for iTrial = 1:nTrials
        data = [ ...
            normalize(traj(iExp).contra.x(iTrial, :), 'zscore', 'robust'); ...
            normalize(traj(iExp).contra.y(iTrial, :), 'zscore', 'robust'); ...
            normalize(traj(iExp).contra.z(iTrial, :), 'zscore', 'robust'); ...
            normalize(traj(iExp).ipsi.x(iTrial, :), 'zscore', 'robust'); ...
            normalize(traj(iExp).ipsi.y(iTrial, :), 'zscore', 'robust'); ...
            normalize(traj(iExp).ipsi.z(iTrial, :), 'zscore', 'robust'); ...
            ];
        trueStartIndex = find(all(abs(data) < 2.5 | isnan(data), 1), 1, 'last');
        if ~isempty(trueStartIndex)
            trueStartTime(iTrial) = t(trueStartIndex) - t(end);
        end
    end
    [~, targetIndex] = ismember(traj(iExp).target, targetNames);
    targetIndex = targetIndex';
    for iTarget = 1:4
        sel = targetIndex==iTarget;
        traj(iExp).eta(iTarget) = pa.exp(iExp).eu.getETA('count', 'press', window=[-4, 0], resolution=0.1, normalize=[-4, -0.5], alignTo='stop', includeInvalid=true, ...
            trials=pa.exp(iExp).eu(1).Trials.Press(sel), correction=trueStartTime(sel));
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


%% Plot ETA
% close all
figure();
for iTarget = 1:4
    ax = subplot(1, 4, iTarget);
    if iTarget == 1
        [~, order] = EphysUnit.plotETA(ax, trajCombined.eta(iTarget), clim=[-1.5, 1.5], sortWindow=[-2, 0.5], signWindow=[-0.5, 0], xlim=[-2, 0]);
    else
        [~, order] = EphysUnit.plotETA(ax, trajCombined.eta(iTarget), clim=[-1.5, 1.5], order=order, xlim=[-2, 0]);
    end
    title(ax, trajCombined.eta(iTarget).target);
    if iTarget < 4
        delete(ax.Colorbar);
    else
        ax.Colorbar.Position = [0.916666666507011,0.10952380952381,0.009523809523809,0.814285714285715];
    end
end
