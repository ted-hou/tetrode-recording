 %% 
read_reachDir_4tgt;
read_reachDir_2tgt;

%% Lever-4-pos
p.fontSize = 9;
p.view = [0, 90];
p.trajDotMultiplier = 1.5;
p.trajDotPower = 1.5;

clear layout
layout.w = 7;
layout.h = 6;
layout.top.h = 2;
layout.bottom.h = 3;
layout.top.left.w = 3;
layout.top.right.w = 4;
layout.top2.h = 2;
layout.bottom2.h = 3;
layout.top2.left.w = 3;
layout.top2.right.w = 4;

close all
fig = figure(Units='inches', Position=[0, 0, layout.w, layout.h]);
layout.tl = tiledlayout(fig, layout.top.h + layout.bottom.h + layout.top2.h + layout.bottom2.h, 1);
layout.top.tl = tiledlayout(layout.tl, 1, layout.top.left.w + layout.top.right.w, TileSpacing='loose', Padding='loose');
l = layout.top.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [layout.top.h, 1];

layout.top.right.tl = tiledlayout(layout.top.tl, 1, 2);
l = layout.top.right.tl; l.Layout.Tile = 1 + layout.top.left.w; l.Layout.TileSpan = [1, layout.top.right.w];

layout.bottom.tl = tiledlayout(layout.tl, 1, 4, TileSpacing='compact');
l = layout.bottom.tl; l.Layout.Tile = 1 + layout.top.h; l.Layout.TileSpan = [layout.bottom.h, 1];

layout.top2.tl = tiledlayout(layout.tl, 1, layout.top2.left.w + layout.top2.right.w, TileSpacing='loose', Padding='loose');
l = layout.top2.tl; l.Layout.Tile = 1 + layout.top.h + layout.bottom.h; l.Layout.TileSpan = [layout.top2.h, 1];

layout.top2.right.tl = tiledlayout(layout.top2.tl, 1, 2);
l = layout.top2.right.tl; l.Layout.Tile = 1 + layout.top2.left.w; l.Layout.TileSpan = [1, layout.top2.right.w];

layout.bottom2.tl = tiledlayout(layout.tl, 1, 2, TileSpacing='compact');
l = layout.bottom2.tl; l.Layout.Tile = 1 + layout.top.h + layout.bottom.h + layout.top2.h; l.Layout.TileSpan = [layout.bottom2.h, 1];

% 7a. (top left, powerpoint) 4tgt mouse/target diagram
% 7b. (top-middle) 4tgt trajectories (contra, ipsi)
ax = nexttile(layout.top.right.tl);
title(ax, 'contra paw')
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
    plot3(ax, x, y, z, LineWidth=1.5, Color=getColor(iTarget, 4, 0.8), DisplayName=targetNames(iTarget));
    scatter3(ax, x, y, z, p.trajDotMultiplier*(selFrames-selFrames(1)+1).^p.trajDotPower, getColor(iTarget, 4, 0.8), Marker='o', DisplayName=targetNames(iTarget));
    h(iTarget) = plot(NaN, NaN, LineStyle='-', LineWidth=1.5, Marker='o', Color=getColor(iTarget, 4, 0.8), DisplayName=targetNames(iTarget));    
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
xticklabels(ax, ["lat", "med"])
yticklabels(ax, ["back", "front"])
zticklabels(ax, ["up", "down"])
set(ax, XMinorGrid='on', YMinorGrid='on', ZMinorGrid='on', Box='off')
ax.XAxis.MinorTickValues=ax.XLim(1) + xrange*[0.375, 0.625];
ax.YAxis.MinorTickValues=ax.YLim(1) + yrange*[0.375, 0.625];
ax.ZAxis.MinorTickValues=ax.ZLim(1) + zrange*[0.375, 0.625];
legend(h, Orientation='horizontal', Location='northoutside', Parent=layout.top.right.tl)

% ipsi
ax = nexttile(layout.top.right.tl);
title(ax, 'ipsi paw')
axis(ax, 'image');
hold(ax, 'on')
nTargets = 4;
nFrames = length(trajCombined.t);
targetNames = ["contra-out", "contra-front", "contra-in", "ipsi-front"];
for iTarget = 1:nTargets
    selFrames = nFrames - nt + 1:nFrames;
    sel = trajCombined.target == targetNames(iTarget);
    x = mean(trajCombined.ipsi.x(sel, selFrames), 1, 'omitnan');
    y = mean(trajCombined.ipsi.y(sel, selFrames), 1, 'omitnan');
    z = mean(trajCombined.ipsi.z(sel, selFrames), 1, 'omitnan');
    plot3(ax, x, y, z, LineWidth=1.5, Color=getColor(iTarget, 4, 0.8));
    scatter3(ax, x, y, z, p.trajDotMultiplier*(selFrames-selFrames(1)+1).^p.trajDotPower, getColor(iTarget, 4, 0.8), Marker='o');
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

% 7c. Plot population ETAs, 4 pos side by side
p.etaWindow = [-2, 0.5];
p.etaSortWindow = [-2, 0];
p.etaSignWindow = [-0.5, 0];
for iTarget = 1:4
    ax = nexttile(layout.bottom.tl);
    if iTarget == 1
        [~, order] = EphysUnit.plotETA(ax, trajCombined.eta(iTarget), event='reach onset', ...
            clim=[-1.5, 1.5], xlim=p.etaWindow, sortWindow=p.etaSortWindow, signWindow=p.etaSignWindow, ...
            sortThreshold=0.3, negativeSortThreshold=0.15);
        yticks(ax, unique([1, 50:50:size(trajCombined.eta(1).X, 1), size(trajCombined.eta(1).X, 1)]))
    else
        EphysUnit.plotETA(ax, trajCombined.eta(iTarget), event='reach onset', order=order, clim=[-1.5, 1.5], xlim=p.etaWindow);
        yticks(ax, []);
    end
    ylabel(ax, '');
    xlabel(ax, '');    
    title(ax, targetNames{iTarget})
    if iTarget < 4
        colorbar(ax, 'off')
    else
        ax.Colorbar.Layout.Tile = 'east';
    end
    fontsize(ax, p.fontSize, 'points')
    fontname(ax, 'Arial')
end
xlabel(layout.bottom.tl, 'Time to reach onset (s)', FontSize=p.fontSize, FontName='Arial')
ylabel(layout.bottom.tl, 'Unit', FontSize=p.fontSize, FontName='Arial')

% 7e. Trajectories, 2tgts
ax = nexttile(layout.top2.right.tl);
title(ax, 'contra paw')
axis(ax, 'image');
hold(ax, 'on')
nTargets = 2;
nFrames = length(trajCombined.t);
targetNames = ["contra-out", "contra-in"];
h = gobjects(nTargets, 1);
for iTarget = 1:nTargets
    switch iTarget
        case 1
            iColor = 1;
        case 2
            iColor = 3;
        otherwise
            iColor = iTarget;
    end
    selFrames = nFrames - nt + 1:nFrames;
    sel = trajCombined2tgt.target == targetNames(iTarget);
    x = mean(trajCombined2tgt.contra.x(sel, selFrames), 1, 'omitnan');
    y = mean(trajCombined2tgt.contra.y(sel, selFrames), 1, 'omitnan');
    z = mean(trajCombined2tgt.contra.z(sel, selFrames), 1, 'omitnan');
    plot3(ax, x, y, z, LineWidth=1.5, Color=getColor(iColor, 4, 0.8), DisplayName=targetNames(iTarget));
    scatter3(ax, x, y, z, p.trajDotMultiplier*(selFrames-selFrames(1)+1).^p.trajDotPower, getColor(iColor, 4, 0.8), Marker='o', DisplayName=targetNames(iTarget));
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
xticklabels(ax, ["lat", "med"])
yticklabels(ax, ["back", "front"])
zticklabels(ax, ["up", "down"])
set(ax, XMinorGrid='on', YMinorGrid='on', ZMinorGrid='on', Box='off')
ax.XAxis.MinorTickValues=ax.XLim(1) + xrange*[0.375, 0.625];
ax.YAxis.MinorTickValues=ax.YLim(1) + yrange*[0.375, 0.625];
ax.ZAxis.MinorTickValues=ax.ZLim(1) + zrange*[0.375, 0.625];

% ipsi
ax = nexttile(layout.top2.right.tl);
title(ax, 'ipsi paw')
axis(ax, 'image');
hold(ax, 'on')
nTargets = 2;
nFrames = length(trajCombined2tgt.t);
targetNames = ["contra-out", "contra-in"];
for iTarget = 1:nTargets
    switch iTarget
        case 1
            iColor = 1;
        case 2
            iColor = 3;
        otherwise
            iColor = iTarget;
    end
    selFrames = nFrames - nt + 1:nFrames;
    sel = trajCombined2tgt.target == targetNames(iTarget);
    x = mean(trajCombined2tgt.ipsi.x(sel, selFrames), 1, 'omitnan');
    y = mean(trajCombined2tgt.ipsi.y(sel, selFrames), 1, 'omitnan');
    z = mean(trajCombined2tgt.ipsi.z(sel, selFrames), 1, 'omitnan');
    plot3(ax, x, y, z, LineWidth=1.5, Color=getColor(iColor, 4, 0.8));
    scatter3(ax, x, y, z, p.trajDotMultiplier*(selFrames-selFrames(1)+1).^p.trajDotPower, getColor(iColor, 4, 0.8), Marker='o');
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

% 7f. ETA, 2tgt2
for iTarget = 1:2
    ax = nexttile(layout.bottom2.tl);
    if iTarget == 1
        [~, order] = EphysUnit.plotETA(ax, trajCombined2tgt.eta(iTarget), event='reach onset', ...
            clim=[-1.5, 1.5], xlim=p.etaWindow, sortWindow=p.etaSortWindow, signWindow=p.etaSignWindow, ...
            sortThreshold=0.3, negativeSortThreshold=0.15);
        yticks(ax, unique([1, 50:50:size(trajCombined2tgt.eta(1).X, 1), size(trajCombined2tgt.eta(1).X, 1)]))
    else
        EphysUnit.plotETA(ax, trajCombined2tgt.eta(iTarget), event='reach onset', order=order, clim=[-1.5, 1.5], xlim=p.etaWindow);
        yticks(ax, []);
    end
    ylabel(ax, '');
    xlabel(ax, '');    
    title(ax, targetNames{iTarget})
    if iTarget < 2
        colorbar(ax, 'off')
    else
        ax.Colorbar.Layout.Tile = 'east';
    end
    fontsize(ax, p.fontSize, 'points')
    fontname(ax, 'Arial')
end
xlabel(layout.bottom.tl, 'Time to reach onset (s)', FontSize=p.fontSize, FontName='Arial')
ylabel(layout.bottom.tl, 'Unit', FontSize=p.fontSize, FontName='Arial')


copygraphics(fig, ContentType='vector')

clear fig ax iTarget order

figure, histogram(cat(1, trueStartTime4tgts{:})), xlabel('Reach duration (s)'), ylabel('Num trials')
figure, histogram(cat(1, trueStartTime2tgts{:})), xlabel('Reach duration (s)'), ylabel('Num trials')

%% Supplement. Example trajectories from individual sessions

close all

fig = figure(Units='inches', Position=[2, 2, 10, p.firstRowHeight]);

for iExp = 1:length(expReachDir)
    ax = subplot(1, length(expReachDir), iExp);
    hold(ax, 'on')
    h = gobjects(5, 1);

    x0 = mean(trajectories(iExp).jaw.XAll, 'all', 'omitnan');
    y0 = mean(trajectories(iExp).jaw.YAll, 'all', 'omitnan');

    for iTarget = 1:4
        rawPos = posOrder(iExp, iTarget);
        switch expReachDir(iExp).animalName
            case 'desmond29'
                x = trajectories(iExp).handContra.Resampled.X(rawPos, :) - x0;
            case {'desmond28', 'desmond30'}
                x = -trajectories(iExp).handContra.Resampled.X(rawPos, :) + x0;
        end
        y = trajectories(iExp).handContra.Resampled.Y(rawPos, :) - y0;
        n = trajectories(iExp).handContra.Resampled.n(rawPos);
        plot(ax, x, y, ...
            Color=getColor(iTarget, 4, 0.8), LineStyle='-');
        h(iTarget) = plot(ax, x(end), y(end), Color=getColor(iTarget, 4, 0.8), Marker='.', MarkerSize=25, ...
            DisplayName=sprintf('%s (n=%i)', posNames{iTarget}, n));
    end

    h(5) = plot(ax, 0, 0, Marker='o', LineStyle='none', MarkerSize=5, Color='black', DisplayName='jaw');

    axis(ax, 'image');
    ax.YDir = 'reverse';
    hold(ax, 'off')
    h = h(:);
    legend(ax, h(:), Interpreter='none', Location='northoutside')
    ax.XLim(2) = 0;
    ax.YLim(1) = 0;    
    xticks(ax, ax.XLim)
    xticklabels(ax, {'out', 'in'})
    yticks(ax, ax.YLim)
    yticklabels(ax, {'up', 'down'})
    xlim(ax, ax.XLim + [-10, 10])
    ylim(ax, ax.YLim + [-10, 10])
    % title(sprintf('%s', exp(iExp).name), Interpreter='none')
    fontsize(ax, p.fontSize, 'points')
    fontname(ax, 'Arial')
end

%% Supplement, True start time of either paw

ax = axes(figure(Units='inches', Position=[0, 0, 3, 2]));
histogram(cat(1, tstReachDir{:}), Normalization='probability')
xlabel('True start time (s)')
ylabel('Probability')
fontname(ax, 'Arial')
fontsize(ax, p.fontSize, 'points')


%% Supplements Draw per-session trajectories (4tgts, each session)
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
