 %% 
read_reachDir_4tgt;
read_reachDir_2tgt;

%% Lever-2-pos
p.fontSize = 9;
p.view = [0, 90];
p.etaWindow = [-4, 0.5];
p.etaSortWindow = [-2, 0.5];
p.etaSignWindow = [-0.2, 0.2];
p.minNumTrials = 4;

hLetters = gobjects(1, 3);

clear layout
layout.w = 5;
layout.h = 6;
layout.top.h = 4;
layout.bottom.h = 8;
layout.top.left.w = 4;
layout.top.right.w = 2;

close all
fig = figure(Units='inches', Position=[0, 0, layout.w, layout.h]);
layout.tl = tiledlayout(fig, layout.top.h + layout.bottom.h, 1, TileSpacing='loose');
layout.top.tl = tiledlayout(layout.tl, 1, layout.top.left.w + layout.top.right.w, TileSpacing='loose');
l = layout.top.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [layout.top.h, 1];

% Topleft (traj 1x2)
layout.top.left.tl = tiledlayout(layout.top.tl, 1, 2, TileSpacing='compact', Padding='compact');
l = layout.top.left.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [1, layout.top.left.w];

% Topright (scatter 1x1)
layout.top.right.tl = tiledlayout(layout.top.tl, 1, 1, TileSpacing='compact');
l = layout.top.right.tl; l.Layout.Tile = 1 + layout.top.left.w; l.Layout.TileSpan = [1, layout.top.right.w];

% Bottom (ETA 1x2)
layout.bottom.tl = tiledlayout(layout.tl, 1, 8, TileSpacing='loose');
l = layout.bottom.tl; l.Layout.Tile = 1 + layout.top.h; l.Layout.TileSpan = [layout.bottom.h, 1];

layout.bottom.left.tl = tiledlayout(layout.bottom.tl, 1, 2, TileSpacing='compact');
l = layout.bottom.left.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [1, 7];

% 7a. Trajectories, 2tgts (a: good trials, b: trials where ipsi paw moved ~30%, based on LM separation, see fig generated in read_reachDir_2tgt)
PAWNAME = ["contra", "ipsi"];
SELTRIALS = {~trajCombined2tgt.usedIpsiPaw, ~trajCombined2tgt.usedIpsiPaw};
XL = {[-35, 20], [-20, 35]};
XTICKLABELS = {["L", "M"], ["M", "L"]};
TITLES = ["contra paw", "ipsi paw"];
DOTFACTOR = [1, 1];
DOTPOWER = [1.25, 1.25];
AX = gobjects(1, 2);
hDummy = gobjects(1, 2);
for iCol = 1:2
    ax = nexttile(layout.top.left.tl);
    AX(iCol) = ax;
    title(ax, TITLES(iCol))
    axis(ax, 'equal');
    hold(ax, 'on')
    nTargets = 2;
    nFrames = length(trajCombined.t);
    targetNames = ["contra-out", "contra-in"];
    targetNamesDisp = ["Lateral reach", "Medial reach"];
    for iTarget = 1:nTargets
        switch iTarget
            case 1
                iColor = 1;
            case 2
                iColor = 3;
        end
        selFrames = nFrames - nt + 1:nFrames;
        selTrials = trajCombined2tgt.target == targetNames(iTarget) & SELTRIALS{iCol};
        x = mean(trajCombined2tgt.(PAWNAME(iCol)).x(selTrials, selFrames), 1, 'omitnan');
        y = mean(trajCombined2tgt.(PAWNAME(iCol)).y(selTrials, selFrames), 1, 'omitnan');
        z = mean(trajCombined2tgt.(PAWNAME(iCol)).z(selTrials, selFrames), 1, 'omitnan');
        plot3(ax, x, y, z, LineWidth=1.5, Color=getColor(iColor, 4, 0.8), DisplayName=targetNamesDisp(iTarget));
        scatter3(ax, x, y, z, DOTFACTOR(iCol)*(selFrames-selFrames(1)+1).^DOTPOWER(iCol), getColor(iColor, 4, 0.8), Marker='o', DisplayName=targetNamesDisp(iTarget));
        hDummy(iTarget) = plot3(ax, NaN, NaN, NaN, '-o', LineWidth=1.5, Color=getColor(iColor, 4, 0.8), DisplayName=targetNamesDisp(iTarget));
    end
    ax.XAxis.Direction = 'normal';
    ax.YAxis.Direction = 'normal';
    ax.ZAxis.Direction = 'reverse';
    ax.View = p.view;
    
    set(ax, XLim=XL{iCol}, YLim=[-15, 35], ZLim=[-100, 100])
    xrange = diff(ax.XLim);
    yrange = diff(ax.YLim);
    zrange = diff(ax.ZLim);
    
    xticks(ax, ax.XLim + xrange*[0.125, -0.125])
    yticks(ax, ax.YLim + yrange*[0.125, -0.125])
    zticks(ax, ax.ZLim + zrange*[0.125, -0.125])
    xticklabels(ax, XTICKLABELS{iCol})
    yticklabels(ax, ["P", "A"])
    zticklabels(ax, ["D", "V"])
    set(ax, XMinorGrid='on', YMinorGrid='on', ZMinorGrid='on', Box='off')
    ax.XAxis.MinorTickValues=ax.XLim(1) + xrange*[0.375, 0.625];
    ax.YAxis.MinorTickValues=ax.YLim(1) + yrange*[0.375, 0.625];
    ax.ZAxis.MinorTickValues=ax.ZLim(1) + zrange*[0.375, 0.625];
end

ax = AX; clear AX;
yticklabels(ax(2), [])
% hLegend = legend(hDummy, Location='layout', Orientation='horizontal');
% hLegend.Layout.Tile = 'south';


hLetters(1) = text(ax(1), 0, 0, 'a', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax(1).Units = 'inches';
hLetters(1).HorizontalAlignment = 'right';
hLetters(1).VerticalAlignment = 'top';
hLetters(1).Position = [-0.2, ax(1).Position(4)+0.3, 0];


clear PAWNAME SELTRIALS XL XTICKLABELS DOTFACTOR DOTPOWER AX hDummy iRow iCol ax nTargets nFrames iTarget selFrames selTrials x y z hLegend


% 7b. Scatter META
% Select units with enough trials per condition
ETA = trajCombined2tgt.eta;
N = horzcat(ETA.N);
selUnits = all(N >= p.minNumTrials, 2);

metaWindow = [-0.2, 0.2];
ax = nexttile(layout.top.right.tl);
hold(ax, 'on')
ex = trajCombined2tgt.eta(1);
ey = trajCombined2tgt.eta(2);
metaX = mean(ex.X(selUnits, ex.t > metaWindow(1) & ex.t < metaWindow(2)), 2);
metaY = mean(ey.X(selUnits, ey.t > metaWindow(1) & ey.t < metaWindow(2)), 2);
scatter(ax, metaX, metaY, 5, 'k')
plot(ax, [-1, 3], [-1, 3], 'k:')
plot(ax, [-1, 3], [0, 0], 'k:')
plot(ax, [0, 0], [-1, 3], 'k:')
xlabel(ax, "Lateral (a.u.)")
ylabel(ax, "Medial (a.u.)")
title(ax, 'SNr response')
axis(ax, 'equal')
xlim(ax, [-1, 2])
ylim(ax, [-1, 2])

hLetters(2) = text(ax, 0, 0, 'b', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetters(2).HorizontalAlignment = 'right';
hLetters(2).VerticalAlignment = 'top';
hLetters(2).Position = [-0.2, ax.Position(4)+0.4, 0];

clear metaWindow ETAX ETAY iRow iCol ax ex ey metaX metaY

% 7c ETA HEATMAPS
TARGETNAME = ["contra-out", "contra-in"];
TARGETNAMEDISP = ["Lateral reach", "Medial reach"];
SELTRIALS = {~trajCombined2tgt.usedIpsiPaw, ~trajCombined2tgt.usedIpsiPaw};
FIELDNAME = ["eta", "eta"];
ITARGET = [1, 2];
AX = gobjects(1, 2);

targetNamesDisp = ["lat", "med"];
nUnits = nnz(selUnits);%size(trajCombined2tgt.eta(1).X, 1);
for iCol = 1:2
    ax = nexttile(layout.bottom.left.tl);
    AX(iCol) = ax;
    iTarget = ITARGET(iCol);
    if iCol == 1
        [~, order] = EphysUnit.plotETA(ax, trajCombined2tgt.(FIELDNAME(iCol))(iTarget), selUnits, event='reach onset', ...
            clim=[-2, 2], xlim=p.etaWindow, sortWindow=p.etaSortWindow, signWindow=p.etaSignWindow, ...
            sortThreshold=0.25, negativeSortThreshold=0.25);
        yticks(ax, unique([1, 50:50:nUnits, nUnits]))
    else
        EphysUnit.plotETA(ax, trajCombined2tgt.(FIELDNAME(iCol))(iTarget), selUnits, event='reach onset', order=order, clim=[-1, 1], xlim=p.etaWindow);
        yticks(ax, []);
    end
    hold(ax, 'on')
    plot(ax, [0, 0], [0, size(trajCombined2tgt.(FIELDNAME(iCol))(1).X, 1)], 'k')
    ylim(ax, [0, nUnits])
    xlim(ax, p.etaWindow)
    ylabel(ax, '');
    xlabel(ax, '');    
    title(ax, TARGETNAMEDISP(iCol))
    if iCol < 2
        colorbar(ax, 'off')
    else
        ax.Colorbar.Layout.Tile = 'east';
    end
    fontsize(ax, p.fontSize, 'points')
    fontname(ax, 'Arial')
    ax.YAxis.Direction = 'reverse';
end
ax = AX;
xlabel(layout.bottom.tl, 'Time from reach onset (s)', FontSize=p.fontSize, FontName='Arial')
ylabel(layout.bottom.tl, 'Unit', FontSize=p.fontSize, FontName='Arial')

hLetters(3) = text(ax(1), 0, 0, 'c', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax(1).Units = 'inches';
hLetters(3).HorizontalAlignment = 'right';
hLetters(3).VerticalAlignment = 'top';
hLetters(3).Position = [-0.2, ax(1).Position(4)+0.3, 0];

% Set fontsize
fontsize(fig, p.fontSize, 'points')
fontname(fig, 'Arial')

fontsize(hLetters, 16, 'points')
fontname(hLetters, 'Arial')

copygraphics(fig, ContentType='vector')
%% S7
boot_amplitude_difference;
%%
close all
% 7g. 4tgt trajectories (contra, ipsi)
DOTFACTOR = 1;
DOTPOWER = 1.25;

clear layout
layout.w = 4.5;
layout.h = 5;
layout.left.w = 2;
layout.right.w = 1;

fig = figure(Units='inches', Position=[1 1 layout.w, layout.h]);

layout.tl = tiledlayout(fig, 1, layout.left.w + layout.right.w);

layout.left.tl = tiledlayout(layout.tl, 3, 2);
l = layout.left.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [1, layout.left.w];

layout.right.tl = tiledlayout(layout.tl, 4, 1);
l = layout.right.tl; l.Layout.Tile = 1 + layout.left.w; l.Layout.TileSpan = [1, layout.right.w];

ax = nexttile(layout.left.tl);
title(ax, 'contra paw')
axis(ax, 'image');
hold(ax, 'on')
nTargets = 4;
nFrames = length(trajCombined.t);
targetNames = ["contra-out", "contra-front", "contra-in", "ipsi-front"];
pawNames = ["contra", "contra", "contra", "ipsi"];
nTrials = zeros(1, 4);
for iTarget = 1:nTargets
    selFrames = nFrames - nt + 1:nFrames;
    selTrials = trajCombined.target == targetNames(iTarget) & trajCombined.paw == pawNames(iTarget);
    nTrials(iTarget) = nnz(selTrials);
    x = mean(trajCombined.contra.x(selTrials, selFrames), 1, 'omitnan');
    y = mean(trajCombined.contra.y(selTrials, selFrames), 1, 'omitnan');
    z = mean(trajCombined.contra.z(selTrials, selFrames), 1, 'omitnan');
    plot3(ax, x, y, z, LineWidth=1.5, Color=getColor(iTarget, 4, 0.8), DisplayName=targetNames(iTarget));
    scatter3(ax, x, y, z, DOTFACTOR*(selFrames-selFrames(1)+1).^DOTPOWER, getColor(iTarget, 4, 0.8), Marker='o', DisplayName=targetNames(iTarget));
end
ax.ZAxis.Direction = 'reverse';
ax.View = p.view;

xl = [-25, 10];
yl = [-5, 40];
zl = [-100, 100];
set(ax, XLim=xl, YLim=yl, ZLim=zl)
xrange = diff(ax.XLim);
yrange = diff(ax.YLim);
zrange = diff(ax.ZLim);

xticks(ax, ax.XLim + xrange*[0.125, -0.125])
yticks(ax, ax.YLim + yrange*[0.125, -0.125])
zticks(ax, ax.ZLim + zrange*[0.125, -0.125])
xticklabels(ax, ["L", "M`"])
yticklabels(ax, ["P", "A"])
zticklabels(ax, ["D", "V"])
set(ax, XMinorGrid='on', YMinorGrid='on', ZMinorGrid='on', Box='off')
ax.XAxis.MinorTickValues=ax.XLim(1) + xrange*[0.375, 0.625];
ax.YAxis.MinorTickValues=ax.YLim(1) + yrange*[0.375, 0.625];
ax.ZAxis.MinorTickValues=ax.ZLim(1) + zrange*[0.375, 0.625];
fontsize(ax, p.fontSize, 'points')
% legend(h, Orientation='horizontal', Location='northoutside', Parent=layout.middle.right.tl)

% ipsi
ax = nexttile(layout.left.tl);
title(ax, 'ipsi paw')
axis(ax, 'image');
hold(ax, 'on')
nTargets = 4;
nFrames = length(trajCombined.t);
targetNames = ["contra-out", "contra-front", "contra-in", "ipsi-front"];
h = gobjects(4, 1);
for iTarget = 1:nTargets
    selFrames = nFrames - nt + 1:nFrames;
    selTrials = trajCombined.target == targetNames(iTarget);
    x = mean(trajCombined.ipsi.x(selTrials, selFrames), 1, 'omitnan');
    y = mean(trajCombined.ipsi.y(selTrials, selFrames), 1, 'omitnan');
    z = mean(trajCombined.ipsi.z(selTrials, selFrames), 1, 'omitnan');
    plot3(ax, x, y, z, LineWidth=1.5, Color=getColor(iTarget, 4, 0.8));
    scatter3(ax, x, y, z, DOTFACTOR*(selFrames-selFrames(1)+1).^DOTPOWER, getColor(iTarget, 4, 0.8), Marker='o');
    h(iTarget) = plot(NaN, NaN, LineStyle='-', LineWidth=1.5, Marker='o', Color=getColor(iTarget, 4, 0.8), DisplayName=targetNames(iTarget));    
end
ax.ZAxis.Direction = 'reverse';
ax.View = p.view;

set(ax, XLim=flip(-xl), YLim=yl, ZLim=zl)

xticks(ax, ax.XLim + xrange*[0.125, -0.125])
yticks(ax, ax.YLim + yrange*[0.125, -0.125])
zticks(ax, ax.ZLim + zrange*[0.125, -0.125])
ax.XAxis.MinorTickValues=ax.XLim(1) + xrange*[0.375, 0.625];
ax.YAxis.MinorTickValues=ax.YLim(1) + yrange*[0.375, 0.625];
ax.ZAxis.MinorTickValues=ax.ZLim(1) + zrange*[0.375, 0.625];
xticklabels(ax, ["M", "L"])
yticklabels(ax, ["P", "A"])
zticklabels(ax, ["D", "V"])
set(ax, XMinorGrid='on', YMinorGrid='on', ZMinorGrid='on', Box='off')
fontsize(ax, p.fontSize, 'points')

% 7c. Plot population ETAs, 4 pos side by side
N = arrayfun(@(eta) eta.N, trajCombined.eta, 'UniformOutput', false);
IPAW = [1, 1, 1, 3];
assert(minNumTrials == 4)
selUnit = N{1, 1} >= minNumTrials & N{2, 1} >= minNumTrials & N{3, 1} >= minNumTrials & N{4, 3} >= minNumTrials;
for iTarget = [2 4 1 3]
    ax = nexttile(layout.left.tl);
    if iTarget == 2 || iTarget == 1
        [~, order] = EphysUnit.plotETA(ax, trajCombined.eta(iTarget, IPAW(iTarget)), selUnit, event='reach onset', ...
            clim=[-2, 2], xlim=p.etaWindow, sortWindow=p.etaSortWindow, signWindow=p.etaSignWindow, ...
            sortThreshold=0.25, negativeSortThreshold=0.25);
        yt = 0:30:nnz(selUnit);
        yt(1) = 1;
        if round(yt(end)./30) == round(nnz(selUnit)./30)
            yt(end) = nnz(selUnit);
        else
            yt(end + 1) = nnz(selUnit);
        end
        yticks(ax, yt)
        ylabel(ax, 'Unit')
    else
        EphysUnit.plotETA(ax, trajCombined.eta(iTarget, IPAW(iTarget)), selUnit, event='reach onset', order=order, clim=[-2, 2], xlim=p.etaWindow);
        yticks(ax, []);
        ylabel(ax, '')
    end
    hold(ax, 'on')
    plot(ax, [0, 0], [0, nnz(selUnit)+1], 'k--');
    ylim(ax, [0, nnz(selUnit)+1])
    xlabel(ax, '');    
    title(ax, targetNames{iTarget})
    colorbar(ax, 'off')
%     if iTarget < 4
%         colorbar(ax, 'off')
%     else
%         ax.Colorbar.Layout.Tile = 'east';
%     end
    fontsize(ax, p.fontSize, 'points')
    fontname(ax, 'Arial')
end
xlabel(layout.left.tl, 'Time from reach onset (s)', FontSize=p.fontSize, FontName='Arial')
fontsize(ax, p.fontSize, 'points')

% Scatter
sz = 5;
ax = nexttile(layout.right.tl);
sel = c.hasPress & c.hasLick;
subselResp = sel & c.isSelective.pressVsLick;
subselNone = sel & ~c.isSelective.pressVsLick;
scatter(ax, meta.lick(subselNone), meta.press(subselNone), sz, 'black'), hold(ax, 'on')
scatter(ax, meta.lick(subselResp), meta.press(subselResp), sz, 'red')
plot(ax, [0 0], [-2 4], 'k:')
plot(ax, [-2 4], [0 0], 'k:')
plot(ax, [-2 4], [-2 4], 'k:')
axis(ax, 'equal')
xlim(ax, [-2, 4])
ylim(ax, [-2, 4])
xlabel(ax, 'Peri-lick')
ylabel(ax, 'Peri-reach')
fontsize(ax, p.fontSize, 'points')
fprintf('%i total, %i responsive.\n', nnz(subselNone) + nnz(subselResp), nnz(subselResp))

ax = nexttile(layout.right.tl);
assert(minNumTrials == 4)
minNumTrialsDisp = 4;
N = arrayfun(@(eta) eta.N, trajCombined.eta, 'UniformOutput', false);
sel = N{1, 1} >= minNumTrials & N{2, 1} >= minNumTrialsDisp & N{3, 1} >= minNumTrials & N{4, 3} >= minNumTrialsDisp;
subselResp = sel(:) & c.isSelective.contraFrontVsIpsiFront4tgt(:);
subselNone = sel(:) & ~c.isSelective.contraFrontVsIpsiFront4tgt(:);
ss = [N{2, 1}, N{4, 3}];
ss = max(ss, [], 2) ./ min(ss, [], 2);
ss = 10./ss;
scatter(ax, trajCombined.meta{2, 1}(subselNone), trajCombined.meta{4, 3}(subselNone), sz, 'black'), hold(ax, 'on')
scatter(ax, trajCombined.meta{2, 1}(subselResp), trajCombined.meta{4, 3}(subselResp), sz, 'red')
plot(ax, [0 0], [-2 4], 'k:')
plot(ax, [-2 4], [0 0], 'k:')
plot(ax, [-2 4], [-2 4], 'k:')
axis(ax, 'equal')
xlim(ax, [-2, 4])
ylim(ax, [-2, 4])
xlabel(ax, 'Contra-front')
ylabel(ax, 'Ipsi-front')
fontsize(ax, p.fontSize, 'points')
fprintf('%i total, %i responsive.\n', nnz(subselNone) + nnz(subselResp), nnz(subselResp))

ax = nexttile(layout.right.tl);
% assert(minNumTrials == 4)
minNumTrialsDisp = 4;
N = arrayfun(@(eta) eta.N, trajCombined.eta, 'UniformOutput', false);
sel = N{1, 1} >= minNumTrialsDisp & N{2, 1} >= minNumTrials & N{3, 1} >= minNumTrialsDisp & N{4, 3} >= minNumTrials;
subselResp = sel(:) & c.isSelective.contraOutVsContraIn4tgt(:);
subselNone = sel(:) & ~c.isSelective.contraOutVsContraIn4tgt(:);
ss = [N{1, 1}, N{3, 1}];
ss = max(ss, [], 2) ./ min(ss, [], 2);
ss = 10./ss;
scatter(ax, trajCombined.meta{1, 1}(subselNone), trajCombined.meta{3, 1}(subselNone), sz, 'black'), hold(ax, 'on')
scatter(ax, trajCombined.meta{1, 1}(subselResp), trajCombined.meta{3, 1}(subselResp), sz, 'red')
plot(ax, [0 0], [-2 4], 'k:')
plot(ax, [-2 4], [0 0], 'k:')
plot(ax, [-2 4], [-2 4], 'k:')
axis(ax, 'equal')
xlim(ax, [-2, 4])
ylim(ax, [-2, 4])
xlabel(ax, 'Contra-out')
ylabel(ax, 'Contra-in')
fontsize(ax, p.fontSize, 'points')
fprintf('%i total, %i responsive.\n', nnz(subselNone) + nnz(subselResp), nnz(subselResp))


ax = nexttile(layout.right.tl);
ETA = trajCombined2tgt.eta;
N = horzcat(ETA.N);
minNumTrialsDisp = 4;
sel = all(N >= minNumTrialsDisp, 2);
subselResp = sel(:) & c.isSelective.contraOutVsContraIn2tgt(:);
subselNone = sel(:) & ~c.isSelective.contraOutVsContraIn2tgt(:);
ex = trajCombined2tgt.eta(1);
ey = trajCombined2tgt.eta(2);
metaX = mean(ex.X(:, ex.t > -0.2 & ex.t < 0.1), 2);
metaY = mean(ey.X(:, ey.t > -0.2 & ey.t < 0.1), 2);
ss = N;
ss = max(ss, [], 2) ./ min(ss, [], 2);
ss = 10./ss;
scatter(ax, metaX(subselNone), metaY(subselNone), sz, 'black'), hold(ax, 'on')
scatter(ax, metaX(subselResp), metaY(subselResp), sz, 'red')
plot(ax, [0 0], [-2 4], 'k:')
plot(ax, [-2 4], [0 0], 'k:')
plot(ax, [-2 4], [-2 4], 'k:')
axis(ax, 'equal')
xlim(ax, [-2, 4])
ylim(ax, [-2, 4])
xlabel(ax, 'Lateral (2tgt)')
ylabel(ax, 'Medial (2tgt)')
fontsize(ax, p.fontSize, 'points')
fprintf('%i total, %i responsive.\n', nnz(subselNone) + nnz(subselResp), nnz(subselResp))

lgd = legend(h, Orientation='horizontal', NumColumns=2);
lgd.Layout.Tile = 'north';

copygraphics(fig, ContentType='vector')
