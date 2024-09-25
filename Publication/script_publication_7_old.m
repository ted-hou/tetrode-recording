 %% 
read_reachDir_4tgt;
read_reachDir_2tgt;

%% Lever-4-pos
p.fontSize = 9;
p.view = [0, 90];
p.etaWindow = [-2, 0.5];
p.etaSortWindow = [-2, 0.5];
p.etaSignWindow = [-0.2, 0.2];
p.minNumTrials = 2;

clear layout
layout.w = 7;
layout.h = 6;
layout.top.h = 4;
layout.bottom.h = 6;
layout.top.left.w = 2;
layout.top.right.w = 2;

close all
fig = figure(Units='inches', Position=[0, 0, layout.w, layout.h]);
layout.tl = tiledlayout(fig, layout.top.h + layout.bottom.h, 1, TileSpacing='loose');
layout.top.tl = tiledlayout(layout.tl, 1, layout.top.left.w + layout.top.right.w, TileSpacing='loose');
l = layout.top.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [layout.top.h, 1];

% Topleft (traj 2x2)
layout.top.left.tl = tiledlayout(layout.top.tl, 2, 2, TileSpacing='compact', Padding='compact');
l = layout.top.left.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [1, layout.top.left.w];

% Topright (scatter 2x2)
layout.top.right.tl = tiledlayout(layout.top.tl, 2, 2, TileSpacing='compact');
l = layout.top.right.tl; l.Layout.Tile = 1 + layout.top.left.w; l.Layout.TileSpan = [1, layout.top.right.w];

% Bottom (ETA 1x4)
layout.bottom.tl = tiledlayout(layout.tl, 1, 4, TileSpacing='loose');
l = layout.bottom.tl; l.Layout.Tile = 1 + layout.top.h; l.Layout.TileSpan = [layout.bottom.h, 1];

% 7a-b. Trajectories, 2tgts (a: good trials, b: trials where ipsi paw moved ~30%, based on LM separation, see fig generated in read_reachDir_2tgt)
PAWNAME = ["contra", "ipsi"];
SELTRIALS = {~trajCombined2tgt.usedIpsiPaw, ~trajCombined2tgt.usedIpsiPaw; trajCombined2tgt.usedIpsiPaw, trajCombined2tgt.usedIpsiPaw};
XL = {[-35, 20], [-20, 35]};
XTICKLABELS = {["L", "M"], ["M", "L"]};
DOTFACTOR = [1, 1];
DOTPOWER = [1.25, 1.25];
AX = gobjects(2, 2);
hDummy = gobjects(1, 2);
for iRow = 1:2
    for iCol = 1:2
        ax = nexttile(layout.top.left.tl);
        AX(iRow, iCol) = ax;
        title(ax, sprintf('%s paw', PAWNAME(iCol)))
        axis(ax, 'equal');
        hold(ax, 'on')
        nTargets = 2;
        nFrames = length(trajCombined.t);
        targetNames = ["contra-out", "contra-in"];
        targetNamesDisp = ["lat", "med"];
        for iTarget = 1:nTargets
            switch iTarget
                case 1
                    iColor = 1;
                case 2
                    iColor = 3;
            end
            selFrames = nFrames - nt + 1:nFrames;
            selTrials = trajCombined2tgt.target == targetNames(iTarget) & SELTRIALS{iRow, iCol};
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
end
ax = AX; clear AX;
xticklabels(ax(1, :), [])
yticklabels(ax(:, 2), [])
title(ax(2, :), [])
hLegend = legend(hDummy, Location='layout', Orientation='horizontal');
hLegend.Layout.Tile = 'north';

clear PAWNAME SELTRIALS XL XTICKLABELS DOTFACTOR DOTPOWER AX hDummy iRow iCol ax nTargets nFrames iTarget selFrames selTrials x y z hLegend


% 7c-f. Scatter META
% Select units with enough trials per condition
ETA = [trajCombined2tgt.eta, trajCombined2tgt.etaIpsiPaw];
N = horzcat(ETA.N);
selUnits = all(N >= p.minNumTrials, 2);

metaWindow = [-0.2, 0.2];
ETAX = {trajCombined2tgt.eta(1), trajCombined2tgt.etaIpsiPaw(1); trajCombined2tgt.eta(1), trajCombined2tgt.eta(2)};
ETAY = {trajCombined2tgt.eta(2), trajCombined2tgt.etaIpsiPaw(2); trajCombined2tgt.etaIpsiPaw(1), trajCombined2tgt.etaIpsiPaw(2)};
XNAME = ["lat", "lat w/ ipsi"; "lat", "med"];
YNAME = ["med", "med w/ ipsi"; "lat w/ ipsi", "med w/ ipsi"];

AX = gobjects(2, 2);
for iRow = 1:2
    for iCol = 1:2
        ax = nexttile(layout.top.right.tl);
        AX(iRow, iCol) = ax;
        hold(ax, 'on')
        ex = ETAX{iRow, iCol};
        ey = ETAY{iRow, iCol};
        metaX = mean(ex.X(selUnits, ex.t > metaWindow(1) & ex.t < metaWindow(2)), 2);
        metaY = mean(ey.X(selUnits, ey.t > metaWindow(1) & ey.t < metaWindow(2)), 2);
        scatter(ax, metaX, metaY, 5, 'k')
        plot(ax, [-1, 3], [-1, 3], 'k:')
        plot(ax, [-1, 3], [0, 0], 'k:')
        plot(ax, [0, 0], [-1, 3], 'k:')
        xlabel(ax, XNAME(iRow, iCol))
        ylabel(ax, YNAME(iRow, iCol))
    end
end
ax = AX;
axis(ax, 'equal')
xlim(ax, [-1, 2])
ylim(ax, [-1, 2])
xticks(ax(1, :), [])
yticks(ax(:, 2), [])
clear metaWindow ETAX ETAY iRow iCol ax ex ey metaX metaY

% 7g-h ETA HEATMAPS
TARGETNAME = ["contra-out", "contra-in", "contra-out", "contra-in"];
TARGETNAMEDISP = ["lat", "med", "lat w/ ipsi", "med w/ ipsi"];
SELTRIALS = {~trajCombined2tgt.usedIpsiPaw, ~trajCombined2tgt.usedIpsiPaw, trajCombined2tgt.usedIpsiPaw, trajCombined2tgt.usedIpsiPaw};
FIELDNAME = ["eta", "eta", "etaIpsiPaw", "etaIpsiPaw"];
ITARGET = [1, 2, 1, 2];

targetNamesDisp = ["lat", "med"];
nUnits = nnz(selUnits);%size(trajCombined2tgt.eta(1).X, 1);
for iCol = 1:4
    ax = nexttile(layout.bottom.tl);
    iTarget = ITARGET(iCol);
    if iCol == 1
        [~, order] = EphysUnit.plotETA(ax, trajCombined2tgt.(FIELDNAME(iCol))(iTarget), selUnits, event='reach onset', ...
            clim=[-2, 2], xlim=p.etaWindow, sortWindow=p.etaSortWindow, signWindow=p.etaSignWindow, ...
            sortThreshold=0.25, negativeSortThreshold=0.25);
        yticks(ax, unique([1, 50:50:nUnits, nUnits]))
    else
        EphysUnit.plotETA(ax, trajCombined2tgt.(FIELDNAME(iCol))(iTarget), selUnits, event='reach onset', order=order, clim=[-2, 2], xlim=p.etaWindow);
        yticks(ax, []);
    end
    hold(ax, 'on')
    plot(ax, [0, 0], [0, size(trajCombined2tgt.(FIELDNAME(iCol))(1).X, 1)], 'k')
    ylim(ax, [0, nUnits])
    xlim(ax, [-2, 0.5])
    ylabel(ax, '');
    xlabel(ax, '');    
    title(ax, TARGETNAMEDISP(iCol))
    if iCol < 4
        colorbar(ax, 'off')
    else
        ax.Colorbar.Layout.Tile = 'east';
    end
    fontsize(ax, p.fontSize, 'points')
    fontname(ax, 'Arial')
    ax.YAxis.Direction = 'reverse';
end
xlabel(layout.bottom.tl, 'Time from reach onset (s)', FontSize=p.fontSize, FontName='Arial')
ylabel(layout.bottom.tl, 'Unit', FontSize=p.fontSize, FontName='Arial')

%%
ax = nexttile(layout.top.tl, 1 + layout.top.left.w + layout.top.middle.w, [1, layout.top.right.w]);

hold(ax, 'on')
tst = cat(1, trueStartTime2tgts{:});
histogram(ax, tst(trajCombined2tgt.target == 'contra-out' & ~trajCombined2tgt.usedIpsiPaw), -0.5:0.05:0, DisplayName='lat', Normalization='probability')
histogram(ax, tst(trajCombined2tgt.target == 'contra-in' & ~trajCombined2tgt.usedIpsiPaw), -0.5:0.05:0, DisplayName='med', Normalization='probability')
xlabel('Reach onset latency (s)'), ylabel('Probability')
legend(ax, Location='northwest')
hold(ax, 'off')

copygraphics(fig, ContentType='vector')

%%
% 7g. 4tgt trajectories (contra, ipsi)
ax = nexttile(layout.middle.right.tl);
title(ax, 'contra paw')
axis(ax, 'image');
hold(ax, 'on')
nTargets = 4;
nFrames = length(trajCombined.t);
targetNames = ["contra-out", "contra-front", "contra-in", "ipsi-front"];
for iTarget = 1:nTargets
    selFrames = nFrames - nt + 1:nFrames;
    selTrials = trajCombined.target == targetNames(iTarget);
    x = mean(trajCombined.contra.x(selTrials, selFrames), 1, 'omitnan');
    y = mean(trajCombined.contra.y(selTrials, selFrames), 1, 'omitnan');
    z = mean(trajCombined.contra.z(selTrials, selFrames), 1, 'omitnan');
    plot3(ax, x, y, z, LineWidth=1.5, Color=getColor(iTarget, 4, 0.8), DisplayName=targetNames(iTarget));
    scatter3(ax, x, y, z, p.trajDotMultiplier*(selFrames-selFrames(1)+1).^p.trajDotPower, getColor(iTarget, 4, 0.8), Marker='o', DisplayName=targetNames(iTarget));
end
ax.ZAxis.Direction = 'reverse';
ax.YAxis.Direction = 'reverse';
ax.View = p.view;

xl = [-25, 10];
yl = [-5, 30];
zl = [-100, 100];
set(ax, XLim=xl, YLim=yl, ZLim=zl)
xrange = diff(ax.XLim);
yrange = diff(ax.YLim);
zrange = diff(ax.ZLim);

xticks(ax, ax.XLim + xrange*[0.125, -0.125])
yticks(ax, ax.YLim + yrange*[0.125, -0.125])
zticks(ax, ax.ZLim + zrange*[0.125, -0.125])
xticklabels(ax, ["lat", "lat"])
yticklabels(ax, ["back", "front"])
zticklabels(ax, ["up", "down"])
set(ax, XMinorGrid='on', YMinorGrid='on', ZMinorGrid='on', Box='off')
ax.XAxis.MinorTickValues=ax.XLim(1) + xrange*[0.375, 0.625];
ax.YAxis.MinorTickValues=ax.YLim(1) + yrange*[0.375, 0.625];
ax.ZAxis.MinorTickValues=ax.ZLim(1) + zrange*[0.375, 0.625];
% legend(h, Orientation='horizontal', Location='northoutside', Parent=layout.middle.right.tl)

% ipsi
ax = nexttile(layout.middle.right.tl);
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
    scatter3(ax, x, y, z, p.trajDotMultiplier*(selFrames-selFrames(1)+1).^p.trajDotPower, getColor(iTarget, 4, 0.8), Marker='o');
    h(iTarget) = plot(NaN, NaN, LineStyle='-', LineWidth=1.5, Marker='o', Color=getColor(iTarget, 4, 0.8), DisplayName=targetNames(iTarget));    
end
hLegend = legend(h, NumColumns=2, Location='southoutside');
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
    colorbar(ax, 'off')
%     if iTarget < 4
%         colorbar(ax, 'off')
%     else
%         ax.Colorbar.Layout.Tile = 'east';
%     end
    fontsize(ax, p.fontSize, 'points')
    fontname(ax, 'Arial')
end
xlabel(layout.bottom.tl, 'Time from reach onset (s)', FontSize=p.fontSize, FontName='Arial')
ylabel(layout.bottom.tl, 'Unit', FontSize=p.fontSize, FontName='Arial')


% 7c. Histogram of movement latency
ax = nexttile(layout.top.tl, 1 + layout.top.left.w + layout.top.middle.w, [1, layout.top.right.w]);

hold(ax, 'on')
histogram(ax, cat(1, trueStartTime2tgts{:}), DisplayName='2tgt', Normalization='probability')
histogram(ax, cat(1, trueStartTime4tgts{:}), DisplayName='4tgt', Normalization='probability')
xlabel('Reach onset latency (s)'), ylabel('Probability')
legend(ax, Location='northwest')
hold(ax, 'off')


copygraphics(fig, ContentType='vector')