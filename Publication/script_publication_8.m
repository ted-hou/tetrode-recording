 %% 
% read_reachDir_4tgt;
% read_reachDir_2tgt;
% euReachDir2tgt = EphysUnit.load('C:\SERVER\Units\acute_3cam_reach_direction_2tgts\SingleUnits_NonDuplicate', waveforms=false, spikecounts=false, spikerates=false);
% euReachDir4tgt = EphysUnit.load('C:\SERVER\Units\acute_3cam_reach_direction\SingleUnits_NonDuplicate', waveforms=false, spikecounts=false, spikerates=false);

load('C:\SERVER\Units\traj_reachDir_2tgt.mat')
load('C:\SERVER\Units\traj_reachDir_4tgt.mat')
load('C:\SERVER\Units\meta_Lite_NonDuplicate_NonDrift.mat')
load('C:\SERVER\Units\boot_20241118_Figure7.mat')
read_lda
nt = 16;
% boot_amplitude_difference;
%% Fig7. Lever-2-pos
p.fontSize = 9;
p.view = [0, 90];
p.etaWindow = [-3.5, 0.5];
p.etaSortWindow = [-2.5, 0.5];
p.etaSignWindow = [-0.1, 0.2];
p.minNumTrials = 4;

hLetters = gobjects(1, 4);

clear layout
layout.w = 5;
layout.h = 6;
layout.top.h = 4;
layout.bottom.h = 4;
layout.top.left.w = 4;
layout.top.right.w = 2;
layout.top.right.top.h = 5;
layout.top.right.bottom.h = 2;

close all
fig = figure(Units='inches', Position=[0, 0, layout.w, layout.h]);
layout.tl = tiledlayout(fig, layout.top.h + layout.bottom.h, 1, TileSpacing='loose');
layout.top.tl = tiledlayout(layout.tl, 1, layout.top.left.w + layout.top.right.w, TileSpacing='loose');
l = layout.top.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [layout.top.h, 1];

% Topleft (traj 1x2)
layout.top.left.tl = tiledlayout(layout.top.tl, 1, 2, TileSpacing='compact', Padding='compact');
l = layout.top.left.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [1, layout.top.left.w];

% Topright (scatter top, decoder bottom)
layout.top.right.tl = tiledlayout(layout.top.tl, layout.top.right.top.h + layout.top.right.bottom.h, 1, TileSpacing='compact');
l = layout.top.right.tl; l.Layout.Tile = 1 + layout.top.left.w; l.Layout.TileSpan = [1, layout.top.right.w];

% Bottom (ETA 1x2)
layout.bottom.tl = tiledlayout(layout.tl, 1, 8, TileSpacing='loose');
l = layout.bottom.tl; l.Layout.Tile = 1 + layout.top.h; l.Layout.TileSpan = [layout.bottom.h, 1];

layout.bottom.left.tl = tiledlayout(layout.bottom.tl, 1, 3, TileSpacing='compact');
l = layout.bottom.left.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [1, 7];

% 8a. Trajectories, 2tgts (a: good trials, b: trials where ipsi paw moved ~30%, based on LM separation, see fig generated in read_reachDir_2tgt)
PAWNAME = ["contra", "ipsi"];
SELTRIALS = {~trajCombined2tgt.usedIpsiPaw, ~trajCombined2tgt.usedIpsiPaw};
XL = {[-35, 20], [-20, 35]};
XTICKLABELS = {["L", "M"], ["M", "L"]};
TITLES = ["Contra paw", "Ipsi paw"];
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
    nFrames = length(trajCombined4tgt.t);
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

axl = AX(1);
ax = AX; clear AX;
yticklabels(ax(2), [])
% hLegend = legend(hDummy, Location='layout', Orientation='horizontal');
% hLegend.Layout.Tile = 'south';


hLetters(1) = text(ax(1), 0, 0, 'a', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax(1).Units = 'inches';
hLetters(1).HorizontalAlignment = 'right';
hLetters(1).VerticalAlignment = 'top';
hLetters(1).Position = [-0.2, ax(1).Position(4)+0.3, 0];


clear PAWNAME SELTRIALS XL XTICKLABELS DOTFACTOR DOTPOWER AX iRow iCol ax nTargets nFrames iTarget selFrames selTrials x y z hLegend


% 8b. Scatter META
% Select units with enough trials per condition
ETA = trajCombined2tgt.eta;
N = horzcat(ETA.N);
selUnits = all(N >= p.minNumTrials, 2);

metaWindow = [-0.1, 0.2];
ax = nexttile(layout.top.right.tl, 1, [layout.top.right.top.h, 1]);
hold(ax, 'on')
ex = trajCombined2tgt.eta(1);
ey = trajCombined2tgt.eta(2);
metaX = mean(ex.X(:, ex.t > metaWindow(1) & ex.t < metaWindow(2)), 2);
metaY = mean(ey.X(:, ey.t > metaWindow(1) & ey.t < metaWindow(2)), 2);
scatter(ax, metaX(selUnits), metaY(selUnits), 5, 'k', MarkerEdgeAlpha=0.5)
scatter(ax, metaX(selUnits & any(cat(2, c.isPressResponsive2tgt{:}), 2)), metaY(selUnits & any(cat(2, c.isPressResponsive2tgt{:}), 2)), 5, 'k', 'filled')
% scatter(ax, metaX(selUnits & c.isSelective.contraOutVsContraIn2tgt), metaY(selUnits & c.isSelective.contraOutVsContraIn2tgt), 5, 'r', 'filled')
size(selUnits & all(cat(2, c.isPressResponsive2tgt{:}), 2))
plot(ax, [-1, 3], [-1, 3], 'k:')
plot(ax, [-1, 3], [0, 0], 'k:')
plot(ax, [0, 0], [-1, 3], 'k:')
xlabel(ax, "Lateral (a.u.)")
ylabel(ax, "Medial (a.u.)")
title(ax, 'SNr response')
axis(ax, 'equal')
xlim(ax, [-1, 2])
ylim(ax, [-1, 2])

hLetters(2) = text(ax, 0, 0, 'c', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetters(2).HorizontalAlignment = 'right';
hLetters(2).VerticalAlignment = 'top';
hLetters(2).Position = [-0.2, ax.Position(4)+0.4, 0];

clear metaWindow ETAX ETAY iRow iCol ax ex ey

% 8c Decoder (LDA) for reach2tgt
t = likelihood2tgt(1).t;
latTrialLat = arrayfun(@(llh) llh.contraOut(llh.trueLabel=="lateral", :), likelihood2tgt, UniformOutput=false);
latTrialMed = arrayfun(@(llh) llh.contraIn(llh.trueLabel=="lateral", :), likelihood2tgt, UniformOutput=false);
medTrialLat = arrayfun(@(llh) llh.contraOut(llh.trueLabel=="medial", :), likelihood2tgt, UniformOutput=false);
medTrialMed = arrayfun(@(llh) llh.contraIn(llh.trueLabel=="medial", :), likelihood2tgt, UniformOutput=false);
latTrialLat = cat(1, latTrialLat{:});
latTrialMed = cat(1, latTrialMed{:});
medTrialLat = cat(1, medTrialLat{:});
medTrialMed = cat(1, medTrialMed{:});

DATA = { ...
    latTrialLat, latTrialMed; ...
    medTrialLat, medTrialMed ...
    };
COLOR = {getColor(1, 4, 0.8), getColor(3, 4, 0.8)};
LABEL = ["lateral", "medial"];

ax = nexttile(layout.top.right.tl, 1 + layout.top.right.top.h, [layout.top.right.bottom.h, 1]);
hold(ax, 'on')
h = gobjects(3, 1);
for iMove = 1:2
    mu = mean(DATA{iMove, 1} - DATA{iMove, 2}, 1, 'omitnan');
    h(iMove) = plot(ax, t, mu, Color=COLOR{iMove}, LineWidth=1.5, DisplayName=sprintf('%s', LABEL(iMove)));
end
% h(3) = plot(ax, t, dfBootStats2tgt.all.mu, 'black', LineStyle='--', LineWidth=1.5, DisplayName='shuffle');
patch(ax, [t, flip(t)], [dfBootStats2tgt.all.ci(1, :), flip(dfBootStats2tgt.all.ci(2, :))], 'black', FaceAlpha=0.1, EdgeAlpha=0.5, DisplayName='99% CI');
patch(ax, [pLDA.responseWindow2tgt, flip(pLDA.responseWindow2tgt)], [-1, -1, 1, 1], 'yellow', FaceAlpha=0.1, EdgeAlpha=0.5, DisplayName='training');
ylim(ax, [-0.2, 0.2])
xline(ax, 0, ':')
yline(ax, 0, ':')
xticks(ax, -4:2:2)
xlim(ax, [-4, 2])
hold(ax, 'off')

xlabel(ax, 'Time to reach onset (s)')
ylabel(ax, 'p(lat) - p(med)')

hLetters(3) = text(ax, 0, 0, 'd', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetters(3).HorizontalAlignment = 'right';
hLetters(3).VerticalAlignment = 'top';
hLetters(3).Position = [-0.2, ax.Position(4)+0.4, 0];

% lgd = legend(ax, h(1:2), Orientation='horizontal', Location='northoutside', FontSize=p.fontSize-1, AutoUpdate=false);

selUnits2tgt = ismember(string({euReachDir2tgt.ExpName}), expNames2tgt(goodExpIndices2tgt));
nAnimals2tgt = length(unique(euReachDir2tgt(selUnits2tgt).getAnimalName()));

fprintf('\nReach (2tgt) lateral vs. medial decoder:\n');
fprintf('\t1) Selected %i sessions (%i animals) with >=%i units (mean=%g, total=%i).\n', length(goodExpIndices2tgt), nAnimals2tgt, pLDA.minNumUnits, mean(nUnits2tgt), sum(nUnits2tgt))
fprintf('\t\t Units per session (sorted): %s\n', num2str(sort(nUnits2tgt)));
fprintf('\t2) %i lateral trials, %i medial trials on aggregate.\n', size(latTrialLat, 1), size(medTrialLat, 1))
fprintf('\t3) Shuffled trial labels %i times to retrain LDA, shuffled mean and 99%%CI is calculated by averaging across all %i trials for each shuffle, then averaging across all shuffles.\n', pLDA.nBoot, size(latTrialLat, 1)+size(medTrialLat, 1))
fprintf('\t4) Training window was chosen at [%g, %g] s.\n', pLDA.responseWindow2tgt(1), pLDA.responseWindow2tgt(2))

% 8d ETA HEATMAPS
TARGETNAME = ["contra-out", "contra-in", "contra-in"];
TARGETNAMEDISP = {{'Lateral reach', 'Sort A'}, {'Medial reach', 'Sort A'}, {'Medial reach', 'Sort B'}};
SELTRIALS = {~trajCombined2tgt.usedIpsiPaw, ~trajCombined2tgt.usedIpsiPaw, ~trajCombined2tgt.usedIpsiPaw};
FIELDNAME = ["eta", "eta", "eta"];
ITARGET = [1, 2, 2];
AX = gobjects(1, 3);

targetNamesDisp = ["lat", "med", "med"];
nUnits = nnz(selUnits);%size(trajCombined2tgt.eta(1).X, 1);

% groupVar = NaN(nUnits, 1);
% groupVar(metaX(selUnits) < 0 & metaY(selUnits) > 0) = 0;
% groupVar(metaX(selUnits) > 0 & metaY(selUnits) < 0) = 1;
% groupVar(metaX(selUnits) < 0 & metaY(selUnits) < 0) = 2;
% groupVar(metaX(selUnits) > 0 & metaY(selUnits) > 0) = 3;
CL = [-1.5, 1.5];

for iCol = 1:3
    ax = nexttile(layout.bottom.left.tl);
    AX(iCol) = ax;
    iTarget = ITARGET(iCol);
    if iCol == 1
        [~, order] = EphysUnit.plotETA(ax, trajCombined2tgt.(FIELDNAME(iCol))(iTarget), selUnits, event='reach onset', ...
            clim=CL, sortWindow=p.etaSortWindow, signWindow=p.etaSignWindow, ...
            sortThreshold=0.25, negativeSortThreshold=0.25);
        yticks(ax, unique([1, 50:50:nUnits, nUnits]))
    elseif iCol == 2
        EphysUnit.plotETA(ax, trajCombined2tgt.(FIELDNAME(iCol))(iTarget), selUnits, event='reach onset', order=order, clim=CL);
        yticks(ax, []);
    elseif iCol == 3
        EphysUnit.plotETA(ax, trajCombined2tgt.(FIELDNAME(iCol))(iTarget), selUnits, event='reach onset', ...
            clim=CL, sortWindow=p.etaSortWindow, signWindow=p.etaSignWindow, ...
            sortThreshold=0.25, negativeSortThreshold=0.25);
        yticks(ax, unique([1, 50:50:nUnits, nUnits]))
    end
    hold(ax, 'on')
    xline(ax, 0, 'k--')
%     yline(ax, nnz(groupVar <= 1), 'k--')
    ylim(ax, [0, nUnits])
    xlim(ax, [-2.5, 0.5])
    ylabel(ax, '');
    xlabel(ax, '');    
    title(ax, TARGETNAMEDISP{iCol})
    if iCol ~= 2
        colorbar(ax, 'off')
    else
        ax.Colorbar.Layout.Tile = 'east';
    end
    fontsize(ax, p.fontSize, 'points')
    fontname(ax, 'Arial')
    ax.YAxis.Direction = 'reverse';
    applyCustomColormap(ax, [-1.5, 1.5], hlim=[0.375, 0, 0, -0.375], llim=[0.125, 0.5, 0.5, 0.25], hpwr=.5, lpwr=1, h0=0.33);
end

ax = AX;
xlabel(layout.bottom.tl, 'Time from reach onset (s)', FontSize=p.fontSize, FontName='Arial')
ylabel(layout.bottom.tl, 'Unit', FontSize=p.fontSize, FontName='Arial')

hLetters(4) = text(ax(1), 0, 0, 'b', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax(1).Units = 'inches';
hLetters(4).HorizontalAlignment = 'right';
hLetters(4).VerticalAlignment = 'top';
hLetters(4).Position = [-0.2, ax(1).Position(4)+0.3, 0];


% Set fontsize
fontsize(fig, p.fontSize, 'points')
fontname(fig, 'Arial')

fontsize(hLetters, 16, 'points')
fontname(hLetters, 'Arial')

lgd = legend(hDummy, Orientation='horizontal');
lgd.Layout.Tile = 'south';

copygraphics(fig, ContentType='vector', BackgroundColor='none')

%% Fig S8. Lever-4-pos and scatter META comparisons for any A vs. B movement
close all
% 8g. 4tgt trajectories (contra, ipsi)
DOTFACTOR = 1;
DOTPOWER = 1.25;

clear layout
layout.w = 3.9;
layout.h = 5;
layout.left.w = 2;
layout.right.w = 1;

fig = figure(Units='inches', Position=[1 1 layout.w, layout.h]);

layout.tl = tiledlayout(fig, 1, layout.left.w + layout.right.w);

layout.left.tl = tiledlayout(layout.tl, 3, 2);
l = layout.left.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [1, layout.left.w];

layout.right.tl = tiledlayout(layout.tl, 2, 1);
l = layout.right.tl; l.Layout.Tile = 1 + layout.left.w; l.Layout.TileSpan = [1, layout.right.w];

ax = nexttile(layout.left.tl);
title(ax, 'Contra paw')
axis(ax, 'image');
hold(ax, 'on')
nTargets = 4;
nFrames = length(trajCombined4tgt.t);
targetNames = ["contra-out", "contra-front", "contra-in", "ipsi-front"];
pawNames = ["contra", "contra", "contra", "ipsi"];
nTrials = zeros(1, 4);
for iTarget = 1:nTargets
    selFrames = nFrames - nt + 1:nFrames;
    selTrials = trajCombined4tgt.target == targetNames(iTarget) & trajCombined4tgt.paw == pawNames(iTarget);
    nTrials(iTarget) = nnz(selTrials);
    x = mean(trajCombined4tgt.contra.x(selTrials, selFrames), 1, 'omitnan');
    y = mean(trajCombined4tgt.contra.y(selTrials, selFrames), 1, 'omitnan');
    z = mean(trajCombined4tgt.contra.z(selTrials, selFrames), 1, 'omitnan');
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
title(ax, 'Ipsi paw')
axis(ax, 'image');
hold(ax, 'on')
nTargets = 4;
nFrames = length(trajCombined4tgt.t);
targetNames = ["contra-out", "contra-front", "contra-in", "ipsi-front"];
h = gobjects(4, 1);
for iTarget = 1:nTargets
    selFrames = nFrames - nt + 1:nFrames;
    selTrials = trajCombined4tgt.target == targetNames(iTarget);
    x = mean(trajCombined4tgt.ipsi.x(selTrials, selFrames), 1, 'omitnan');
    y = mean(trajCombined4tgt.ipsi.y(selTrials, selFrames), 1, 'omitnan');
    z = mean(trajCombined4tgt.ipsi.z(selTrials, selFrames), 1, 'omitnan');
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

% 8c. Plot population ETAs, 4 pos side by side
N = arrayfun(@(eta) eta.N, trajCombined4tgt.eta, 'UniformOutput', false);
IPAW = [1, 1, 1, 3];
minNumTrials = 4;
assert(minNumTrials == 4);
selUnit = N{1, 1} >= minNumTrials & N{2, 1} >= minNumTrials & N{3, 1} >= minNumTrials & N{4, 3} >= minNumTrials;
for iTarget = [2 4 1 3]
    ax = nexttile(layout.left.tl);
    if iTarget == 2 || iTarget == 1
        [~, order] = EphysUnit.plotETA(ax, trajCombined4tgt.eta(iTarget, IPAW(iTarget)), selUnit, event='reach onset', ...
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
        EphysUnit.plotETA(ax, trajCombined4tgt.eta(iTarget, IPAW(iTarget)), selUnit, event='reach onset', order=order, clim=[-2, 2], xlim=p.etaWindow);
        yticks(ax, []);
        ylabel(ax, '')
    end
    hold(ax, 'on')
    plot(ax, [0, 0], [0, nnz(selUnit)+1], 'k--');
    xlim(ax, [-2.5, 0.5])
    ylim(ax, [0, nnz(selUnit)+1])
    xlabel(ax, '');    
    title(ax, targetNames{iTarget})
    colorbar(ax, 'off')
    applyCustomColormap(ax, [-1.5, 1.5], hlim=[0.375, 0, 0, -0.375], llim=[0.125, 0.5, 0.5, 0.25], hpwr=.5, lpwr=1, h0=0.33);

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
assert(minNumTrials == 4)
minNumTrialsDisp = 4;
N = arrayfun(@(eta) eta.N, trajCombined4tgt.eta, 'UniformOutput', false);
sel = N{1, 1} >= minNumTrials & N{2, 1} >= minNumTrialsDisp & N{3, 1} >= minNumTrials & N{4, 3} >= minNumTrialsDisp;
subselSign = sel(:) & (c.isPressResponsive4tgt{2}(:) | c.isPressResponsive4tgt{4}(:));
subselAmp = sel(:) & c.isSelective.contraFrontVsIpsiFront4tgt(:);
subselNone = sel(:) & ~subselSign & ~subselAmp;
ss = [N{2, 1}, N{4, 3}];
ss = max(ss, [], 2) ./ min(ss, [], 2);
ss = 10./ss;
scatter(ax, trajCombined4tgt.meta{2, 1}(subselNone), trajCombined4tgt.meta{4, 3}(subselNone), sz, 'black', MarkerEdgeAlpha=0.5), hold(ax, 'on')
scatter(ax, trajCombined4tgt.meta{2, 1}(subselAmp), trajCombined4tgt.meta{4, 3}(subselAmp), sz, 'red')
scatter(ax, trajCombined4tgt.meta{2, 1}(subselSign), trajCombined4tgt.meta{4, 3}(subselSign), sz, 'black', 'filled')
plot(ax, [0 0], [-2 4], 'k:')
plot(ax, [-2 4], [0 0], 'k:')
plot(ax, [-2 4], [-2 4], 'k:')
axis(ax, 'equal')
xlim(ax, [-2, 4])
ylim(ax, [-2, 4])
xlabel(ax, 'Contra-front (4tgt)')
ylabel(ax, 'Ipsi-front (4tgt)')
fontsize(ax, p.fontSize, 'points')
fprintf('4tgt contra-front vs. ipsi-front: %i total, %i sign-change, %i amplitude change.\n', nnz(sel), nnz(subselSign), nnz(subselAmp))

ax = nexttile(layout.right.tl);
% assert(minNumTrials == 4)
minNumTrialsDisp = 4;
N = arrayfun(@(eta) eta.N, trajCombined4tgt.eta, 'UniformOutput', false);
sel = N{1, 1} >= minNumTrialsDisp & N{2, 1} >= minNumTrials & N{3, 1} >= minNumTrialsDisp & N{4, 3} >= minNumTrials;
subselSign = sel(:) & (c.isPressResponsive4tgt{1}(:) | c.isPressResponsive4tgt{3}(:));
subselAmp = sel(:) & c.isSelective.contraOutVsContraIn4tgt(:);
subselNone = sel(:) & ~subselSign & ~subselAmp;
ss = [N{1, 1}, N{3, 1}];
ss = max(ss, [], 2) ./ min(ss, [], 2);
ss = 10./ss;
scatter(ax, trajCombined4tgt.meta{1, 1}(subselNone), trajCombined4tgt.meta{3, 1}(subselNone), sz, 'black', MarkerEdgeAlpha=0.5), hold(ax, 'on')
scatter(ax, trajCombined4tgt.meta{1, 1}(subselAmp), trajCombined4tgt.meta{3, 1}(subselAmp), sz, 'red')
scatter(ax, trajCombined4tgt.meta{1, 1}(subselSign), trajCombined4tgt.meta{3, 1}(subselSign), sz, 'black', 'filled')
plot(ax, [0 0], [-2 4], 'k:')
plot(ax, [-2 4], [0 0], 'k:')
plot(ax, [-2 4], [-2 4], 'k:')
axis(ax, 'equal')
xlim(ax, [-2, 4])
ylim(ax, [-2, 4])
xlabel(ax, 'Contra-out (4tgt)')
ylabel(ax, 'Contra-in (4tgt)')
fontsize(ax, p.fontSize, 'points')
fprintf('4tgt contra-out vs. contra-in: %i total, %i sign-change, %i amplitude change.\n', nnz(sel), nnz(subselSign), nnz(subselAmp))


% ax = nexttile(layout.right.tl);
% ETA = trajCombined2tgt.eta;
% N = horzcat(ETA.N);
% sel = all(N >= p.minNumTrials, 2);
% subselAmp = sel(:) & c.isSelective.contraOutVsContraIn2tgt(:);
% subselSign = sel(:) & (c.isPressResponsive2tgt{1}(:) | c.isPressResponsive2tgt{2}(:));
% subselNone = sel(:) & ~subselSign & ~subselAmp;
% ex = trajCombined2tgt.eta(1);
% ey = trajCombined2tgt.eta(2);
% metaX = mean(ex.X(:, ex.t > -0.1 & ex.t < 0.2), 2);
% metaY = mean(ey.X(:, ey.t > -0.1 & ey.t < 0.2), 2);
% ss = N;
% ss = max(ss, [], 2) ./ min(ss, [], 2);
% ss = 10./ss;
% scatter(ax, metaX(subselNone), metaY(subselNone), sz, 'black', MarkerEdgeAlpha=0.5), hold(ax, 'on')
% scatter(ax, metaX(subselAmp), metaY(subselAmp), sz, 'red')
% scatter(ax, metaX(subselSign), metaY(subselSign), sz, 'black', 'filled')
% yline(ax, 0, 'k:')
% xline(ax, 0, 'k:')
% plot(ax, [-2 4], [-2 4], 'k:')
% axis(ax, 'equal')
% xlim(ax, [-2, 4])
% ylim(ax, [-2, 4])
% xlabel(ax, 'Lateral (2tgt)')
% ylabel(ax, 'Medial (2tgt)')
% fontsize(ax, p.fontSize, 'points')
% fprintf('2tgt: %i total, %i sign-change, %i amplitude change.\n', nnz(sel), nnz(subselSign), nnz(subselAmp))

% ax = nexttile(layout.right.tl);
% sel = c.hasPress & c.hasLick;
% subselSign = sel & (c.isPressResponsive & c.isLickResponsive);
% subselAmp = sel & c.isPressVsLickSelective;
% subselNone = sel & ~subselSign & ~subselAmp;
% scatter(ax, meta.lick(subselNone), meta.press(subselNone), sz, 'black', MarkerEdgeAlpha=0.5), hold(ax, 'on')
% scatter(ax, meta.lick(subselAmp), meta.press(subselAmp), sz, 'red')
% scatter(ax, meta.lick(subselSign), meta.press(subselSign), sz, 'black', 'filled')
% plot(ax, [0 0], [-2 4], 'k:')
% plot(ax, [-2 4], [0 0], 'k:')
% plot(ax, [-2 4], [-2 4], 'k:')
% axis(ax, 'equal')
% xlim(ax, [-2, 4])
% ylim(ax, [-2, 4])
% xlabel(ax, 'Peri-lick')
% ylabel(ax, 'Peri-reach')
% fontsize(ax, p.fontSize, 'points')
% mdl = fitlm(meta.lick(sel), meta.press(sel));
% fprintf('press vs. lick: %i total, %i sign-change, %i amplitude change (LM slope p<%g).\n', nnz(sel), nnz(subselSign), nnz(subselAmp), mdl.Coefficients.pValue(2))



lgd = legend(h, Orientation='horizontal', NumColumns=2);
lgd.Layout.Tile = 'north';

copygraphics(fig, ContentType='vector', BackgroundColor='none')
