%%
% 8a. Spontaneous reach task, trial structure
Done
% 8b. Distribution of Inter-touch-intervals
% 8c. Example units (raster vs ETA)

% 8d. Heatmap PETH for all units

%%
% read_spontaneous; % This takes a while because of boostrapping, also does
% drift/multiunit/duplicate removal
clear
euSpontaneous = EphysUnit.load('C:\SERVER\Units\acute_spontaneous_reach\SNr_SingleUnit_NonDuplicate_NonDrift');
load('C:\SERVER\Units\acute_spontaneous_reach\meta\SNr_SingleUnit_NonDuplicate_NonDrift.mat')

%%
p.fontSize = 9;

clear layout
layout.w = 7;
layout.h = 4.5;
layout.left.w = 4;
layout.right.w = 3;
layout.left.top.h = 3;
layout.left.bottom.h = 6;
layout.right.top.h = 3;
layout.right.bottom.h = 5;

close all
fig = figure(Units='inches', Position=[0, 0, layout.w, layout.h]);
layout.tl = tiledlayout(fig, 1, layout.left.w + layout.right.w, TileSpacing='loose', Padding='loose');
layout.left.tl = tiledlayout(layout.tl, layout.left.top.h + layout.left.bottom.h, 1, TileSpacing='loose', Padding='loose');
l = layout.left.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [1, layout.left.w];

layout.right.tl = tiledlayout(layout.tl, layout.right.top.h + layout.right.bottom.h, 1, TileSpacing='loose', Padding='loose');
l = layout.right.tl; l.Layout.Tile = 1 + layout.left.w; l.Layout.TileSpan = [1, layout.right.w];

layout.left.bottom.tl = tiledlayout(layout.left.tl, 2, 2, TileSpacing='compact', Padding='loose');
l = layout.left.bottom.tl; l.Layout.Tile = 1 + layout.left.top.h; l.Layout.TileSpan = [layout.left.bottom.h, 1];

% 8b. Distribution of Inter-touch-intervals
ax = nexttile(layout.right.tl, 1, [layout.right.top.h, 1]);
interTouchIntervals = cell(length(expSpontaneous), 1);
for iExp = 1:length(expSpontaneous)
    touchTimes = [expSpontaneous(iExp).eu(1).Trials.Press.Stop];
    interTouchIntervals{iExp} = diff(touchTimes);
end
interTouchIntervals = cat(2, interTouchIntervals{:});
% edges = 0:2:max(ceil(interTouchIntervals/2)*2);
edges = 0.5:0.5:15;
histogram(ax, interTouchIntervals, edges, Normalization='probability', ...
    EdgeAlpha=1, FaceColor='black')
xlabel(ax, 'Inter-reach interval (s)')
ylabel(ax, 'Probability')
yticks(ax, 0:0.1:0.2)
xticks(ax, [0.5, 5, 10, 15])
fontsize(ax, p.fontSize, 'points')
fontname(ax, 'Arial')
% copygraphics(fig)

% 8c. Example unit Raster (maybe ETA overlayed)
% Let's try to find two examples by how much spike rate increases relative
% to baseline
t = etaSpontaneousRaw.t;
X = etaSpontaneousRaw.X;
xBase = mean(X(:, t < -2 & t > -4), 2, 'omitnan');
xMove = mean(X(:, t < -0 & t > -0.5), 2, 'omitnan');
xDiff = xMove - xBase;
[~, IDIFF] = sort(xDiff, 'ascend');
[~, IMOVE] = sort(xMove, 'ascend');
[~, IDIFFMOVE] = sort(xMove + xDiff, 'ascend');
[~, IZ] = sort(bootSpontaneous.press.muDiffObs, 'ascend');
iMin = IDIFF(2);
iMax = IDIFF(end-4);

exampleUnitNames = { ...
%         euSpontaneous(iMax).getName, ...
%         euSpontaneous(iMin).getName, ...
        'desmond31_20230804_Channel49_Unit1', ... % Up at -1, 40sp/s
        'daisy18_20230802_Channel66_Unit1', ... % Down at -1.5, wierd for 4-8s trials
    };
YLIMS = {[20, 80], [0, 80]};
NSKIP = [5, 5];

% PETH Rasster
for i = 1:length(exampleUnitNames)
    iEu = find(strcmpi(euSpontaneous.getName(), exampleUnitNames{i}));
    iExp = find(strcmpi(euSpontaneous(iEu).ExpName, {expSpontaneous.name}));
    ax = nexttile(layout.left.bottom.tl);
    nnz(isnan(onset(iExp).contra))
    trials = euSpontaneous(iEu).getTrials('press');
    selTrials = onset(iExp).contra >= onsetThreshold;
    thisRd = euSpontaneous(iEu).getRasterData('press', window=[-0, 2], sort=false, MinTrialDuration=6, ...
        correction=onset(iExp).contra(selTrials), trials=trials(selTrials));
    hold(ax, 'on')
    yyaxis(ax, 'right')
    EphysUnit.plotRaster(ax, thisRd, xlim=[-4, 2], sz=1, iti=false, ...
        maxTrials=50, maxTrialsMethod='randomsample', everyNth=NSKIP(i));
    ylabel(ax, 'Trial')
    yticks(ax, [1, 25, 50])
    ax.YAxis(2).Direction = 'reverse';
    yyaxis(ax, 'left')
    plot(ax, etaSpontaneousRaw.t, etaSpontaneousRaw.X(iEu, :)./0.1, LineWidth=1.5, Color='black')
    set(ax.YAxis, FontSize=p.fontSize, Color=[0.15, 0.15, 0.15]);
    ylabel(ax, 'Spike rate (sp/s)')
    yticks(ax, YLIMS{i}(1):20:YLIMS{i}(end))
    ylim(ax, YLIMS{i})
    title(ax, '')
    legend(ax, 'off')
    hold(ax, 'on')
    plot(ax, [0, 0], [0, 100], 'k--')
    hold(ax, 'off')
    fontsize(ax, p.fontSize, 'points')
    xlabel(ax, '')
    xticks(ax, [-4, -2, 0, 1])
    xlim(ax, [-4, 0.5])
end

SEL = {cSpontaneous.isPressUp, cSpontaneous.isPressDown};
AX = gobjects(1, 2);
for i = 1:2
    ax = nexttile(layout.left.bottom.tl);
    AX(i) = ax;
    plot(ax, etaSpontaneousRaw.t, mean(etaSpontaneousRaw.X(SEL{i}, :)./0.1, 1, 'omitnan'), 'k', LineWidth=1.5)
    hold(ax, 'on')
    iSel = find(SEL{i});
    for iTrial = iSel(:)'
        plot(ax, etaSpontaneousRaw.t, etaSpontaneousRaw.X(iTrial, :)./0.1, Color=[0 0 0 0.1])
    end
    plot(ax, [0, 0], [0, 100], 'k--')
    ylim(ax, [0, 80])
    delete(legend(ax))
    set(ax, FontSize=p.fontSize, FontName='Arial')
    ylabel(ax, 'Spike rate (sp/s)', FontSize=p.fontSize, FontName='Arial')
    xticks(ax, [-4, -2, 0, 1])
    xlim(ax, [-4, 0.5])
end
xlabel(layout.left.bottom.tl, 'Time to reach onset (s)', FontSize=p.fontSize, FontName='Arial')

% 8d. Plot ETA (touch time)
ax = nexttile(layout.right.tl, 1 + layout.right.top.h, [layout.right.bottom.h, 1]);
EphysUnit.plotETA(ax, etaSpontaneous, xlim=[-4,0.5], clim=[-1.5, 1.5], sortWindow=[-2, 0.5], signWindow=[-0.3, 0.2], sortThreshold=0.25, negativeSortThreshold=0.25);
title(ax, '')
yticks(ax, [1, 20:20:length(euSpontaneous), length(euSpontaneous)])
xlabel(ax, 'Time to reach onset (s)')
fontsize(ax, p.fontSize, 'points')
fontname(ax, 'Arial')

copygraphics(fig, ContentType='vector')
