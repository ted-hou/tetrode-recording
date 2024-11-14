%%% Figure 5. Lick vs. Reach ephys

%% Load all units
load_ephysunits;
find_osci_lick_circ;

% Plot ETA from peri-move to cyclick
eta.lickBoutEnd = eu.getETA('count', 'lickboutend', window=[0, 2], resolution=[2*pi/5, 0.025], normalize='none', ...
    maxInterval=0.2, minInterval=0.05, minBoutCycles=6);

eta.correctLickBout = eu.getETA('count', 'lick+lickbout', window=[-4, 0], resolution=[0.025, 2*pi/5], normalize='none', minTrialDuration=4, ...
    maxInterval=0.2, minInterval=0.05, minBoutCycles=2, maxBoutCycles=6);


eta.correctPressBout = eu.getETA('count', 'press+lickbout', window=[-4, 0], resolution=[0.025, 2*pi/32, 2*pi/5], normalize='none', minTrialDuration=4, ...
    maxInterval=0.2, minInterval=0.05, minBoutCycles=2, maxBoutCycles=6);


eta.incorrectLickBout = eu.getETA('count', 'lick+lickbout', window=[-4, 0], resolution=[0.025, 2*pi/5], normalize='none', minTrialDuration=2, maxTrialDuration=4, ...
    maxInterval=0.2, minInterval=0.05, minBoutCycles=2, maxBoutCycles=6);

eta.pressCueRaw = eu.getETA('count', 'press', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize='none', includeInvalid=true, resolution=0.025);
eta.lickCueRaw = eu.getETA('count', 'lick', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize='none', includeInvalid=true, resolution=0.025);


etaFine.press = eu.getETA('count', 'press', [-4, 2], minTrialDuration=2, normalize=[-4, -2], resolution=0.025);
etaFine.lick = eu.getETA('count', 'lick', [-4, 2], minTrialDuration=2, normalize=[-4, -2], resolution=0.025);

etaFine.correctPress = eu.getETA('count', 'press', [-4, 2], minTrialDuration=4, normalize=[-4, -2], resolution=0.025);
etaFine.correctLick = eu.getETA('count', 'lick', [-4, 2], minTrialDuration=4, normalize=[-4, -2], resolution=0.025);
etaFine.incorrectPress = eu.getETA('count', 'press', [-4, 2], minTrialDuration=2, maxTrialDuration=4, normalize=[-4, -2], resolution=0.025);
etaFine.incorrectLick = eu.getETA('count', 'lick', [-4, 2], minTrialDuration=2, maxTrialDuration=4, normalize=[-4, -2], resolution=0.025);

etaFine.lickCorrectRaw = eu.getETA('count', 'lick', p.etaWindow, minTrialDuration=4, normalize='none', resolution=0.025);
etaFine.lickIncorrectRaw = eu.getETA('count', 'lick', p.etaWindow, minTrialDuration=2, maxTrialDuration=4, normalize='none', resolution=0.025);
etaFine.p.minTrialDuration = p.minTrialDuration;
etaFine.p.norm = p.etaNorm;
etaFine.p.resolution = 0.025;

%%
% [boot.pressVsLick.h, boot.pressVsLick.p, boot.pressVsLick.ci, boot.pressVsLick.obs] = bootstrapAmplitude(eu(c.hasLick & c.hasPress), 'press', 'lick', ...
%     responseWindowA=p.metaWindowPress, responseWindowB=p.metaWindowLick, allowedTrialDuration=[2, Inf], withReplacement=false, alpha=0.01);
% c.isPressVsLickSelective = false(1, length(eu));
% c.isPressVsLickSelective(c.hasLick & c.hasPress) = boot.pressVsLick.h;
save('C:\SERVER\Units\meta_Lite_NonDuplicate_NonDrift.mat', 'ai', 'boot', 'c', 'eta', 'etaSmooth', 'euPos', 'meta', 'msr', 'onset', 'p', 'bouts', 'trialsCircLick', 'trialsCircLickBaseline', 'bootCLick', 'lickHist', 'durations')

%% plot params
nEgUnits = 3;
p.fontSize = 9;
p.etaSortWindow = [-3, 0.1];
p.etaSignWindow = [-0.2, 0.1];
p.etaLatencyThresholdPos = 0.25;
p.etaLatencyThresholdNeg = 0.25;

layout.w = 7;
layout.h = 8;
layout.top.h = 3;
layout.middle.h = 3;
layout.middle.left.w = 1;
layout.middle.right.w = 1;
layout.bottom.h = 3;
layout.bottom.left.w = 3;
layout.bottom.right.w = 3;
layout.bottom.right.ww = [(2+0.8+6/8)*100, (2+6/8)*100];

close all
fig = figure(Units='inches', Position=[0, 0, layout.w, layout.h]);
layout.tl = tiledlayout(fig, layout.top.h + layout.middle.h + layout.bottom.h, 1, TileSpacing='loose');

layout.top.tl = tiledlayout(layout.tl, 2, nEgUnits, TileSpacing='tight', Padding='tight');
l = layout.top.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [layout.top.h, 1];


layout.middle.tl = tiledlayout(layout.tl, 1, layout.middle.left.w + layout.middle.right.w, TileSpacing='compact', Padding='compact');
l = layout.middle.tl; l.Layout.Tile = 1 + layout.top.h; l.Layout.TileSpan = [layout.middle.h, 1];

layout.middle.left.tl = tiledlayout(layout.middle.tl, 1, 2, TileSpacing='compact', Padding='compact');
l = layout.middle.left.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [1, layout.middle.left.w];

layout.middle.right.tl = tiledlayout(layout.middle.tl, 1, 2, TileSpacing='compact', Padding='compact');
l = layout.middle.right.tl; l.Layout.Tile = 1 + layout.middle.left.w; l.Layout.TileSpan = [1, layout.middle.right.w];

layout.bottom.tl = tiledlayout(layout.tl, 1, layout.bottom.left.w + layout.bottom.right.w, TileSpacing='compact', Padding='compact');
l = layout.bottom.tl; l.Layout.Tile = 1 + layout.top.h + layout.middle.h; l.Layout.TileSpan = [layout.bottom.h, 1];

layout.bottom.left.tl = tiledlayout(layout.bottom.tl, 3, 3, TileSpacing='compact', Padding='compact');
l = layout.bottom.left.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [1, layout.bottom.left.w];

layout.bottom.right.tl = tiledlayout(layout.bottom.tl, 2, sum(layout.bottom.right.ww), TileSpacing='compact', Padding='compact');
l = layout.bottom.right.tl; l.Layout.Tile = 1 + layout.bottom.left.w; l.Layout.TileSpan = [1, layout.bottom.right.w];

% 5a. Double rasters (reach vs lick) 4 example units
unitNames = { ...
    'desmond24_20220510_Channel44_Unit1'; ...
    'Daisy2_20180420_Channel14_Unit1'; ...
    'daisy13_20220106_Electrode97_Unit1'; ...
%     'daisy8_20210709_Channel7_Unit1'; ...
    };
files = cellfun(@(name) sprintf('C:\\SERVER\\Units\\Lite_NonDuplicate\\%s.mat', name), unitNames, UniformOutput=false);
euEg = EphysUnit.load(files);
assert(nEgUnits == length(unitNames))

AX = gobjects(2, nEgUnits);
for i = 1:2
    for iEu = 1:nEgUnits
        AX(i, iEu) = nexttile(layout.top.tl);
    end
end
TRIALTYPE = {'press'; 'lick'};
TITLE = ["Reach trials"; "Lick trials"];
XLABEL = ["Time to bar-contact(s)"; "Time to spout-contact (s)"];
MINTRIALDURATION = [...
        2, 2, 2; ...
        2, 2, 2; ...
    ];
MAXTRIALDURATION = [...
        Inf, Inf, Inf; ...
        Inf, Inf, Inf; ...
    ];
EVERYNTH = [4, 4, 4];
YLIM = {[0, 30], [20, 140], [20, 80]};
YTICKS = {0:15:30, 20:60:140, 20:30:80};

for iEu = 1:nEgUnits
    for i = 1:2
        ax = AX(i, iEu);
        theseTrials = euEg(iEu).getTrials(TRIALTYPE{i});
        theseTrials = theseTrials(theseTrials.duration() >= MINTRIALDURATION(i, iEu) & theseTrials.duration() <= MAXTRIALDURATION(i, iEu));
        thisRD = euEg(iEu).getRasterData(TRIALTYPE{i}, window=[0, 2], sort=true, trials=theseTrials);
        thisETA = euEg(iEu).getETA('count', TRIALTYPE{i}, p.etaWindow, normalize='none', trials=theseTrials, includeInvalid=false);
        yyaxis(ax, 'right')
        EphysUnit.plotRaster(ax, thisRD, xlim=[-4, 2], iti=false, sz=1, maxTrials=40, maxTrialsMethod='uniformsample', ...
            sz=1, everyNth=EVERYNTH(iEu), timingCriterion=4);
        ylabel(ax, '')
        yticks(ax, [])
        ax.YAxis(2).Direction = 'reverse';
        yyaxis(ax, 'left')
        plot(ax, thisETA.t, thisETA.X./0.1, LineWidth=1.5, Color='black')
        hold(ax, 'on')
        set(ax.YAxis, FontSize=p.fontSize, Color=[0.15, 0.15, 0.15]);
        ylabel(ax, 'Spike rate (sp/s)')
        delete(ax.Legend)
        title(ax, TITLE(i))
        xlabel(ax, XLABEL{i})
        plot(ax, [0, 0], [0, 100], 'k--', LineWidth=1)
        ylim(ax, YLIM{iEu})
        yticks(ax, YTICKS{iEu})
        fontsize(ax, p.fontSize, 'points')
    end
end

ax = AX(1, 1);
hLetter = text(ax, 0, 0, 'a', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.3, ax.Position(4) + 0.2, 0];

% xlabel(AX, '')
ylabel(AX, '')
% xlabel(layout.top.tl, 'Time to bar/spout contact (s)', FontSize=p.fontSize)
ylabel(layout.top.tl, 'Spike rate (sp/s)', FontSize=p.fontSize)
clear thisRD ax iEu AX

% 5b/c. PETH Reach vs. Lick vs. Osci Lick
% Reach, Lick, Osci Lick, sorted amongst themselves
% close all
ax = gobjects(1, 3);
ax(1) = nexttile(layout.middle.left.tl);
ax(2) = nexttile(layout.middle.left.tl);
ax(3) = nexttile(layout.middle.right.tl);
ax(4) = nexttile(layout.middle.right.tl);


sel = c.hasPress & c.hasLick;

c.isLickUnresponsiveButUp = ~c.isLickResponsive & meta.lick > 0;
c.isLickUnresponsiveButDown = ~c.isLickResponsive & meta.lick < 0;
c.isPressUnresponsiveButUp = ~c.isPressResponsive & meta.press > 0;
c.isPressUnresponsiveButDown = ~c.isPressResponsive & meta.press < 0;

groupVar = NaN(length(sel), 1);
groupVar(c.isPressUnresponsiveButDown & c.isLickUp) = 0;
groupVar(c.isPressUnresponsiveButDown & c.isLickUnresponsiveButUp) = 0;
groupVar(c.isPressDown & c.isLickUp) = 0;
groupVar(c.isPressDown & c.isLickUnresponsiveButUp) = 0;
groupVar(c.isPressUp & c.isLickUnresponsiveButDown) = 1;
groupVar(c.isPressUp & c.isLickDown) = 1;
groupVar(c.isPressUnresponsiveButUp & c.isLickUnresponsiveButDown) = 1;
groupVar(c.isPressUnresponsiveButUp & c.isLickDown) = 1;

groupVar(c.isPressUnresponsiveButDown & c.isLickUnresponsiveButDown) = 2;
groupVar(c.isPressUnresponsiveButDown & c.isLickDown) = 2;
groupVar(c.isPressDown & c.isLickUnresponsiveButDown) = 2;
groupVar(c.isPressDown & c.isLickDown) = 2;
groupVar(c.isPressUp & c.isLickUp) = 3;
groupVar(c.isPressUp & c.isLickUnresponsiveButUp) = 3;
groupVar(c.isPressUnresponsiveButUp & c.isLickUp) = 3;
groupVar(c.isPressUnresponsiveButUp & c.isLickUnresponsiveButUp) = 3;
groupVar = groupVar(sel);
EphysUnit.plotETA(ax(1), etaFine.press, sel, ...
    clim=[-1.5, 1.5], xlim=[-2.5, 0.5], sortWindow=[-3, 0], signWindow=[-0.3, 0], ...
    sortThreshold=0.25, negativeSortThreshold=0.25, hidecolorbar=true);
EphysUnit.plotETA(ax(2), etaFine.lick, sel, ...
    clim=[-1.5, 1.5], xlim=[-2.5, 0.5], sortWindow=[-3, 0], signWindow=[-0.3, 0], ...
    sortThreshold=0.25, negativeSortThreshold=0.25, hidecolorbar=true);

[~, order] = EphysUnit.plotETA(ax(3), etaFine.press, sel, sortGroup=groupVar, ...
    clim=[-1.5, 1.5], xlim=[-2.5, 0.5], sortWindow=[-3, 0], signWindow=[-0.3, 0], ...
    sortThreshold=0.25, negativeSortThreshold=0.25, hidecolorbar=true);
EphysUnit.plotETA(ax(4), etaFine.lick, sel, order=order, sortGroup=groupVar, ...
    clim=[-1.5, 1.5], xlim=[-2.5, 0.5], hidecolorbar=true);

N = histcounts(groupVar, -0.5:2:3.5);
yline(ax(3), cumsum(N(1:end-1)) + 1, 'k--');
yline(ax(4), cumsum(N(1:end-1)) + 1, 'k--');

xline(ax(1), 0, 'k--')
xline(ax(2), 0, 'k--')
xline(ax(3), 0, 'k--')
xline(ax(4), 0, 'k--')
ylim(ax(1:4), [0, nnz(sel)+1])
yt = 0:100:nnz(sel);
yt(1) = 1;
if round(yt(end)./100) == round(nnz(sel)./100)
    yt(end) = nnz(sel);
else
    yt(end + 1) = nnz(sel);
end
yt = unique(yt);
yticks(ax(1), yt)
yticks(ax(2:end), [])

title(ax([1, 3]), 'Reach')
title(ax([2, 4]), 'Lick')
ylabel(layout.middle.tl, 'Unit', FontSize=p.fontSize)
ylabel(ax, '')
xlabel(layout.middle.tl, 'Time to bar/spout contact (s)', FontSize=p.fontSize)
xlabel(ax, '')
fontsize(ax, p.fontSize, 'points')
fontname(ax, 'Arial')
axc = ax(4);

hLetter = text(ax(1), 0, 0, 'b', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax(1).Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.3, ax(1).Position(4) + 0.1, 0];

hLetter = text(ax(3), 0, 0, 'c', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax(1).Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.1, ax(1).Position(4) + 0.1, 0];

clear yt

%5d Mean ETA grouped by press/reach up/down
SEL_HASTRIAL = c.hasPress & c.hasLick;
SEL_I = {c.isPressUp; ~c.isPressResponsive; c.isPressDown};
SEL_J = {c.isLickDown; ~c.isLickResponsive; c.isLickUp};

h = gobjects(2, 1);
for i = 1:length(SEL_I)
    for j = 1:length(SEL_J)
        ax = nexttile(layout.bottom.left.tl);
        sel = SEL_HASTRIAL & SEL_I{i} & SEL_J{j};
        assert(numel(sel) == numel(eu))
        hold(ax, 'on')
        h(1) = plot(ax, etaFine.press.t, mean(etaFine.press.X(sel, :), 1, 'omitnan'), 'r', LineWidth=1.5, DisplayName='reach');
        h(2) = plot(ax, etaFine.lick.t, mean(etaFine.lick.X(sel, :), 1, 'omitnan'), 'b', LineWidth=1.5, DisplayName='lick');
        text(ax, -2.2, -1.8, sprintf('n=%i', nnz(sel)), fontSize=p.fontSize, VerticalAlignment='bottom')
        hold(ax, 'off')
        xlim(ax, [-2.5, 0.5])
        ylim(ax, [-2, 2])
        fontsize(ax, p.fontSize, 'points')
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        xticks(ax, [-4, -2, 0])
        yticks(ax, [-2, 0, 2])
        if i < length(SEL_I)
            ax.XTickLabel = [];
        end
        if j > 1
            ax.YTickLabel = [];
        end
        fontsize(ax, p.fontSize, 'points')
        if i == 1 && j == 1
            hLetter = text(ax, 0, 0, 'd', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
            ax.Units = 'inches';
            hLetter.HorizontalAlignment = 'right';
            hLetter.VerticalAlignment = 'top';
            hLetter.Position = [-0.3, ax.Position(4) + 0.3, 0];
        end
    end
end
lgd = legend(h, Orientation='horizontal');
lgd.Layout.Tile = 'north';
xlabel(layout.bottom.left.tl, 'Time to bar/spout contact (s)', FontSize=p.fontSize)
ylabel(layout.bottom.left.tl, 'Normalized spike rate (a.u.)', FontSize=p.fontSize)

% 5e. Correct trial PETHs
yl = [0, 120];
selCommon = c.hasPress & c.hasLick;
TASK = ["press", "lick"];
ETACUE = {eta.pressCueRaw, eta.lickCueRaw};
ETAPERIMOVE = {eta.correctPressBout, eta.correctLickBout};
ETABOUTEND = {eta.lickBoutEnd, eta.lickBoutEnd};
SEL = {selCommon & c.isPressUp, selCommon & c.isPressUp; ...
    selCommon & c.isPressDown, selCommon & c.isPressDown};
SUBSEL = {c.isLickUp, c.isLickDown};
TITLE = ["Self-timed Reach", "Self-timed Lick"];
LABEL1 = ["Reach Inc", "Reach Inc"; "Reach Dec", "Reach Dec"];
LABEL2 = ["lick inc", "lick dec"];
COLORS = {[236, 131, 5]./255, [136, 194, 115]./255};

nRows = size(SEL, 1);
nCols = size(SEL, 2);

AX = gobjects(nRows, nCols);
hLine = gobjects(2, 2, 2);
for iRow = 1:nRows
    for iCol = 1:nCols
        sel = SEL{iRow, iCol};
        % Peri-move
        ax = nexttile(layout.bottom.right.tl); AX(iRow, iCol) = ax;
        ax.Layout.TileSpan = [1, layout.bottom.right.ww(iCol)];
        hold(ax, 'on')
        t = ETAPERIMOVE{iCol}.t;
        switch TASK(iCol)
            case 'press'
                t(t>0 & t<=2*pi) = t(t>0 & t<=2*pi) / (2*pi) * 0.8;
                t(t>2*pi) = (t(t>2*pi) - 2*pi) ./ (2*pi) / 8 + 0.8;
            case 'lick'
                t(t>0) = t(t>0) ./ (2*pi) / 8;
        end
        tt = [t, flip(t)];

        for iSubSel = 1:2
            subSel = sel & SUBSEL{iSubSel};
            X = ETAPERIMOVE{iCol}.X(subSel, :);
            err = std(X, 0, 1, 'omitnan') ./ sqrt(size(X, 1));
            xx = [mean(X, 1, 'omitnan') + err, flip(mean(X, 1, 'omitnan') - err)];

            hLine(iRow, iCol, iSubSel) = plot(ax, t, mean(X, 1, 'omitnan'), Color=COLORS{iSubSel}, LineWidth=1.5, DisplayName=sprintf('%s', LABEL2(iSubSel)));
            patch(ax, tt(~isnan(xx)), xx(~isnan(xx)), COLORS{iSubSel}, FaceAlpha=0.25, EdgeColor='none')
            text(ax, -1.8, 120 - (iSubSel-1)*25, sprintf('n=%i', nnz(subSel)), FontSize=p.fontSize, Color=COLORS{iSubSel}, VerticalAlignment='top')
        end

        switch TASK(iCol)
            case 'press'
                xline(ax, [0.8,  0.8+(1:5)/8], LineStyle=':')
                xlim(ax, [-2, 0.8+5/8])
                xticks(ax, [-2, -1, 0, 0.8,  0.8+(1:5)/8])
                xticklabels(ax, {'-2', '-1', '0', '2\pi', '', '', '', '', '12\pi'})
            case 'lick'        
                xline(ax, (1:6)/8, LineStyle=':')       
                xlim(ax, [-2, 6/8])
                xticks(ax, [-2:1:0, (1:6)/8])
                xticklabels(ax, {'-2', '-1', '0', '', '', '', '', '', '12\pi'}) 
        end
        ax.XAxis.TickLabelRotation = 0;
        xline(ax, 0, '-')
    end
end

xlabel(layout.bottom.right.tl, 'Time (s) or lick phase', FontSize=p.fontSize)
ylabel(layout.bottom.right.tl, 'Spike rate (sp/s)', FontSize=p.fontSize)
ylim(AX(:), yl)
yticks(AX(:, 2), [])
title(AX(1, 1), 'Reach')
title(AX(1, 2), 'Lick')
fontsize(AX(:), p.fontSize, 'points')

hLetter = text(AX(1, 1), 0, 0, 'e', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.3, AX(1, 1).Position(4) + 0.75, 0];

lgd = legend(squeeze(hLine(1, 1, :)), Orientation='horizontal');
lgd.Layout.Tile = 'north';

h = colorbar(axc); 
% h.Location = 'eastoutside';
h.Layout.Tile = 'east';
h.Label.String = 'Normalized spike rate (a.u.)';
copygraphics(fig, ContentType='vector')


%% Supplement S5

%% Supplement S5
minBoutCycles = 4;
maxBoutCycles = 16;
close all
p.fontSize = 9;
clear layout
layout.w = 7.5;
layout.h = 4;
layout.top.h = 12;

fig = figure(Units='inches', Position=[0.2 0.2 layout.w layout.h]); 
layout.tl = tiledlayout(fig, layout.top.h, 1, TileSpacing='compact', Padding='compact');

yl = [0, 120];
selCommon = c.hasPress & c.hasLick;
TASK = ["press", "lick"];
ETACUE = {eta.pressCueRaw, eta.lickCueRaw};
ETAPERIMOVE = {eta.correctPressBout, eta.correctLickBout};
ETABOUTEND = {eta.lickBoutEnd, eta.lickBoutEnd};
SEL = {selCommon & c.isPressUp, selCommon & c.isPressUp; ...
    selCommon & c.isPressDown, selCommon & c.isPressDown};
SUBSEL = {c.isLickUp, c.isLickDown};
TITLE = ["Self-timed Reach", "Self-timed Lick"];
LABEL1 = ["Reach Inc", "Reach Inc"; "Reach Dec", "Reach Dec"];
LABEL2 = ["lick inc", "lick dec"];
COLORS = 'rb';

nRows = size(SEL, 1);
nCols = size(SEL, 2);
layout.wp = [0, sum([1.8, 2+0.8+6/8, 1+6/8]*100), sum([1.8, 2+6/8, 1+6/8]*100)];

layout.top.tl = tiledlayout(layout.tl, nRows, sum(layout.wp), TileSpacing='compact', Padding='none');
l = layout.top.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [layout.top.h, 1];

% Plot ETA aligned to cue (i.e. flinchiness) vs. aligned to lever touch
% close all
AX = gobjects(nRows, nCols, 3);
hLine = gobjects(2, 2, 2);
hPatch = gobjects(2, 2, 2);
hPatchDummy = gobjects(2, 2, 2);
TL = gobjects(2, 2);
for iRow = 1:nRows
    for iCol = 1:nCols
        switch TASK(iCol)
            case "press"
                layout.ww = [1.8, 2+0.8+6/8, 1+6/8]*100;
            case "lick"
                layout.ww = [1.8, 2+6/8, 1+6/8]*100;
        end
        layout.gap = [2, 2]*10;

        tl = tiledlayout(layout.top.tl, 1, sum(layout.ww) + sum(layout.gap), TileSpacing='none', Padding='tight');
        tl.Layout.Tile = (iRow-1)*sum(layout.wp) + 1 + sum(layout.wp(1:iCol));
        tl.Layout.TileSpan = [1, layout.wp(iCol + 1)];
        TL(iRow, iCol) = tl;

        bgAx = axes(tl, XTick=[], YTick=[], Box='off');
        bgAx.Layout.Tile = 1;
        bgAx.Layout.TileSpan = [1, sum(layout.ww)];
        if iRow == 1
            title(bgAx, TITLE(iCol));
        end

        % 1. Peri-cue
        ax = axes(tl); AX(iRow, iCol, 1) = ax;
        ax.Layout.Tile = 1;
        ax.Layout.TileSpan = [1, layout.ww(1)];
        sel = SEL{iRow, iCol};

        hold(ax, 'on')
        t = ETACUE{iCol}.t;
        tt = [t, flip(t)];

        for iSubSel = 1:2
            subSel = sel & SUBSEL{iSubSel};
            X = ETACUE{iCol}.X(subSel, :)./0.025;
            err = std(X, 0, 1, 'omitnan') ./ sqrt(size(X, 1));
            xx = [mean(X, 1, 'omitnan') + err, flip(mean(X, 1, 'omitnan') - err)];

            plot(ax, t, mean(X, 1, 'omitnan'), COLORS(iSubSel), LineWidth=1.5);
            patch(ax, tt(~isnan(xx)), xx(~isnan(xx)), COLORS(iSubSel), FaceAlpha=0.25, EdgeColor='none') 
        end

        switch TASK(iCol)
            case "press"
                name = 'bar deploy';
            case "lick"
                name = 'spout deploy';
        end
        hPatch(iRow, iCol, 1) = patch(ax, XData=[-0.8, 0, 0, -0.8], YData=[yl(1), yl(1), yl(2), yl(2)], FaceColor=[0.15, 0.15, 0.15], FaceAlpha=0.1, EdgeColor='none', DisplayName=name);
        hPatch(iRow, iCol, 2) = patch(ax, XData=[0, 0.1, 0.1, 0], YData=[yl(1), yl(1), yl(2), yl(2)], FaceColor=[0.8, 0.2, 0.2], FaceAlpha=0.3, EdgeColor='none', DisplayName='tone');

        ax.Box = 'off';
        xlim(ax, [-0.95, 1])
        xticks(ax, [-0.8, 0, 1])
        xline(ax, 0, '-')
        if iRow == 2
            xlabel(ax, 'start cue')
        end
        ax.XAxis.TickLabelRotation = 0;
        
        % 2. Peri-move
        ax = axes(tl); AX(iRow, iCol, 2) = ax;
        ax.Layout.Tile = sum(layout.ww(1)) + 1 + layout.gap(1);
        ax.Layout.TileSpan = [1, layout.ww(2)];
        hold(ax, 'on')
        t = ETAPERIMOVE{iCol}.t;
        switch TASK(iCol)
            case 'press'
                t(t>0 & t<=2*pi) = t(t>0 & t<=2*pi) / (2*pi) * 0.8;
                t(t>2*pi) = (t(t>2*pi) - 2*pi) ./ (2*pi) / 8 + 0.8;
            case 'lick'
                t(t>0) = t(t>0) ./ (2*pi) / 8;
        end
        tt = [t, flip(t)];

        for iSubSel = 1:2
            subSel = sel & SUBSEL{iSubSel};
            X = ETAPERIMOVE{iCol}.X(subSel, :);
            err = std(X, 0, 1, 'omitnan') ./ sqrt(size(X, 1));
            xx = [mean(X, 1, 'omitnan') + err, flip(mean(X, 1, 'omitnan') - err)];

            hLine(iRow, iCol, iSubSel) = plot(ax, t, mean(X, 1, 'omitnan'), COLORS(iSubSel), LineWidth=1.5, DisplayName=sprintf('%s (n=%i)', LABEL2(iSubSel), nnz(subSel)));
            patch(ax, tt(~isnan(xx)), xx(~isnan(xx)), COLORS(iSubSel), FaceAlpha=0.25, EdgeColor='none')
        end
        hPatchDummy(iRow, iCol, 1) = patch(ax, XData=NaN, YData=NaN, FaceColor=[0.15, 0.15, 0.15], FaceAlpha=0.1, EdgeColor='none', DisplayName=name);
        hPatchDummy(iRow, iCol, 2) = patch(ax, XData=NaN, YData=NaN, FaceColor=[0.8, 0.2, 0.2], FaceAlpha=0.3, EdgeColor='none', DisplayName='tone');
        

        switch TASK(iCol)
            case 'press'
                xline(ax, [0.8,  0.8+(1:5)/8], LineStyle=':')
                xlim(ax, [-2, 0.8+5/8])
                xticks(ax, [-2, -1, 0, 0.8,  0.8+(1:5)/8])
                if iRow == 1
                    xticklabels(ax, {'-2', '-1', '0', '2\pi', '', '', '', '', '12\pi'})
                else
                    xticklabels(ax, {'-2', '-1', '        0      \newlinebar contact', '2\pi', '', '', '', '', '12\pi'})
                end
            case 'lick'        
                xline(ax, (1:6)/8, LineStyle=':')       
                xlim(ax, [-2, 6/8])
                xticks(ax, [-2:1:0, (1:6)/8])
                if iRow == 1
                    xticklabels(ax, {'-2', '-1', '0', '', '', '', '', '', '12\pi'}) 
                else
                    xticklabels(ax, {'-2', '-1', '          0      \newlinespout contact', '', '', '', '', '', '12\pi'}) 
                end
        end
        ax.XAxis.TickLabelRotation = 0;
        xline(ax, 0, '-')
        ax.Box = 'off';
        ax.YAxis.Visible = 'off';
    
        % 3. End of bout
        ax = axes(tl); AX(iRow, iCol, 3) = ax;
        ax.Layout.Tile = sum(layout.ww(1:2)) + 1 + sum(layout.gap(1:2));
        ax.Layout.TileSpan = [1, layout.ww(3)];
        hold(ax, 'on')
        
        t = eta.lickBoutEnd.t;
        t(t<0) = t(t<0) ./ (2*pi) / 8;
        tt = [t, flip(t)];

        for iSubSel = 1:2
            subSel = sel & SUBSEL{iSubSel};
            X = ETABOUTEND{iCol}.X(subSel, :);
            err = std(X, 0, 1, 'omitnan') ./ sqrt(size(X, 1));
            xx = [mean(X, 1, 'omitnan') + err, flip(mean(X, 1, 'omitnan') - err)];

            plot(ax, t, mean(X, 1, 'omitnan'), COLORS(iSubSel), LineWidth=1.5);
            patch(ax, tt(~isnan(xx)), xx(~isnan(xx)), COLORS(iSubSel), FaceAlpha=0.25, EdgeColor='none')
        end

        yline(ax, 0, 'k')
        xlim(ax, [-6/8, 1])
        xticks(ax, [(-6:-1)/8, 0:1])
        xline(ax, (-6:-1)/8, LineStyle=':')       
        xticklabels(ax, {'', '', '', '-6\pi', '', '', '0', '1'})
        ax.XAxis.TickLabelRotation = 0;
        if iRow == 2
            xlabel(ax, 'last lick')
        end
        ax.Box = 'off';
        ax.YAxis.Visible = 'off';
        xline(ax, 0, '-')
    
        ylim(AX(iRow, iCol, :), yl)
        linkaxes(AX(iRow, iCol, :), 'y')
    end
end

xlabel(layout.top.tl, 'Time (s) or lick phase', FontSize=p.fontSize)
ylabel(layout.top.tl, 'Spike rate (sp/s)', FontSize=p.fontSize)
fontsize(AX(:), p.fontSize, 'points')


hLineDummy = gobjects(2, 1);
lgdPositions = {[0.08,0.925,0.18,0.033], [0.58,0.925,0.18,0.033]; ...
    [0.08,0.79,0.18,0.033], [0.58,0.79,0.18,0.033]};
for iRow = 1:2
    for iCol = 1:2
        axDummy = axes(fig, Position=TL(iRow, iCol).Position);
        hold(axDummy, 'on')
        subSel = SEL{iRow, 1} & SUBSEL{1};
        hLineDummy(1) = plot(axDummy, NaN, NaN, 'red', LineWidth=1.5, DisplayName=sprintf('%s (n=%i)', LABEL2(1), nnz(subSel)));
        subSel = SEL{iRow, 1} & SUBSEL{2};
        hLineDummy(2) = plot(axDummy, NaN, NaN, 'blue', LineWidth=1.5, DisplayName=sprintf('%s (n=%i)', LABEL2(2), nnz(subSel)));
        axDummy.Visible = 'off';
        legend(axDummy, hLineDummy, Location='northwest', fontSize=8, Position=lgdPositions{iRow, iCol})
    end
end


copygraphics(fig, ContentType='vector')
