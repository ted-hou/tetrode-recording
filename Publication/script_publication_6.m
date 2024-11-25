%%% Figure 6. Lick vs. Reach ephys

%% Load all units
load_ephysunits;
find_osci_lick_circ;


etaFine.press = eu.getETA('count', 'press', [-4, 2], minTrialDuration=2, normalize=[-4, -2], resolution=0.025);
etaFine.lick = eu.getETA('count', 'lick', [-4, 2], minTrialDuration=2, normalize=[-4, -2], resolution=0.025);

etaFine.correctPress = eu.getETA('count', 'press', [-4, 2], minTrialDuration=4, normalize=etaFine.press.stats, resolution=0.025);
etaFine.correctLick = eu.getETA('count', 'lick', [-4, 2], minTrialDuration=4, normalize=etaFine.lick.stats, resolution=0.025);
etaFine.incorrectPress = eu.getETA('count', 'press', [-4, 2], minTrialDuration=2, maxTrialDuration=4, normalize=etaFine.press.stats, resolution=0.025);
etaFine.incorrectLick = eu.getETA('count', 'lick', [-4, 2], minTrialDuration=2, maxTrialDuration=4, normalize=etaFine.lick.stats, resolution=0.025);

etaFine.lickCorrectRaw = eu.getETA('count', 'lick', p.etaWindow, minTrialDuration=4, normalize='none', resolution=0.025);
etaFine.lickIncorrectRaw = eu.getETA('count', 'lick', p.etaWindow, minTrialDuration=2, maxTrialDuration=4, normalize='none', resolution=0.025);
etaFine.p.minTrialDuration = p.minTrialDuration;
etaFine.p.norm = p.etaNorm;
etaFine.p.resolution = 0.025;

% Plot ETA from peri-move to cyclick
eta.lickBoutEnd = eu.getETA('count', 'lickboutend', window=[0, 2], resolution=[2*pi/5, 0.025], normalize='none', ...
    maxInterval=0.2, minInterval=0.05, minBoutCycles=6);

eta.correctLickBout = eu.getETA('count', 'lick+lickbout', window=[-4, 0], resolution=[0.025, 2*pi/5], normalize='none', minTrialDuration=4, ...
    maxInterval=0.2, minInterval=0.05, minBoutCycles=2, maxBoutCycles=6);
eta.correctLickBoutNorm = eu.getETA('count', 'lick+lickbout', window=[-4, 0], resolution=[0.025, 2*pi/5], normalize=[-4, -2], minTrialDuration=4, ...
    maxInterval=0.2, minInterval=0.05, minBoutCycles=2, maxBoutCycles=6);

eta.correctPressBout = eu.getETA('count', 'press+lickbout', window=[-4, 0], resolution=[0.025, 2*pi/32, 2*pi/5], normalize='none', minTrialDuration=4, ...
    maxInterval=0.2, minInterval=0.05, minBoutCycles=2, maxBoutCycles=6);
eta.correctPressBoutNorm = eu.getETA('count', 'press+lickbout', window=[-4, 0], resolution=[0.025, 2*pi/32, 2*pi/5], normalize=[-4, -2], minTrialDuration=4, ...
    maxInterval=0.2, minInterval=0.05, minBoutCycles=2, maxBoutCycles=6);


% eta.incorrectLickBout = eu.getETA('count', 'lick+lickbout', window=[-4, 0], resolution=[0.025, 2*pi/5], normalize='none', minTrialDuration=2, maxTrialDuration=4, ...
%     maxInterval=0.2, minInterval=0.05, minBoutCycles=2, maxBoutCycles=6);

eta.pressCueRaw = eu.getETA('count', 'press', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize='none', includeInvalid=true, resolution=0.025);
eta.lickCueRaw = eu.getETA('count', 'lick', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize='none', includeInvalid=true, resolution=0.025);

eta.pressCueNorm = eu.getETA('count', 'press', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize=etaFine.press.stats, includeInvalid=true, resolution=0.025);
eta.lickCueNorm = eu.getETA('count', 'lick', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize=etaFine.lick.stats, includeInvalid=true, resolution=0.025);

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

% close all
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

% 6a. Double rasters (reach vs lick) 3 example units

% selI = find(c.hasPress & c.hasLick & c.isPressResponsive & c.isLickResponsive);
% [~, downI] = sort(meta.press(selI) + meta.lick(selI), 'ascend');
% [~, upI] = sort(meta.press(selI) + meta.lick(selI), 'descend');
% [~, oppositeI] = sort(meta.press(selI) .* meta.lick(selI), 'ascend');
% downI = selI(downI);
% upI = selI(upI);
% oppositeI = selI(oppositeI);

unitNames = { ...
    'desmond25_20220430_Channel24_Unit1' % eu(downI(7)).getName(); ...
    'Daisy2_20180422_Channel21_Unit1' % eu(upI(4)).getName(); ...
    'daisy14_20220506_Channel14_Unit1' % eu(oppositeI(10)).getName(); ...
    };
[~, locb] = ismember(unitNames, eu.getName());
euEg = eu(locb);
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
EVERYNTH = [5, 5, 5];
% YLIM = {[20, 80], [20, 140], [20, 80]};
% YTICKS = {20:30:80, 20:60:140, 20:30:80};

for iEu = 1:nEgUnits
    for i = 1:2
        ax = AX(i, iEu);
        theseTrials = euEg(iEu).getTrials(TRIALTYPE{i});
        theseTrials = theseTrials(theseTrials.duration() >= MINTRIALDURATION(i, iEu) & theseTrials.duration() <= MAXTRIALDURATION(i, iEu));
        thisRD = euEg(iEu).getRasterData(TRIALTYPE{i}, window=[0, 2], sort=true, trials=theseTrials);
        thisETA = euEg(iEu).getETA('count', TRIALTYPE{i}, p.etaWindow, normalize='none', trials=theseTrials, includeInvalid=false);
        yyaxis(ax, 'right')
        EphysUnit.plotRaster(ax, thisRD, xlim=[-4, 2], iti=false, sz=1, maxTrials=40, maxTrialsMethod='uniformsample', ...
            everyNth=EVERYNTH(iEu), timingCriterion=4);
        hRaster = ax.Children(3);
        hRaster.MarkerFaceAlpha = 0.5;
        ylabel(ax, '')
        yticks(ax, [])
        ax.YAxis(2).Direction = 'reverse';
        yyaxis(ax, 'left')
        plot(ax, thisETA.t, thisETA.X./0.1, LineWidth=1.5, Color=[0.2, 0.2, 0.8, 1.0])
        hold(ax, 'on')
        set(ax.YAxis, FontSize=p.fontSize, Color=[0.15, 0.15, 0.15]);
        ylabel(ax, 'Spike rate (sp/s)')
        delete(ax.Legend)
        title(ax, TITLE(i))
        xlabel(ax, XLABEL{i})
        xline(ax, 0, 'k--', LineWidth=1)
        ylim(ax, 'auto')
        % ylim(ax, YLIM{iEu})
        % yticks(ax, YTICKS{iEu})
        fontsize(ax, p.fontSize, 'points')
    end
    yl1 = AX(1, iEu).YLim;
    yl2 = AX(2, iEu).YLim;
    yl = [min(yl1(1), yl2(1)), max(yl1(2), yl2(2))];
    ylim(AX(:, iEu), yl);
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

% 6b/c. PETH Reach vs. Lick vs. Osci Lick
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

% 6d Mean ETA grouped by press/reach up/down
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

% 6e. Correct trial PETHs
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
h.Layout.Tile = 'east';
h.Label.String = 'Normalized spike rate (a.u.)';
copygraphics(fig, ContentType='vector', BackgroundColor='none')


%% Supplement S6

close all
p.fontSize = 9;
p.lineWidth = 1.5;
nBoutsDisp = 6;
W = [0, 0.5 + nBoutsDisp/8, 0.5 + 0.8+(nBoutsDisp-1)/8]*40;

clear layout l
layout.w = 7;
layout.h = 8;
layout.left.w = 3;
layout.right.w = 6;
layout.left.h = [3, 9, 4, 12];
layout.right.h = [3, 6];

fig = figure(Units='inches', Position=[1, 1, layout.w, layout.h], DefaultAxesFontSize=p.fontSize);
layout.tl = tiledlayout(fig, 1, layout.left.w + layout.right.w, TileSpacing='compact', Padding='compact');

layout.left.tl = tiledlayout(layout.tl, sum(layout.left.h), 1, TileSpacing='compact', Padding='compact');
l = layout.left.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [1, layout.left.w];

layout.right.tl = tiledlayout(layout.tl, sum(layout.right.h), 1, TileSpacing='compact', Padding='compact');
l = layout.right.tl; l.Layout.Tile = 1 + layout.left.w; l.Layout.TileSpan = [1, layout.right.w];

layout.left.top.tl = tiledlayout(layout.left.tl, sum(layout.left.h(1:3)), 1, TileSpacing='compact', Padding='compact');
l = layout.left.top.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [sum(layout.left.h(1:3)), 1];

layout.left.bottom.tl = tiledlayout(layout.left.tl, 4, 1, TileSpacing='compact', Padding='compact');
l = layout.left.bottom.tl; l.Layout.Tile = 1 + sum(layout.left.h(1:3)); l.Layout.TileSpan = [sum(layout.left.h(4)), 1];

layout.right.top.tl = tiledlayout(layout.right.tl, 1, 2, TileSpacing='compact', Padding='compact');
l = layout.right.top.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [layout.right.h(1), 1];

layout.right.bottom.tl = tiledlayout(layout.right.tl, 4, sum(W), TileSpacing='compact', Padding='compact');
l = layout.right.bottom.tl; l.Layout.Tile = 1 + sum(layout.right.h(1)); l.Layout.TileSpan = [layout.right.h(2), 1];


%S6a. Peri-lick lick prob histogram
% S6f. Correct lick bout lick histogram
ax = nexttile(layout.left.top.tl, [layout.left.h(1), 1]);
counts = lickHist.correctLickOsci.count;
smoothedCounts = smoothdata(counts, 'gaussian', 25);
[pks, locs] = findpeaks(smoothedCounts, lickHist.correctLickOsci.t, MinPeakProminence=0.5);
histogram(ax, BinEdges=lickHist.correctLickOsci.edges, BinCounts=counts, Normalization='probability', EdgeColor='none', FaceColor='black', FaceAlpha=0.9)
hold(ax, 'on')
xline(ax, locs, LineStyle=':', Color=[0, 0.5, 0.5, 1], LineWidth=2.5)
ylabel(ax, ' lick\newlineprob')
xlabel(ax, 'Time from first lick (s)')
xticks(ax, -1:0.5:2)
xlim(ax, [0, 1/8*nBoutsDisp])
yticks(ax, [])
ax.Box = 'off';
fontsize(ax, p.fontSize, 'points')

hLetter = text(ax, 0, 0, 'a', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.35, ax.Position(4) + 0.15, 0];


% S6b ETA Heatmap osci lick
ax = nexttile(layout.left.top.tl, [layout.left.h(2), 1]);

eta.circLick.Z = eta.circLick.X.*exp(eta.circLick.t*1i);
eta.circLick.Z(:, 1) = mean(eta.circLick.X(:, [2, 30]), 2).*exp(eta.circLick.t(1)*1i);
meanZ = mean(eta.circLick.Z, 2);

eta.lickBoutNorm = eta.lickBout;
eta.lickBoutNorm.X = normalize(eta.lickBout.X, 2, 'zscore', 'robust');

maxBoutCycles = 4;
sel = c.hasPress & c.hasLick & c.isLick;
phase = angle(meanZ(sel));
amp = abs(meanZ(sel));
phase(phase < 0) = phase(phase < 0) + 2*pi;
amp = amp(:);
phase = phase(:);
[sortedPhase, I] = sort(phase);
[~, ~] = EphysUnit.plotETA(ax, eta.lickBoutNorm, sel, order=I, ...
    clim=[-2, 2], xlim=[0, 2*pi*maxBoutCycles], hidecolorbar=false);
xticks (ax, (0:2:8).*pi);
xticklabels(ax, [{'0'}, arrayfun(@(x) sprintf('%i\\pi', x), 2:2:8, UniformOutput=false)]);
title(ax, 'Lick-entrained')
ylabel(ax, 'Unit')
xlabel(ax, 'Lick phase')
xlim(ax, [0, 2*pi*maxBoutCycles])
ylim(ax, [0, nnz(sel)+1])
yt = 0:100:nnz(sel);
yt(1) = 1;
if round(yt(end)./100) == round(nnz(sel)./100)
    yt(end) = nnz(sel);
else
    yt(end + 1) = nnz(sel);
end
yt = unique(yt);
yticks(ax, yt)
fontsize(ax, p.fontSize, 'points')
fontname(ax, 'Arial')
axc = ax;

hLetter = text(ax, 0, 0, 'b', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.35, ax.Position(4) + 0.25, 0];


% Phase calculations
nBoutsDisp = 6;
sel = find(c.isLick & c.hasPress & c.hasLick);
rectifiedPhase = phase;
rectifiedPhase(phase < 0) = rectifiedPhase(phase < 0) + 2*pi;
ampThreshold = 0;
isInPhase = phase < 0.25*pi | phase >= 1.75*pi;
isAntiPhase = phase >= 0.75*pi & phase < 1.25*pi;
isFirstQuarterPhase = phase >= 0.25*pi & phase < 0.75*pi;
isThirdQuarterPhase = phase >= 1.25*pi & phase < 1.75*pi;
isHighAmp = amp >= quantile(amp, ampThreshold);
idInPhase = sel(isInPhase & isHighAmp);
idInPhaseLickUp = sel(isInPhase & isHighAmp & c.isLickUp(sel)');
idInPhaseLickFlat = sel(isInPhase & isHighAmp & ~c.isLickResponsive(sel)');
idInPhaseLickDown = sel(isInPhase & isHighAmp & c.isLickDown(sel)');
idInPhasePressUp = sel(isInPhase & isHighAmp & c.isPressUp(sel)');
idInPhasePressFlat = sel(isInPhase & isHighAmp & ~c.isPressResponsive(sel)');
idInPhasePressDown = sel(isInPhase & isHighAmp & c.isPressDown(sel)');
idAntiPhase = sel(isAntiPhase & isHighAmp);
idAntiPhaseLickUp = sel(isAntiPhase & isHighAmp & c.isLickUp(sel)');
idAntiPhaseLickFlat = sel(isAntiPhase & isHighAmp & ~c.isLickResponsive(sel)');
idAntiPhaseLickDown = sel(isAntiPhase & isHighAmp & c.isLickDown(sel)');
idAntiPhasePressUp = sel(isAntiPhase & isHighAmp & c.isPressUp(sel)');
idAntiPhasePressFlat = sel(isAntiPhase & isHighAmp & ~c.isPressResponsive(sel)');
idAntiPhasePressDown = sel(isAntiPhase & isHighAmp & c.isPressDown(sel)');
idFirstQuarterPhase = sel(isFirstQuarterPhase & isHighAmp);
idFirstQuarterPhaseLickUp = sel(isFirstQuarterPhase & isHighAmp & c.isLickUp(sel)');
idFirstQuarterPhaseLickFlat = sel(isFirstQuarterPhase & isHighAmp & ~c.isLickResponsive(sel)');
idFirstQuarterPhaseLickDown = sel(isFirstQuarterPhase & isHighAmp & c.isLickDown(sel)');
idFirstQuarterPhasePressUp = sel(isFirstQuarterPhase & isHighAmp & c.isPressUp(sel)');
idFirstQuarterPhasePressFlat = sel(isFirstQuarterPhase & isHighAmp & ~c.isPressResponsive(sel)');
idFirstQuarterPhasePressDown = sel(isFirstQuarterPhase & isHighAmp & c.isPressDown(sel)');
idThirdQuarterPhase = sel(isThirdQuarterPhase & isHighAmp);
idThirdQuarterPhaseLickUp = sel(isThirdQuarterPhase & isHighAmp & c.isLickUp(sel)');
idThirdQuarterPhaseLickFlat = sel(isThirdQuarterPhase & isHighAmp & ~c.isLickResponsive(sel)');
idThirdQuarterPhaseLickDown = sel(isThirdQuarterPhase & isHighAmp & c.isLickDown(sel)');
idThirdQuarterPhasePressUp = sel(isThirdQuarterPhase & isHighAmp & c.isPressUp(sel)');
idThirdQuarterPhasePressFlat = sel(isThirdQuarterPhase & isHighAmp & ~c.isPressResponsive(sel)');
idThirdQuarterPhasePressDown = sel(isThirdQuarterPhase & isHighAmp & c.isPressDown(sel)');
colors = getColor(1:4, 4, 0.6);

euSel = eu(sel);
[~, expEuIndices] = unique({euSel.ExpName});
FsLick = 500;
lickHistEdges = 0.01:1/FsLick:2;
lickHistCenters = 0.5*(lickHistEdges(2:end) + lickHistEdges(1:end-1));

lickHistCounts = zeros(size(lickHistCenters));
lickHistNLicks = 0;
for iExp = 1:length(expEuIndices)
    iEu = expEuIndices(iExp);
    trials = euSel(iEu).getTrials('lick');
    trials = trials(trials.duration() >= 4);
    firstLickTimes = [trials.Stop];
    allLickTimes = euSel(iEu).EventTimes.Lick;
    for iLick = 1:length(firstLickTimes)
        edgesGlobal = firstLickTimes(iLick) + lickHistEdges;
        n = histcounts(allLickTimes, edgesGlobal);
        lickHistCounts = lickHistCounts + n;
        lickHistNLicks = lickHistNLicks + 1;
    end
end

lickHist.correctLickOsci.count = lickHistCounts;
lickHist.correctLickOsci.pdf = lickHistCounts ./ lickHistNLicks;
lickHist.correctLickOsci.t = lickHistCenters;
lickHist.correctLickOsci.edges = lickHistEdges;


% S5c. Osci lick phase distribution histogram
ax = nexttile(layout.left.top.tl, [layout.left.h(3), 1]);
hold(ax, 'on')
edges = -2*pi:2*pi/50:4*pi;
histogram(ax, rectifiedPhase(isInPhase), edges, FaceColor=colors(1, :), FaceAlpha=1, EdgeAlpha=0.5);
histogram(ax, rectifiedPhase(isFirstQuarterPhase), edges, FaceColor=colors(2, :), FaceAlpha=1, EdgeAlpha=0.5);
histogram(ax, rectifiedPhase(isAntiPhase), edges, FaceColor=colors(4, :), FaceAlpha=1, EdgeAlpha=0.5);
histogram(ax, rectifiedPhase(isThirdQuarterPhase), edges, FaceColor=colors(3, :), FaceAlpha=1, EdgeAlpha=0.5);
xticks(ax, 0:pi:2*pi)
xlim(ax, pi*[0, 2])
xticklabels(ax, {'0', '\pi', '2\pi'});
xlabel('Lick phase')
ylabel('# units')
fontsize(ax, p.fontSize, 'points')

hLetter = text(ax, 0, 0, 'c', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.35, ax.Position(4) + 0.25, 0];

ID = {idInPhase; idFirstQuarterPhase; idAntiPhase; idThirdQuarterPhase};
IDSplit = { ...
    {idInPhaseLickUp, idInPhaseLickFlat, idInPhaseLickDown}, {idInPhasePressUp, idInPhasePressFlat, idInPhasePressDown}; ...
    {idFirstQuarterPhaseLickUp, idFirstQuarterPhaseLickFlat, idFirstQuarterPhaseLickDown}, {idFirstQuarterPhasePressUp, idFirstQuarterPhasePressFlat, idFirstQuarterPhasePressDown}; ...
    {idAntiPhaseLickUp, idAntiPhaseLickFlat, idAntiPhaseLickDown}, {idAntiPhasePressUp, idAntiPhasePressFlat, idAntiPhasePressDown}; ...
    {idThirdQuarterPhaseLickUp, idThirdQuarterPhaseLickFlat, idThirdQuarterPhaseLickDown}, {idThirdQuarterPhasePressUp, idThirdQuarterPhasePressFlat, idThirdQuarterPhasePressDown}; ...
    };
ETAMOVEBOUT = {eta.correctLickBoutNorm, eta.correctPressBoutNorm};
TASK = ["lick", "press"];
TASKTITLE = ["Self-timed lick", "Self-timed reach"];
PHASENAME = ["2\pi", "1/2\pi", "\pi", "3/2\pi"];
ICOLOR = [1, 2, 4, 3];

% S6d (left) and S6g (right)
AX = gobjects(4, 3);
for iAx = 1:4
    iEu = ID{iAx};
    for iTask = 1:2
        % First lick
        ax = nexttile(layout.right.bottom.tl, (iAx-1)*sum(W) + 1 + sum(W(1:iTask)), [1, W(iTask + 1)]); 
        AX(iAx, iTask + 1) = ax;
        hold(ax, 'on')
        t = ETAMOVEBOUT{iTask}.t;
        switch TASK(iTask)
            case "press"
                t(t>0 & t<=2*pi) = t(t>0 & t<=2*pi) / (2*pi) * 0.8;
                t(t>2*pi) = (t(t>2*pi) - 2*pi) ./ (2*pi) / 8 + 0.8;
            case "lick"
                t(t>0) = t(t>0) ./ (2*pi) / 8;
        end

        for iDir = [1, 3]
            iEuDir = IDSplit{iAx, iTask}{iDir};
            X = ETAMOVEBOUT{iTask}.X(iEuDir, :);
            X = smoothdata(X, 2, 'gaussian', 5);
            plot(ax, t, mean(X, 1, 'omitnan'), Color=colors(ICOLOR(iAx), :), LineWidth=1.5)
            plot(ax, t, X, Color=[0.15, 0.15, 0.15, 0.05], LineWidth=0.5)
        end

%         X = ETAMOVEBOUT{iTask}.X(iEu, :);
%         X = smoothdata(X, 2, 'gaussian', 5);
%         plot(ax, t, mean(X, 1, 'omitnan'), Color=colors(ICOLOR(iAx), :), LineWidth=1.5)
%         plot(ax, t, X, Color=[0.15, 0.15, 0.15, 0.1], LineWidth=0.5)

        ylim(ax, [-1.25, 1.75])
        yticks(ax, [-1, 0, 1])
        
        switch TASK(iTask)
            case "press"
                xline(ax, [0.8,  0.8+(1:(nBoutsDisp-1))/8], LineStyle=':')
                xlim(ax, [-0.5, 0.8+(nBoutsDisp-1)/8])
                xticks(ax, [-2, -0.5, 0, 0.8,  0.8+(1:(nBoutsDisp-1))/8])
                xticklabels(ax, {'-2', '-0.5', '0', '2\pi', '', '', '', '', '12\pi'})
            case "lick"    
                xline(ax, (0:nBoutsDisp)/8, LineStyle=':')       
                xlim(ax, [-0.5, nBoutsDisp/8])
                xticks(ax, [-2, -0.5, 0, (1:nBoutsDisp)/8])
                xticklabels(ax, {'-2', '-0.5', '0', '', '', '', '', '', '12\pi'}) 
        end
        ax.XAxis.TickLabelRotation = 0;

        text(ax, 0.025, 1, sprintf('inc(n=%i)', nnz(IDSplit{iAx, iTask}{1})), Unit='normalized', HorizontalAlignment='left', VerticalAlignment='top', Interpreter='none', FontSize=p.fontSize-1)
        text(ax, 0.025, 0.025, sprintf('dec(n=%i)', nnz(IDSplit{iAx, iTask}{3})), Unit='normalized', HorizontalAlignment='left', VerticalAlignment='bottom', Interpreter='none', FontSize=p.fontSize-1)

        if iAx == 1
            title(ax, TASKTITLE(iTask))
        end

        ax.XGrid = 'on';
        hold(ax, 'off')
        fontsize(ax, p.fontSize, 'points')
    end

    % 6d. Any bout
    ax = nexttile(layout.left.bottom.tl); AX(iAx, 1) = ax;
    X = eta.lickBoutNorm.X(iEu, :);
    X = smoothdata(X, 2, 'gaussian', 5);
    plot(ax, eta.lickBoutNorm.t, mean(X, 1, 'omitnan'), Color=colors(ICOLOR(iAx), :), LineWidth=1.5)
    hold(ax, 'on')
    plot(ax, eta.lickBoutNorm.t, X, Color=[0.15, 0.15, 0.15, 4./nnz(iEu)], LineWidth=0.5)
    text(ax, 0.05, 0.025, sprintf('n=%i', nnz(iEu)), Unit='normalized', HorizontalAlignment='left', VerticalAlignment='bottom', Interpreter='none', FontSize=p.fontSize)
    xticks(ax, 0:2*pi:8*pi)
    xticklabels(ax, [{'0'}, arrayfun(@(x) sprintf('%i\\pi', x), 2:2:8, UniformOutput=false)]);
    xlim(ax, [0, 8*pi])
    ylim(ax, [-4, 4])
    yticks(ax, [-3, 0, 3])
    ax.XGrid = 'on';
end
xlabel(layout.right.bottom.tl, 'Time from bar/spout contact (s) & lick phase', FontSize=p.fontSize);
ylabel(layout.right.bottom.tl, 'Normalized spike rate (a.u.)', FontSize=p.fontSize)
xlabel(layout.left.bottom.tl, 'Lick phase', FontSize=p.fontSize)
ylabel(layout.left.bottom.tl, 'Normalized spike rate (a.u.)', FontSize=p.fontSize)

ax = AX(1, 1);
hLetter = text(ax, 0, 0, 'd', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.35, ax.Position(4) + 0.25, 0];

ax = AX(1, 2);
hLetter = text(ax, 0, 0, 'g', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.3, ax.Position(4) + 1.0, 0];

% S5f. Reach vs. Lick (scatter)
AX = gobjects(2, 1);
sz = 7;
ax = nexttile(layout.right.top.tl); AX(1) = ax;
hold(ax, 'on')
sel = c.hasLick & c.hasPress;
x = meta.lick;
y = meta.press;
subselResp = sel & c.isLick;
subselNone = sel & ~c.isLick;
h = gobjects(2, 1);
h(2) = scatter(ax, x(subselNone), y(subselNone), sz, 'black', 'filled', Marker='o', MarkerFaceAlpha=0.5, MarkerEdgeAlpha=0.5, DisplayName=sprintf('others (%i)', nnz(subselNone)));   
h(1) = scatter(ax, x(subselResp), y(subselResp), sz, [0, 0.5, 0.5], 'filled', Marker='o', MarkerFaceAlpha=0.5, MarkerEdgeAlpha=0.5, DisplayName=sprintf('lick-entrained (%i)', nnz(subselResp)));

plot(ax, [-10, 10], [0, 0], 'k:');
plot(ax, [0, 0], [-10, 10], 'k:');
plot(ax, [-10, 10], [-10, 10], 'k:')


% S5e right, Additional plot, scatter press vs lick META, color by lick entrainment phase: 
sz = 7;
ax = nexttile(layout.right.top.tl); AX(2) = ax;
hold(ax, 'on')
x = meta.lick;
y = meta.press;
sel = c.hasLick & c.hasPress & ~c.isLick;
h = gobjects(4, 1);
for i = 1:4
    sel = ID{i};
    h(i) = scatter(ax, x(sel), y(sel), sz, colors(ICOLOR(i), :), 'filled', Marker='o', MarkerFaceAlpha=0.75, MarkerEdgeAlpha=1, DisplayName=PHASENAME(i));
end

plot(ax, [-10, 10], [0, 0], 'k:');
plot(ax, [0, 0], [-10, 10], 'k:');
plot(ax, [-10, 10], [-10, 10], 'k:')

axis(AX, 'equal')
xlim(AX, [-2, 5])
ylim(AX, [-2, 5])

fontsize(AX, p.fontSize, 'points')
fontname(AX, 'Arial')

ylabel(layout.right.top.tl, 'Peri-reach activity (a.u.)', FontSize=p.fontSize)
xlabel(AX, 'Peri-lick activity (a.u.)', FontSize=p.fontSize)

ax = AX(1);
hLetter = text(ax, 0, 0, 'e', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.3, ax.Position(4) + 0.3, 0];

ax = AX(2);
hLetter = text(ax, 0, 0, 'f', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.25, ax.Position(4) + 0.3, 0];

% colorbar(axc, Location='eastoutside')

copygraphics(fig, ContentType='vector', BackgroundColor='none')
