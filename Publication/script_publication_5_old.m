%%% Figure 6. Lick vs. Reach ephys

%% Load all units
load_ephysunits;
find_osci_lick_circ;
load('C:\SERVER\Units\lda_pressVsLick_20241216.mat')

%%
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


eta.correctLickBoutInSeconds = eu.getETA('count', 'lick+lickbout_in_seconds', window=[-4, 3], resolution=[0.025, 2*pi/5], normalize='none', minTrialDuration=4, ...
    maxInterval=0.2, minInterval=0.05, minBoutCycles=2, maxBoutCycles=6);
eta.correctLickBoutInSecondsNorm = eu.getETA('count', 'lick+lickbout_in_seconds', window=[-4, 3], resolution=[0.025, 2*pi/5], normalize=[-4, -2], minTrialDuration=4, ...
    maxInterval=0.2, minInterval=0.05, minBoutCycles=2, maxBoutCycles=6);

eta.correctPressBoutInSeconds = eu.getETA('count', 'press+lickbout_in_seconds', window=[-4, 3], resolution=[0.025, 2*pi/32, 2*pi/5], normalize='none', minTrialDuration=4, ...
    maxInterval=0.2, minInterval=0.05, minBoutCycles=2, maxBoutCycles=6);
eta.correctPressBoutInSecondsNorm = eu.getETA('count', 'press+lickbout_in_seconds', window=[-4, 3], resolution=[0.025, 2*pi/32, 2*pi/5], normalize=[-4, -2], minTrialDuration=4, ...
    maxInterval=0.2, minInterval=0.05, minBoutCycles=2, maxBoutCycles=6);

% eta.incorrectLickBout = eu.getETA('count', 'lick+lickbout', window=[-4, 0], resolution=[0.025, 2*pi/5], normalize='none', minTrialDuration=2, maxTrialDuration=4, ...
%     maxInterval=0.2, minInterval=0.05, minBoutCycles=2, maxBoutCycles=6);

% eta.pressCueRaw = eu.getETA('count', 'press', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize='none', includeInvalid=true, resolution=0.025);
% eta.lickCueRaw = eu.getETA('count', 'lick', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize='none', includeInvalid=true, resolution=0.025);
% 
% eta.pressCueNorm = eu.getETA('count', 'press', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize=etaFine.press.stats, includeInvalid=true, resolution=0.025);
% eta.lickCueNorm = eu.getETA('count', 'lick', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize=etaFine.lick.stats, includeInvalid=true, resolution=0.025);

%%
% [boot.pressVsLick.h, boot.pressVsLick.p, boot.pressVsLick.ci, boot.pressVsLick.obs] = bootstrapAmplitude(eu(c.hasLick & c.hasPress), 'press', 'lick', ...
%     responseWindowA=p.metaWindowPress, responseWindowB=p.metaWindowLick, allowedTrialDuration=[2, Inf], withReplacement=false, alpha=0.01);
% c.isPressVsLickSelective = false(1, length(eu));
% c.isPressVsLickSelective(c.hasLick & c.hasPress) = boot.pressVsLick.h;
save('C:\SERVER\Units\meta_Lite_NonDuplicate_NonDrift.mat', 'ai', 'boot', 'c', 'eta', 'etaSmooth', 'euPos', 'meta', 'msr', 'onset', 'p', 'bouts', 'trialsCircLick', 'trialsCircLickBaseline', 'bootCLick', 'lickHist', 'durations')

%% Fig 6
close all

nEgUnits = 3;
p.fontSize = 9;
p.etaSortWindow = [-3, 0.1];
p.etaSignWindow = [-0.2, 0.1];
p.etaLatencyThresholdPos = 0.25;
p.etaLatencyThresholdNeg = 0.25;

layout.w = 7;
layout.h = 7;

layout.top.h = 9;
layout.middle.h = 9;
layout.bottom.h = 7;

layout.middle.left.w = 4;
layout.middle.right.w = 4;
layout.middle.rightMargin.w = 1;

layout.bottom.left.w = 6;
layout.bottom.middle.w = 8;
layout.bottom.right.w = 6;

% close all
fig = figure(Units='inches', Position=[0, 0, layout.w, layout.h]);
layout.tl = tiledlayout(fig, layout.top.h + layout.middle.h + layout.bottom.h, 1, TileSpacing='loose');

layout.top.tl = tiledlayout(layout.tl, 2, nEgUnits, TileSpacing='tight', Padding='tight');
l = layout.top.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [layout.top.h, 1];


layout.middle.tl = tiledlayout(layout.tl, 1, layout.middle.left.w + layout.middle.right.w + layout.middle.rightMargin.w, TileSpacing='tight', Padding='none');
l = layout.middle.tl; l.Layout.Tile = 1 + layout.top.h; l.Layout.TileSpan = [layout.middle.h, 1];

layout.middle.left.tl = tiledlayout(layout.middle.tl, 1, 2, TileSpacing='compact', Padding='compact');
l = layout.middle.left.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [1, layout.middle.left.w];

layout.middle.right.tl = tiledlayout(layout.middle.tl, 1, 2, TileSpacing='compact', Padding='compact');
l = layout.middle.right.tl; l.Layout.Tile = 1 + layout.middle.left.w; l.Layout.TileSpan = [1, layout.middle.right.w];

layout.bottom.tl = tiledlayout(layout.tl, 1, layout.bottom.left.w + layout.bottom.middle.w + layout.bottom.right.w, TileSpacing='tight', Padding='tight');
l = layout.bottom.tl; l.Layout.Tile = 1 + layout.top.h + layout.middle.h; l.Layout.TileSpan = [layout.bottom.h, 1];

layout.bottom.left.tl = tiledlayout(layout.bottom.tl, 1, 1, TileSpacing='tight', Padding='tight');
l = layout.bottom.left.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [1, layout.bottom.left.w];

layout.bottom.middle.tl = tiledlayout(layout.bottom.tl, 1, 1, TileSpacing='tight', Padding='tight');
l = layout.bottom.middle.tl; l.Layout.Tile = 1 + layout.bottom.left.w; l.Layout.TileSpan = [1, layout.bottom.middle.w];

layout.bottom.right.tl = tiledlayout(layout.bottom.tl, 1, 1, TileSpacing='tight', Padding='tight');
l = layout.bottom.right.tl; l.Layout.Tile = 1 + layout.bottom.left.w + layout.bottom.middle.w; l.Layout.TileSpan = [1, layout.bottom.right.w];

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
        thisETA = euEg(iEu).getETA('count', TRIALTYPE{i}, [-4, 2], normalize='none', trials=theseTrials, includeInvalid=false);
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


% sel = c.hasPress & c.hasLick;
sel = 1:length(eu);

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

for iAx = 1:4
    applyCustomColormap(ax(iAx), [-1.5, 1.5], hlim=[0.375, 0, 0, -0.375], llim=[0.125, 0.5, 0.5, 0.25], hpwr=.5, lpwr=1, h0=0.33);
end

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

% SEL for 6d and 6e
lims = [-3, 5];
selCommon = c.hasPress & c.hasLick;
SEL = {selCommon & c.isPressDown & c.isLickUp, selCommon & c.isPressDown & c.isLickDown, selCommon & c.isPressUp & c.isLickUp, selCommon & c.isPressUp & c.isLickDown};
LABELS = ["reach-dec lick-inc", "reach-dec lick-dec", "reach-inc lick-inc", "reach-inc lick-dec"];
COLORS = cellfun(@(hsl) hsl2rgb(hsl), {[0/36, 1, 0.5], [27/36, 1, 0.33], [12/36, 0.9, 0.33], [24/36, 1, 0.5]}, UniformOutput=false);
TEXTPOS = {lims([2, 1]) + [-0.1, 0.1], lims([1, 1]) + [0.1, 0.1], lims([2, 2]) + [-0.1, -0.1], lims([1, 2]) + [0.1, -0.1]};
HORZALIGN = ["right", "left", "right", "left"];
VERTALIGN = ["bottom", "bottom", "top", "top"];

% 6d Scatter META for lick vs reach
sz = 4;
ax = nexttile(layout.bottom.left.tl);
selCommon = c.hasPress & c.hasLick;
subselSign = selCommon & (c.isPressResponsive & c.isLickResponsive);
subselAmp = selCommon & c.isPressVsLickSelective;
subselNone = selCommon & ~subselSign;
scatter(ax, meta.lick(subselNone), meta.press(subselNone), sz, 'black', MarkerEdgeAlpha=0.125), hold(ax, 'on')
for iScat = 1:length(SEL)
    scatter(ax, meta.lick(SEL{iScat}), meta.press(SEL{iScat}), sz, COLORS{iScat}, 'filled', MarkerFaceAlpha=0.67)
end
% scatter(ax, meta.lick(subselSign), meta.press(subselSign), sz, 'black', 'filled', MarkerEdgeAlpha=0, MarkerFaceAlpha=1)
% scatter(ax, meta.lick(subselAmp), meta.press(subselAmp), sz, 'red', MarkerEdgeAlpha=0.5, LineWidth=0.7)
xline(ax, 0, 'k:')
yline(ax, 0, 'k:')
plot(ax, lims, lims, 'k:')
axis(ax, 'equal')
xlim(ax, lims)
ylim(ax, lims)
xlabel(ax, 'Lick (a.u.)')
ylabel(ax, 'Reach (a.u.)')
fontsize(ax, p.fontSize, 'points')

for iScat = 1:length(SEL)
    text(ax, TEXTPOS{iScat}(1), TEXTPOS{iScat}(2), sprintf('n=%i', nnz(SEL{iScat})), Color=COLORS{iScat}, ...
        HorizontalAlignment=HORZALIGN(iScat), VerticalAlignment=VERTALIGN(iScat), FontSize=p.fontSize-1)
end

mdl = fitlm(meta.lick(selCommon), meta.press(selCommon));
fprintf('press vs. lick: %i total, %i sign-change, %i amplitude change (LM slope p<%g).\n', nnz(selCommon), nnz(subselSign), nnz(subselAmp), mdl.Coefficients.pValue(2))

hLetter = text(ax, 0, 0, 'd', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.4, ax.Position(4) + 0.25, 0];
ax.Box = 'on';

% 6e. Correct trial PETHs
yl = [-1.2, 2.2];
ETA = eta.correctPressBoutNorm;
% ETA = eta.correctPressBoutInSecondsNorm;
% COLORS = {'b', 'r'};


% Peri-move
ax = nexttile(layout.bottom.middle.tl);
ax.Box = 'on';
hold(ax, 'on')
t = ETA.t;
t(t>0 & t<=2*pi) = t(t>0 & t<=2*pi) / (2*pi) * 0.8;
t(t>2*pi) = (t(t>2*pi) - 2*pi) ./ (2*pi) / 8 + 0.8;
tt = [t, flip(t)];

for iLine = 1:4
    sel = SEL{iLine};
    X = ETA.X(sel, :);
    err = std(X, 0, 1, 'omitnan') ./ sqrt(size(X, 1));
    xx = [mean(X, 1, 'omitnan') + err, flip(mean(X, 1, 'omitnan') - err)];

    plot(ax, t, mean(X, 1, 'omitnan'), Color=[COLORS{iLine}, 0.7], LineWidth=1.5);
%     plot(ax, t, X, Color=[COLORS{iLine}, LINEALPHA(iLine)], LineWidth=0.5)
    patch(ax, tt(~isnan(xx)), xx(~isnan(xx)), COLORS{iLine}, FaceAlpha=0.25, EdgeColor='none')
end

xline(ax, [0.8,  0.8+(1:5)/8], LineStyle=':')
xlim(ax, [-2, 0.8+5/8])
ylim(ax, yl)
xticks(ax, [-2, -1, 0, 0.8,  0.8+(1:5)/8])
xticklabels(ax, {'-2', '-1', '0', '1^s^t', '', '', '', '', '6^t^h'})
yticks(ax, 0:40:80)

ax.XAxis.TickLabelRotation = 0;
xline(ax, 0, '-')

% xlabel(layout.bottom.middle.tl, 'Time to bar contact (s)    lick phase', FontSize=p.fontSize)
text(ax, -1, yl(1)-0.75, 'Time to bar contact (s)', HorizontalAlignment='center', VerticalAlignment='top', Color=[0.15, 0.15, 0.15])
text(ax, 0.8+3/8, yl(1)-0.75, 'Lick', HorizontalAlignment='center', VerticalAlignment='top', Color=[0.15, 0.15, 0.15])
ylabel(ax, {'Normalized','spike rate (a.u.)'}, FontSize=p.fontSize)
fontsize(ax, p.fontSize, 'points')

hLetter = text(ax, 0, 0, 'e', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.2, ax.Position(4) + 0.3, 0];

% lgd = legend(hLine, Orientation='horizontal');
% lgd.Layout.Tile = 'north';

% 6f Decoder LDA for lick vs reach
t = likelihood(1).t;
pressTrialPress = arrayfun(@(llh) llh.press(llh.trueLabel=="press", :), likelihood, UniformOutput=false); % Press likelihood for true-press trials
pressTrialLick = arrayfun(@(llh) llh.lick(llh.trueLabel=="press", :), likelihood, UniformOutput=false); % Lick likelihood for true-press trials
lickTrialPress = arrayfun(@(llh) llh.press(llh.trueLabel=="lick", :), likelihood, UniformOutput=false);
lickTrialLick = arrayfun(@(llh) llh.lick(llh.trueLabel=="lick", :), likelihood, UniformOutput=false);
pressTrialPress = cat(1, pressTrialPress{:});
pressTrialLick = cat(1, pressTrialLick{:});
lickTrialPress = cat(1, lickTrialPress{:});
lickTrialLick = cat(1, lickTrialLick{:});   

DATA = { ...
    pressTrialPress, pressTrialLick; ...
    lickTrialPress, lickTrialLick ...
    };
COLOR = ["red", "blue"];
LABEL = ["reach", "lick"];

ax = nexttile(layout.bottom.right.tl);
ax.Box = 'on';
hold(ax, 'on')
h = gobjects(3, 1);
for iMove = 1:2
    mu = mean(DATA{iMove, 1} - DATA{iMove, 2}, 1, 'omitnan');
    h(iMove) = plot(ax, t, mu, COLOR(iMove), LineWidth=1.5, DisplayName=sprintf('%s', LABEL(iMove)));
end
% h(3) = plot(ax, t, dfBootStats.all.mu, 'black', LineStyle='--', LineWidth=1.5, DisplayName='shuffle');
patch(ax, [t, flip(t)], [dfBootStats.all.ci(1, :), flip(dfBootStats.all.ci(2, :))], 'black', FaceAlpha=0.1, EdgeAlpha=0.5, DisplayName='99% CI');
patch(ax, [pLDA.responseWindow, flip(pLDA.responseWindow)], [-1, -1, 1, 1], 'yellow', FaceAlpha=0.1, EdgeAlpha=0.5, DisplayName='training');
ylim(ax, [-1, 1])
xlim(ax, [-2.5, 1.5])
xline(ax, 0, ':')
yline(ax, 0, ':')
hold(ax, 'off')

xlabel(ax, 'Time to contact (s)')
ylabel(ax, 'p(reach) - p(lick)')

fontsize(ax, p.fontSize, 'points')

% Manual legends
text(ax, -2.4, 0.8, sprintf('(%i trials)\nreach', size(pressTrialPress, 1)), Color='r', HorizontalAlignment='left', VerticalAlignment='top', FontSize=p.fontSize-1)
text(ax, -2.4, -0.8, sprintf('lick\n(%i trials)', size(lickTrialLick, 1)), Color='b', HorizontalAlignment='left', VerticalAlignment='bottom', FontSize=p.fontSize-1)
% lgd = legend(ax, h(1:2), Orientation='horizontal', Location='northoutside', FontSize=p.fontSize-1, AutoUpdate=false);

nAnimals = length(unique(eu(ismember(string({eu.ExpName}), goodExpNames)).getAnimalName()));

nUnits = cellfun(@(x) size(x, 2), {resp.press}, UniformOutput=true);
fprintf('\nReach vs. lick decoder:\n');
fprintf('\t1) Selected %i sessions (%i animals) with >=%i units (mean=%g, total=%i).\n', length(goodExpNames), nAnimals, pLDA.minNumUnits, mean(nUnits), sum(nUnits))
fprintf('\t\t Units per session (sorted): %s\n', num2str(sort(nUnits)));

hLetter = text(ax, 0, 0, 'f', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.35, ax.Position(4) + 0.25, 0];

% Colorbar for 6b
h = colorbar(axc); 
h.Layout.Tile = 'east';
h.Label.String = 'Normalized spike rate (a.u.)';

% lgd.Position = [lgd.Position(1), lgd.Position(2) + 0.05, lgd.Position(3:4)];
  
copygraphics(fig, ContentType='vector', BackgroundColor='none')
