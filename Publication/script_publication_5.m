%%% Figure 6. Lick vs. Reach ephys + Osci lick and somatotopy

%% Load all units
load_ephysunits;
% find_osci_lick_circ;
% [boot.pressVsLick.h, boot.pressVsLick.p, boot.pressVsLick.ci, boot.pressVsLick.obs] = bootstrapAmplitude(eu(c.hasLick & c.hasPress), 'press', 'lick', ...
%     responseWindowA=p.metaWindowPress, responseWindowB=p.metaWindowLick, allowedTrialDuration=[2, Inf], withReplacement=false, alpha=0.01);
% c.isPressVsLickSelective = false(1, length(eu));
% c.isPressVsLickSelective(c.hasLick & c.hasPress) = boot.pressVsLick.h;
% save('C:\SERVER\Units\meta_Lite_NonDuplicate_NonDrift.mat', 'ai', 'boot', 'c', 'eta', 'etaSmooth', 'euPos', 'meta', 'msr', 'onset', 'p', 'bouts', 'trialsCircLick', 'trialsCircLickBaseline', 'bootCLick', 'lickHist', 'durations')

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
layout.bottom.h = 5;
layout.bottom.left.w = 4;
layout.bottom.right.w = 2;
layout.bottom.left.top.h = 3;
layout.bottom.left.bottom.h = 3;
layout.bottom.right.top.h = 8;
layout.bottom.right.middle.h = 3;
layout.bottom.right.bottom.h = 8;

close all
fig = figure(Units='inches', Position=[0, 0, layout.w, layout.h]);
layout.tl = tiledlayout(fig, layout.top.h + layout.bottom.h, 1, TileSpacing='loose');

layout.top.tl = tiledlayout(layout.tl, 2, nEgUnits, TileSpacing='tight', Padding='tight');
l = layout.top.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [layout.top.h, 1];

layout.bottom.tl = tiledlayout(layout.tl, 1, layout.bottom.left.w + layout.bottom.right.w, TileSpacing='compact', Padding='compact');
l = layout.bottom.tl; l.Layout.Tile = 1 + layout.top.h; l.Layout.TileSpan = [layout.bottom.h, 1];

layout.bottom.left.tl = tiledlayout(layout.bottom.tl, layout.bottom.left.top.h + layout.bottom.left.bottom.h, 1, TileSpacing='compact', Padding='compact');
l = layout.bottom.left.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [1, layout.bottom.left.w];

layout.bottom.left.top.tl = tiledlayout(layout.bottom.left.tl, 1, 2, TileSpacing='compact', Padding='compact');
l = layout.bottom.left.top.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [layout.bottom.left.top.h, 1];

layout.bottom.left.bottom.tl = tiledlayout(layout.bottom.left.tl, 3, 3, TileSpacing='tight', Padding='compact');
l = layout.bottom.left.bottom.tl; l.Layout.Tile = 1 + layout.bottom.left.top.h; l.Layout.TileSpan = [layout.bottom.left.bottom.h, 1];

layout.bottom.right.tl = tiledlayout(layout.bottom.tl, layout.bottom.right.top.h + layout.bottom.right.middle.h + layout.bottom.right.bottom.h, 1, TileSpacing='compact', Padding='compact');
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
ax(1) = nexttile(layout.bottom.left.top.tl);
ax(2) = nexttile(layout.bottom.left.top.tl);
ax(3) = nexttile(layout.bottom.right.tl, [layout.bottom.right.top.h, 1]);
EphysUnit.plotETA(ax(1), eta.press, c.hasPress & c.hasLick, ...
    clim=[-2, 2], xlim=[-4, 0.5], sortWindow=p.etaSortWindow, signWindow=p.etaSignWindow, ...
    sortThreshold=p.etaLatencyThresholdPos, negativeSortThreshold=p.etaLatencyThresholdNeg, hidecolorbar=true);
EphysUnit.plotETA(ax(2), eta.lick, c.hasPress & c.hasLick, ...
    clim=[-2, 2], xlim=[-4, 0.5], sortWindow=p.etaSortWindow, signWindow=p.etaSignWindow, ...
    sortThreshold=p.etaLatencyThresholdPos, negativeSortThreshold=p.etaLatencyThresholdNeg, hidecolorbar=true);
hold(ax(1:2), 'on')
plot(ax(1), [0, 0], [0, nnz(c.hasPress & c.hasLick)], 'k--')
plot(ax(2), [0, 0], [0, nnz(c.hasPress & c.hasLick)], 'k--')
ylim(ax(1:2), [0, nnz(c.hasPress & c.hasLick)+1])
yt = 0:100:nnz(c.hasPress & c.hasLick);
yt(1) = 1;
if round(yt(end)./100) == round(nnz(c.hasPress & c.hasLick)./100)
    yt(end) = nnz(c.hasPress & c.hasLick);
else
    yt(end + 1) = nnz(c.hasPress & c.hasLick);
end
yt = unique(yt);
yticks(ax(1:2), yt)

% Visualize onset
% hold(ax(1:2), 'on')
% plot(ax(1), onset.press(order.press), 1:length(order.press))
% plot(ax(2), onset.lick(order.lick), 1:length(order.lick))

% meanLickInterval = mean(cat(2, durations{:}));
% meanBinWidth = meanLickInterval / length(eta.circLick.t);
% eta.circLick.Z = eta.circLick.X.*exp(eta.circLick.t*1i)./meanBinWidth;
% eta.circLick.Z(:, 1) = mean(eta.circLick.X(:, [2, 30]), 2).*exp(eta.circLick.t(1)*1i)./meanBinWidth;
eta.circLick.Z = eta.circLick.X.*exp(eta.circLick.t*1i);
eta.circLick.Z(:, 1) = mean(eta.circLick.X(:, [2, 30]), 2).*exp(eta.circLick.t(1)*1i);
meanZ = mean(eta.circLick.Z, 2);
eta.lickBoutNorm = eta.lickBout;
eta.lickBoutNorm.X = normalize(eta.lickBout.X, 2, 'zscore', 'robust');
% mu = {eta.lick.stats.mean}';
% sel = cellfun(@isempty, mu);
% mu(sel) = {NaN};
% mu = cellfun(@(x) x, mu);
% sd = {eta.lick.stats.sd}';
% sel = cellfun(@isempty, sd);
% sd(sel) = {NaN};
% sd = cellfun(@(x) x, sd);
% eta.lickBoutNorm.X = (eta.lickBout.X./(1/(9*2*15)) - mu./0.1)./(sd./0.1);

maxBoutCycles = 4;
sel = c.hasPress & c.hasLick & c.isLick; 
phase = angle(meanZ(sel));
phase(phase < 0) = phase(phase < 0) + 2*pi;
[sortedPhase, I] = sort(phase);
[~, ~] = EphysUnit.plotETA(ax(3), eta.lickBoutNorm, sel, order=I, ...
    clim=[-2, 2], xlim=[0, 2*pi*maxBoutCycles], hidecolorbar=true);
xticks (ax(3), (0:2:8).*pi);
xticklabels(ax(3), [{'0'}, arrayfun(@(x) sprintf('%i\\pi', x), 2:2:8, UniformOutput=false)]);
title(ax(1), 'Pre-reach')
title(ax(2), 'Pre-lick')
title(ax(3), 'Lick-entrained')
ylabel(layout.bottom.left.top.tl, 'Unit', FontSize=p.fontSize)
ylabel(ax(3), 'Unit')
ylabel(ax(1:2), '')
xlabel(ax(1), 'Time to bar contact (s)')
xlabel(ax(2), 'Time to spout contact (s)')
xlabel(ax(3), 'Lick phase')
xlim(ax(3), [0, 2*pi*maxBoutCycles])
ylim(ax(3), [0, nnz(sel)+1])
yt = 0:100:nnz(sel);
yt(1) = 1;
if round(yt(end)./100) == round(nnz(sel)./100)
    yt(end) = nnz(sel);
else
    yt(end + 1) = nnz(sel);
end
yt = unique(yt);
yticks(ax(3), yt)
fontsize(ax, p.fontSize, 'points')
fontname(ax, 'Arial')
axc = ax(3);

hLetter = text(ax(1), 0, 0, 'b', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax(1).Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.3, ax(1).Position(4) + 0.1, 0];

hLetter = text(ax(3), 0, 0, 'd', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax(3).Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.2, ax(3).Position(4) + 0.2, 0];

clear yt

%5d. Peri-lick lick prob histogram
ax = nexttile(layout.bottom.right.tl, [layout.bottom.right.middle.h, 1]);
histogram(ax, BinEdges=lickHist.firstLick.edges, BinCounts=lickHist.firstLick.count, Normalization='probability', EdgeColor='none', FaceColor='black', FaceAlpha=1)
xlim(ax, [0, 1/9*4]); % Assuming 9 Hz, 5 licks (4 cycles)
xticks(ax, 0:0.2:0.5)
yticks(ax, [])
xlabel(ax, 'Time from first lick (s)')
ylabel(ax, ' lick\newlineprob', HorizontalAlignment='center')
fontsize(ax, p.fontSize, 'points')
hold(ax, 'off')
hLetter = text(ax, 0, 0, 'e', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.25, ax.Position(4) + 0.25, 0];

%5d Mean ETA grouped by press/reach up/down
SEL_HASTRIAL = c.hasPress & c.hasLick;
SEL_I = {c.isPressUp; ~c.isPressResponsive; c.isPressDown};
SEL_J = {c.isLickDown; ~c.isLickResponsive; c.isLickUp};

h = gobjects(2, 1);
for i = 1:length(SEL_I)
    for j = 1:length(SEL_J)
        ax = nexttile(layout.bottom.left.bottom.tl);
        sel = SEL_HASTRIAL & SEL_I{i} & SEL_J{j};
        assert(numel(sel) == numel(eu))
        hold(ax, 'on')
        h(1) = plot(ax, eta.press.t, mean(eta.press.X(sel, :), 1, 'omitnan'), 'r', LineWidth=1.5, DisplayName='reach');
        h(2) = plot(ax, eta.lick.t, mean(eta.lick.X(sel, :), 1, 'omitnan'), 'b', LineWidth=1.5, DisplayName='lick');
        text(ax, -3.8, 2, sprintf('n = %i', nnz(sel)), fontSize=p.fontSize, VerticalAlignment='top')
        hold(ax, 'off')
        xlim(ax, [-4, 0.5])
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
            hLetter = text(ax, 0, 0, 'c', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
            ax.Units = 'inches';
            hLetter.HorizontalAlignment = 'right';
            hLetter.VerticalAlignment = 'top';
            hLetter.Position = [-0.3, ax.Position(4) + 0.3, 0];
        end
    end
end
lgd = legend(h, Orientation='horizontal');
lgd.Layout.Tile = 'north';
xlabel(layout.bottom.left.bottom.tl, 'Time to bar/spout contact (s)', FontSize=p.fontSize)
ylabel(layout.bottom.left.bottom.tl, 'Normalized spike rate (a.u.)', FontSize=p.fontSize)

% 5e. Reach vs. Lick (scatter), Reach vs. Lick (osci subset, scatter), Pie-chart
sz = 7;
AX = gobjects(1, 2);
% 1. Reach vs. Lick scatter
ax = nexttile(layout.bottom.right.tl, [layout.bottom.right.bottom.h, 1]);
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

xlabel(ax, 'Peri-lick activity (a.u.)')
ylabel(ax, 'Peri-reach activity (a.u.)')
% axis(ax, 'equal')
xlim(ax, [-2, 5])
ylim(ax, [-2, 5])
% legend(h, Location='north')
fontsize(ax, p.fontSize, 'points')

hLetter = text(ax, 0, 0, 'f', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.25, ax.Position(4) + 0.25, 0];


h = colorbar(axc); 
h.Location = 'eastoutside';
% h.Layout.Tile = 'south';
h.Label.String = 'normalized spike rate (a.u.)';
copygraphics(fig, ContentType='vector')


%% Supplement S5
% All sign-changing units
%% Plot ETA from peri-move to cyclick
% eta.correctLickBout = eu.getETA('count', 'lick+lickbout', window=[-4, 0], resolution=[0.01, 2*pi/12], normalize='none', minTrialDuration=4, ...
%     maxInterval=0.2, minInterval=0.05, minBoutCycles=minBoutCycles, maxBoutCycles=maxBoutCycles);
% eta.correctLickBoutNorm = eu.getETA('count', 'lick+lickbout', window=[-4, 0], resolution=[0.01, 2*pi/12], normalize=[-4, -2], minTrialDuration=4, ...
%     maxInterval=0.2, minInterval=0.05, minBoutCycles=minBoutCycles, maxBoutCycles=maxBoutCycles);
%% Supplement S5
minBoutCycles = 4;
maxBoutCycles = 16;
close all
p.fontSize = 9;
clear layout
layout.w = 6;
layout.h = 10;
layout.top.h = 8;
layout.middle.h = 2;
layout.bottom.h = 8;
layout.left.w = 4;
layout.right.w = 2;
layout.psbottom.h = 3;

fig = figure(Units='inches', Position=[0.2 0.2 layout.w layout.h]); 
layout.tl = tiledlayout(fig, layout.top.h + layout.middle.h + layout.bottom.h + layout.psbottom.h, 1, TileSpacing='loose', Padding='loose');

nRows = 6;
nCols = 6;
layout.top.tl = tiledlayout(layout.tl, nRows, nCols, TileSpacing='tight', Padding='loose');
l = layout.top.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [layout.top.h, 1];

layout.middle.tl = tiledlayout(layout.tl, 1, layout.left.w + layout.right.w, TileSpacing='loose', Padding='loose');
l = layout.middle.tl; l.Layout.Tile = 1 + layout.top.h; l.Layout.TileSpan = [layout.middle.h, 1];

layout.bottom.tl = tiledlayout(layout.tl, 1, layout.left.w + layout.right.w, TileSpacing='tight', Padding='tight');
l = layout.bottom.tl; l.Layout.Tile = 1 + layout.top.h + layout.middle.h; l.Layout.TileSpan = [layout.bottom.h, 1];

layout.bottom.left.tl = tiledlayout(layout.bottom.tl, 4, 1, TileSpacing='compact', Padding='tight');
l = layout.bottom.left.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [1, layout.left.w];

layout.bottom.right.tl = tiledlayout(layout.bottom.tl, 4, 1, TileSpacing='compact', Padding='tight');
l = layout.bottom.right.tl; l.Layout.Tile = 1 + layout.left.w; l.Layout.TileSpan = [1, layout.right.w];

layout.psbottom.tl = tiledlayout(layout.tl, 1, 1, TileSpacing='tight', Padding='tight');
l = layout.psbottom.tl; l.Layout.Tile = 1 + layout.top.h + layout.middle.h + layout.bottom.h; l.Layout.TileSpan = [layout.psbottom.h, 1];


signChangingUnitIndices = [find(c.isPressUp & c.isLickDown), find(c.isPressDown & c.isLickUp)];
assert(length(signChangingUnitIndices) <= nRows*nCols)
iAx = 0;
for i = signChangingUnitIndices
    ax = nexttile(layout.top.tl);
    iAx = iAx + 1;
    hold(ax, 'on')
    plot(ax, eta.pressRaw.t, eta.pressRaw.X(i, :)./0.1, 'r', LineWidth=1.5, DisplayName='reach')
    plot(ax, eta.lickRaw.t, eta.lickRaw.X(i, :)./0.1, 'b', LineWidth=1.5, DisplayName='lick')
%     plot(ax, [0, 0], [0, 100], 'k:')
    hold(ax, 'off')
    xlim(ax, [-4, 0.5])
    ymedian = median([mean(eta.pressRaw.X(i, :), 1, 'omitnan')./0.1, mean(eta.lickRaw.X(i, :), 1, 'omitnan')./0.1]);
    ylim(ax, ymedian + [-40, 40])
    yticks(ax, round(ymedian/10)*10 + [-30, 0, 30])
    ax.YGrid = 'on';
    ax.XGrid = 'on';
    fontsize(ax, p.fontSize, 'points')
end
xlabel(layout.top.tl, 'Time to bar/spout contact (s)', FontSize=p.fontSize)
ylabel(layout.top.tl, 'Spike rate (sp/s)', FontSize=p.fontSize)
lgd = legend(ax, Orientation='horizontal');
lgd.Layout.Tile = 'north';

nBoutsDisp = 6;
sel = find(c.isLick & c.hasPress & c.hasLick);
amp = abs(meanZ(sel));
phase = angle(meanZ(sel));
rectifiedPhase = phase;
rectifiedPhase(phase < 0) = rectifiedPhase(phase < 0) + 2*pi;
ampThreshold = 0;
isInPhase = abs(phase) <= 0.25*pi;
isAntiPhase = abs(phase) >= 0.75*pi;
isFirstQuarterPhase = phase > 0.25*pi & phase < 0.75*pi;
isThirdQuarterPhase = phase > -0.75*pi & phase < -0.25*pi;
isHighAmp = amp >= quantile(amp, ampThreshold);
idInPhase = sel(isInPhase & isHighAmp);
idAntiPhase = sel(isAntiPhase & isHighAmp);
idFirstQuarterPhase = sel(isFirstQuarterPhase & isHighAmp);
idThirdQuarterPhase = sel(isThirdQuarterPhase & isHighAmp);
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

% clear trials expEuIndices FsLick lickHistEdges lickHistCenters lickHistNLicks lickHistCounts iExp iEu firstLickTimes allLickTimes iLick edgesGlobal n

% Correct lick bout
ax = nexttile(layout.middle.tl, [1, layout.left.w]);
counts = lickHist.correctLickOsci.count;
smoothedCounts = smoothdata(counts, 'gaussian', 25);
[pks, locs] = findpeaks(smoothedCounts, lickHist.correctLickOsci.t, MinPeakProminence=0.5);
locs = [0, locs];
histogram(ax, BinEdges=lickHist.correctLickOsci.edges, BinCounts=counts, Normalization='probability', EdgeColor='none', FaceColor='black', FaceAlpha=1)
hold(ax, 'on')
for iLoc = 1:length(locs)
    plot(ax, locs(iLoc)*[1, 1], ax.YLim, 'k:')
end
ylabel(ax, ' lick\newlineprob')
xlabel(ax, 'Time from first lick (s)')
xticks(ax, -1:0.5:2)
xlim(ax, [-0.5, 1/8*nBoutsDisp])
yticks(ax, [])
fontsize(ax, p.fontSize, 'points')

% Osci lick phase distribution histogram
ax = nexttile(layout.middle.tl, [1, layout.right.w]);
hold(ax, 'on')
edges = 0:2*pi/30:2*pi;
histogram(ax, rectifiedPhase(isInPhase), edges, FaceColor=colors(1, :));
histogram(ax, rectifiedPhase(isFirstQuarterPhase), edges, FaceColor=colors(2, :));
histogram(ax, rectifiedPhase(isAntiPhase), edges, FaceColor=colors(4, :));
histogram(ax, rectifiedPhase(isThirdQuarterPhase), edges, FaceColor=colors(3, :));
xticks(ax, 0:pi:2*pi)
xlim(ax, pi*[0, 2])
xticklabels(ax, {'0', '\pi', '2\pi'});
xlabel('Lick phase')
ylabel('# units')

ID = {idInPhase; idFirstQuarterPhase; idAntiPhase; idThirdQuarterPhase};
PHASENAME = ["2\pi", "1/2\pi", "\pi", "3/2\pi"];
ICOLOR = [1, 2, 4, 3];

for iAx = 1:4
    iEu = ID{iAx};
    % First lick
    ax = nexttile(layout.bottom.left.tl);
    t = eta.correctLickBout.t;
    t(t>0) = t(t>0) ./ (2*pi) / 8;
    X = eta.correctLickBout.X(iEu, :);
    X = smoothdata(X, 2, 'gaussian', 5);
    plot(ax, t, mean(X, 1, 'omitnan'), Color=colors(ICOLOR(iAx), :), LineWidth=1.5)
    hold(ax, 'on')
    plot(ax, t, X, Color=[0.15, 0.15, 0.15, 0.1], LineWidth=0.5)
    xlim(ax, [-0.5, nBoutsDisp/8])
    ylim(ax, [0, 120])
    text(ax, 0.05, 0.95, sprintf('n=%i', nnz(iEu)), Unit='normalized', HorizontalAlignment='left', VerticalAlignment='top', Interpreter='none', FontSize=p.fontSize)
    xticks(ax, [-1:0.5:0, (1:nBoutsDisp)/8])
    xticklabels(ax, [{'-1', '-0.5', '0'}, arrayfun(@(x) sprintf('%i\\pi', x), 2*(1:nBoutsDisp), UniformOutput=false)]);
    ax.XGrid = 'on';
    hold(ax, 'off')
    fontsize(ax, p.fontSize, 'points')
    if iAx == 4
        text(ax, 0, -0.42, 'Time before first lick (s)', Units='normalized', VerticalAlignment='top', FontSize=p.fontSize)
        text(ax, 0.55, -0.42, 'Phase after first lick', Units='normalized', VerticalAlignment='top', FontSize=p.fontSize)
    end

    % Any bout
    ax = nexttile(layout.bottom.right.tl);
    X = eta.lickBoutNorm.X(iEu, :);
    X = smoothdata(X, 2, 'gaussian', 3);
    plot(ax, eta.lickBoutNorm.t, mean(X, 1, 'omitnan'), Color=colors(ICOLOR(iAx), :), LineWidth=1.5)
    hold(ax, 'on')
    plot(ax, eta.lickBoutNorm.t, X, Color=[0.15, 0.15, 0.15, 0.025], LineWidth=0.5)
    xticks(ax, 0:2*pi:8*pi)
    xticklabels(ax, [{'0'}, arrayfun(@(x) sprintf('%i\\pi', x), 2:2:8, UniformOutput=false)]);
    xlim(ax, [0, 8*pi])
    ylim(ax, [-4, 4])
    ax.XGrid = 'on';
end
% h = xlabel(layout.bottom.left.tl, 'Time from first lick (s)       Lick Phase', FontSize=p.fontSize, HorizontalAlignment='left');
ylabel(layout.bottom.left.tl, 'Spike rate (sp/s)', FontSize=p.fontSize)
xlabel(layout.bottom.right.tl, 'Lick phase', FontSize=p.fontSize)
ylabel(layout.bottom.right.tl, 'Normalized spike rate (a.u.)', FontSize=p.fontSize)

% Additional plot, scatter press vs lick META, color by lick entrainment phase: 
sz = 3;
% 1. Reach vs. Lick scatter
ax = nexttile(layout.psbottom.tl, [1, 1]);
hold(ax, 'on')
x = meta.lick;
y = meta.press;
sel = c.hasLick & c.hasPress & ~c.isLick;
h = gobjects(4, 1);
for i = 1:4
    sel = ID{i};
    h(i) = scatter(ax, x(sel), y(sel), sz, colors(ICOLOR(i), :), 'filled', Marker='o', MarkerFaceAlpha=0.5, MarkerEdgeAlpha=0.75, DisplayName=PHASENAME(i));
end

plot(ax, [-10, 10], [0, 0], 'k:');
plot(ax, [0, 0], [-10, 10], 'k:');
plot(ax, [-10, 10], [-10, 10], 'k:')

xlabel(ax, '    Peri-lick\newlineresponse (a.u.)', HorizontalAlignment='center')
ylabel(ax, '    Peri-reach\newlineresponse (a.u.)', HorizontalAlignment='center')
axis(ax, 'equal')
xlim(ax, [-2, 5])
ylim(ax, [-2, 5])
legend(ax, h, Location='eastoutside')
fontsize(ax, p.fontSize, 'points')

copygraphics(fig, ContentType='vector')


% %% Single unit PETHs (pre and post first lick)
% IS = {isInPhase; isFirstQuarterPhase; isAntiPhase; isThirdQuarterPhase};
% ID = {idInPhase; idFirstQuarterPhase; idAntiPhase; idThirdQuarterPhase};
% ICOLOR = [1, 2, 4, 3];
% IEU = cell(1, 4);
% nEg = 15;
% for i = 1:4
%     isThisPhase = IS{i};
%     [SORTEDAMP{i}, I] = sort(amp(isThisPhase), 'descend');
%     IEU{i} = ID{i}(I);
% end
% 
% fig = figure(Units='inches', Position=[1 1 8 4]);
% tl = tiledlayout(fig, nEg, 2, TileSpacing='tight', Padding='compact');
% for iEg = 1:nEg
%     for i = [1, 3]
%         ax = nexttile(tl);
%         iEu = IEU{i}(iEg);
%         t = eta.correctLickBout.t;
%         t(t>0) = t(t>0) ./ (2*pi) / 8;
%         X = eta.correctLickBout.X(iEu, :);
%         X = smoothdata(X, 2, 'gaussian', 5);
%         plot(ax, t, X, LineWidth=1.5, Color=colors(ICOLOR(i), :))
%         xlim(ax, [-0.5, nBoutsDisp/8])
% %         ylim(ax, [0, 200])
%         text(ax, 0.05, 0.95, eu(iEu).getName(), Unit='normalized', HorizontalAlignment='left', VerticalAlignment='bottom', Interpreter='none', FontSize=p.fontSize)
%         xticks(ax, [-1:0.5:0, (1:nBoutsDisp)/8])
%         xticklabels(ax, [{'-1', '-0.5', '0'}, arrayfun(@(x) sprintf('%i\\pi', x), 2*(1:nBoutsDisp), UniformOutput=false)]);
%         ax.XGrid = 'on';
%         fontsize(ax, p.fontSize, 'points')
%     end
% end
% 
% 
% %% circlicks, bouts, trials analysis
% nCirclicks = cellfun(@length, trials);
% nBouts = cellfun(@(x) size(x, 1), bouts);
% nCorrectLickTrials = arrayfun(@(eu) nnz(eu.Trials.Lick.duration()>=4), eu(:));
% nCorrectPressTrials = arrayfun(@(eu) nnz(eu.Trials.Press.duration()>=4), eu(:));
% nCorrectTrials = nCorrectLickTrials + nCorrectPressTrials;
% tbl = table(nCirclicks, nBouts, nCorrectTrials, nCorrectLickTrials, nCorrectPressTrials);