eu = EphysUnit.load('C:\SERVER\Units\PressVsLick_ArtifactsRemoved_Full\FixedEventsAndTrials');
load('C:\SERVER\Units\meta_PressVsLick_ArtifactsRemoved_Full.mat');
c = cc;

%% Calculate peri-lick lick frequency histograms for any first lick in trial
[~, expEuIndices] = unique({eu.ExpName});
FsLick = 500;
lickHistEdges = 0.01:1/FsLick:2;
lickHistCenters = 0.5*(lickHistEdges(2:end) + lickHistEdges(1:end-1));

lickHistCounts = zeros(size(lickHistCenters));
lickHistNLicks = 0;
for iExp = 1:length(expEuIndices)
    iEu = expEuIndices(iExp);
    firstLickTimes = [eu(iEu).makeTrials('firstlick').Stop];
    allLickTimes = eu(iEu).EventTimes.Lick;
    for iLick = 1:length(firstLickTimes)
        edgesGlobal = firstLickTimes(iLick) + lickHistEdges;
        n = histcounts(allLickTimes, edgesGlobal);
        lickHistCounts = lickHistCounts + n;
        lickHistNLicks = lickHistNLicks + 1;
    end
end
lickHist.firstLick.count = lickHistCounts;
lickHist.firstLick.pdf = lickHistCounts ./ lickHistNLicks;
lickHist.firstLick.t = lickHistCenters;
lickHist.firstLick.edges = lickHistEdges;

clear expEuIndices FsLick lickHistEdges lickHistCenters lickHistNLicks lickHistCounts iExp iEu firstLickTimes allLickTimes iLick edgesGlobal n

% Calculate peri-correct-lick lick frequency histograms for osci cells
sel = c.isLick;
[~, expEuIndices] = unique({eu(sel).ExpName});
FsLick = 500;
lickHistEdges = 0.01:1/FsLick:2;
lickHistCenters = 0.5*(lickHistEdges(2:end) + lickHistEdges(1:end-1));

lickHistCounts = zeros(size(lickHistCenters));
lickHistNLicks = 0;
for iExp = 1:length(expEuIndices)
    iEu = expEuIndices(iExp);
    trials = eu(iEu).getTrials('lick');
    trials = trials(trials.duration() >= 4);
    firstLickTimes = [trials.Stop];
    allLickTimes = eu(iEu).EventTimes.Lick;
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

clear sel trials expEuIndices FsLick lickHistEdges lickHistCenters lickHistNLicks lickHistCounts iExp iEu firstLickTimes allLickTimes iLick edgesGlobal n
%%
eta.circLickNaive = eu.getETA('count', 'circlick_naive', window=[0, 2*pi], resolution=2*pi/30, normalize='none',  minInterval=0.05, maxInterval=0.20, artifacts=artifactParams);
eta.lickBoutNaive = eu.getETA('count', 'lickbout_naive', window=[0, 2*pi*4], resolution=2*pi/30, normalize='none',  minInterval=0.05, maxInterval=0.20, ...
    minBoutCycles=2, maxBoutCycles=4, artifacts=artifactParams);

eta.correctLickBout = eu.getETA('count', 'lick+lickbout', window=[-4, 0], resolution=[0.025, 2*pi/5], normalize='none', minTrialDuration=4, ...
    maxInterval=0.2, minInterval=0.05, minBoutCycles=2, maxBoutCycles=6);
eta.correctLickBoutNorm = eu.getETA('count', 'lick+lickbout', window=[-4, 0], resolution=[0.025, 2*pi/5], normalize=[-4, -2], minTrialDuration=4, ...
    maxInterval=0.2, minInterval=0.05, minBoutCycles=2, maxBoutCycles=6);

eta.correctPressBout = eu.getETA('count', 'press+lickbout', window=[-4, 0], resolution=[0.025, 2*pi/32, 2*pi/5], normalize='none', minTrialDuration=4, ...
    maxInterval=0.2, minInterval=0.05, minBoutCycles=2, maxBoutCycles=6);
eta.correctPressBoutNorm = eu.getETA('count', 'press+lickbout', window=[-4, 0], resolution=[0.025, 2*pi/32, 2*pi/5], normalize=[-4, -2], minTrialDuration=4, ...
    maxInterval=0.2, minInterval=0.05, minBoutCycles=2, maxBoutCycles=6);

%% Fig6

close all
p.fontSize = 9;
p.lineWidth = 1.5;
nBoutsDisp = 6;
W = [0, 0.5 + nBoutsDisp/8, 0.5 + 0.8+(nBoutsDisp-1)/8]*40;

clear layout l
layout.w = 7;
layout.h = 6.5;
layout.left.w = 3;
layout.right.w = 6;
layout.left.h = [4, 9, 4, 12];
layout.right.h = [3, 6];

fig = figure(Units='inches', Position=[1, 1, layout.w, layout.h], DefaultAxesFontSize=p.fontSize);
layout.tl = tiledlayout(fig, 1, layout.left.w + layout.right.w, TileSpacing='loose', Padding='loose');

layout.left.tl = tiledlayout(layout.tl, sum(layout.left.h), 1, TileSpacing='loose', Padding='loose');
l = layout.left.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [1, layout.left.w];

layout.right.tl = tiledlayout(layout.tl, sum(layout.right.h), 1, TileSpacing='compact', Padding='loose');
l = layout.right.tl; l.Layout.Tile = 1 + layout.left.w; l.Layout.TileSpan = [1, layout.right.w];

layout.left.top.tl = tiledlayout(layout.left.tl, sum(layout.left.h(1:3)), 1, TileSpacing='compact', Padding='compact');
l = layout.left.top.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [sum(layout.left.h(1:3)), 1];

layout.left.bottom.tl = tiledlayout(layout.left.tl, 4, 1, TileSpacing='compact', Padding='compact');
l = layout.left.bottom.tl; l.Layout.Tile = 1 + sum(layout.left.h(1:3)); l.Layout.TileSpan = [sum(layout.left.h(4)), 1];

layout.right.top.tl = tiledlayout(layout.right.tl, 1, 2, TileSpacing='compact', Padding='compact');
l = layout.right.top.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [layout.right.h(1), 1];

layout.right.bottom.tl = tiledlayout(layout.right.tl, 4, sum(W), TileSpacing='compact', Padding='compact');
l = layout.right.bottom.tl; l.Layout.Tile = 1 + sum(layout.right.h(1)); l.Layout.TileSpan = [layout.right.h(2), 1];


% 7a. Peri-lick lick prob histogram
% 7f. Correct lick bout lick histogram
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
hLetter.Position = [-0.35, ax.Position(4) + 0.1, 0];


% 7b ETA Heatmap osci lick
tlp = tiledlayout(layout.left.top.tl, 1, 6+1, TileSpacing='tight');
tlp.Layout.Tile = 1 + layout.left.h(1); tlp.Layout.TileSpan = [layout.left.h(2), 1];

tl = tiledlayout(tlp, 1, 1);
tl.Layout.Tile = 1; tl.Layout.TileSpan = [1, 6];
ax = nexttile(tl);

% eta.circLickNaive.Z = eta.circLickNaive.X.*exp(eta.circLickNaive.t*1i);
% eta.circLickNaive.Z(:, 1) = mean(eta.circLickNaive.X(:, [2, 30]), 2).*exp(eta.circLickNaive.t(1)*1i);
% meanZ = mean(eta.circLickNaive.Z, 2);

eta.lickBoutNaiveNorm = eta.lickBoutNaive;
eta.lickBoutNaiveNorm.X = normalize(eta.lickBoutNaive.X, 2, 'zscore', 'robust');

% Special plot: let's use pre-reach baseline for z-scoring osci-lick spike
% rates
% eta.lickBoutNaiveNorm.X(c.hasPress, :) = (eta.lickBoutNaive.X(c.hasPress, :) - vertcat(eta.press.stats(c.hasPress).mean)./0.1) ./ (vertcat(eta.press.stats(c.hasPress).sd)./0.1);
% eta.lickBoutNaiveNorm.X(c.hasPress, :) = (eta.lickBoutNaive.X(c.hasPress, :) - vertcat(eta.press.stats(c.hasPress).mean)./0.1) ./ (vertcat(eta.press.stats(c.hasPress).sd)./0.1);


maxBoutCycles = 4;
sel = c.isLick;%c.hasPress & c.hasLick & c.isLick;
phase = angle(circlick.Z(sel));
amp = abs(circlick.Z(sel));
phase(phase < 0) = phase(phase < 0) + 2*pi;
amp = amp(:);
phase = phase(:);
[sortedPhase, I] = sort(phase);
[~, ~] = EphysUnit.plotETA(ax, eta.lickBoutNaiveNorm, sel, order=I, ...
    clim=[-5, 5], xlim=[0, 2*pi*maxBoutCycles], hidecolorbar=false);
applyCustomColormap(ax, [-5, 5], hlim=[0.375, 0, 0, -0.375], llim=[0.125, 0.5, 0.5, 0.25], hpwr=.5, lpwr=1, h0=0.33);
xticks(ax, (0:2:8).*pi);
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
% sel = find(c.isLick & c.hasPress & c.hasLick);
sel = find(c.isLick);
rectifiedPhase = phase;
rectifiedPhase(phase < 0) = rectifiedPhase(phase < 0) + 2*pi;
ampThreshold = 0;
isInPhase = phase < 0.25*pi | phase >= 1.75*pi;
isAntiPhase = phase >= 0.75*pi & phase < 1.25*pi;
isFirstQuarterPhase = phase >= 0.25*pi & phase < 0.75*pi;
isThirdQuarterPhase = phase >= 1.25*pi & phase < 1.75*pi;
isHighAmp = amp >= quantile(amp, ampThreshold);
idInPhase = sel(isInPhase & isHighAmp);
idInPhaseLickUp = sel(isInPhase & isHighAmp & c.isLickUp(sel));
idInPhaseLickFlat = sel(isInPhase & isHighAmp & ~c.isLickResponsive(sel));
idInPhaseLickDown = sel(isInPhase & isHighAmp & c.isLickDown(sel));
idInPhasePressUp = sel(isInPhase & isHighAmp & c.isPressUp(sel));
idInPhasePressFlat = sel(isInPhase & isHighAmp & ~c.isPressResponsive(sel));
idInPhasePressDown = sel(isInPhase & isHighAmp & c.isPressDown(sel));
idAntiPhase = sel(isAntiPhase & isHighAmp);
idAntiPhaseLickUp = sel(isAntiPhase & isHighAmp & c.isLickUp(sel));
idAntiPhaseLickFlat = sel(isAntiPhase & isHighAmp & ~c.isLickResponsive(sel));
idAntiPhaseLickDown = sel(isAntiPhase & isHighAmp & c.isLickDown(sel));
idAntiPhasePressUp = sel(isAntiPhase & isHighAmp & c.isPressUp(sel));
idAntiPhasePressFlat = sel(isAntiPhase & isHighAmp & ~c.isPressResponsive(sel));
idAntiPhasePressDown = sel(isAntiPhase & isHighAmp & c.isPressDown(sel));
idFirstQuarterPhase = sel(isFirstQuarterPhase & isHighAmp);
idFirstQuarterPhaseLickUp = sel(isFirstQuarterPhase & isHighAmp & c.isLickUp(sel));
idFirstQuarterPhaseLickFlat = sel(isFirstQuarterPhase & isHighAmp & ~c.isLickResponsive(sel));
idFirstQuarterPhaseLickDown = sel(isFirstQuarterPhase & isHighAmp & c.isLickDown(sel));
idFirstQuarterPhasePressUp = sel(isFirstQuarterPhase & isHighAmp & c.isPressUp(sel));
idFirstQuarterPhasePressFlat = sel(isFirstQuarterPhase & isHighAmp & ~c.isPressResponsive(sel));
idFirstQuarterPhasePressDown = sel(isFirstQuarterPhase & isHighAmp & c.isPressDown(sel));
idThirdQuarterPhase = sel(isThirdQuarterPhase & isHighAmp);
idThirdQuarterPhaseLickUp = sel(isThirdQuarterPhase & isHighAmp & c.isLickUp(sel));
idThirdQuarterPhaseLickFlat = sel(isThirdQuarterPhase & isHighAmp & ~c.isLickResponsive(sel));
idThirdQuarterPhaseLickDown = sel(isThirdQuarterPhase & isHighAmp & c.isLickDown(sel));
idThirdQuarterPhasePressUp = sel(isThirdQuarterPhase & isHighAmp & c.isPressUp(sel));
idThirdQuarterPhasePressFlat = sel(isThirdQuarterPhase & isHighAmp & ~c.isPressResponsive(sel));
idThirdQuarterPhasePressDown = sel(isThirdQuarterPhase & isHighAmp & c.isPressDown(sel));
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
edges = -0:2*pi/32:2*pi;
histogram(ax, rectifiedPhase(isInPhase), edges, FaceColor=colors(1, :), FaceAlpha=1, EdgeAlpha=0.5);
histogram(ax, rectifiedPhase(isFirstQuarterPhase), edges, FaceColor=colors(2, :), FaceAlpha=1, EdgeAlpha=0.5);
histogram(ax, rectifiedPhase(isAntiPhase), edges, FaceColor=colors(4, :), FaceAlpha=1, EdgeAlpha=0.5);
histogram(ax, rectifiedPhase(isThirdQuarterPhase), edges, FaceColor=colors(3, :), FaceAlpha=1, EdgeAlpha=0.5);
xticks(ax, 0:pi:2*pi)
yticks(ax, [0, 20])
xlim(ax, pi*[0, 2])
xticklabels(ax, {'0', '\pi', '2\pi'});
xlabel('Lick phase')
ylabel('# units')
fontsize(ax, p.fontSize, 'points')

hLetter = text(ax, 0, 0, 'c', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.4, ax.Position(4) + 0.25, 0];

ID = {idInPhase; idFirstQuarterPhase; idAntiPhase; idThirdQuarterPhase};
IDSplit = { ...
    {idInPhaseLickUp, idInPhaseLickFlat, idInPhaseLickDown}, {idInPhasePressUp, idInPhasePressFlat, idInPhasePressDown}; ...
    {idFirstQuarterPhaseLickUp, idFirstQuarterPhaseLickFlat, idFirstQuarterPhaseLickDown}, {idFirstQuarterPhasePressUp, idFirstQuarterPhasePressFlat, idFirstQuarterPhasePressDown}; ...
    {idAntiPhaseLickUp, idAntiPhaseLickFlat, idAntiPhaseLickDown}, {idAntiPhasePressUp, idAntiPhasePressFlat, idAntiPhasePressDown}; ...
    {idThirdQuarterPhaseLickUp, idThirdQuarterPhaseLickFlat, idThirdQuarterPhaseLickDown}, {idThirdQuarterPhasePressUp, idThirdQuarterPhasePressFlat, idThirdQuarterPhasePressDown}; ...
    };
ETAMOVEBOUT = {eta.correctLickBoutNorm, eta.correctPressBoutNorm};
TASKS = ["lick", "press"];
TASKTITLE = ["Self-timed lick", "Self-timed reach"];
PHASENAME = ["2\pi", "1/2\pi", "\pi", "3/2\pi"];
ICOLOR = [1, 2, 4, 3];

% 7d (left) and 7g (right)
AX = gobjects(4, 3);
for iAx = 1:4
    iEu = ID{iAx};
    for iTask = 1:2
        % First lick
        ax = nexttile(layout.right.bottom.tl, (iAx-1)*sum(W) + 1 + sum(W(1:iTask)), [1, W(iTask + 1)]); 
        AX(iAx, iTask + 1) = ax;
        hold(ax, 'on')
        t = ETAMOVEBOUT{iTask}.t;
        switch TASKS(iTask)
            case "press"
                t(t>0 & t<=2*pi) = t(t>0 & t<=2*pi) / (2*pi) * 0.8;
                t(t>2*pi) = (t(t>2*pi) - 2*pi) ./ (2*pi) / 8 + 0.8;
            case "lick"
                t(t>0) = t(t>0) ./ (2*pi) / 8;
        end

        for iDir = [1, 3]
            iEuDir = IDSplit{iAx, iTask}{iDir};
            X = ETAMOVEBOUT{iTask}.X(iEuDir, :);
            if isempty(X)
                continue
            end
            X = smoothdata(X, 2, 'gaussian', 5);
            plot(ax, t, mean(X, 1, 'omitnan'), Color=colors(ICOLOR(iAx), :), LineWidth=1.5)
            plot(ax, t, X, Color=[0.15, 0.15, 0.15, 0.05], LineWidth=0.5)
        end

        ylim(ax, [-1.25, 1.75])
        yticks(ax, [-1, 0, 1])
        
        switch TASKS(iTask)
            case "press"
                xline(ax, [0, 0.8,  0.8+(1:(nBoutsDisp-1))/8], LineStyle=':')
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

        if iAx == 1
            title(ax, TASKTITLE(iTask))
        end

        ax.XGrid = 'on';
        hold(ax, 'off')
        fontsize(ax, p.fontSize, 'points')
        text(ax, 0.025, 1, sprintf('inc(n=%i)', nnz(IDSplit{iAx, iTask}{1})), Unit='normalized', HorizontalAlignment='left', VerticalAlignment='top', Interpreter='none', FontSize=p.fontSize-1)
        text(ax, 0.025, 0.025, sprintf('dec(n=%i)', nnz(IDSplit{iAx, iTask}{3})), Unit='normalized', HorizontalAlignment='left', VerticalAlignment='bottom', Interpreter='none', FontSize=p.fontSize-1)
        yline(ax, 0, '--')
    end

    % 6d. Any bout
    ax = nexttile(layout.left.bottom.tl); AX(iAx, 1) = ax;
    X = eta.lickBoutNaiveNorm.X(iEu, :);
    X = smoothdata(X, 2, 'gaussian', 5);
    plot(ax, eta.lickBoutNaiveNorm.t, mean(X, 1, 'omitnan'), Color=colors(ICOLOR(iAx), :), LineWidth=1.5)
    hold(ax, 'on')
    plot(ax, eta.lickBoutNaiveNorm.t, X, Color=[0.15, 0.15, 0.15, 1./nnz(iEu)], LineWidth=0.5)
    xticks(ax, 0:2*pi:8*pi)
    xticklabels(ax, [{'0'}, arrayfun(@(x) sprintf('%i\\pi', x), 2:2:8, UniformOutput=false)]);
    xlim(ax, [0, 8*pi])
    ylim(ax, [-4, 4])
    yticks(ax, [-3, 3])
    ax.XGrid = 'on';
    fontsize(ax, p.fontSize, 'points')
    text(ax, 0.05, -0.025, sprintf('n=%i', nnz(iEu)), Unit='normalized', HorizontalAlignment='left', VerticalAlignment='bottom', Interpreter='none', FontSize=p.fontSize-1.5)
    yline(ax, 0, '--')
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
hLetter.Position = [-0.4, ax.Position(4) + 0.25, 0];

ax = AX(1, 2);
hLetter = text(ax, 0, 0, 'g', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.1, ax.Position(4) + 0.8, 0];

% 7f. Reach vs. Lick (scatter)
AX = gobjects(2, 1);
sz = 7;
ax = nexttile(layout.right.top.tl); AX(1) = ax;
hold(ax, 'on')
% sel = true(size(eu));
x = meta.lickNorm;
y = meta.pressNorm;
subselResp = c.isLick;
subselNone = ~c.isLick;
h = gobjects(2, 1);
h(2) = scatter(ax, x(subselNone), y(subselNone), sz, 'black', 'filled', Marker='o', MarkerFaceAlpha=0.5, MarkerEdgeAlpha=0.5, DisplayName=sprintf('others (%i)', nnz(subselNone)));   
h(1) = scatter(ax, x(subselResp), y(subselResp), sz, [0, 0.5, 0.5], 'filled', Marker='o', MarkerFaceAlpha=0.5, MarkerEdgeAlpha=0.5, DisplayName=sprintf('lick-entrained (%i)', nnz(subselResp)));

mdl = fitlm(meta.lickNorm(subselResp), meta.pressNorm(subselResp));
fprintf('press vs. lick (osci): LM slope p<%g.\n', mdl.Coefficients.pValue(2))

plot(ax, [-10, 10], [0, 0], 'k:');
plot(ax, [0, 0], [-10, 10], 'k:');
plot(ax, [-10, 10], [-10, 10], 'k:')


% 7e right, Additional plot, scatter press vs lick META, color by lick entrainment phase: 
sz = 7;
ax = nexttile(layout.right.top.tl); AX(2) = ax;
hold(ax, 'on')
x = meta.lickNorm;
y = meta.pressNorm;
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

ylabel(AX, 'Peri-reach activity (a.u.)', FontSize=p.fontSize)
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

hCb = colorbar(axc);
hCb.Label.String = 'Norm spike rate (a.u.)';
hCb.Label.Position(1) = 0;
hCb.Label.VerticalAlignment = 'bottom';
% hCb.Layout.Tile = 'east';

copygraphics(fig, ContentType='vector', BackgroundColor='none')
