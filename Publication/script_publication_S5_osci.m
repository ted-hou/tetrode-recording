p.fontSize = 9;
p.lineWidth = 1.5;

clear layout
layout.w = 7;
layout.h = 8;


fig = figure(Units='inches', Position=[1, 1, layout.w, layout.h], DefaultAxesFontSize=p.fontSize);
layout.tl = tiledlayout(fig, );

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


% 5f

eta.circLick.Z = eta.circLick.X.*exp(eta.circLick.t*1i);
eta.circLick.Z(:, 1) = mean(eta.circLick.X(:, [2, 30]), 2).*exp(eta.circLick.t(1)*1i);
meanZ = mean(eta.circLick.Z, 2);
eta.lickBoutNorm = eta.lickBout;
eta.lickBoutNorm.X = normalize(eta.lickBout.X, 2, 'zscore', 'robust');

maxBoutCycles = 4;
sel = c.hasPress & c.hasLick & c.isLick; 
phase = angle(meanZ(sel));
phase(phase < 0) = phase(phase < 0) + 2*pi;
[sortedPhase, I] = sort(phase);
[~, ~] = EphysUnit.plotETA(ax(4), eta.lickBoutNorm, sel, order=I, ...
    clim=[-2, 2], xlim=[0, 2*pi*maxBoutCycles], hidecolorbar=true);
xticks (ax(4), (0:2:8).*pi);
xticklabels(ax(4), [{'0'}, arrayfun(@(x) sprintf('%i\\pi', x), 2:2:8, UniformOutput=false)]);
title(ax(4), 'Lick-entrained')
ylabel(ax(4), 'Unit')
xlabel(ax(4), 'Lick phase')
xlim(ax(4), [0, 2*pi*maxBoutCycles])
ylim(ax(4), [0, nnz(sel)+1])
yt = 0:100:nnz(sel);
yt(1) = 1;
if round(yt(end)./100) == round(nnz(sel)./100)
    yt(end) = nnz(sel);
else
    yt(end + 1) = nnz(sel);
end
yt = unique(yt);
yticks(ax(4), yt)
fontsize(ax, p.fontSize, 'points')
fontname(ax, 'Arial')
axc = ax(4);


%%



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

ttll = tiledlayout(figure, 4, 2);
for iAx = 1:4
    iEu = ID{iAx};
    % First lick
    % ax = nexttile(layout.bottom.left.tl);
    ax = nexttile(ttll);
    t = eta.correctLickBout.t;
    t(t>0) = t(t>0) ./ (2*pi) / 8;
    X = eta.correctLickBoutNorm.X(iEu, :);
    X = smoothdata(X, 2, 'gaussian', 5);
    plot(ax, t, mean(X, 1, 'omitnan'), Color=colors(ICOLOR(iAx), :), LineWidth=1.5)
    hold(ax, 'on')
    plot(ax, t, X, Color=[0.15, 0.15, 0.15, 0.025], LineWidth=0.5)
    xlim(ax, [-1, nBoutsDisp/8])
    ylim(ax, [-0.5, 1.5])
    % yticks(ax, [40, 100])
    text(ax, 0.05, 0.95, sprintf('n=%i', nnz(iEu)), Unit='normalized', HorizontalAlignment='left', VerticalAlignment='top', Interpreter='none', FontSize=p.fontSize)
    xticks(ax, [-1:0.5:0, (1:nBoutsDisp)/8])
    xticklabels(ax, [{'-1', '-0.5', '0'}, arrayfun(@(x) sprintf('%i\\pi', x), 2*(1:nBoutsDisp), UniformOutput=false)]);
    
    ax.XGrid = 'on';
    hold(ax, 'off')
    fontsize(ax, p.fontSize, 'points')
    if iAx == 4
        text(ax, 0, -0.5, 'Time before first lick (s)', Units='normalized', VerticalAlignment='top', FontSize=p.fontSize)
        text(ax, 0.55, -0.55, 'Lick phase', Units='normalized', VerticalAlignment='top', FontSize=p.fontSize)
    end

    % Any bout
    ax = nexttile(ttll);
    % ax = nexttile(layout.bottom.right.tl);
    X = eta.lickBoutNorm.X(iEu, :);
    X = smoothdata(X, 2, 'gaussian', 5);
    plot(ax, eta.lickBoutNorm.t, mean(X, 1, 'omitnan'), Color=colors(ICOLOR(iAx), :), LineWidth=1.5)
    hold(ax, 'on')
    plot(ax, eta.lickBoutNorm.t, X, Color=[0.15, 0.15, 0.15, 0.025], LineWidth=0.5)
    xticks(ax, 0:2*pi:8*pi)
    xticklabels(ax, [{'0'}, arrayfun(@(x) sprintf('%i\\pi', x), 2:2:8, UniformOutput=false)]);
    xlim(ax, [0, 8*pi])
    ylim(ax, [-4, 4])
    yticks(ax, [-3, 0, 3])
    ax.XGrid = 'on';
end
% h = xlabel(layout.bottom.left.tl, 'Time from first lick (s)       Lick Phase', FontSize=p.fontSize, HorizontalAlignment='left');
ylabel(layout.bottom.left.tl, 'Spike rate (sp/s)', FontSize=p.fontSize)
xlabel(layout.bottom.right.tl, 'Lick phase', FontSize=p.fontSize)
ylabel(layout.bottom.right.tl, 'Normalized spike rate (a.u.)', FontSize=p.fontSize)

% Additional plot, scatter press vs lick META, color by lick entrainment phase: 
sz = 7;
% 1. Reach vs. Lick scatter
ax = nexttile(layout.psbottom.tl, 1 + layout.left.w, [1, layout.right.w]);
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

xlabel(ax, 'peri-lick response (a.u.)', HorizontalAlignment='center')
ylabel(ax, 'peri-reach response (a.u.)', HorizontalAlignment='center')
axis(ax, 'equal')
xlim(ax, [-2, 5])
ylim(ax, [-2, 5])
% lgd = legend(ax, h);
% lgd.Layout.Tile = 'north';

fontsize(ax, p.fontSize, 'points')



