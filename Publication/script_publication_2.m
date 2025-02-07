% Fig 3. GLM, neural modulation onset vs. movement onset
%% 
load_ephysunits
% load('C:\SERVER\bootMoveResponse_20240830.mat')
% boot_baselineVsUpDown
read_DLC_data

% fit_GLM
load('C:\SERVER\acute_glm_20241218.mat') % Load previously saved fit_GLM output


% boot_bta
load('C:\SERVER\boot_bta_20250121.mat') % load boot_bta results


%% TST calculated from contra paw spd (2 samples conseq above 0.25 zscore)
% tst = [-0.233333333333334	-0.266666666666680	-0.0999999999999943	-0.200000000000045	-0.533333333333303	-0.166666666666629	-0.166666666666629	-0.100000000000023	-0.133333333333326	-0.100000000000023	-0.133333333333326	-0.133333333333326	-0.133333333333326	-0.100000000000023	-0.166666666666629	-0.566666666666720	-0.133333333333326	-0.133333333333326	-0.200000000000045	-0.200000000000045	-0.166666666666629	-0.133333333333439	-0.133333333333439	-0.433333333333394	-0.166666666666742	-0.133333333333439	-0.133333333333212	-0.166666666666515	-0.0999999999999091	-0.166666666666515	-0.133333333333212	-0.0999999999999091	-0.133333333333212	-0.133333333333212	-0.133333333333212	-0.133333333333212	-0.199999999999818	-0.133333333333212	-0.133333333333212	-0.0999999999999091	-0.0999999999999091	-0.166666666666515	-0.133333333333212	-0.333333333333485	-0.166666666666515	-0.199999999999818	0	-0.0999999999999091	0	-0.166666666666970	-0.133333333333212	-0.166666666666970	-0.166666666666970	-0.199999999999818	-0.233333333333576	-0.166666666666970	-0.199999999999818	-0.166666666666970	-0.133333333333212	0	-0.133333333333212	-0.166666666666657	-0.0999999999999943	0	-0.300000000000011	-0.199999999999989	-0.266666666666652	-0.166666666666686	-0.133333333333326	0	-0.100000000000023	0	0	-0.100000000000023	-0.166666666666629	-0.333333333333371	-0.166666666666742	-0.166666666666515	-0.166666666666515	-0.133333333333212	-0.266666666666879	-0.233333333333121	-0.400000000000091	0	0	-0.233333333333121	-0.133333333333212	-0.233333333333121	-0.300000000000182	-0.133333333333212	-0.0999999999999091	0	-0.300000000000182	-0.400000000000091	-0.166666666666970	-0.100000000000364	-0.233333333333576	-0.199999999999818	-0.133333333333212	-0.133333333333212	-0.166666666666970	-0.199999999999989	-0.0999999999999943	-0.133333333333326	-0.100000000000023	-0.133333333333326	-0.133333333333326	-0.300000000000011	-0.166666666666686	-0.0666666666667197	-0.133333333333326	-0.100000000000023	-0.100000000000023	-0.133333333333326	-0.0666666666667197	-0.133333333333326	0	-0.0666666666666060	-0.0666666666666060	-0.133333333333212	-0.166666666666515	-0.133333333333212	-0.0999999999999091	-0.266666666666879	-0.0666666666666060	-0.0999999999999091	-0.166666666666515	-0.133333333333212	-0.233333333333121	-0.133333333333212	-0.0999999999999091	-0.0666666666666060	-0.0666666666666060	-0.166666666666515	-0.0999999999999091	-0.0999999999999091	-0.199999999999818	-0.133333333333212	-0.133333333333212	-0.133333333333212	0	0	-0.166666666666686	0	-0.166666666666686	-0.100000000000023	-0.133333333333326	0	-0.100000000000023	-0.133333333333326	0	-0.0666666666667197	-0.233333333333349	-0.133333333333326	-0.166666666666742	-0.133333333333439	-0.200000000000045	-0.233333333333349	-0.166666666666515	-0.266666666666879	-0.133333333333212	-0.0999999999999091	-0.0333333333333030	-0.0333333333333030	-0.0666666666666060	-0.233333333333121	-0.133333333333212	-0.233333333333576	-0.0333333333337578	0	0	0	-0.433333333333337	-0.199999999999989	-0.333333333333314	0	-0.133333333333326	-0.600000000000023	0	0	-0.433333333333280	-0.433333333333280	-0.0666666666667197	-0.100000000000023	0	-0.166666666666515	-0.233333333333121	-0.166666666666515	-0.0666666666666060	-0.199999999999818	0	0	-0.233333333333576	-0.133333333333212	-0.199999999999818	-0.100000000000364	-0.266666666666680	-0.0999999999999943	-0.166666666666657	0	-0.166666666666629	-0.199999999999818	-0.0999999999999091	-0.233333333333121	-0.266666666666879	-0.166666666666515	0	-0.199999999999818	-0.199999999999818	0	-0.166666666666515	-0.166666666666515	-0.199999999999818	-0.233333333333121	-0.199999999999818	-0.199999999999818	0	-0.133333333333212	-0.300000000000182	0	0	-0.333333333333030	-0.266666666666424	-0.300000000000182	0	0	-0.266666666666424	-0.133333333333212	-0.233333333333576	-0.199999999999818	-0.0666666666666060	-0.600000000000001	-0.133333333333326	-0.166666666666629	0	0	0	-0.633333333333212	-0.233333333333121	-0.533333333333303	0	0	0	0	-0.199999999999818	-0.133333333333212	-0.300000000000182	-0.333333333333030	-0.666666666666970	0	0	-0.766666666666424	0	-0.199999999999818	-0.500000000000000	-0.233333333333576];

% p.etaSortWindow = [-3, 0];
% p.etaSignWindow = [-0.3, 0];
% p.etaLatencyThresholdPos = 0.25;
% p.etaLatencyThresholdNeg = 0.25;

% close all
figTemp = figure(Units='normalized', Position=[0.1, 0.1, 0.5, 0.5]);
axTemp = subplot(1, 2, 1);
EphysUnit.plotETA(axTemp, eta.press, c.hasPress, xlim=[-4,0.5], clim=[-1.5, 1.5], ...
    order=onset.pressOrder(c.hasPress)); 
hold(axTemp, 'on')
pressLatencySorted = onset.press(c.hasPress); pressLatencySorted = pressLatencySorted(onset.pressOrder(c.hasPress));
plot(axTemp, pressLatencySorted, 1:nnz(c.hasPress))
xline(axTemp, 0, 'k--')

axTemp = subplot(1, 2, 2);
EphysUnit.plotETA(axTemp, eta.lick, c.hasLick, xlim=[-4,0.5], clim=[-1.5, 1.5], ...
    order=onset.lickOrder(c.hasLick)); 
hold(axTemp, 'on')
lickLatencySorted = onset.lick(c.hasLick); lickLatencySorted = lickLatencySorted(onset.lickOrder(c.hasLick));
plot(axTemp, lickLatencySorted, 1:nnz(c.hasLick))
xline(axTemp, 0, 'k--')

latency.press = onset.press;
latency.lick = onset.lick;
contraOnset4tgt = {trajCombined4tgt.onset.contra};

latency.contraPaw = [fAll.press.onset; cat(1, trajCombined2tgt.onset{:}); cat(2, contraOnset4tgt{:})'];
latency.contraPawN = struct( ...
    trials=nnz(~isnan(latency.contraPaw)), ...
    sessions=length(fCorrect) + length(traj2tgt) + length(traj4tgt), ...
    animals=length(unique(euAcute.getAnimalName())) + length(unique(euReachDir2tgt.getAnimalName())) + length(unique(euReachDir4Tgt.getAnimalName())) ...
    );
latency.pRankSum.pressVsContraPaw = ranksum(latency.contraPaw , latency.press(c.isPressResponsive), tail='right');
latency.pRankSum.pressUpVsContraPaw = ranksum(latency.contraPaw , latency.press(c.isPressUp), tail='right');
latency.pRankSum.pressDownVsContraPaw = ranksum(latency.contraPaw , latency.press(c.isPressDown), tail='right');
latency.pRankSum.pressVsLick = ranksum(latency.lick, latency.press, tail='right');
latency.pRankSum.pressVsContraPawPlus200 = ranksum(latency.contraPaw -0.2, latency.press(c.isPressResponsive), tail='right');
latency.pRankSum.pressUpVsContraPawPlus200 = ranksum(latency.contraPaw -0.2, latency.press(c.isPressUp), tail='right');
latency.pRankSum.pressDownVsContraPawPlus200 = ranksum(latency.contraPaw -0.2, latency.press(c.isPressDown), tail='right');
latency.pRankSum.pressUpVsPressDown = ranksum(latency.press(c.isPressUp), latency.press(c.isPressDown), tail='both');

fprintf(1, 'Median pre-press spiking onset latency = %.1f ms \n', median(latency.press(c.isPressResponsive)*1000, 'omitnan'))
fprintf(1, 'Median pre-press spiking onset latency (excited) = %.1f ms \n', median(latency.press(c.isPressUp)*1000, 'omitnan'))
fprintf(1, 'Median pre-press spiking onset latency (suppressed) = %.1f ms \n', median(latency.press(c.isPressDown)*1000, 'omitnan'))
fprintf(1, 'Median pre-lick spiking onset latency = %.1f ms \n', median(latency.lick(c.isLickResponsive)*1000, 'omitnan'))
fprintf(1, 'Median contralateral paw movement onset latency = %.3f ms \n', median(latency.contraPaw*1000, 'omitnan'))
fprintf(1, '90%% of forearm movements were initiated within %.3f ms prior to bar contact \n', quantile(latency.contraPaw*1000, 0.1))
fprintf(1, 'Press spiking precedes paw: One-tailed ranksum test p = %g\n', latency.pRankSum.pressVsContraPaw)
fprintf(1, 'Press spiking precedes paw (-200ms): One-tailed ranksum test p = %g\n', latency.pRankSum.pressVsContraPawPlus200)
fprintf(1, 'Press up spiking precedes paw (-200ms): One-tailed ranksum test p = %g\n', latency.pRankSum.pressUpVsContraPawPlus200)
fprintf(1, 'Press down spiking precedes paw (-200ms): One-tailed ranksum test p = %g\n', latency.pRankSum.pressDownVsContraPawPlus200)
fprintf(1, 'Excited press spiking precedes paw: One-tailed ranksum test p = %g\n', latency.pRankSum.pressUpVsContraPaw)
fprintf(1, 'Inhibited press spiking precedes paw: One-tailed ranksum test p = %g\n', latency.pRankSum.pressDownVsContraPaw)
fprintf(1, 'Press spiking precedes lick spiking: One-tailed ranksum test p = %g\n', latency.pRankSum.pressVsLick)
fprintf(1, 'Press excided vs. press supressed: two-tailed ranksum test p = %g\n', latency.pRankSum.pressUpVsPressDown)


%% Figure 3
p.fontSize = 9;

variantColors = hsl2rgb([linspace(0.3, 0.9, 5)', linspace(0.75, 0.25, 5)', 0.5*ones(5, 1)]);

nBTABins = length(p.binnedTrialEdges) - 1;
btaColors = hsl2rgb([linspace(0.7, 0, nBTABins)'.^1.3, linspace(0.8, 0.6, nBTABins)', 0.5*ones(nBTABins, 1)]);

clear layout
layout.w = 7;
layout.h = 8;
layout.top.h = 11;
layout.top.left.w = 3;
layout.top.left.top.h = 3;
layout.top.left.middle.h = 2;
layout.top.left.bottom.h = 2;
layout.top.right.w = 5;
layout.top.right.top.h = 3;
layout.top.right.middle.h = 3;
layout.top.right.bottom.h = 4;
layout.bottom.h = 5;

p.lineWidth = 1.5;
close all

% Create figure and layout panels
fig = figure(Units='inches', Position=[0, 0, layout.w, layout.h]);
layout.tl = tiledlayout(fig, layout.top.h + layout.bottom.h, 1, TileSpacing='compact');

layout.top.tl = tiledlayout(layout.tl, 1, layout.top.left.w + layout.top.right.w, TileSpacing='loose');
l = layout.top.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [layout.top.h, 1];

layout.top.left.tl = tiledlayout(layout.top.tl, layout.top.left.top.h + layout.top.left.middle.h + layout.top.left.bottom.h, 1, TileSpacing='tight');
l = layout.top.left.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [1, layout.top.left.w]; 

layout.top.right.tl = tiledlayout(layout.top.tl, layout.top.right.top.h + layout.top.right.middle.h + layout.top.right.bottom.h, 1, TileSpacing='tight', Padding='tight');
l = layout.top.right.tl; l.Layout.Tile = layout.top.left.w + 1; l.Layout.TileSpan = [1, layout.top.right.w];

layout.top.right.top.tl = tiledlayout(layout.top.right.tl, 1, 3, TileSpacing='tight', Padding='tight');
l = layout.top.right.top.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [layout.top.right.top.h, 1];

layout.top.right.middle.tl = tiledlayout(layout.top.right.tl, 1, 2, TileSpacing='tight', Padding='tight');
l = layout.top.right.middle.tl; l.Layout.Tile = 1 + layout.top.right.top.h; l.Layout.TileSpan = [layout.top.right.middle.h, 1];

layout.top.right.bottom.tl = tiledlayout(layout.top.right.tl, 1, 2, TileSpacing='tight', Padding='tight');
l = layout.top.right.bottom.tl; l.Layout.Tile = 1 + layout.top.right.top.h + layout.top.right.middle.h; l.Layout.TileSpan = [layout.top.right.bottom.h, 1];

layout.bottom.tl = tiledlayout(layout.tl, 2, 2, TileSpacing='tight');
l = layout.bottom.tl; l.Layout.Tile = 1 + layout.top.h; l.Layout.TileSpan = [layout.bottom.h, 1];

% 3a PETH of increase vs. decrease populations with arrows pointing to mean
% onset time
ax = gobjects(3, 1);
ax(1) = nexttile(layout.top.left.tl, 1, [layout.top.left.top.h, 1]);
hold(ax(1), 'on')
SEL = {c.isPressUp & c.hasPress, c.isPressDown & c.hasPress};
LABEL = [sprintf("inc (n=%i)", nnz(c.isPressUp)), sprintf("dec (n=%i)", nnz(c.isPressDown))];
LABELSHORT = ["inc", "dec"];
% ERR = {boot.pressBaseline.muErrUp.ci, boot.pressBaseline.muErrDown.ci};
COLOR = 'rb';
minY = NaN;
maxY = NaN;
h = gobjects(2, 1);
% Plot mean spike rate and SEM (across trials then across units)
for i = 1:2
    selUnits = SEL{i};
    sr = mean(eta.pressRaw.X(selUnits, :), 1)./0.1;
    err = std(eta.pressRaw.X(selUnits, :)./0.1, 0, 1)./sqrt(nnz(selUnits)); % SEM
%     err = ERR{i}; % 95CI from permtest
    minY = min(minY, min(sr));
    maxY = max(maxY, max(sr));
    h(i) = plot(ax(1), eta.pressRaw.t, sr, LineWidth=1, Color=COLOR(i), DisplayName=LABEL(i));
    patch(ax(1), [eta.pressRaw.t, flip(eta.pressRaw.t)], [sr - err, flip(sr + err)], COLOR(i), FaceAlpha=0.1, EdgeColor=COLOR(i), EdgeAlpha=0.2);
end

% Plot arrows pointing to neural onset
for i = 1:2
    selUnits = SEL{i};
    sr = mean(eta.pressRaw.X(selUnits, :), 1)./0.1;
    t0 = median(latency.press(selUnits));
    sr0 = interp1(eta.pressRaw.t, sr, t0, 'linear');
    plot(ax(1), [t0, t0], [sr0, sr0 + (maxY - minY)*0.15], LineStyle='-', Color=COLOR(i))
    scatter(ax(1), t0, sr0 + (maxY - minY)*0.04, 25, COLOR(i), 'v', 'filled')
end
ylim(ax(1), [10, 80])
hold(ax(1), 'off')
xlabel(ax(1), 'Time to bar-contact (s)')
ylabel(ax(1), 'Spike rate (sp/s)')
hLgd = legend(ax(1), h, Location='northwest', AutoUpdate=false);
% hLgd.Position(1) = 0.24;
% hLgd.Position(2) = 0.85;

% 4b. Histogram of neural onset times
ax(2) = nexttile(layout.top.left.tl, [layout.top.left.middle.h, 1]);
edges = -4:0.1:0;
hold(ax(2), 'on')
hHist = gobjects(2, 1);
hHist(1) = histogram(ax(2), latency.press(c.isPressUp), edges, Normalization='count', DisplayName=LABELSHORT(1), FaceColor=COLOR(1), FaceAlpha=0.6);
hHist(2) = histogram(ax(2), latency.press(c.isPressDown), edges, Normalization='count', DisplayName=LABELSHORT(2), FaceColor=COLOR(2), FaceAlpha=0.6);
xlabel(ax(2), 'Time to bar-contact (s)')
ylabel(ax(2), 'No. units')
title(ax(2), 'SNr response onset')
hLgd = legend(ax(2), hHist, Location='northwest', AutoUpdate=false);
hLgd.Position(1) = 0.22;
hLgd.Position(2) = 0.46;

% 4c. Histogram of movement onset times
ax(3) = nexttile(layout.top.left.tl, [layout.top.left.bottom.h, 1]);
histogram(latency.contraPaw, edges, Normalization='count', FaceColor='black');
hLgd = legend(sprintf('%g trials\n%g sessions\n%g animals', latency.contraPawN.trials, latency.contraPawN.sessions, latency.contraPawN.animals), Location='northwest', AutoUpdate=false);
hLgd.Position(1) = 0.22;
hLgd.Position(2) = 0.17;
title('Forepaw movement onset')
xlabel('Time to bar-contact (s)')
ylabel('No. trials')

fontsize(ax, p.fontSize, 'points')
fontname(ax, 'Arial')
xlim(ax, [-4, 0.5])
for iAx = 1:3
    xline(ax(iAx), 0, 'k--')
end
clear iAx

hLetters = gobjects(1, 7);
hLetters(1) = text(ax(1), 0, 0, 'a', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax(1).Units = 'inches';
hLetters(1).HorizontalAlignment = 'right';
hLetters(1).VerticalAlignment = 'top';
hLetters(1).Position = [-0.5, ax(1).Position(4)+0.3, 0];

hLetters(2) = text(ax(2), 0, 0, 'b', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax(2).Units = 'inches';
hLetters(2).HorizontalAlignment = 'right';
hLetters(2).VerticalAlignment = 'top';
hLetters(2).Position = [-0.5, ax(2).Position(4)+0.2, 0];

hLetters(3) = text(ax(3), 0, 0, 'c', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax(3).Units = 'inches';
hLetters(3).HorizontalAlignment = 'right';
hLetters(3).VerticalAlignment = 'top';
hLetters(3).Position = [-0.5, ax(3).Position(4)+0.2, 0];

% 4d. GLM: w/ vs. w/o ramp predictor

% 2.3 Plot R^2 distribution for all units, compare different models.

% fig = figure(Units='inches', Position=[0 0 4 2], DefaultAxesFontSize=12);

paramNames = {'Ordinary', 'Adjusted', 'AdjGeneralized', 'LLR', 'Deviance'};
variantNamesShort = {'const', '+C', '+V', '+R', '+T'};
variantNamesMedium = {'const', '+Cue', '+Vel', '+Ramp', '+Time'};
lineStyles = {'-', '-', '-', '-', ':'};
% variantColors = repmat(linspace(0.9, 0, 5)', 1, 3);
% variantColors = hsl2rgb([zeros(5, 1), linspace(0.2, 1, 5)', 0.5*ones(5, 1)]);

% CDF of R^2
clear ax
ax = nexttile(layout.top.right.top.tl);
hold(ax, 'on')
edges = 0:0.05:1;
centers = (edges(1:end-1) + edges(2:end))*0.5;
for iVariant = 2:nVariants
    N = histcounts(R2(:, iVariant), edges, Normalization='probability');
    plot(ax, edges, [0, cumsum(N)], Color=variantColors(iVariant, :), LineStyle=lineStyles{iVariant}, LineWidth=1.5, DisplayName=variantNamesMedium{iVariant})
end
xlabel('R^2')
ylabel('CDF')
h = legend(ax, Orientation='horizontal');
h.Layout.Tile = 'north';
hold(ax, 'off')
ax.FontSize = p.fontSize;
fontsize(ax, p.fontSize, 'points');
fontname(ax, 'Arial')

hLetters(4) = text(ax, 0, 0, 'd', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetters(4).HorizontalAlignment = 'right';
hLetters(4).VerticalAlignment = 'top';
hLetters(4).Position = [-0.3, ax.Position(4), 0];

% Swarm/boxplot of DeltaR^2
ax = nexttile(layout.top.right.top.tl);
hold(ax, 'on')
dR2 = diff(R2, 1, 2);

x = repmat(1:nVariants-1, [size(dR2, 1), 1]);
x = x(:);
y = dR2(:);
colors = reshape(repmat(variantColors(2:nVariants, :), 1, size(dR2, 1))', 3, [])';
swarmchart(ax, x, y, 1.5, colors, 'filled')

boxplot(ax, dR2, Symbol='.', OutlierSize=0.000001, Color=variantColors(2:nVariants, :), Whisker=0)

xticks(ax, 1:nVariants-1)
xticklabels(ax, variantNamesShort(2:end))
%     xtickangle(ax, 315)
ylabel('\DeltaR^2')
ylim(ax, [0, max(y)+0.01])
xlim(ax, [0,nVariants])
fontsize(ax, p.fontSize, 'points');
fontname(ax, 'Arial')

% Swarm/boxplot of DeltaAIC
ax = nexttile(layout.top.right.top.tl);
hold(ax, 'on')
dAIC = diff(modelCriterion.AIC, 1, 2);

x = repmat(1:nVariants-1, [size(dAIC, 1), 1]);
x = x(:);
y = dAIC(:);
colors = reshape(repmat(variantColors(2:nVariants, :), 1, size(dAIC, 1))', 3, [])';
swarmchart(ax, x, y, 1.5, colors, 'filled')

boxplot(ax, dAIC, Symbol='.', OutlierSize=0.000001, Color=variantColors(2:nVariants, :), Whisker=0)

xticks(ax, 1:nVariants-1)
xticklabels(ax, variantNamesShort(2:end))
%     xtickangle(ax, 315)
ylabel('\DeltaAIC')
ylim(ax, [min(y)-0.01, 0])
xlim(ax, [0,nVariants])
fontsize(ax, p.fontSize, 'points');
fontname(ax, 'Arial')


% Plot fitted vs. observed for 2 example units and population average
SEL = { ...
    84, ...
    39, ...
%     cAcute.isPressDown, ...
%     cAcute.isPressUp, ...
    };
TITLE = { ...
    '', ...
    '', ...
%     sprintf('Example unit (R^2=%.2f)', R2(84, end)), ...
%     sprintf('Example unit (R^2=%.2f)', R2(39, end)), ...
%     sprintf('Population average (N=%d)', nnz(SEL{3})), ...
%     sprintf('Population average (N=%d)', nnz(SEL{4})), ...
    };
LOCATION = { ...
    'southwest', ...
    'northwest', ...
%     'southwest', ...
%     'northwest', ...
    };
SHOW_LEGEND = { ...
    false, ...
    true, ...
%     false, ...
%     false, ...
    };


for i = 1:length(SEL)
    ax = nexttile(layout.top.right.middle.tl);
    hold(ax, 'on')
    clear h
    selUnits = SEL{i}';
    h = gobjects(nVariants - 1, 1);
    for iVariant = 2:nVariants-1
        h(iVariant - 1) = plot(ax, tHat, mean(msrHatAcute(:, iVariant, selUnits), 3, 'omitnan'), ...
            Color=variantColors(iVariant, :), LineWidth=1.5, ...
            DisplayName=variantNamesMedium{iVariant});
    end
    h(end) = plot(ax, tHat, mean(msrObs(:, :, selUnits), 3, 'omitnan'), Color='black', LineStyle='--', LineWidth=2, DisplayName='Obs');
    if SHOW_LEGEND{i}
        l = legend(ax, h, Location=LOCATION{i}, Orientation='horizontal');
        l.Layout.Tile = 'north';
    end
    xlim(ax, [-2, 0])
    title(ax, TITLE{i});
    fontsize(ax, p.fontSize, 'points');
    fontname(ax, 'Arial')
    if i == 1
        hLetters(5) = text(ax, 0, 0, 'e', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
        ax.Units = 'inches';
        hLetters(5).HorizontalAlignment = 'right';
        hLetters(5).VerticalAlignment = 'top';
        hLetters(5).Position = [-0.3, ax.Position(4), 0];
    end
end
xlabel(layout.top.right.middle.tl, 'Time to bar-contact (s)', FontSize=p.fontSize)
ylabel(layout.top.right.middle.tl, 'sp/s', FontSize=p.fontSize)

% Plot fitted vs observed ramp onset times and peak SR
% fig = figure(Units='inches', Position=[0, 0, 4, 1.5]);
SEL = { ...
    cAcute.isPressResponsive, ...
    cAcute.isPressResponsive, ...
    };
XDATA = { ...
    peakAcute, ...
    tOnsetAcute, ...
    };

YDATA = { ...
    peakHatAcute, ...
    tOnsetHatAcute, ...
    };
SHOW_LEGEND = { ...
    false, ...
    true, ...
    };
XLIM = { ...
    'auto', ...
    [-2, 0], ...
    };
UNITS = { ...
    'sp/s', ...
    's', ...
    };

TITLE = {'Peak spike rate (sp/s)', 'Onset time (s)'};

for i = 1:length(SEL)
    ax = nexttile(layout.top.right.bottom.tl);
    hold(ax, 'on')
    clear h
    selUnits = SEL{i}';
    x = XDATA{i}(selUnits);
    y = YDATA{i}(selUnits, :);
    h = gobjects(nVariants - 2, 1);
    for iVariant = 2:nVariants-1
        h(iVariant - 1) = scatter(ax, x, y(:, iVariant), 5, variantColors(iVariant, :), ...
            'filled', DisplayName=variantNamesMedium{iVariant});
    end
    if SHOW_LEGEND{i}
        l = legend(ax, h, AutoUpdate=false, Orientation='horizontal');
        l.Layout.Tile = 'north';
    end
    title(ax, TITLE{i})
    xlim(ax, XLIM{i})
    xl = ax.XLim;
    ylim(ax, xl);
    plot(ax, xl, xl, 'k:')
    fontsize(ax, p.fontSize, 'points');
    fontname(ax, 'Arial')
    if i == 1
        hLetters(6) = text(ax, 0, 0, 'f', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
        ax.Units = 'inches';
        hLetters(6).HorizontalAlignment = 'right';
        hLetters(6).VerticalAlignment = 'top';
        hLetters(6).Position = [-0.3, ax.Position(4), 0];
    end
end
xlabel(layout.top.right.bottom.tl, 'Observed', FontSize=p.fontSize)
ylabel(layout.top.right.bottom.tl, 'Fitted', FontSize=p.fontSize)


% 4g top PETH of two example units (turn on vs. turn off)
unitNames = { ... 
    'daisy13_20220106_Electrode39_Unit1'; ... % Down
    'daisy9_20211013_Electrode23_Unit1'; ... % Up
    };
files = cellfun(@(name) sprintf('C:\\SERVER\\Units\\Lite_NonDuplicate\\%s.mat', name), unitNames, UniformOutput=false);
euEg = EphysUnit.load(files);

YLIM = {[0, 40], [0, 150]};
for iEu = 1:length(euEg)
    ax = nexttile(layout.bottom.tl);
    clear btaEg
    [btaEg.X, btaEg.T, btaEg.N, btaEg.S, btaEg.B] = euEg(iEu).getBinnedTrialAverage('count', p.binnedTrialEdges, 'press', ...
        alignTo='stop', window=[-4, 0.5], resolution=0.1, normalize=false, startBlankWindow=[0, 0.5]);
    btaEg.X = btaEg.X ./ 0.1;
    btaEg.S = btaEg.S ./ 0.1;
    EphysUnit.plotBinnedTrialAverage(ax, btaEg, [-4, 0.5], nsigmas=1, sem=true, showTrialNum=false, colors=btaColors)
    fontsize(ax, p.fontSize, 'points');
    fontname(ax, 'Arial');
    delete(ax.Legend)
    if iEu == 1
        hLetters(7) = text(ax, 0, 0, 'g', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
        ax.Units = 'inches';
        hLetters(7).HorizontalAlignment = 'right';
        hLetters(7).VerticalAlignment = 'top';
        hLetters(7).Position = [-0.5, ax.Position(4) + 0.3, 0];
    end
    xline(ax, 0, 'k--')
    ylim(ax, YLIM{iEu})
    xlim(ax, [-2.5, 0.5])
end


% 4g bottom Ramps are the same regardless of trial length (BTA bootstrap)
clear ax
ax = gobjects(2, 1);
ax(1) = nexttile(layout.bottom.tl);
ax(2) = nexttile(layout.bottom.tl);
EphysUnit.plotBinnedTrialAverage(ax(1), bta.pressDownRaw, [-4, 0.5], nsigmas=1, sem=true, numFormat='%i', colors=btaColors);
EphysUnit.plotBinnedTrialAverage(ax(2), bta.pressUpRaw, [-4, 0.5], nsigmas=1, sem=true, numFormat='%i', colors=btaColors);
delete(ax(2).Legend)
fontsize(ax, p.fontSize, 'points')
fontname(ax, 'Arial')
xlabel(layout.bottom.tl, 'Time to bar-contact (s)', FontSize=p.fontSize)
ylabel(layout.bottom.tl, 'Spike rate (sp/s)', FontSize=p.fontSize)
ax(1).Legend.Orientation = 'horizontal';
ax(1).Legend.Layout.Tile = 'north';
ax(1).Legend.Title.String = 'Time from cue to movement';
ax(1).Legend.AutoUpdate = false;
xline(ax(1), 0, 'k--')
xline(ax(2), 0, 'k--')
ylim(ax, [20, 80])
xlim(ax, [-2.5, 0.5])

copygraphics(fig, ContentType='vector', BackgroundColor='none')

%% Fig S3 top
close all

clear layout
layout.w = 7;
layout.h = 4;
layout.top.h = 2;
layout.bottom.h = 2;
layout.bottom.left.w = 1;
layout.bottom.right.w = 2;

fig = figure(Units='inches', Position=[1 1 layout.w layout.h]);
layout.tl = tiledlayout(fig, layout.top.h + layout.bottom.h, 1, TileSpacing='loose', Padding='compact');

layout.top.tl = tiledlayout(layout.tl, 1, 2*(4+7), TileSpacing='compact');
l = layout.top.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [layout.top.h, 1];

layout.bottom.tl = tiledlayout(layout.tl, 1, layout.bottom.left.w + layout.bottom.right.w, TileSpacing='tight', Padding='tight');
l = layout.bottom.tl; l.Layout.Tile = 1 + layout.top.h; l.Layout.TileSpan = [layout.bottom.h, 1];

layout.bottom.left.tl = tiledlayout(layout.bottom.tl, 1, 1, TileSpacing='tight', Padding='tight');
l = layout.bottom.left.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [1, layout.bottom.left.w];

layout.bottom.right.tl = tiledlayout(layout.bottom.tl, 1, 2, TileSpacing='tight', Padding='tight');
l = layout.bottom.right.tl; l.Layout.Tile = 1 + layout.bottom.left.w; l.Layout.TileSpan = [1, layout.bottom.right.w];

% S3a
% YL = {[30, 90], [10, 70]};
YL = {[-2, 2], [-2, 2]};
SEL = {c.hasPress & c.isPressUp, c.hasPress & c.isPressDown};
TITLE = [sprintf("Increase Units (n=%i)", nnz(SEL{1})), sprintf("Decrease Units (n=%i)", nnz(SEL{2}))];

% Plot ETA aligned to cue (i.e. flinchiness) vs. aligned to lever touch
% close all
AX = gobjects(2, 2);
for iRow = 1:2
    yl = YL{iRow};
    bgAx = axes(layout.top.tl, XTick=[], YTick=[], Box='off');
    bgAx.Layout.Tile = (iRow-1) * (4+7) + 1;
    bgAx.Layout.TileSpan = [1, 4 + 7];
    title(bgAx, TITLE(iRow))
    
    ax = axes(layout.top.tl); AX(iRow, 1) = ax;
    ax.Layout.Tile = (iRow-1) * (4+7) + 1;
    ax.Layout.TileSpan = [1, 4];
    sel = SEL{iRow};
    hold(ax, 'on')
    tt = eta.pressCue.t;
    mu = mean(eta.pressCue.X(sel, :), 1, 'omitnan');
    err = std(eta.pressCue.X(sel, :), 0, 1, 'omitnan');
    plot(ax, eta.pressCue.t, mu, 'k', LineWidth=1.5)
    patch(ax, [tt, flip(tt)], [mu+err, flip(mu-err)], [0.15, 0.15, 0.15], FaceAlpha=0.1, EdgeColor='none');
%     plot(ax, eta.pressCueRaw.t, eta.pressCueRaw.X(sel, :)./0.1, Color=[0.15, 0.15, 0.15, 0.05])
    
    if iRow == 1
        h = gobjects(2, 1);
        h(1) = patch(ax, XData=[-0.8, 0, 0, -0.8], YData=[yl(1), yl(1), yl(2), yl(2)], FaceColor=[0.15, 0.15, 0.15], FaceAlpha=0.1, EdgeColor='none', DisplayName='bar deploy');
        h(2) = patch(ax, XData=[0, 0.1, 0.1, 0], YData=[yl(1), yl(1), yl(2), yl(2)], FaceColor=[0.8, 0.2, 0.2], FaceAlpha=0.3, EdgeColor='none', DisplayName='tone');
    else
        patch(ax, XData=[-0.8, 0, 0, -0.8], YData=[yl(1), yl(1), yl(2), yl(2)], FaceColor=[0.15, 0.15, 0.15], FaceAlpha=0.1, EdgeColor='none', DisplayName='bar deploy');
        patch(ax, XData=[0, 0.1, 0.1, 0], YData=[yl(1), yl(1), yl(2), yl(2)], FaceColor=[0.8, 0.2, 0.2], FaceAlpha=0.3, EdgeColor='none', DisplayName='tone');
    end
    
    ax.Box = 'off';
    xlim(ax, [-0.95, 1])
    xticks(ax, [-0.8, 0, 1])
    xline(ax, 0, '-')
    xline(ax, 1, ':')
    xlabel(ax, 'Time to start cue (s)')
    
    ax = axes(layout.top.tl); AX(iRow, 2) = ax;
    ax.Layout.Tile = (iRow-1) * (4+7) + 4 + 1;
    ax.Layout.TileSpan = [1, 7];
    hold(ax, 'on')
    mu = mean(eta.press.X(sel, :), 1, 'omitnan');
    err = std(eta.press.X(sel, :), 0, 1, 'omitnan');
    plot(ax, eta.press.t, mu, 'k', LineWidth=1.5)
    patch(ax, [eta.press.t, flip(eta.press.t)], [mu+err, flip(mu-err)], [0.15, 0.15, 0.15], FaceAlpha=0.1, EdgeColor='none');
%     plot(ax, eta.pressRaw.t, eta.pressRaw.X(sel, :)./0.1, Color=[0.15, 0.15, 0.15, 0.05])
    xlim(ax, [-3, 0.5])
    xticks(ax, [-3, -2, -1, 0])
    xline(ax, -3, ':')
    xline(ax, 0, '-')
    xlabel(ax, 'Time to bar contact (s)')
    ax.Box = 'off';
    ax.YAxis.Visible = 'off';

    ylim(AX(iRow, :), yl)
    linkaxes(AX(iRow, :), 'y')

    fontsize(AX(iRow, :), p.fontSize, 'points')
end
ylabel(layout.top.tl, 'Spike rate (sp/s)', FontSize=p.fontSize)

% S3b
ax = nexttile(layout.bottom.left.tl);
hold(ax, 'on')

sel = c.hasPress;
% isOutlier = (meta.press < quantile(meta.press(sel), 0.01) | meta.press > quantile(meta.press(sel), 0.99)) | (meta.pressCue < quantile(meta.pressCue(sel), 0.01) | meta.pressCue > quantile(meta.pressCue(sel), 0.99));

% x = meta.pressRaw./0.1 - msr;
% y = meta.pressCueRaw./0.1 - msr;

x = meta.press;
y = meta.pressCue;

hScat = scatter(ax, x(sel & c.isPressResponsive), y(sel & c.isPressResponsive), 4, 'k', 'filled', MarkerFaceAlpha=0.9, DisplayName=sprintf('Responsive (%i units)', nnz(sel&c.isPressResponsive)));
scatter(ax, x(sel & ~c.isPressResponsive), y(sel & ~c.isPressResponsive), 3, 'k', MarkerEdgeAlpha=0.1, DisplayName=sprintf('Unresponsive (n=%i)', nnz(sel&~c.isPressResponsive)));
lgd = legend(ax, hScat, Orientation='horizontal', AutoUpdate=false);
lgd.Layout.Tile = 'north';
xline(ax, 0, ':')
yline(ax, 0, ':')

mdl = fitlm(x(sel), y(sel));
xl = ax.XLim;
plot(ax, xl, mdl.predict(xl(:)), 'k-', LineWidth=1.5)
axis(ax, 'equal')
xlim(ax, [-2, 4])
ylim(ax, [-2, 4])

xlabel(layout.bottom.left.tl, 'Peri-reach (a.u.)', FontSize=p.fontSize)
ylabel('Peri-cue (a.u.)')
title('Cue vs. reach response')
% ylabel({'Peri-cue', 'response (a.u.)'})

fontsize(ax, p.fontSize, 'points')

% S3c/d (correct/incorrect press trial video-tracked paw/spine/lick)
SELTRIALS = {~fAll.press.correct, fAll.press.correct};
RESULTNAMES = {'incorrect', 'correct'};
RESULTNAMESDISP = {'Incorrect', 'Correct'};
FTNAMES = {'spine_yVel', 'handContra_xVel', 'tongue'};
FTDISPNAMES = {'Spine velocity', 'Contra forepaw velocity', 'Lick probability'};
% COLORS = {'red', 'black', 'blue'};
COLORS = arrayfun(@(i) getColor(i, 3, 0.7), [2, 1, 3], UniformOutput=false);
YYAXIS = {'left', 'left', 'right'};
for iResult = 1:2
    ax = nexttile(layout.bottom.right.tl);
    colororder(ax, [[0.15, 0.15, 0.15]; getColor(3, 3, 0.7)])
    hold(ax, 'on');
    selTrials = SELTRIALS{iResult};
    t = fAll.press.t;
    n = nnz(selTrials);
    hFt = gobjects(1, length(FTNAMES));
    for iFt = 1:length(FTNAMES)
        yyaxis(ax, YYAXIS{iFt})
        X = fAll.press.(FTNAMES{iFt});
        X = X(selTrials, :);
        mu = mean(X, 1, 'omitnan');
        sd = std(X, 0, 1, 'omitnan');
        col = COLORS{iFt};
        hFt(iFt) = plot(ax, t, mu, Color=col, LineStyle='-', LineWidth=1.5, DisplayName=FTDISPNAMES{iFt});
        sel = ~isnan(mu + sd);
        patch(ax, [t(sel), flip(t(sel))], [mu(sel)-sd(sel), flip(mu(sel)+sd(sel))], 'r', ...
            LineStyle='-', FaceAlpha=0.075, FaceColor=col, EdgeAlpha=0.075, EdgeColor=col)
        xline(ax, 0, 'k--')
    end
    xlim(ax, [-2, 2])

    yyaxis(ax, 'left')
    ylabel(ax, 'Velocity (a.u.)')
    ylim(ax, [-5, 5])

    yyaxis(ax, 'right')
    ylabel(ax, 'Lick probability')
    ylim(ax, [-1.5, 1.5])
    title(ax, RESULTNAMESDISP{iResult})
end
xlabel(layout.bottom.right.tl, 'Time to bar contact (s)', FontSize=p.fontSize)

lgd = legend(ax, hFt, Orientation='horizontal');
lgd.Layout.Tile = 'north';

lgd = legend(AX(1, 2), h, Location='northwest', AutoUpdate='off', FontSize=p.fontSize);
lgd.Position(1) = lgd.Position(1) - 0.15;

copygraphics(fig, ContentType='vector', BackgroundColor='none')

%% S3 bottom
p.fontSize = 9;

nBTABins = length(p.binnedTrialEdgesFine) - 1;
btaColors = hsl2rgb([linspace(0.8, 0, nBTABins)'.^1.3, linspace(0.8, 0.6, nBTABins)', 0.5*ones(nBTABins, 1)]);

clear layout
layout.w = 7;
layout.h = 4;

p.lineWidth = 1.5;

fig = figure(Units='inches', Position=[1, 1, layout.w, layout.h], DefaultAxesFontSize=p.fontSize);
layout.tl = tiledlayout(fig, 2, 2, TileSpacing='compact', TileIndexing='columnmajor');

AX = gobjects(2, 2);

% S3a 
ax = nexttile(layout.tl);
AX(1, 1) = ax;
EphysUnit.plotBinnedTrialAverage(ax, btaUp, [-4, 0.5], nsigmas=1, sem=true, showTrialNum=false, numFormat='%i', colors=btaColors, lineWidth=1.2);
hLgd = ax.Legend;
hLgd.Layout.Tile = 'north';
hLgd.Orientation = 'horizontal';
hLgd.Title.String = 'Time from cue to movement';
hLgd.AutoUpdate = false;
xline(ax, 0, 'k--')
title(ax, sprintf('Significant increase units (n=%i)', nnz(c.isPressBTADifferentUp & c.isPressUp)))


% S3b
ax = nexttile(layout.tl);
AX(2, 1) = ax;
EphysUnit.plotBinnedTrialAverage(ax, btaDown, [-4, 0.5], nsigmas=1, sem=true, showTrialNum=false, numFormat='%i', colors=btaColors, lineWidth=1.2);
xline(ax, 0, 'k--')
title(ax, sprintf('Significant decrease units (n=%i)', nnz(c.isPressBTADifferentDown & c.isPressDown)))
delete(ax.Legend);

% S3c
ax = nexttile(layout.tl);
AX(1, 2) = ax;
EphysUnit.plotBinnedTrialAverage(ax, btaNulUp, [-4, 0.5], nsigmas=1, sem=true, showTrialNum=false, numFormat='%i', colors=btaColors, lineWidth=1.2);
xline(ax, 0, 'k--')
title(ax, sprintf('Null increase units (n=%i)', nnz(~c.isPressBTADifferentUp & c.isPressUp)))
delete(ax.Legend);

% S3d
ax = nexttile(layout.tl);
AX(2, 2) = ax;
EphysUnit.plotBinnedTrialAverage(ax, btaNulDown, [-4, 0.5], nsigmas=1, sem=true, showTrialNum=false, numFormat='%i', colors=btaColors, lineWidth=1.2);
xline(ax, 0, 'k--')
title(ax, sprintf('Null decrease units (n=%i)', nnz(~c.isPressBTADifferentDown & c.isPressDown)))
delete(ax.Legend);

ax = AX;
xlim(ax, [-2.5, 0.5]);
ylim(ax(1, :), [30, 100])
ylim(ax(2, :), [10, 80])
xlabel(ax, '')
ylabel(ax, '')
xlabel(layout.tl, 'Time to bar contact (s)', FontSize=p.fontSize)
ylabel(layout.tl, 'Spike rate (sp/s)', FontSize=p.fontSize)

for i = 1:4
    patch(ax(i), [-2, -0.2, -0.2, -2], [0, 0, 100, 100], [0.15, 0.15, 0.15], FaceAlpha=0.15, EdgeColor='none')
end

copygraphics(fig, ContentType='vector', BackgroundColor='none')

% %% S3 bottom (lick)
% p.fontSize = 9;
% 
% nBTABins = length(p.binnedTrialEdgesFine) - 1;
% btaColors = hsl2rgb([linspace(0.8, 0, nBTABins)'.^1.3, linspace(0.8, 0.6, nBTABins)', 0.5*ones(nBTABins, 1)]);
% 
% clear layout
% layout.w = 7;
% layout.h = 4;
% 
% p.lineWidth = 1.5;
% 
% fig = figure(Units='inches', Position=[1, 1, layout.w, layout.h], DefaultAxesFontSize=p.fontSize);
% layout.tl = tiledlayout(fig, 2, 2, TileSpacing='compact', TileIndexing='columnmajor');
% 
% AX = gobjects(2, 2);
% 
% % S3a 
% ax = nexttile(layout.tl);
% AX(1, 1) = ax;
% EphysUnit.plotBinnedTrialAverage(ax, btaLickUp, [-4, 0.5], nsigmas=1, sem=true, showTrialNum=false, numFormat='%i', colors=btaColors, lineWidth=1.2);
% hLgd = ax.Legend;
% hLgd.Layout.Tile = 'north';
% hLgd.Orientation = 'horizontal';
% hLgd.Title.String = 'Time from cue to movement';
% hLgd.AutoUpdate = false;
% xline(ax, 0, 'k--')
% title(ax, sprintf('Significant increase units (n=%i)', nnz(c.isLickBTADifferentUp & c.isLickUp & c.hasPress)))
% 
% 
% % S3b
% ax = nexttile(layout.tl);
% AX(2, 1) = ax;
% EphysUnit.plotBinnedTrialAverage(ax, btaLickDown, [-4, 0.5], nsigmas=1, sem=true, showTrialNum=false, numFormat='%i', colors=btaColors, lineWidth=1.2);
% xline(ax, 0, 'k--')
% title(ax, sprintf('Significant decrease units (n=%i)', nnz(c.isLickBTADifferentDown & c.isLickDown & c.hasPress)))
% delete(ax.Legend);
% 
% % S3c
% ax = nexttile(layout.tl);
% AX(1, 2) = ax;
% EphysUnit.plotBinnedTrialAverage(ax, btaLickNulUp, [-4, 0.5], nsigmas=1, sem=true, showTrialNum=false, numFormat='%i', colors=btaColors, lineWidth=1.2);
% xline(ax, 0, 'k--')
% title(ax, sprintf('Null increase units (n=%i)', nnz(~c.isLickBTADifferentUp & c.isLickUp & c.hasPress)))
% delete(ax.Legend);
% 
% % S3d
% ax = nexttile(layout.tl);
% AX(2, 2) = ax;
% EphysUnit.plotBinnedTrialAverage(ax, btaLickNulDown, [-4, 0.5], nsigmas=1, sem=true, showTrialNum=false, numFormat='%i', colors=btaColors, lineWidth=1.2);
% xline(ax, 0, 'k--')
% title(ax, sprintf('Null decrease units (n=%i)', nnz(~c.isLickBTADifferentDown & c.isLickDown & c.hasPress)))
% delete(ax.Legend);
% 
% ax = AX;
% xlim(ax, [-2.5, 0.5]);
% ylim(ax(1, 2), [30, 100])
% ylim(ax(2, 2), [10, 80])
% ylim(ax(1, 1), [30, 200])
% ylim(ax(2, 1), [0, 100])
% xlabel(ax, '')
% ylabel(ax, '')
% xlabel(layout.tl, 'Time to spout contact (s)', FontSize=p.fontSize)
% ylabel(layout.tl, 'Spike rate (sp/s)', FontSize=p.fontSize)
% 
% for i = 1:4
%     patch(ax(i), [-2, -0, -0, -2], [0, 0, 500, 500], [0.15, 0.15, 0.15], FaceAlpha=0.15, EdgeColor='none')
% end
% 
% copygraphics(fig, ContentType='vector', BackgroundColor='none')
