%% 
read_DLC_data

%% 4a. DeepLabCut camera screenshot
%% 4b. Histograms: movement initiation vs. activity initiation (turn on) vs. activity initiation (turn off)
% To do: recalculate 'tst' (true movement onset time).
tst = [-0.233333333333334	-0.266666666666680	-0.0999999999999943	-0.200000000000045	-0.533333333333303	-0.166666666666629	-0.166666666666629	-0.100000000000023	-0.133333333333326	-0.100000000000023	-0.133333333333326	-0.133333333333326	-0.133333333333326	-0.100000000000023	-0.166666666666629	-0.566666666666720	-0.133333333333326	-0.133333333333326	-0.200000000000045	-0.200000000000045	-0.166666666666629	-0.133333333333439	-0.133333333333439	-0.433333333333394	-0.166666666666742	-0.133333333333439	-0.133333333333212	-0.166666666666515	-0.0999999999999091	-0.166666666666515	-0.133333333333212	-0.0999999999999091	-0.133333333333212	-0.133333333333212	-0.133333333333212	-0.133333333333212	-0.199999999999818	-0.133333333333212	-0.133333333333212	-0.0999999999999091	-0.0999999999999091	-0.166666666666515	-0.133333333333212	-0.333333333333485	-0.166666666666515	-0.199999999999818	0	-0.0999999999999091	0	-0.166666666666970	-0.133333333333212	-0.166666666666970	-0.166666666666970	-0.199999999999818	-0.233333333333576	-0.166666666666970	-0.199999999999818	-0.166666666666970	-0.133333333333212	0	-0.133333333333212	-0.166666666666657	-0.0999999999999943	0	-0.300000000000011	-0.199999999999989	-0.266666666666652	-0.166666666666686	-0.133333333333326	0	-0.100000000000023	0	0	-0.100000000000023	-0.166666666666629	-0.333333333333371	-0.166666666666742	-0.166666666666515	-0.166666666666515	-0.133333333333212	-0.266666666666879	-0.233333333333121	-0.400000000000091	0	0	-0.233333333333121	-0.133333333333212	-0.233333333333121	-0.300000000000182	-0.133333333333212	-0.0999999999999091	0	-0.300000000000182	-0.400000000000091	-0.166666666666970	-0.100000000000364	-0.233333333333576	-0.199999999999818	-0.133333333333212	-0.133333333333212	-0.166666666666970	-0.199999999999989	-0.0999999999999943	-0.133333333333326	-0.100000000000023	-0.133333333333326	-0.133333333333326	-0.300000000000011	-0.166666666666686	-0.0666666666667197	-0.133333333333326	-0.100000000000023	-0.100000000000023	-0.133333333333326	-0.0666666666667197	-0.133333333333326	0	-0.0666666666666060	-0.0666666666666060	-0.133333333333212	-0.166666666666515	-0.133333333333212	-0.0999999999999091	-0.266666666666879	-0.0666666666666060	-0.0999999999999091	-0.166666666666515	-0.133333333333212	-0.233333333333121	-0.133333333333212	-0.0999999999999091	-0.0666666666666060	-0.0666666666666060	-0.166666666666515	-0.0999999999999091	-0.0999999999999091	-0.199999999999818	-0.133333333333212	-0.133333333333212	-0.133333333333212	0	0	-0.166666666666686	0	-0.166666666666686	-0.100000000000023	-0.133333333333326	0	-0.100000000000023	-0.133333333333326	0	-0.0666666666667197	-0.233333333333349	-0.133333333333326	-0.166666666666742	-0.133333333333439	-0.200000000000045	-0.233333333333349	-0.166666666666515	-0.266666666666879	-0.133333333333212	-0.0999999999999091	-0.0333333333333030	-0.0333333333333030	-0.0666666666666060	-0.233333333333121	-0.133333333333212	-0.233333333333576	-0.0333333333337578	0	0	0	-0.433333333333337	-0.199999999999989	-0.333333333333314	0	-0.133333333333326	-0.600000000000023	0	0	-0.433333333333280	-0.433333333333280	-0.0666666666667197	-0.100000000000023	0	-0.166666666666515	-0.233333333333121	-0.166666666666515	-0.0666666666666060	-0.199999999999818	0	0	-0.233333333333576	-0.133333333333212	-0.199999999999818	-0.100000000000364	-0.266666666666680	-0.0999999999999943	-0.166666666666657	0	-0.166666666666629	-0.199999999999818	-0.0999999999999091	-0.233333333333121	-0.266666666666879	-0.166666666666515	0	-0.199999999999818	-0.199999999999818	0	-0.166666666666515	-0.166666666666515	-0.199999999999818	-0.233333333333121	-0.199999999999818	-0.199999999999818	0	-0.133333333333212	-0.300000000000182	0	0	-0.333333333333030	-0.266666666666424	-0.300000000000182	0	0	-0.266666666666424	-0.133333333333212	-0.233333333333576	-0.199999999999818	-0.0666666666666060	-0.600000000000001	-0.133333333333326	-0.166666666666629	0	0	0	-0.633333333333212	-0.233333333333121	-0.533333333333303	0	0	0	0	-0.199999999999818	-0.133333333333212	-0.300000000000182	-0.333333333333030	-0.666666666666970	0	0	-0.766666666666424	0	-0.199999999999818	-0.500000000000000	-0.233333333333576];

p.etaSortWindow = [-3, 0];
p.etaSignWindow = [-0.3, 0];
p.etaLatencyThresholdPos = 0.5;
p.etaLatencyThresholdNeg = 0.5;

close all
pressLatency = NaN(size(eu));
lickLatency = NaN(size(eu));
[ax, ~, ~, pressLatency(c.hasPress)] = EphysUnit.plotETA(eta.press, c.hasPress, xlim=[-4,0], clim=[-2, 2], ...
    sortWindow=p.etaSortWindow, signWindow=p.etaSignWindow, ...
    sortThreshold=p.etaLatencyThresholdPos, negativeSortThreshold=p.etaLatencyThresholdNeg); 
title('Pre-reach PETH')
xlabel('Time to touch (s)')
ax.Parent.Position(3) = 0.25;

[ax, ~, ~, lickLatency(c.hasLick)] = EphysUnit.plotETA(eta.lick, c.hasLick, xlim=[-4,0], clim=[-2, 2], ...
    sortWindow=p.etaSortWindow, signWindow=p.etaSignWindow, ...
    sortThreshold=p.etaLatencyThresholdPos, negativeSortThreshold=p.etaLatencyThresholdNeg); 
title('Pre-lick PETH')
xlabel('Time to spout-contact (s)')
ax.Parent.Position(3) = 0.25;

latency.press = pressLatency;
latency.lick = lickLatency;
latency.contraPaw = tst;
latency.pRankSum.pressVsContraPaw = ranksum(tst, pressLatency(c.isPressResponsive), tail='right');
latency.pRankSum.pressUpVsContraPaw = ranksum(tst, pressLatency(c.isPressUp), tail='right');
latency.pRankSum.pressDownVsContraPaw = ranksum(tst, pressLatency(c.isPressDown), tail='right');
latency.pRankSum.pressVsLick = ranksum(lickLatency, pressLatency, tail='right');
latency.pRankSum.pressVsContraPawPlus200 = ranksum(tst-0.2, pressLatency(c.isPressResponsive), tail='right');
latency.pRankSum.pressUpVsPressDown = ranksum(pressLatency(c.isPressUp), pressLatency(c.isPressDown), tail='both');

fprintf(1, 'Median pre-press spiking onset latency = %.1f ms \n', median(latency.press(c.isPressResponsive)*1000, 'omitnan'))
fprintf(1, 'Median pre-press spiking onset latency (excited) = %.1f ms \n', median(latency.press(c.isPressUp)*1000, 'omitnan'))
fprintf(1, 'Median pre-press spiking onset latency (suppressed) = %.1f ms \n', median(latency.press(c.isPressDown)*1000, 'omitnan'))
fprintf(1, 'Median pre-lick spiking onset latency = %.1f ms \n', median(latency.lick(c.isLickResponsive)*1000, 'omitnan'))
fprintf(1, 'Median contralateral paw movement onset latency = %.3f ms \n', median(latency.contraPaw*1000, 'omitnan'))
fprintf(1, 'Press spiking precedes paw: One-tailed ranksum test p = %g\n', latency.pRankSum.pressVsContraPaw)
fprintf(1, 'Press spiking precedes paw (-200ms): One-tailed ranksum test p = %g\n', latency.pRankSum.pressVsContraPawPlus200)
fprintf(1, 'Excited press spiking precedes paw: One-tailed ranksum test p = %g\n', latency.pRankSum.pressUpVsContraPaw)
fprintf(1, 'Inhibited press spiking precedes paw: One-tailed ranksum test p = %g\n', latency.pRankSum.pressDownVsContraPaw)
fprintf(1, 'Press spiking precedes lick spiking: One-tailed ranksum test p = %g\n', latency.pRankSum.pressVsLick)
fprintf(1, 'Press excided vs. press supressed: two-tailed ranksum test p = %g\n', latency.pRankSum.pressUpVsPressDown)


fig = figure(DefaultAxesFontSize=p.fontSize, Units='inches', Position=[0 0 3 5]);

ax(1) = subplot(4, 1, 1);
histogram(latency.press(c.isPressResponsive), (-2:0.1:0), Normalization='probability', FaceColor='black')
title('SNr response onset')
legend(sprintf('%d units', nnz(~isnan(pressLatency(c.isPressResponsive & c.hasPress)))), Location='northwest')
ylabel('Probability')

ax(2) = subplot(4, 1, 2);
histogram(latency.press(c.isPressUp), (-2:0.1:0), Normalization='probability', FaceColor='red')
title('SNr response onset (excited)')
legend(sprintf('%d units', nnz(~isnan(pressLatency(c.isPressUp & c.hasPress)))), Location='northwest')
ylabel('Probability')

ax(3) = subplot(4, 1, 3);
histogram(latency.press(c.isPressDown), (-2:0.1:0), Normalization='probability', FaceColor='blue')
title('SNr response onset (suppressed)')
legend(sprintf('%d units', nnz(~isnan(pressLatency(c.isPressDown & c.hasPress)))), Location='northwest')
ylabel('Probability')

ax(4) = subplot(4, 1, 4);
histogram(latency.contraPaw, (-2:0.1:0), Normalization='probability', FaceColor='black')
legend(sprintf('%g trials, %g sessions', nnz(~isnan(tst)), 7), Location='northwest')
title('Forepaw movement onset')
xlabel('Time to touch (s)')
ylabel('Probability')
yticks(ax, [])
ylabel(ax, [])
fontsize(fig, p.fontSize, 'points')
fontname(fig, 'Arial')

%% 4c. GLM: w/ vs. w/o ramp predictor
fit_GLM
%%

%% 2.3 Plot R^2 distribution for all units, compare different models.
edges = 0:0.05:1;
centers = (edges(1:end-1) + edges(2:end))*0.5;

fig = figure(Units='inches', Position=[0 0 4 2], DefaultAxesFontSize=12);

paramNames = {'Ordinary', 'Adjusted', 'AdjGeneralized', 'LLR', 'Deviance'};
np = length(paramNames);

clear ax
np = 1;
for ip = 1:np
    ax = subplot(np, 2, 2*(ip-1)+1);
    hold(ax, 'on')
    for iVariant = 2:nVariants
        N = histcounts(R2(:, iVariant), edges, Normalization='probability');
        plot(ax, edges, [0, cumsum(N)], Color=getColor(iVariant-1, nVariants-1), LineWidth=1.5, DisplayName=variantNames{iVariant})
    end
%     xlabel(sprintf('R^2 %s', paramNames{ip}))
    xlabel('R^2')
    ylabel('Cumulative probability')
    legend(ax, Location='southeast');
    hold(ax, 'off')
    ax.FontSize = p.fontSize;

    ax = subplot(np, 2, 2*(ip-1)+2);
    hold(ax, 'on')
    dR2 = diff(R2, 1, 2);

    x = repmat(1:nVariants-1, [size(dR2, 1), 1]);
    x = x(:);
    y = dR2(:);
    swarmchart(ax, x, y, 1, 'filled', 'k')
    
    boxplot(ax, dR2, Symbol='.', OutlierSize=0.000001, Color='k', Whisker=0)
    
    xticks(ax, 1:nVariants-1)
    xticklabels(ax, variantNames(2:end))
    xtickangle(ax, 315)
    ylabel('\DeltaR^2')
    ylim(ax, [0, max(y)+0.01])
    xlim(ax, [0,nVariants])
end
fontsize(fig, p.fontSize, 'points');
fontname(fig, 'Arial')
%% Plot fitted vs. observed for 2 example units and population average
SEL = { ...
    84, ...
    39, ...
%     cAcute.isPressDown, ...
%     cAcute.isPressUp, ...
    };
TITLE = { ...
    sprintf('Example unit (R^2=%.2f)', R2(84, end)), ...
    sprintf('Example unit (R^2=%.2f)', R2(39, end)), ...
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


% fig = figure(Units='inches', Position=[0, 0, 4, 2.5]);
fig = figure(Units='inches', Position=[0, 0, 4, 1.5]);
for i = 1:length(SEL)
%     ax = subplot(2, 2, i);
    ax = subplot(1, 2, i);
    hold(ax, 'on')
    clear h
    sel = SEL{i}';
    h(1) = plot(ax, t, mean(msrObs(:, :, sel), 3, 'omitnan'), 'k:', LineWidth=3, DisplayName='Observed');
    for iVariant = 2:nVariants-1
        h(iVariant) = plot(ax, t, mean(msrHatAcute(:, iVariant, sel), 3, 'omitnan'), ...
            Color=getColor(iVariant-1, nVariants-1), LineWidth=1.5, ...
            DisplayName=variantNames{iVariant});
    end
    if SHOW_LEGEND{i}
        legend(ax, h, Location=LOCATION{i})
    end
    xlim(ax, [-2, 0])
    xlabel(ax, sprintf('Time to %s (s)', 'touch'))
    if i == 1
        ylabel(ax, 'Spike rate (sp/s)')
    end
    title(ax, TITLE{i});
    fontsize(ax, p.fontSize, 'points');
    fontname(ax, 'Arial')
end

%% Plot fitted vs observed ramp onset times and peak SR
fig = figure(Units='inches', Position=[0, 0, 4, 1.5]);
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

TITLE = {'Peak spike rate', 'Onset time'};

for i = 1:length(SEL)
    ax = subplot(1, 2, i);
    hold(ax, 'on')
    clear h
    sel = SEL{i}';
    x = XDATA{i}(sel);
    y = YDATA{i}(sel, :);
    h = gobjects(nVariants - 2, 1);
    for iVariant = 2:nVariants-1
        h(iVariant - 1) = scatter(ax, x, y(:, iVariant), 5, getColor(iVariant-1, nVariants-1), ...
            'filled', DisplayName=variantNames{iVariant});
    end
    if SHOW_LEGEND{i}
        legend(ax, h, Location='northwest', AutoUpdate=false)
    end
    title(ax, TITLE{i})
    xlim(ax, XLIM{i})
    xl = ax.XLim;
    ylim(ax, xl);
    plot(ax, xl, xl, 'k:')
    xlabel(ax, sprintf('Observed (%s)', UNITS{i}));
    ylabel(ax, sprintf('Predicted (%s)', UNITS{i}));
    fontsize(ax, p.fontSize, 'points');
    fontname(ax, 'Arial')
end




%% 4d. Ramps are the same regardless of trial length (BTA bootstrap)
[bta.pressUpRaw.X, bta.pressUpRaw.T, bta.pressUpRaw.N, bta.pressUpRaw.S, bta.pressUpRaw.B] = eu(c.isPressUp).getBinnedTrialAverage('rate', p.binnedTrialEdges, 'press', 'window', [-10, 1], 'normalize', false);
[bta.pressDownRaw.X, bta.pressDownRaw.T, bta.pressDownRaw.N, bta.pressDownRaw.S, bta.pressDownRaw.B] = eu(c.isPressDown).getBinnedTrialAverage('rate', p.binnedTrialEdges, 'press', 'window', [-10, 1], 'normalize', false);
%%
figure(Units='inches', Position=[0, 0, 7, 2.5]);
clear ax
ax(1) = subplot(2, 1, 1);
ax(2) = subplot(2, 1, 2);
EphysUnit.plotBinnedTrialAverage(ax(1), bta.pressDownRaw, [-8, 0], nsigmas=1, sem=true);
EphysUnit.plotBinnedTrialAverage(ax(2), bta.pressUpRaw, [-8, 0], nsigmas=1, sem=true);
title(ax(1), sprintf('Reach-inhibited (populaton-average, N=%i)', nnz(c.isPressDown)))
title(ax(2), sprintf('Reach-excited (populaton-average, N=%i)', nnz(c.isPressUp)))
xlabel(ax(2), 'Time to touchbar-contact (s)')
ylabel(ax, 'Spike rate (sp/s)')
h = legend(ax(1)); h.Location='southwest';
fontsize(ax, p.fontSize, 'points')
fontname(ax, 'Arial')
clear ax