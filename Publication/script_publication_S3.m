%% Load ephysUnits
load_ephysunits;

%% Bootstrap modulation direction
boot_response_dir;

%% Plot S2 (baseline vs. dir)

close all

sz = 5;

fig = figure(Units='inches', Position=[0, 0, 6.5, 6], DefaultAxesFontSize=p.fontSize);
ax = subplot(2, 2, 2);

nBins = 15;
edges = linspace(0, max(msr, [], 'all', 'omitnan'), nBins);
NAll = histcounts(msr(c.hasPress), edges, 'Normalization', 'probability');
NUp = histcounts(msr(c.isPressUp), edges, 'Normalization', 'probability');
NDown = histcounts(msr(c.isPressDown), edges, 'Normalization', 'probability');
centers = 0.5*(edges(2:end) + edges(1:end-1));
hold(ax, 'on')
plot(ax, centers, cumsum(NDown), 'LineWidth', 2, 'Color', 'blue', 'DisplayName', sprintf('suppressed (N=%g)', nnz(c.isPressDown)))
plot(ax, centers, cumsum(NUp), 'LineWidth', 2, 'Color', 'red', 'DisplayName', sprintf('excited (N=%g)', nnz(c.isPressUp)));
hold(ax, 'off')
legend(ax, Location='southeast');
ylim(ax, [0, 1])
xlabel(ax, 'Baseline spike rate (sp/s)')
ylabel(ax, 'Cumulative probability')
[ks.h, ks.p] = kstest2(msr(c.isPressDown), msr(c.isPressUp), Tail='larger');
title(ax, 'Suppressed vs. Excited (CDF)')
fprintf(1, 'Median: suppressed=%.1f, excited=%.1f, one-tailed rank-sum (suppressed > excited) p=%g\n', median(msr(c.isPressDown)), median(msr(c.isPressUp)), ranksum(msr(c.isPressDown), msr(c.isPressUp), tail='right'))
fprintf(1, 'One-tailed K-S test (suppressed > excited), p=%g\n', ks.p)



ax = subplot(2, 2, 1);
histogram(ax, msr(c.isPressDown), edges, FaceColor='blue')
xlabel(ax, 'Baseline spike rate (sp/s)'), ylabel(ax, 'Count')
legend(ax, sprintf('N=%g', nnz(c.isPressDown)), Location='northeast')
title(ax, 'Suppressed population')

ax = subplot(2, 2, 3);
histogram(ax, msr(c.isPressUp), edges, FaceColor='red');
xlabel(ax, 'Baseline spike rate (sp/s)'), ylabel(ax, 'Count')
legend(ax, sprintf('N=%g', nnz(c.isPressUp)), Location='northeast')
title(ax, 'Excited population')

ax = subplot(2, 2, 4);
hold(ax, 'on')
isResp = c.isPressResponsive;
isNonResp = c.hasPress & ~c.isPressResponsive;
clear h
sel = c.isPressResponsive;
h(1) = scatter(ax, msr(sel), meta.pressRaw(sel)./0.1 - msr(sel), sz, 'black', 'filled', DisplayName=sprintf('N=%g', nnz(sel)));
xl = ax.XLim;
plot(ax, ax.XLim, [0, 0], 'k--', LineWidth=1.5);
hold(ax, 'off')
xlabel(ax, 'Baseline spike rate (sp/s)'), ylabel(ax, 'Pre-move response (\Deltasp/s)')
legend(ax, h, Location='best')
title(ax, 'Responsive population')


fontsize(fig, p.fontSize, 'points')
fontname(fig, 'Arial')
