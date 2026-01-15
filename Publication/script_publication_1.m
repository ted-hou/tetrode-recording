%% Only run once (make bs.mat and save to C:\SERVER_PRIVATE)
% make_behavior_sessions

%% Load bs.mat, do some processing
load_behavior_sessions

%% Load all units
load_ephysunits;
% boot_response_dir;
% load('C:\SERVER\Units\boot_20241024_perimovement_0.3_0.mat')
etaFine.press = eu.getETA('count', 'press', [-4, 2], minTrialDuration=2, normalize=[-4, -2], resolution=0.025);
% etaFine.lick = eu.getETA('count', 'lick', [-4, 2], minTrialDuration=2, normalize=[-4, -2], resolution=0.025);
%% Report baseline spike rate tests, stats
bsr = NaN(length(eta.press.stats), 1);
bsr(c.hasPress) = [eta.press.stats(c.hasPress).mean]./0.1;
assert(nnz(~isnan(bsr)) == nnz(c.hasPress))
[~, test.ks2BaselineSpikeRateUpVsDown.p] = kstest2(bsr(c.hasPress & c.isPressUp), bsr(c.hasPress & c.isPressDown), Tail='larger');
test.ranksumBaselineSpikeRateUpVsDown.p = ranksum(bsr(c.hasPress & c.isPressUp), bsr(c.hasPress & c.isPressDown), tail='left');

% ROC (using bsr to classify inc/dec)
roc = rocmetrics(categorical(c.isPressDown(c.isPressResponsive), [false, true], ["inc", "dec"]), bsr(c.isPressResponsive), "dec");
ax = axes(figure());
plot(ax, roc);
clear ax

fprintf('%i SNr units, median baseline spike rate [-4, -2] = %.2f, median absolute deviation = %.2f.\n', nnz(c.hasPress), median(bsr(c.hasPress)), mad(bsr(c.hasPress), 1));
fprintf(['%i/%i (%i%%) is reach-modulated, ' ...
    '%i(%i%%) increase, ' ...
    '%i(%i%%) decrease.\n'], ...
    nnz(c.isPressResponsive), nnz(c.hasPress), round(100*nnz(c.isPressResponsive)/nnz(c.hasPress)), ...
    nnz(c.isPressUp), round(100*nnz(c.isPressUp)/nnz(c.isPressResponsive)), ...
    nnz(c.isPressDown), round(100*nnz(c.isPressDown)/nnz(c.isPressResponsive)))
fprintf(['Baseline spike rate (median+-mad): decrease=%.2f+-%.2f, increase=%.2f+_%.2f;\n' ...
    'one-tailed ranksum p<%.4f, one-tailed KS p<%.4f, AUC=%.2f.\n'], median(bsr(c.isPressDown)), mad(bsr(c.isPressDown), 1), median(bsr(c.isPressUp)), mad(bsr(c.isPressUp), 1), test.ranksumBaselineSpikeRateUpVsDown.p, test.ks2BaselineSpikeRateUpVsDown.p, roc.AUC)



%% Fig 1.


p.fontSize = 9;
p.lineWidth = 1.5;

clear layout
layout.w = 7;
layout.h = 8;
layout.top.h = 5;
layout.middle.h = 6;
layout.bottom.h = 10;

layout.top.left.w = 9;
layout.top.middle.w = 9;
layout.top.right.w = 9;
layout.top.rightcb.w = 1;
layout.middle.top.h = 9;
layout.middle.bottom.h = 7;
layout.bottom.left.w = 3;
layout.bottom.right.w = 2;
layout.bottom.right.hh = [1, 4, 4, 4];

% close all
fig = figure(Units='inches', Position=[0, 0, layout.w, layout.h]);

layout.tl = tiledlayout(fig, layout.top.h + layout.middle.h + layout.bottom.h, 1, TileSpacing='loose');

layout.top.tl = tiledlayout(layout.tl, 1, layout.top.left.w + layout.top.middle.w + layout.top.right.w + layout.top.rightcb.w, TileSpacing='compact', Padding='compact');
l = layout.top.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [layout.top.h, 1];

layout.middle.tl = tiledlayout(layout.tl, layout.middle.top.h + layout.middle.bottom.h, 1, TileSpacing='compact', Padding='loose');
l = layout.middle.tl; l.Layout.Tile = 1 + layout.top.h; l.Layout.TileSpan = [layout.middle.h, 1];

layout.middle.top.tl = tiledlayout(layout.middle.tl, 1, 2, TileSpacing='compact');
l = layout.middle.top.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [layout.middle.top.h, 1];

layout.middle.bottom.tl = tiledlayout(layout.middle.tl, 1, 2, TileSpacing='compact');
l = layout.middle.bottom.tl; l.Layout.Tile = 1 + layout.middle.top.h; l.Layout.TileSpan = [layout.middle.bottom.h, 1];

layout.bottom.tl = tiledlayout(layout.tl, 1, layout.bottom.left.w + layout.bottom.right.w, TileSpacing='loose');
l = layout.bottom.tl; l.Layout.Tile = 1 + layout.top.h + layout.middle.h; l.Layout.TileSpan = [layout.bottom.h, 1];

layout.bottom.right.tlp = tiledlayout(layout.bottom.tl, sum(layout.bottom.right.hh), 1, TileSpacing='tight', Padding='tight');
l = layout.bottom.right.tlp; l.Layout.Tile = 1 + layout.bottom.left.w; l.Layout.TileSpan = [1, layout.bottom.right.w];
layout.bottom.right.tl = tiledlayout(layout.bottom.right.tlp, 3, 1, TileSpacing='tight', Padding='tight');
l = layout.bottom.right.tl; l.Layout.Tile = 1 + layout.bottom.right.hh(1); l.Layout.TileSpan = [sum(layout.bottom.right.hh(2:end)), 1];

% 1c'. Behavior
edges = 0:0.5:10;

% Plot aggregate histograms as line plots
ax = nexttile(layout.top.tl, 1 + layout.top.left.w + layout.top.middle.w, [1, layout.top.right.w]);
centers = 0.5*(edges(2:end) + edges(1:end-1));
hold(ax, 'on')
ndayshown = 18;
for id = 1:ndayshown
    N = histcounts(ptcat{id}, edges, Normalization='probability');
    plot(ax, centers, N, Color=hsl2rgb([0.7*(id-1)/(ndayshown-1), 1, 0.4]), LineWidth=1.5, DisplayName=sprintf('Day %g', daysPress(id)))
    xlabel(ax, 'Bar contact time (s)')
    ylabel(ax, 'Probability')
end
hold(ax, 'off')

% Legends (we're gonna be in this city)
assert(ndayshown==18)
ticks = [1, 18];
cmap = arrayfun(@(id) hsl2rgb([0.7*(id-1)/(ndayshown-1), 1, 0.4]), 1:ndayshown, 'UniformOutput', false);
cmap = cat(1, cmap{:});
colormap(ax, cmap)
hCb = colorbar(ax);
hCb.Ticks = (1:ndayshown) / ndayshown;
ticklabels = arrayfun(@(x) sprintf('%i', x), 1:ndayshown, UniformOutput=false);
for i = 1:ndayshown
    if ~ismember(i, ticks)
        ticklabels{i} = '';
    end
end
hCb.TickLabels = ticklabels;
hCb.Label.String = 'session';
hCb.Label.Position = [1.103333312471708,0.505319625773328,0];

% title(ax, sprintf('Reach task performance (%g animals)', nAnimalsPress), FontSize=p.fontSize)
fontsize(fig, p.fontSize, 'points')
xticks(ax, [0, 4, 10])
yticks(ax, ax.YLim(2))
ax.YLabel.Position = [-0.466049359260518,0.153243392254861,-1];


% 1d. Raster/PETH of two example units (turn on vs. turn off)
% Load example units
unitNames = { ... 
    'daisy13_20220106_Electrode39_Unit1'; ... % Down
    'daisy9_20211013_Electrode23_Unit1'; ... % Up
    };
% files = cellfun(@(name) sprintf('C:\\SERVER\\Units\\Lite_NonDuplicate\\%s.mat', name), unitNames, UniformOutput=false);
% euEg = EphysUnit.load(files);
euEg = eu(ismember(eu.getName(), unitNames));
for iEu = 1:length(euEg)
    ax = nexttile(layout.middle.top.tl);
    thisRd = euEg(iEu).getRasterData('press', window=[0, 0.5], sort=true);
    EphysUnit.plotRaster(ax, thisRd, xlim=[-4, 0.5], sz=1);
    xline(ax, 0, 'k--')
    xlim(ax, [-4, 0.5])
    switch iEu
        case 1
%             hLgd = legend(ax, {'spike', 'tone'});
            delete(legend(ax))
        case 2
            ylabel(ax, '');
            delete(legend(ax))
    end
%     xlabel(ax, 'Time to touch (s)')
    title(ax, '')
    xlabel(ax, '')
    fontsize(ax, p.fontSize, 'points');
    fontname(ax, 'Arial');
    if iEu == 1
        h = text(ax, 0, 0, 'd', FontSize=16, FontName='Arial', FontWeight='bold');
        ax.Units = 'inches'; h.Units = 'inches';
        h.HorizontalAlignment = 'right';
        h.VerticalAlignment = 'top';
        h.Position = [-0.4, ax.Position(4), 0];
    end
end
% PETH of two example units
for iEu = 1:length(euEg)
    ax = nexttile(layout.middle.bottom.tl);
    thisETA = euEg(iEu).getETA('count', 'press', [-4, 0.5], minTrialDuration=2, normalize='none');
    thisETA.X = thisETA.X./0.100;
    plot(ax, thisETA.t, thisETA.X, LineWidth=p.lineWidth, Color='black')
    xline(ax, 0, 'k--')
    xlim(ax, [-4, 0.5])
    switch iEu
        case 1
%             ylim(ax, [5, 30])
            ylabel(ax, 'Spike rate (sp/s)')
        case 2
%             ylim(ax, [20, 75])
    end
    fontsize(ax, p.fontSize, 'points');
    fontname(ax, 'Arial');
    ax.Box = 'off';
end
xlabel(layout.middle.bottom.tl, 'Time to bar-contact (s)', FontSize=p.fontSize)
delete(legend(ax))

% 1e. PETH of all units (heatmap)
% fig = figure(Units='inches', Position=[0, 0, 4, 5]);
ax = nexttile(layout.bottom.tl, [1, layout.bottom.left.w]);
EphysUnit.plotETA(ax, etaFine.press, c.hasPress, xlim=[-4,0.5], clim=[-1.5, 1.5], sortWindow=[-3, 0], signWindow=[-0.3, 0], sortThreshold=0.25, negativeSortThreshold=0.25, ...
    order=onset.pressOrder(c.hasPress)); 
hold(ax, 'on')
xline(ax, 0, 'k--')
% applyCustomColormap(ax, [-1.5, 1.5], hlim=[0.375, 0, 0, -0.375], llim=[0.125, 0.5, 0.5, 0.25], hpwr=.5, lpwr=1, h0=0.33);
applyCustomColormap(ax, [-1.5, 1.5], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);

yt = 0:100:nnz(c.hasPress);
yt(1) = 1;
if round(yt(end)./100) == round(nnz(c.hasPress)./100)
    yt(end) = nnz(c.hasPress);
else
    yt(end + 1) = nnz(c.hasPress);
end
yt = unique(yt);
yticks(ax, yt)

% hold(ax, 'on')
% plot(lat(order), 1:nnz(c.hasPress))
ax.Colorbar.Label.Position = [-0.995833372448878, 0.033151078619351, 0];
title(ax, '')
xlabel(ax, 'Time to bar-contact (s)')
fontsize(ax, p.fontSize, 'points');
fontname(ax, 'Arial')
% annotation(fig, 'textbox', Units='inches', Position=[0.30,4.5,0.4,0.36], EdgeColor='none', String='c', FontSize=16, FontName='Arial', FontWeight='bold');
h = text(ax, 0, 0, 'e', FontSize=16, FontName='Arial', FontWeight='bold');
ax.Units = 'inches'; h.Units = 'inches';
h.HorizontalAlignment = 'right';
h.VerticalAlignment = 'top';
h.Position = [-0.4, ax.Position(4)+0.2, 0];

% 1f,g,h Baseline spike rates, pre-move response, normalized pre-move response
% fig = figure(Units='inches', Position=[5, 0, 2.5, 5], DefaultAxesFontSize=p.fontSize);
ax = nexttile(layout.bottom.right.tl);
hold(ax, 'on')
edges = 0:5:150;
hHist1 = gobjects(3, 1);
hHist1(1) = histogram(ax, bsr(c.hasPress), edges, FaceColor='white', DisplayName='all');
hHist1(2) = histogram(ax, bsr(c.isPressUp), edges, FaceColor='red', DisplayName='inc', EdgeColor='none');
hHist1(3) = histogram(ax, bsr(c.isPressDown), edges, FaceColor='blue', DisplayName='dec', EdgeColor='none');
hold(ax, 'off')
xlabel(ax, 'Baseline spike rate (sp/s)'), ylabel(ax, 'Count')
fontsize(ax, p.fontSize, 'points');
fontname(ax, 'Arial')
h = text(ax, 0, 0, 'f', FontSize=16, FontName='Arial', FontWeight='bold');
ax.Units = 'inches'; h.Units = 'inches';
h.HorizontalAlignment = 'right';
h.VerticalAlignment = 'top';
h.Position = [-0.5, ax.Position(4)+0.2, 0];
hLgd = legend(ax, hHist1, Orientation='horizontal');
hLgd.Layout.Tile = 'north';

ax = nexttile(layout.bottom.right.tl);
hold(ax, 'on')
hHist2 = histogram(ax, meta.pressRaw(c.hasPress)./0.1 - meta.pressRawBaseline(c.hasPress)./0.1, 30, FaceColor='white');
histogram(ax, meta.pressRaw(c.isPressUp)./0.1 - meta.pressRawBaseline(c.isPressUp)./0.1, hHist2.BinEdges, FaceColor='red', EdgeColor='none')
histogram(ax, meta.pressRaw(c.isPressDown)./0.1 - meta.pressRawBaseline(c.isPressDown)./0.1, hHist2.BinEdges, FaceColor='blue', EdgeColor='none')
hold(ax, 'off')
xlabel(ax, 'Peri-reach response (\Deltasp/s)'), ylabel(ax, 'Count')
fontsize(ax, p.fontSize, 'points');
fontname(ax, 'Arial')
h = text(ax, 0, 0, 'g', FontSize=16, FontName='Arial', FontWeight='bold');
ax.Units = 'inches'; h.Units = 'inches';
h.HorizontalAlignment = 'right';
h.VerticalAlignment = 'top';
h.Position = [-0.5, ax.Position(4)+0.2, 0];

ax = nexttile(layout.bottom.right.tl);
hold(ax, 'on')
hHist3 = histogram(ax, meta.press(c.hasPress), 40, FaceColor='white');
histogram(ax, meta.press(c.isPressUp), hHist3.BinEdges, FaceColor='red', EdgeColor='none')
histogram(ax, meta.press(c.isPressDown), hHist3.BinEdges, FaceColor='blue', EdgeColor='none')
hold(ax, 'off')
xlabel(ax, {'Normalized peri-reach', 'response (a.u.)'}), ylabel(ax, 'Count')
fontsize(ax, p.fontSize, 'points');
fontname(ax, 'Arial')
h = text(ax, 0, 0, 'h', FontSize=16, FontName='Arial', FontWeight='bold');
ax.Units = 'inches'; h.Units = 'inches';
h.HorizontalAlignment = 'right';
h.VerticalAlignment = 'top';
h.Position = [-0.5, ax.Position(4)+0.2, 0];

copygraphics(fig, ContentType='vector', BackgroundColor='none')



%% Fig S1. Movement time histograms for each animal (self-timed reach)
close all
animalNames = cellfun(@(bs) bs(1).animalName, bs, UniformOutput=false);


edges = 1:1:10;
centers = 0.5*(edges(2:end) + edges(1:end-1));

ncols = 5;
nrows = ceil(nAnimalsPress/ncols);
fig = figure(Units='inches', Position=[0, 0, 6.5, 6], DefaultAxesFontSize=p.fontSize, DefaultAxesFontName='Arial', Name='Reach task training progress');
tl = tiledlayout(fig, nrows, ncols, TileSpacing='compact');

hasPress = nSessionsPress > 0;
hasLick = nSessionsLick > 0;

[~, bestDayPress] = sort(-cellfun(@(t) nnz(t >= 3 & t <= 7) / nnz(t >= 1), pt), 2);
[~, bestDayLick] = sort(-cellfun(@(t) nnz(t >= 3 & t <= 7) / nnz(t >= 1), lt), 2);


for ia = find(hasPress(:)')
    ax = nexttile(tl);
    hold(ax, 'on')
    ptsel = pt(ia, bestDayPress(ia, 1:3));
    ptsel = cat(1, ptsel{:});
    ptsel = ptsel(ptsel >= edges(1));
    N = histcounts(ptsel, edges, Normalization='probability');
    plot(ax, centers, N, 'r', LineWidth=2)
    ylim(ax, [0, max(N) + 0.01])
    
%     ndays = nSessionsPress(ia);
%     for id = 1:ndays
%         N = histcounts(pt{ia, id}, edges, Normalization='probability');
%         plot(ax, centers, N, Color=hsl2rgb([0.7*(id-1)/(ndays-1), 0.1, 0.5]), LineWidth=0.1, DisplayName=sprintf('Day %g', daysPress(id)))
%     end
end