%%% Figure 3
p.fontSize = 9;
p.lineWidth = 1.5;
close all


%% Load all units
load_ephysunits;
% boot_response_dir;
load('C:\SERVER\Units\boot_20241021_perimovement_0.3_0.mat')
% load('C:\SERVER\bootMoveResponse_20240830.mat')
%% Load example units
unitNames = { ... 
    'daisy13_20220106_Electrode39_Unit1'; ... % Down
    'daisy9_20211013_Electrode23_Unit1'; ... % Up
    };
files = cellfun(@(name) sprintf('C:\\SERVER\\Units\\Lite_NonDuplicate\\%s.mat', name), unitNames, UniformOutput=false);
euEg = EphysUnit.load(files);

%% 3a,b. Raster/PETH of two example units (turn on vs. turn off)
clear layout
layout.w = 7;
layout.h = 8;
layout.top.h = 3;
layout.bottom.h = 5;
layout.top.top.h = 9;
layout.top.bottom.h = 7;
layout.bottom.left.w = 3;
layout.bottom.right.w = 2;

close all
fig = figure(Units='inches', Position=[0, 0, layout.w, layout.h]);

layout.tl = tiledlayout(fig, layout.top.h + layout.bottom.h, 1, TileSpacing='loose');
layout.top.tl = tiledlayout(layout.tl, layout.top.top.h + layout.top.bottom.h, 1, TileSpacing='compact');
l = layout.top.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [layout.top.h, 1];

layout.top.top.tl = tiledlayout(layout.top.tl, 1, 2, TileSpacing='compact');
l = layout.top.top.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [layout.top.top.h, 1];

layout.top.bottom.tl = tiledlayout(layout.top.tl, 1, 2, TileSpacing='compact');
l = layout.top.bottom.tl; l.Layout.Tile = 1 + layout.top.top.h; l.Layout.TileSpan = [layout.top.bottom.h, 1];

layout.bottom.tl = tiledlayout(layout.tl, 1, layout.bottom.left.w + layout.bottom.right.w, TileSpacing='loose');
l = layout.bottom.tl; l.Layout.Tile = 1 + layout.top.h; l.Layout.TileSpan = [layout.bottom.h, 1];

layout.bottom.right.tl = tiledlayout(layout.bottom.tl, 3, 1, TileSpacing='compact');
l = layout.bottom.right.tl; l.Layout.Tile = 1 + layout.bottom.left.w; l.Layout.TileSpan = [1, layout.bottom.right.w];

for iEu = 1:length(euEg)
    ax = nexttile(layout.top.top.tl);
    thisRd = euEg(iEu).getRasterData('press', window=[0, 0], sort=true);
    EphysUnit.plotRaster(ax, thisRd, xlim=[-4, 0], sz=1);
    switch iEu
        case 1
            legend(ax, {'spike', 'tone'})
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
        h = text(ax, 0, 0, 'a', FontSize=16, FontName='Arial', FontWeight='bold');
        ax.Units = 'inches'; h.Units = 'inches';
        h.HorizontalAlignment = 'right';
        h.VerticalAlignment = 'top';
        h.Position = [-0.5, ax.Position(4), 0];
    end
end

% 3b. PETH of two example units
for iEu = 1:length(euEg)
    ax = nexttile(layout.top.bottom.tl);
    thisETA = euEg(iEu).getETA('count', 'press', [-4, 0], minTrialDuration=2, normalize='none');
    thisETA.X = thisETA.X./0.100;
    plot(ax, thisETA.t, thisETA.X, LineWidth=p.lineWidth, Color='black')
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
xlabel(layout.top.bottom.tl, 'Time to bar-contact (s)', FontSize=p.fontSize)
delete(legend(ax))


% 3d. PETH of all units (heatmap)
% fig = figure(Units='inches', Position=[0, 0, 4, 5]);
ax = nexttile(layout.bottom.tl, [1, layout.bottom.left.w]);
[~, order, ~, lat] = EphysUnit.plotETA(ax, eta.press, c.hasPress, xlim=[-4,0], clim=[-2, 2], sortWindow=[-3, 0], signWindow=[-0.3, 0], sortThreshold=0.25, negativeSortThreshold=0.25); 

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
title(ax, 'Reach PETH')
xlabel(ax, 'Time to bar-contact (s)')
fontsize(ax, p.fontSize, 'points');
fontname(ax, 'Arial')
% annotation(fig, 'textbox', Units='inches', Position=[0.30,4.5,0.4,0.36], EdgeColor='none', String='c', FontSize=16, FontName='Arial', FontWeight='bold');
h = text(ax, 0, 0, 'b', FontSize=16, FontName='Arial', FontWeight='bold');
ax.Units = 'inches'; h.Units = 'inches';
h.HorizontalAlignment = 'right';
h.VerticalAlignment = 'top';
h.Position = [-0.4, ax.Position(4)+0.4, 0];

% 3e,f,g Baseline spike rates, pre-move response, normalized pre-move response
% fig = figure(Units='inches', Position=[5, 0, 2.5, 5], DefaultAxesFontSize=p.fontSize);
ax = nexttile(layout.bottom.right.tl);
hold(ax, 'on')
edges = 0:5:150;
hHist1 = gobjects(3, 1);
hHist1(1) = histogram(ax, msr(c.hasPress), edges, FaceColor='white', DisplayName='all');
hHist1(2) = histogram(ax, msr(c.isPressUp), edges, FaceColor='red', DisplayName='inc', EdgeColor='none');
hHist1(3) = histogram(ax, msr(c.isPressDown), edges, FaceColor='blue', DisplayName='dec', EdgeColor='none');
hold(ax, 'off')
xlabel(ax, 'Baseline spike rate (sp/s)'), ylabel(ax, 'Count')
fontsize(ax, p.fontSize, 'points');
fontname(ax, 'Arial')
h = text(ax, 0, 0, 'c', FontSize=16, FontName='Arial', FontWeight='bold');
ax.Units = 'inches'; h.Units = 'inches';
h.HorizontalAlignment = 'right';
h.VerticalAlignment = 'top';
h.Position = [-0.5, ax.Position(4)+0.4, 0];
hLgd = legend(ax, hHist1, Orientation='horizontal', Location='layout');
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
h = text(ax, 0, 0, 'd', FontSize=16, FontName='Arial', FontWeight='bold');
ax.Units = 'inches'; h.Units = 'inches';
h.HorizontalAlignment = 'right';
h.VerticalAlignment = 'top';
h.Position = [-0.5, ax.Position(4)+0.1, 0];

ax = nexttile(layout.bottom.right.tl);
hold(ax, 'on')
hHist3 = histogram(ax, meta.press(c.hasPress), 40, FaceColor='white');
histogram(ax, meta.press(c.isPressUp), hHist3.BinEdges, FaceColor='red', EdgeColor='none')
histogram(ax, meta.press(c.isPressDown), hHist3.BinEdges, FaceColor='blue', EdgeColor='none')
hold(ax, 'off')
xlabel(ax, 'Peri-reach response (a.u.)'), ylabel(ax, 'Count')
fontsize(ax, p.fontSize, 'points');
fontname(ax, 'Arial')
h = text(ax, 0, 0, 'e', FontSize=16, FontName='Arial', FontWeight='bold');
ax.Units = 'inches'; h.Units = 'inches';
h.HorizontalAlignment = 'right';
h.VerticalAlignment = 'top';
h.Position = [-0.5, ax.Position(4)+0.1, 0];

copygraphics(fig, ContentType='vector')
