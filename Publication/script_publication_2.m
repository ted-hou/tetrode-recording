%%% Figure 3
p.fontSize = 9;
p.width = 7;
p.height = 8;

p.h1 = 3;
p.h2 = 5;

p.h1_1 = 9;
p.h1_2 = 7;

p.lineWidth = 1.5;
close all


%% Load all units
load_ephysunits;
boot_response_dir;

%% Load example units
unitNames = { ... 
    'daisy13_20220106_Electrode39_Unit1'; ... % Down
    'daisy9_20211013_Electrode23_Unit1'; ... % Up
    };
files = cellfun(@(name) sprintf('C:\\SERVER\\Units\\Lite_NonDuplicate\\%s.mat', name), unitNames, UniformOutput=false);
euEg = EphysUnit.load(files);

%% 3a,b. Raster/PETH of two example units (turn on vs. turn off)
close all
fig = figure(Units='inches', Position=[0, 0, p.width, p.height]);
tlp = tiledlayout(fig, p.h1+p.h2, 1, TileSpacing='loose');
tl = tiledlayout(tlp, p.h1_1 + p.h1_2, 2, TileSpacing='compact');
tl.Layout.Tile = 1;
tl.Layout.TileSpan = [p.h1, 1];
for iEu = 1:length(euEg)
    ax = nexttile(tl, [p.h1_1, 1]);
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
    ax = nexttile(tl, [p.h1_2, 1]);
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
    if iEu == 1
        % h = text(ax, 0, 0, 'b', FontSize=16, FontName='Arial', FontWeight='bold');
        ax.Units = 'inches'; h.Units = 'inches';
        h.HorizontalAlignment = 'right';
        h.VerticalAlignment = 'top';
        h.Position = [-0.5, ax.Position(4), 0];
    end
end
xlabel(tl, 'Time to bar-contact (s)')
delete(legend(ax))


% 3d. PETH of all units (heatmap)
% fig = figure(Units='inches', Position=[0, 0, 4, 5]);
tl = tiledlayout(tlp, 3, 5, TileSpacing='compact');
tl.Layout.Tile = p.h1 + 1;
tl.Layout.TileSpan = [p.h2, 1];
ax = nexttile(tl, [3, 3]);
EphysUnit.plotETA(ax, eta.press, c.hasPress, xlim=[-4,0], clim=[-2, 2], sortWindow=[-3, 0], signWindow=[-0.2, 0], sortThreshold=0.6, negativeSortThreshold=0.3); 
ax.Colorbar.Label.Position = [-0.995833372448878, 0.033151078619351, 0];
title(ax, 'Pre-reach PETH')
xlabel(ax, 'Time to bar-contact (s)')
fontsize(ax, p.fontSize, 'points');
fontname(ax, 'Arial')
% annotation(fig, 'textbox', Units='inches', Position=[0.30,4.5,0.4,0.36], EdgeColor='none', String='c', FontSize=16, FontName='Arial', FontWeight='bold');
h = text(ax, 0, 0, 'b', FontSize=16, FontName='Arial', FontWeight='bold');
ax.Units = 'inches'; h.Units = 'inches';
h.HorizontalAlignment = 'right';
h.VerticalAlignment = 'top';
h.Position = [-0.5, ax.Position(4), 0];

% 3e,f,g Baseline spike rates, pre-move response, normalized pre-move response
% fig = figure(Units='inches', Position=[5, 0, 2.5, 5], DefaultAxesFontSize=p.fontSize);
ax = nexttile(tl, [1, 2]);
histogram(ax, msr(c.hasPress), 0:5:150, FaceColor='white');
xlabel(ax, 'Baseline spike rate (sp/s)'), ylabel(ax, 'Count')
fontsize(ax, p.fontSize, 'points');
fontname(ax, 'Arial')
h = text(ax, 0, 0, 'c', FontSize=16, FontName='Arial', FontWeight='bold');
ax.Units = 'inches'; h.Units = 'inches';
h.HorizontalAlignment = 'right';
h.VerticalAlignment = 'top';
h.Position = [-0.4, ax.Position(4)+0.1, 0];

ax = nexttile(tl, [1, 2]);
histogram(ax, meta.pressRaw(c.hasPress)./0.1 - msr(c.hasPress), 30, FaceColor='white')
xlabel(ax, 'Pre-move response (\Deltasp/s)'), ylabel(ax, 'Count')
fontsize(ax, p.fontSize, 'points');
fontname(ax, 'Arial')
h = text(ax, 0, 0, 'd', FontSize=16, FontName='Arial', FontWeight='bold');
ax.Units = 'inches'; h.Units = 'inches';
h.HorizontalAlignment = 'right';
h.VerticalAlignment = 'top';
h.Position = [-0.4, ax.Position(4)+0.1, 0];

ax = nexttile(tl, [1, 2]);
histogram(ax, meta.press(c.hasPress), 30, FaceColor='white')
xlabel(ax, 'Normalized pre-move response (a.u.)'), ylabel(ax, 'Count')
fontsize(ax, p.fontSize, 'points');
fontname(ax, 'Arial')
h = text(ax, 0, 0, 'e', FontSize=16, FontName='Arial', FontWeight='bold');
ax.Units = 'inches'; h.Units = 'inches';
h.HorizontalAlignment = 'right';
h.VerticalAlignment = 'top';
h.Position = [-0.4, ax.Position(4)+0.1, 0];
