%%% Figure 3b-c
p.fontSize = 8;
p.width = 3.5;
p.height = 3;
close all


%% Load example units
unitNames = { ... 
    'desmond13_20190508_Channel2_Unit1'; ... % Down
    'desmond17_20200310_Channel8_Unit1'; ... % Up
    };
files = cellfun(@(name) sprintf('C:\\SERVER\\Units\\Lite_NonDuplicate\\%s.mat', name), unitNames, UniformOutput=false);
eu = EphysUnit.load(files);

%% 3b. Raster of two example units (turn on vs. turn off)
for iEu = 1:length(eu)
    thisRd = eu(iEu).getRasterData('press', window=[0, 0], sort=true);
    ax = EphysUnit.plotRaster(thisRd, xlim=[-4, 0], sz=1);
    ax.Parent.Units = "inches";
    ax.Parent.Position = [0, 0, p.width, p.height];
    fontsize(ax.Parent, p.fontSize, 'points');
    fontsize(ax, p.fontSize, 'points');
    fontname(ax.Parent, 'Arial');
    xlabel(ax, 'Time to bar-contact (s)')
    title(ax, '')
end
delete(legend(ax))

%% 3c. PETH of two example units (turn on vs. turn off)
for iEu = 1:length(eu)
    clear bta
    [bta.X, bta.T, bta.N, bta.S, bta.B] = eu(iEu).getBinnedTrialAverage('count', linspace(1, 6, 6), 'press', ...
        alignTo='stop', window=[-4, 0], resolution=0.1, normalize=false);
    fig = figure(Units='inches', Position=[0, 0, p.width, p.height*0.5]);
    ax = axes(fig);
    bta.X = bta.X ./ 0.1;
    bta.S = bta.S ./ 0.1;
    EphysUnit.plotBinnedTrialAverage(ax, bta, [-4, 0], nsigmas=0, showTrialNum=false)
    xlabel(ax, 'Time to bar-contact (s)')
    ylabel(ax, 'Spike rate (sp/s)')
    fontsize(fig, p.fontSize, 'points');
    fontname(fig, 'Arial');
end
delete(legend(ax))

%% 3d. PETH of all units (heatmap)
fig = figure(Units='inches', Position=[0, 0, 4, 5]);
ax = axes(fig);
EphysUnit.plotETA(ax, eta.press, c.hasPress, xlim=[-4,0], clim=[-2, 2], sortWindow=[-3, 0], signWindow=[-0.2, 0], sortThreshold=0.6, negativeSortThreshold=0.3); 
title(ax, 'Pre-reach PETH')
xlabel(ax, 'Time to bar-contact (s)')
fontname(fig, 'Arial')
fontsize(fig, p.fontSize, 'points');
fontsize(ax, p.fontSize, 'points');

%% 3e,f,g Baseline spike rates, pre-move response, normalized pre-move response
fig = figure(Units='inches', Position=[5, 0, 2.5, 5], DefaultAxesFontSize=p.fontSize);
ax = subplot(3, 1, 1);
histogram(ax, msr(c.hasPress), 0:5:150);
xlabel(ax, 'Baseline spike rate (sp/s)'), ylabel(ax, 'Count')
% legend(ax, sprintf('N=%g', length(msr(c.hasPress))))

ax = subplot(3, 1, 2);
histogram(ax, meta.pressRaw(c.hasPress)./0.1 - msr(c.hasPress), 30)
xlabel(ax, 'Pre-move response (\Deltasp/s)'), ylabel(ax, 'Count')
% legend(ax, sprintf('N=%g', nnz(c.hasPress)))

ax = subplot(3, 1, 3);
histogram(ax, meta.press(c.hasPress), 30)
xlabel(ax, 'Normalized pre-move response (a.u.)'), ylabel(ax, 'Count')
% legend(ax, sprintf('N=%g', nnz(c.hasPress)))
fontsize(fig, p.fontSize, 'points');
fontname(fig, 'Arial')


%% Load all units
clear, clc
load_ephysunits;