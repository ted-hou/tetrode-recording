%%
% 8a. Spontaneous reach task, trial structure
Done
% 8b. Distribution of Inter-touch-intervals
% 8c. Example units (raster vs ETA)

% 8d. Heatmap PETH for all units

%%
read_spontaneous;

%%
p.fontSize = 9;

clear layout
layout.w = 7;
layout.h = 6;
layout.left.w = 4;
layout.right.w = 3;
layout.left.top.h = 3;
layout.left.bottom.h = 7;
layout.right.top.h = 3;
layout.right.bottom.h = 7;

close all
fig = figure(Units='inches', Position=[0, 0, layout.w, layout.h]);
layout.tl = tiledlayout(fig, 1, layout.left.w + layout.right.w, TileSpacing='loose', Padding='loose');
layout.left.tl = tiledlayout(layout.tl, layout.left.top.h + layout.left.bottom.h, 1, TileSpacing='loose', Padding='loose');
l = layout.left.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [1, layout.left.w];

layout.right.tl = tiledlayout(layout.tl, layout.right.top.h + layout.right.bottom.h, 1, TileSpacing='loose', Padding='loose');
l = layout.right.tl; l.Layout.Tile = 1 + layout.left.w; l.Layout.TileSpan = [1, layout.right.w];

layout.left.bottom.tl = tiledlayout(layout.left.tl, 2, 2, TileSpacing='compact', Padding='loose');
l = layout.left.bottom.tl; l.Layout.Tile = 1 + layout.left.top.h; l.Layout.TileSpan = [layout.left.bottom.h, 1];

% 8b. Distribution of Inter-touch-intervals
ax = nexttile(layout.right.tl, 1, [layout.right.top.h, 1]);
interTouchIntervals = cell(length(expSpontaneous), 1);
for iExp = 1:length(expSpontaneous)
    touchTimes = [expSpontaneous(iExp).eu(1).Trials.Press.Stop];
    interTouchIntervals{iExp} = diff(touchTimes);
end
interTouchIntervals = cat(2, interTouchIntervals{:});
% edges = 0:2:max(ceil(interTouchIntervals/2)*2);
edges = 0:1:15;
histogram(ax, interTouchIntervals, edges, Normalization='probability', ...
    EdgeAlpha=1, FaceColor='black')
xlabel(ax, 'Inter-reach interval (s)')
ylabel(ax, 'Probability')
yticks(ax, 0:0.2:0.4)
fontsize(ax, p.fontSize, 'points')
fontname(ax, 'Arial')
% copygraphics(fig)

% Plot single unit raster/ETAs
% for iEu = 1:length(euSpontaneous)
%     fig = figure(Units='inches', Position=[0, 0, p.width-p.heatmapWidth, p.secondRowHeight]);
% 
%     % Raster
%     ax = subplot(2, 1, 1);
%     thisRd = euSpontaneous(iEu).getRasterData('press', window=[-0, 0], sort=true, MinTrialDuration=4);
%     EphysUnit.plotRaster(ax, thisRd, xlim=[-4, 0], sz=2, iti=false);
%     % ax.Parent.Units = "inches";
%     % ax.Parent.Position = [0, 0, 2, 2];
%     % fontsize(ax.Parent, p.fontSize, 'points');
%     xlabel(ax, '')
%     title(ax, '')
% %     legend(ax, {'spike', 'previous rewarded reach'}, Orientation='horizontal', ...
% %         Position=[0.209661369813197,0.943898803036696,0.638020822235073,0.044270832324401], ...
% %         FontSize=p.fontSize, fontName='Arial')
%     legend(ax, 'off')
%     fontsize(ax, p.fontSize, 'points');
%     fontname(ax, 'Arial');
% 
%     % PETH
%     ax = subplot(2, 1, 2);
%     clear btaEg
%     [btaEg.X, btaEg.T, btaEg.N, btaEg.S, btaEg.B] = euSpontaneous(iEu).getBinnedTrialAverage('count', [4, 8, 16, 32, 64], 'press', ...
%         alignTo='stop', window=[-4, 0], resolution=0.1, normalize=false);
%     btaEg.X = btaEg.X ./ 0.1;
%     btaEg.S = btaEg.S ./ 0.1;
%     EphysUnit.plotBinnedTrialAverage(ax, btaEg, [-4, 0], nsigmas=0, showTrialNum=false);
%     xlabel(ax, 'Time to bar contact (s)')
%     ylabel(ax, 'Spike rate (sp/s)')
%     fontsize(ax, p.fontSize, 'points');
%     fontname(ax, 'Arial');
%     ax.Legend.Orientation = 'horizontal';
%     ax.Legend.NumColumns = 2;
%     ax.Legend.Position = [0.230494698528723,0.466257448903559,0.567708324330549,0.069010414959242];
% 
%     print(fig, sprintf('C:\\SERVER\\Figures\\Spontaneous\\%s', euSpontaneous(iEu).getName()), '-dpng')
%     close(fig)
% end
% 
% clear iEu fig ax thisRd btaEg

% 8c. Example unit ETAs and population average ETAs
exampleUnitNames = { ...
        'desmond31_20230804_Channel85_Unit1', ... % Up at -1, 40sp/s
        'daisy18_20230728_Channel109_Unit1', ... % Down at -1.5, wierd for 4-8s trials
    };

for i = 1:length(exampleUnitNames)
    iEu = find(strcmpi(euSpontaneous.getName(), exampleUnitNames{i}));
    ax = nexttile(layout.left.bottom.tl);
    clear btaEg
    [btaEg.X, btaEg.T, btaEg.N, btaEg.S, btaEg.B] = euSpontaneous(iEu).getBinnedTrialAverage('count', [8, 16, 32, 64], 'press', ...
        alignTo='stop', window=[-4, 0], resolution=0.1, normalize=false);
    btaEg.X = btaEg.X ./ 0.1;
    btaEg.S = btaEg.S ./ 0.1;
    EphysUnit.plotBinnedTrialAverage(ax, btaEg, [-4, 0], nsigmas=0, showTrialNum=false, numFormat='%i')
    ax.Legend.Orientation = 'horizontal';
    ax.Legend.Layout.Tile = 'north';
    if i == 2
        delete(legend(ax))
    end
    if i == 1
        ylim(ax, [30, 80])
    else
        ylim(ax, [10, 60])
    end
    set(ax, FontSize=p.fontSize, FontName='Arial')
end

DATA = {btaSpontaneous.pressUpRaw, btaSpontaneous.pressDownRaw};
for i = 1:2
    ax = nexttile(layout.left.bottom.tl);
    EphysUnit.plotBinnedTrialAverage(ax, DATA{i}, [-4, 0], nsigmas=1, sem=true, showTrialNum=false, numFormat='%i')
    if i == 1
        ylim(ax, [30, 60])
    else
        ylim(ax, [20, 50])
    end
    delete(legend(ax))
    set(ax, FontSize=p.fontSize, FontName='Arial')
end

xlabel(layout.left.bottom.tl, 'Time to bar contact (s)', FontSize=p.fontSize, FontName='Arial')
ylabel(layout.left.bottom.tl, 'Spike rate (sp/s)', FontSize=p.fontSize, FontName='Arial')


% 8d. Plot ETA (touch time)
ax = nexttile(layout.right.tl, 1 + layout.right.top.h, [layout.right.bottom.h, 1]);
EphysUnit.plotETA(ax, etaSpontaneous, xlim=[-4,0], clim=[-1.5, 1.5], sortWindow=[-2, 0], signWindow=[-0.5, 0], sortThreshold=0.25, negativeSortThreshold=0.125);
title(ax, '')
xlabel(ax, 'Time to bar contact (s)')
fontsize(ax, p.fontSize, 'points')
fontname(ax, 'Arial')

copygraphics(fig, ContentType='vector')
