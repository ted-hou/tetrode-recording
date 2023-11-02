%%
% 8a. Spontaneous reach task, trial structure
Done
% 8b. Distribution of Inter-touch-intervals
% 8c. Example units (raster vs ETA)

% 8d. Heatmap PETH for all units


%%
p.fontSize = 8;
p.width = 7;
p.height = 6;
p.firstRowHeight = 2;
p.secondRowHeight = 4;
p.heatmapWidth = 3;

%% 1.1. Load acute EU objects (duplicates already removed)
eu = EphysUnit.load('C:\SERVER\Units\acute_spontaneous_reach'); 

% Remove multiunit detected by ISI test.
p.ISIThreshold = 0.0015;
for iEu = 1:length(eu)
    st = eu(iEu).SpikeTimes;
    isi = [NaN, diff(st)];
    st(isi == 0) = [];
    isi = [NaN, diff(st)];
    eu(iEu).SpikeTimes = st;
    ISI{iEu} = isi;
end

for iEu = 1:length(eu)
    prcLowISI(iEu) = nnz(ISI{iEu} < p.ISIThreshold) ./ length(ISI{iEu});
end
histogram(prcLowISI, 0:0.01:1)
cat.isMultiUnit = prcLowISI > 0.05;
cat.isSingleUnit = prcLowISI <= 0.05;
eu = eu(cat.isSingleUnit);
clearvars -except p eu

eu = eu';

% 1.2. Load Video Tracking Data (vtd) and ArduinoConnection (ac), and group into experiments
clearvars -except eu
exp = CompleteExperiment2(eu);

% 1.3 Align video and ephys timestamps
exp.alignTimestamps();

% Calculate ETA
eta = eu.getETA('count', 'press', window=[-4, 0], alignTo='stop', includeInvalid=true, ...
        normalize=[-4, -2]);

%% 8b. Distribution of Inter-touch-intervals
fig = figure(Units='inches', Position=[0, 0, p.heatmapWidth, p.firstRowHeight]);
ax = axes(fig);
interTouchIntervals = cell(length(exp), 1);
for iExp = 1:length(exp)
    touchTimes = [exp(iExp).eu(1).Trials.Press.Stop];
    interTouchIntervals{iExp} = diff(touchTimes);
end
interTouchIntervals = cat(2, interTouchIntervals{:});
edges = 0:2:max(ceil(interTouchIntervals/2)*2);
edges = 0:1:100;
histogram(ax, interTouchIntervals, edges, Normalization='probability', ...
    EdgeAlpha=1, FaceColor='white')
xlabel(ax, 'Inter-reach interval (s)')
ylabel(ax, 'Probability')
fontsize(ax, p.fontSize, 'points')
fontname(ax, 'Arial')
copygraphics(fig)

%% 8c. Example unit raster/ETAs
close all
fig = figure(Units='inches', Position=[0, 0, p.width-p.heatmapWidth, p.secondRowHeight]);
exampleUnitNames = { ...
    'daisy17_20230804_Channel34_Unit1', ... % Down at -1
    'daisy17_20230804_Channel45_Unit1', ... % Up at -1
    % 'daisy18_20230802_Channel66_Unit1', ... % Down at different times
    % 'desmond31_20230804_Channel21_Unit1', ... % Down at -1
    };

for i = 1:length(exampleUnitNames)
    ax = subplot(2, 2, i);
    iEu = find(strcmpi(eu.getName(), exampleUnitNames{i}));
    thisRd = eu(iEu).getRasterData('press', window=[-0, 0], sort=true);
    EphysUnit.plotRaster(ax, thisRd, xlim=[-4, 0], sz=1, iti=false);
    % ax.Parent.Units = "inches";
    % ax.Parent.Position = [0, 0, 2, 2];
    % fontsize(ax.Parent, p.fontSize, 'points');
    xlabel(ax, 'Time to touch (s)')
    title(ax, '')
    legend(ax, {'spike', 'previous rewarded reach'}, Orientation='horizontal', ...
        Position=[0.209661369813197,0.943898803036696,0.638020822235073,0.044270832324401], ...
        FontSize=p.fontSize, fontName='Arial')
    fontsize(ax, p.fontSize, 'points');
    fontname(ax, 'Arial');
    if i == 2
        delete(ax.Legend)
        ylabel(ax, '')
    end
    xlabel(ax, '')
end
for i = 1:length(exampleUnitNames)
    ax = subplot(2, 2, i + 2);
    iEu = find(strcmpi(eu.getName(), exampleUnitNames{i}));
    clear btaEg
    [btaEg.X, btaEg.T, btaEg.N, btaEg.S, btaEg.B] = eu(iEu).getBinnedTrialAverage('count', [0.5, 2, 7, 20, 100], 'press', ...
        alignTo='stop', window=[-4, 0], resolution=0.1, normalize=false);
    btaEg.X = btaEg.X ./ 0.1;
    btaEg.S = btaEg.S ./ 0.1;
    EphysUnit.plotBinnedTrialAverage(ax, btaEg, [-4, 0], nsigmas=0, showTrialNum=false)
    xlabel(ax, 'Time to bar contact (s)', Position=[0.456694820734458,-2.390726558461782,-1])
    ylabel(ax, 'Spike rate (sp/s)')
    fontsize(ax, p.fontSize, 'points');
    fontname(ax, 'Arial');
    ax.Legend.Orientation = 'horizontal';
    ax.Legend.NumColumns = 2;
    ax.Legend.Position = [0.230494698528723,0.466257448903559,0.567708324330549,0.069010414959242];
    if i == 2
        delete(legend(ax))
        xlabel(ax, '')
        ylabel(ax, '')
    end
end
delete(legend(ax))
copygraphics(fig)

%% 8d. Plot ETA (touch time)
fig = figure(Units='inches', Position=[0, 0, p.heatmapWidth, p.secondRowHeight]);
ax = axes(fig);
EphysUnit.plotETA(ax, eta, xlim=[-4,0], clim=[-2, 2], sortWindow=[-2, 0], signWindow=[-0.5, 0], sortThreshold=0.25, negativeSortThreshold=0.125);
fontsize(ax, p.fontSize, 'points')
fontname(ax, 'Arial')
title(ax, '')
xlabel(ax, 'Time to bar contact (s)')
copygraphics(fig)
