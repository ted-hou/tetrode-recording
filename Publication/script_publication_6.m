
%%
if exist('E:\Data\Units\PressVsLick_ArtifactsRemoved_Full\FixedEventsAndTrials', 'dir')
    euArtiFree = EphysUnit.load('E:\Data\Units\PressVsLick_ArtifactsRemoved_Full\FixedEventsAndTrials');
    load('E:\Data\Units\meta_PressVsLick_ArtifactsRemoved_Full_20260107.mat');
    etaArtiFree = metaArtiFree.eta;
else
    euArtiFree = EphysUnit.load('C:\SERVER\Units\PressVsLick_ArtifactsRemoved_Full\FixedEventsAndTrials');
    load('C:\SERVER\Units\meta_PressVsLick_ArtifactsRemoved_Full_20260107.mat');
    etaArtiFree = metaArtiFree.eta;
end

% %% Calculate ETA for correct reach vs. incorrect reach; correct vs. incorrect retract; correct vs. incorrect release; correct vs. incorrect lick
% % clear eta
% 
% clear artifactParams;
% artifactParams(1) = struct(event='LickOn', length=10, lengthUnit='ms', direction='both');
% artifactParams(2) = struct(event='LickOff', length=10, lengthUnit='ms', direction='both');
% artifactParams(3) = struct(event='PressOn', length=10, lengthUnit='ms', direction='both');
% artifactParams(4) = struct(event='PressOff', length=10, lengthUnit='ms', direction='both');
% 
% etaArtiFree.artifactParams = artifactParams;
% etaArtiFree.baselineWindow = struct(press=[-4, -2], lick=[-4, -2], release=[-2, 0]);
% 
% etaArtiFree.resolution = 0.025;
% 
% etaArtiFree.correctPress = euArtiFree.getETA('count', 'PressCorrect', [-4, 4], alignTo='stop', ...
%     normalize='none', resolution=etaArtiFree.resolution, artifacts=artifactParams);
% etaArtiFree.incorrectPress = euArtiFree.getETA('count', 'PressIncorrect', [-4, 4], alignTo='stop', ...
%     normalize='none', resolution=etaArtiFree.resolution, artifacts=artifactParams);
% 
% 
% etaArtiFree.correctRelease = euArtiFree.getETA('count', 'CueToLeverReleaseCorrect', [-4, 4], alignTo='stop', ...
%     normalize='none', resolution=etaArtiFree.resolution, artifacts=artifactParams);
% 
% etaArtiFree.correctLick = euArtiFree.getETA('count', 'lick', [-4, 4], minTrialDuration=4, ...
%     normalize='none', resolution=etaArtiFree.resolution, artifacts=artifactParams);
% etaArtiFree.incorrectLick = euArtiFree.getETA('count', 'lick', [-4, 4], minTrialDuration=2, maxTrialDuration=4, ...
%     normalize='none', resolution=etaArtiFree.resolution, artifacts=artifactParams);
% 
% etaArtiFree.correctLickLastLickOff = euArtiFree.getETA('count', 'CueToLastLickOffCorrect', [-4, 4], alignTo='stop', ...
%     normalize='none', resolution=etaArtiFree.resolution, artifacts=artifactParams);
% 
% etaArtiFree.correctPressFirstLick = euArtiFree.getETA('count', 'CorrectPressToFirstRewardLick', [-4, 4], alignTo='stop', normalize='none', resolution=etaArtiFree.resolution, artifacts=artifactParams);
% etaArtiFree.correctPressLastLickOff = euArtiFree.getETA('count', 'CorrectPressToLastLickOff', [-4, 4], alignTo='stop', normalize='none', resolution=etaArtiFree.resolution, artifacts=artifactParams);
% 
% % Convert to sp/s
% for field = ["correctPress", "incorrectPress", "correctLick", "incorrectLick", "correctRelease", ...
%         "correctLickLastLickOff", "correctPressFirstLick", "correctPressLastLickOff"]
%     etaArtiFree.(field).X = etaArtiFree.(field).X ./ etaArtiFree.resolution;
% end
% 
% clear field
% etaArtiFree.pressNorm = euArtiFree.getETA('count', 'press', [-4, 4], alignTo='stop', minTrialDuration=2, normalize=etaArtiFree.baselineWindow.press, resolution=etaArtiFree.resolution, artifacts=artifactParams);
% etaArtiFree.lickNorm = euArtiFree.getETA('count', 'lick', [-4, 4], alignTo='stop', minTrialDuration=2, normalize=etaArtiFree.pressNorm.stats, resolution=etaArtiFree.resolution, artifacts=artifactParams);
% etaArtiFree.correctReleaseNorm = euArtiFree.getETA('count', 'CueToLeverReleaseCorrect', [-4, 4], alignTo='stop', normalize=etaArtiFree.pressNorm.stats, resolution=etaArtiFree.resolution, artifacts=artifactParams);
% etaArtiFree.correctReleaseBeforeRetractNorm = euArtiFree.getETA('count', 'CueToLeverReleaseBeforeRetractCorrect', [-4, 4], alignTo='stop', normalize=etaArtiFree.pressNorm.stats, resolution=etaArtiFree.resolution, artifacts=artifactParams);
% 
% etaArtiFree.lickNormToSelf = euArtiFree.getETA('count', 'lick', [-4, 4], alignTo='stop', minTrialDuration=2, normalize=etaArtiFree.baselineWindow.lick, resolution=etaArtiFree.resolution, artifacts=artifactParams);
% etaArtiFree.correctReleaseNormToSelf = euArtiFree.getETA('count', 'CueToLeverReleaseCorrect', [-4, 4], alignTo='stop', normalize=etaArtiFree.baselineWindow.release, resolution=etaArtiFree.resolution, artifacts=artifactParams);
% 
% etaArtiFree.correctPressNorm = euArtiFree.getETA('count', 'press', [-4, 4], alignTo='stop', minTrialDuration=4, normalize=etaArtiFree.pressNorm.stats, resolution=etaArtiFree.resolution, artifacts=artifactParams);
% etaArtiFree.correctLickNorm = euArtiFree.getETA('count', 'lick', [-4, 4], alignTo='stop', minTrialDuration=4, normalize=etaArtiFree.pressNorm.stats, resolution=etaArtiFree.resolution, artifacts=artifactParams);
% 
% etaArtiFree.incorrectPressNorm = euArtiFree.getETA('count', 'press', [-4, 4], alignTo='stop', minTrialDuration=2, maxTrialDuration=4, normalize=etaArtiFree.pressNorm.stats, resolution=etaArtiFree.resolution, artifacts=artifactParams);
% etaArtiFree.incorrectLickNorm = euArtiFree.getETA('count', 'lick', [-4, 4], alignTo='stop', minTrialDuration=2, maxTrialDuration=4, normalize=etaArtiFree.pressNorm.stats, resolution=etaArtiFree.resolution, artifacts=artifactParams);
% 
% etaArtiFree.correctPressFirstLickNorm = euArtiFree.getETA('count', 'CorrectPressToFirstRewardLick', [-4, 4], alignTo='stop', normalize=etaArtiFree.pressNorm.stats, resolution=etaArtiFree.resolution, artifacts=artifactParams);
% etaArtiFree.correctPressLastLickOffNorm = euArtiFree.getETA('count', 'CorrectPressToLastLickOff', [-4, 4], alignTo='stop', normalize=etaArtiFree.pressNorm.stats, resolution=etaArtiFree.resolution, artifacts=artifactParams);
% 
% etaArtiFree.correctLickLastLickOffNorm = euArtiFree.getETA('count', 'CueToLastLickOffCorrect', [-4, 4], alignTo='stop', ...
%     normalize=etaArtiFree.pressNorm.stats, resolution=etaArtiFree.resolution, artifacts=artifactParams);
% 
% % Lick trials, last lickOff before next cue
% 
% 
% % Lick bouts (norm to pre-press [-4, -2])
% etaArtiFree.lickBoutNaive = euArtiFree.getETA('count', 'lickbout_naive', window=[0, 2*pi*4], resolution=2*pi/8, normalize='none',  minInterval=0.05, maxInterval=0.20, ...
%     minBoutCycles=2, maxBoutCycles=4, artifacts=artifactParams);
% etaArtiFree.lickBoutNaiveNorm = etaArtiFree.lickBoutNaive;
% etaArtiFree.lickBoutNaiveNorm.X = (etaArtiFree.lickBoutNaiveNorm.X - vertcat(etaArtiFree.pressNorm.stats.mean)/etaArtiFree.resolution) ./ (vertcat(etaArtiFree.pressNorm.stats.sd)/etaArtiFree.resolution);
% %%
% metaArtiFree.eta = etaArtiFree;
% save('E:\Data\Units\meta_PressVsLick_ArtifactsRemoved_Full_20260107.mat', 'metaArtiFree')
%%
clearvars -except etaArtiFree euArtiFree metaArtiFree
if ~exist('p', 'var') || ~isfield(p, 'fontSize')
    p.fontSize = 9;
end

%% Calculate peri-lick lick frequency histograms for any first lick in trial
[~, expEuIndices] = unique({euArtiFree.ExpName});
FsLick = 500;
lickHistEdges = 0.01:1/FsLick:2;
lickHistCenters = 0.5*(lickHistEdges(2:end) + lickHistEdges(1:end-1));

lickHistCounts = zeros(size(lickHistCenters));
lickHistNLicks = 0;
for iExp = 1:length(expEuIndices)
    iEu = expEuIndices(iExp);
    firstLickTimes = [euArtiFree(iEu).makeTrials('firstlick').Stop];
    allLickTimes = euArtiFree(iEu).EventTimes.Lick;
    for iLick = 1:length(firstLickTimes)
        edgesGlobal = firstLickTimes(iLick) + lickHistEdges;
        n = histcounts(allLickTimes, edgesGlobal);
        lickHistCounts = lickHistCounts + n;
        lickHistNLicks = lickHistNLicks + 1;
    end
end
lickHist.firstLick.count = lickHistCounts;
lickHist.firstLick.pdf = lickHistCounts ./ lickHistNLicks;
lickHist.firstLick.t = lickHistCenters;
lickHist.firstLick.edges = lickHistEdges;

clear expEuIndices FsLick lickHistEdges lickHistCenters lickHistNLicks lickHistCounts iExp iEu firstLickTimes allLickTimes iLick edgesGlobal n

% Calculate peri-correct-lick lick frequency histograms for osci cells
sel = metaArtiFree.cc.isLick;
[~, expEuIndices] = unique({euArtiFree(sel).ExpName});
FsLick = 500;
lickHistEdges = 0.01:1/FsLick:2;
lickHistCenters = 0.5*(lickHistEdges(2:end) + lickHistEdges(1:end-1));

lickHistCounts = zeros(size(lickHistCenters));
lickHistNLicks = 0;
for iExp = 1:length(expEuIndices)
    iEu = expEuIndices(iExp);
    trials = euArtiFree(iEu).getTrials('lick');
    trials = trials(trials.duration() >= 4);
    firstLickTimes = [trials.Stop];
    allLickTimes = euArtiFree(iEu).EventTimes.Lick;
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

clear sel trials expEuIndices FsLick lickHistEdges lickHistCenters lickHistNLicks lickHistCounts iExp iEu firstLickTimes allLickTimes iLick edgesGlobal n
%%
% etaArtiFree.circLickNaive = euArtiFree.getETA('count', 'circlick_naive', window=[0, 2*pi], resolution=2*pi/30, normalize='none',  minInterval=0.05, maxInterval=0.20, artifacts=metaArtiFree.artifactParams);
% etaArtiFree.lickBoutNaive = euArtiFree.getETA('count', 'lickbout_naive', window=[0, 2*pi*4], resolution=2*pi/30, normalize='none',  minInterval=0.05, maxInterval=0.20, ...
%     minBoutCycles=2, maxBoutCycles=4, artifacts=metaArtiFree.artifactParams);
% 
% etaArtiFree.correctLickBout = euArtiFree.getETA('count', 'lick+lickbout', window=[-4, 0], resolution=[0.025, 2*pi/5], normalize='none', minTrialDuration=4, ...
%     maxInterval=0.2, minInterval=0.05, minBoutCycles=2, maxBoutCycles=6);
% etaArtiFree.correctLickBoutNorm = euArtiFree.getETA('count', 'lick+lickbout', window=[-4, 0], resolution=[0.025, 2*pi/5], normalize=[-4, -2], minTrialDuration=4, ...
%     maxInterval=0.2, minInterval=0.05, minBoutCycles=2, maxBoutCycles=6);
% 
% etaArtiFree.correctPressBout = euArtiFree.getETA('count', 'press+lickbout', window=[-4, 0], resolution=[0.025, 2*pi/32, 2*pi/5], normalize='none', minTrialDuration=4, ...
%     maxInterval=0.2, minInterval=0.05, minBoutCycles=2, maxBoutCycles=6);
% etaArtiFree.correctPressBoutNorm = euArtiFree.getETA('count', 'press+lickbout', window=[-4, 0], resolution=[0.025, 2*pi/32, 2*pi/5], normalize=[-4, -2], minTrialDuration=4, ...
%     maxInterval=0.2, minInterval=0.05, minBoutCycles=2, maxBoutCycles=6);
% 
% 
% metaArtiFree.eta = etaArtiFree;
% save('E:\Data\Units\meta_PressVsLick_ArtifactsRemoved_Full_20260107.mat', 'metaArtiFree')

%% Fig 6
XLIM = {[-1, 0.3], [-1, 0.3], [-0.3, 0.3], [-0.3, 0.3], [-0.3, 0.3]};
W = cellfun(@(xl) diff(xl*10), XLIM, UniformOutput=true);
CW = cumsum([0, W]);
SORTWINDOW = {[-0.3, 0.3], [-0.3, 0.3], [-0.1, 0.3], [-0.1, 0.3], [-0.1, 0.3]};
DATAWINDOW = {[-2, 0.5], [-2, 0.5], [-0.5, 0.5], [-0.5, 0.5], [-0.5, 0.5]};
NAME = ["reach", "lick", "lick\nstart", "lick\nend", "retract"];
TRIALTYPE = {'press', 'lick', 'CorrectPressToFirstRewardLick', 'CueToLastLickOffCorrect', 'CueToLeverReleaseCorrect'};
% XTICKS = {[-1, 0], [-1, 0], [0, 0.3], [0, 0.3], [0, 0.3]};
% XTICKLABELS = {["-1", "touch"], ["-1", "lick"], ["lick", "0.3"], ["lick", "0.3"], ["release", "0.3"]};
XTICKS = {[-1, 0], [-1, 0], [0], [0], [0]};
XTICKLABELS = {["-1", "0"], ["-1", "0"], ["0"], ["0"], ["0"]};
nEgUnits = 2;
close all

% Figure layout
fig = figure(Units='inches', Position=[1, 1, 7, 8]);

clear layout
layout.w = [5, 5];
layout.h = [4, 8, 7];
layout.ch = cumsum([0, layout.h]).*sum(layout.w);
layout.tl = tiledlayout(fig, sum(layout.h), sum(layout.w), TileSpacing='tight', Padding='loose');

% 1st row (left) (examples)
layout.child(1).h = nEgUnits;
layout.child(1).w = W;
layout.child(1).cw = cumsum([0, W]);
layout.child(1).tl = tiledlayout(layout.tl, sum(layout.child(1).h), sum(layout.child(1).w), TileSpacing='tight', Padding='tight');
l = layout.child(1).tl; l.Layout.Tile = 1 + layout.ch(1); l.Layout.TileSpan = [layout.h(1), layout.w(1)];

% 2nd row (left) (heatmap)
layout.child(2).h = 1;
layout.child(2).w = W;
layout.child(2).cw = cumsum([0, W]);
layout.child(2).tl = tiledlayout(layout.tl, sum(layout.child(2).h), sum(layout.child(2).w), TileSpacing='tight', Padding='tight');
l = layout.child(2).tl; l.Layout.Tile = 1 + layout.ch(2); l.Layout.TileSpan = [layout.h(2), layout.w(1)];

% 3rd row (osci)
layout.child(3).h = [9, 20];
layout.child(3).w = [3, 3, 3];
layout.child(3).cw = cumsum([0, layout.child(3).w]);
layout.child(3).ch = cumsum([0, layout.child(3).h]);
layout.child(3).tl = tiledlayout(layout.tl, sum(layout.child(3).h), sum(layout.child(3).w), TileSpacing='loose', Padding='compact');
l = layout.child(3).tl; l.Layout.Tile = 1 + layout.ch(3); l.Layout.TileSpan = [layout.h(3), sum(layout.w)];


% 6a. Raster/PETH examples, 5 phases of movement
unitNames = { ...
    'Daisy3_20180611_Channel15_Unit1' % euArtiFree(28).getName(); ...
    'desmond10_20180909_Channel10_Unit1' % euArtiFree(153).getName(); ...
    };
[~, locb] = ismember(unitNames, euArtiFree.getName());
euEg = euArtiFree(locb);
assert(nEgUnits == length(unitNames))

MINTRIALDURATION = [...
        2, 2, 0, 0, 0; ...
        2, 2, 0, 0, 0; ...
    ];
MAXTRIALDURATION = [...
        Inf, Inf, Inf, Inf, Inf; ...
        Inf, Inf, Inf, Inf, Inf; ...
    ];
EVERYNTH = [5, 5, 5, 5, 5];
YLIM = {[40, 160]; [10, 130]};
YTICKS = {[40, 100, 160]; [10, 70, 130]};

AX = gobjects(nEgUnits, 5);
for iEu = 1:nEgUnits
    for iAx = 1:5
        AX(iEu, iAx) = nexttile(layout.child(1).tl, 1 + CW(iAx) + (iEu-1)*sum(W), [1, W(iAx)]);
    end
end
for iEu = 1:nEgUnits
    for iAx = 1:5
        ax = AX(iEu, iAx);
        thisRD = euEg(iEu).getRasterData(TRIALTYPE{iAx}, window=DATAWINDOW{iAx}, sort=true, alignTo='stop', minTrialDuration=MINTRIALDURATION(iEu, iAx), maxTrialDuration=MAXTRIALDURATION(iEu, iAx));
        thisETA = euEg(iEu).getETA('count', TRIALTYPE{iAx}, DATAWINDOW{iAx}, normalize='none', includeInvalid=false, alignTo='stop', minTrialDuration=MINTRIALDURATION(iEu, iAx), maxTrialDuration=MAXTRIALDURATION(iEu, iAx));
        yyaxis(ax, 'right')
        EphysUnit.plotRaster(ax, thisRD, xlim=XLIM{iAx}, iti=false, sz=1, maxTrials=40, maxTrialsMethod='uniformsample', ...
            everyNth=EVERYNTH(iEu), timingCriterion=NaN, onlyPlotSpikes=true);
        hRaster = ax.Children(1);
        hRaster.MarkerFaceAlpha = 0.5;
        ylabel(ax, '')
        yticks(ax, [])
        ax.YAxis(2).Direction = 'reverse';
        yyaxis(ax, 'left')
        plot(ax, thisETA.t, thisETA.X./0.1, LineWidth=1.5, Color=[0.2, 0.2, 0.8, 1.0])
        hold(ax, 'on')
        % set(ax.YAxis, FontSize=p.fontSize, Color=[0.15, 0.15, 0.15]);
        set(ax.YAxis(1), FontSize=p.fontSize, Color=[0.2, 0.2, 0.8]);%0.15, 0.15, 0.15]);
        set(ax.YAxis(2), FontSize=p.fontSize, Color=[0.15, 0.15, 0.15]);%0.15, 0.15, 0.15]);
        ax.YAxis(1).TickLength = [0.025, 0.1];
        ylabel(ax, 'Spike rate (sp/s)')
        delete(ax.Legend)
        xticks(ax, XTICKS{iAx})
        xticklabels(ax, XTICKLABELS{iAx})
        xtickangle(ax, 0)
        if iEu == 1
            title(ax, strsplit(NAME(iAx), '\\n'))
            xticks(ax, [])
        else
            title(ax, '')
        end
        xline(ax, 0, 'k--', LineWidth=1)
        ylim(ax, YLIM{iEu})
        yticks(ax, YTICKS{iEu})
        xlabel(ax, '')
        fontsize(ax, p.fontSize, 'points')
    end
end
% ylim(AX, [-10, 160]);
% yticks(AX(:, 1), [0, 75, 150])
yticklabels(AX(:, 2:end), [])

ax = AX(1, 1);
hLetter = text(ax, 0, 0, 'a', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.3, ax.Position(4) + 0.2, 0];

% xlabel(AX, '')
ylabel(AX, '')
ylabel(layout.child(1).tl, 'Spike rate (sp/s)', FontSize=p.fontSize, Color=[0.2, 0.2, 0.8, 1.0])
xlabel(layout.child(1).tl, 'Time (s)', FontSize=p.fontSize)
clear thisRD ax iEu AX

% 6b. Heatmap, 5 phases of movement
selUnits = 1:length(euArtiFree);
ETASORT = {etaArtiFree.pressNorm, etaArtiFree.lickNorm, etaArtiFree.correctPressFirstLickNorm, etaArtiFree.correctLickLastLickOffNorm, etaArtiFree.correctReleaseNorm};
ETA = {etaArtiFree.pressNorm, etaArtiFree.lickNorm, etaArtiFree.correctPressFirstLickNorm, etaArtiFree.correctLickLastLickOffNorm, etaArtiFree.correctReleaseNorm};

% Combine ETA, PCA, and sort along 1st dimension
etaCombined = struct(X=[], t=[]);
etaCombined.X = cellfun(@(eta) eta.X, ETASORT, UniformOutput=false);
etaCombined.X = cat(2, etaCombined.X{:});
etaCombined.t = cellfun(@(eta) eta.t, ETASORT, UniformOutput=false);
etaCombined.t = cat(2, etaCombined.t{:});
etaCombined.epoch = arrayfun(@(i) i*ones(1, length(ETASORT{i}.t)), 1:length(ETASORT), UniformOutput=false);
etaCombined.epoch = cat(2, etaCombined.epoch{:});
etaCombined.X(etaCombined.X>1.5) = 1.5;
etaCombined.X(etaCombined.X<-1.5) = -1.5;

etaCombined.X = etaCombined.X(selUnits, :);

% Make templates to project onto
clear template
template(length(ETASORT)) = struct(t=[], x=[]);
for iETA = 1:length(ETASORT)
    template(iETA).t = etaCombined.t;
    template(iETA).x = zeros(1, length(etaCombined.t));
    template(iETA).x(1, isin(etaCombined.t, SORTWINDOW{iETA}) & etaCombined.epoch==iETA) = 1;
end

score = zeros(size(etaCombined.X, 1), length(ETASORT));
etaCombined.X(isnan(etaCombined.X)) = 0;
for iETA = 1:length(ETASORT)
    score(:, iETA) = etaCombined.X * template(iETA).x';
end
groupVar = arrayfun(@(i) bitshift(int16(score(:, i)>0), length(ETASORT)-i), 1:size(score, 2), UniformOutput=false);
groupVar = sum(horzcat(groupVar{:}), 2);

% First, sort by number of negative modulations
numNeg = sum(score<0, 2);
numNeg(numNeg == 5) = -1;
numNeg(numNeg > 1) = 2;
[uniqueGroupVars, ia] = unique(groupVar);
[~, I] = sort(numNeg(ia), 'ascend');
groupVar = changem(groupVar, 0:length(uniqueGroupVars)-1, uniqueGroupVars(I));

% Tighten up the groupvars
uniqueGroupVars = unique(groupVar);
groupVar = changem(groupVar, 0:length(uniqueGroupVars)-1, uniqueGroupVars);
groupSize = histcounts(groupVar, 0:length(uniqueGroupVars));
groupSizeCum = cumsum(groupSize);
sortVal = groupVar*10 - score(:, 1)./max(abs(score(:, 1)));
% selMultiNeg = groupVar >= 7;
% sortVal(selMultiNeg) = max(uniqueGroupVars)*10 + (score(selMultiNeg, :)./max(abs(score(:)))*[-2; -1; 0; 1; 2]);
[~, sortOrder] = sort(sortVal, 'ascend');

ax = gobjects(1, length(XLIM));
for i = 1:length(XLIM)
    ax(i) = nexttile(layout.child(2).tl, CW(i)+1, [1, W(i)]);
end
for iAx = 1:length(ETA)
    hidecb = iAx < length(ETA);
    EphysUnit.plotETA(ax(iAx), ETA{iAx}, selUnits, xlim=XLIM{iAx}, clim=[-1.5, 1.5], order=sortOrder, hidecolorbar=hidecb);
    applyCustomColormap(ax(iAx), [-1.5, 1.5], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);    
    if ~hidecb
        ax(iAx).Colorbar.Layout.Tile = 'east';
        % ax(iAx).Colorbar.Location = 'eastoutside';
        ax(iAx).Colorbar.Label.String = 'Normalized spike rate (a.u.)';
        axc2 = ax(iAx);
    end
    if iAx > 1
        yticks(ax(iAx), [])
    else
        yticks(ax(iAx), groupSizeCum([1, 2, 6, end])+0.5)
        yticklabels(ax(iAx), string(groupSizeCum([1, 2, 6, end])))
    end
    title(ax(iAx), strsplit(NAME{iAx}, "\\n"))
    xlabel(ax(iAx), "")
    ylabel(ax(iAx), "")
    xticks(ax(iAx), XTICKS{iAx})
    xticklabels(ax(iAx), XTICKLABELS{iAx})
    xtickangle(ax(iAx), 0)
    xline(ax(iAx), 0, 'k-')
    yline(ax(iAx), groupSizeCum([2, 7])+0.5, 'k--', LineWidth=2) % Thick lines
    yline(ax(iAx), groupSizeCum([1, 3,4,5,6])+0.5, 'k--', LineWidth=0.5) % Small lines
    fontsize(ax, p.fontSize, 'points')
end
xlabel(layout.child(2).tl, "Time (s)", FontSize=p.fontSize)
ylabel(layout.child(2).tl, "Unit", FontSize=p.fontSize)
% clear etaCombined nDims coeff score explained sortOrder fig tl ax iAx ETASORT NAME ZEROLABEL XLIM hidecp

ax = ax(1);
hLetter = text(ax, 0, 0, 'b', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.25, ax.Position(4) + 0.25, 0];

% 6c. Correct lick bout lick histogram
nBoutsDisp = 6;
ax = nexttile(layout.child(3).tl, [layout.child(3).h(1), layout.child(3).w(1)]);
counts = lickHist.correctLickOsci.count;
smoothedCounts = smoothdata(counts, 'gaussian', 25);
[pks, locs] = findpeaks(smoothedCounts, lickHist.correctLickOsci.t, MinPeakProminence=0.5);
histogram(ax, BinEdges=lickHist.correctLickOsci.edges, BinCounts=counts, Normalization='probability', EdgeColor='none', FaceColor='black', FaceAlpha=0.9)
hold(ax, 'on')
xline(ax, locs, LineStyle=':', Color=[0, 0.5, 0.5, 1], LineWidth=2.5)
ylabel(ax, ' lick\newlineprob')
xlabel(ax, 'Time from first lick (s)')
xticks(ax, -1:0.5:2)
xlim(ax, [0, 1/8*nBoutsDisp])
yticks(ax, [])
ax.Box = 'off';
fontsize(ax, p.fontSize, 'points')

hLetter = text(ax, 0, 0, 'c', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.35, ax.Position(4) + 0.1, 0];


% 6d ETA Heatmap osci lick
tlp = tiledlayout(layout.child(3).tl, 1, 6+1, TileSpacing='tight');
tlp.Layout.Tile = 1 + layout.child(3).ch(2)*sum(layout.child(3).w); tlp.Layout.TileSpan = [layout.child(3).h(2), layout.child(3).w(1)];

ax = nexttile(tlp, 1, [1, 6]);

% etaArtiFree.circLickNaive.Z = etaArtiFree.circLickNaive.X.*exp(etaArtiFree.circLickNaive.t*1i);
% etaArtiFree.circLickNaive.Z(:, 1) = mean(etaArtiFree.circLickNaive.X(:, [2, 30]), 2).*exp(etaArtiFree.circLickNaive.t(1)*1i);
% meanZ = mean(etaArtiFree.circLickNaive.Z, 2);

etaArtiFree.lickBoutNaiveNorm = etaArtiFree.lickBoutNaive;
etaArtiFree.lickBoutNaiveNorm.X = normalize(etaArtiFree.lickBoutNaive.X, 2, 'zscore', 'robust');

% Special plot: let's use pre-reach baseline for z-scoring osci-lick spike
% rates
% etaArtiFree.lickBoutNaiveNorm.X(metaArtiFree.cc.hasPress, :) = (etaArtiFree.lickBoutNaive.X(metaArtiFree.cc.hasPress, :) - vertcat(etaArtiFree.press.stats(metaArtiFree.cc.hasPress).mean)./0.1) ./ (vertcat(etaArtiFree.press.stats(metaArtiFree.cc.hasPress).sd)./0.1);
% etaArtiFree.lickBoutNaiveNorm.X(metaArtiFree.cc.hasPress, :) = (etaArtiFree.lickBoutNaive.X(metaArtiFree.cc.hasPress, :) - vertcat(etaArtiFree.press.stats(metaArtiFree.cc.hasPress).mean)./0.1) ./ (vertcat(etaArtiFree.press.stats(metaArtiFree.cc.hasPress).sd)./0.1);

maxBoutCycles = 4;
sel = metaArtiFree.cc.isLick;%metaArtiFree.cc.hasPress & metaArtiFree.cc.hasLick & metaArtiFree.cc.isLick;
phase = angle(metaArtiFree.circlick.Z(sel));
amp = abs(metaArtiFree.circlick.Z(sel));
phase(phase < 0) = phase(phase < 0) + 2*pi;
amp = amp(:);
phase = phase(:);
[sortedPhase, I] = sort(phase);
[~, ~] = EphysUnit.plotETA(ax, etaArtiFree.lickBoutNaiveNorm, sel, order=I, ...
    clim=[-5, 5], xlim=[0, 2*pi*maxBoutCycles], hidecolorbar=false);
applyCustomColormap(ax, [-5, 5], hlim=[0.375, 0, 0, -0.375], llim=[0.125, 0.5, 0.5, 0.25], hpwr=.5, lpwr=1, h0=0.33);
xline(ax, (2:2:6)*pi, 'k--', LineWidth=1)
xticks(ax, (0:2:8).*pi);
xticklabels(ax, [{'0'}, arrayfun(@(x) sprintf('%i\\pi', x), 2:2:8, UniformOutput=false)]);
xtickangle(ax, 0)
title(ax, 'Lick-entrained')
ylabel(ax, 'Unit')
xlabel(ax, 'Lick phase')
xlim(ax, [0, 2*pi*maxBoutCycles])
ylim(ax, [0, nnz(sel)+1])
yt = 0:100:nnz(sel);
yt(1) = 1;
if round(yt(end)./100) == round(nnz(sel)./100)
    yt(end) = nnz(sel);
else
    yt(end + 1) = nnz(sel);
end
yt = unique(yt);
yticks(ax, yt)
fontsize(ax, p.fontSize, 'points')
fontname(ax, 'Arial')
axc = ax;

hLetter = text(ax, 0, 0, 'd', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.35, ax.Position(4) + 0.25, 0];


% Phase calculations
sel = find(metaArtiFree.cc.isLick);
rectifiedPhase = phase;
rectifiedPhase(phase < 0) = rectifiedPhase(phase < 0) + 2*pi;
ampThreshold = 0;
isInPhase = phase < 0.25*pi | phase >= 1.75*pi;
isAntiPhase = phase >= 0.75*pi & phase < 1.25*pi;
isFirstQuarterPhase = phase >= 0.25*pi & phase < 0.75*pi;
isThirdQuarterPhase = phase >= 1.25*pi & phase < 1.75*pi;
isHighAmp = amp >= quantile(amp, ampThreshold);
idInPhase = sel(isInPhase & isHighAmp);
idInPhaseLickUp = sel(isInPhase & isHighAmp & metaArtiFree.cc.isLickUp(sel));
idInPhaseLickFlat = sel(isInPhase & isHighAmp & ~metaArtiFree.cc.isLickResponsive(sel));
idInPhaseLickDown = sel(isInPhase & isHighAmp & metaArtiFree.cc.isLickDown(sel));
idInPhasePressUp = sel(isInPhase & isHighAmp & metaArtiFree.cc.isPressUp(sel));
idInPhasePressFlat = sel(isInPhase & isHighAmp & ~metaArtiFree.cc.isPressResponsive(sel));
idInPhasePressDown = sel(isInPhase & isHighAmp & metaArtiFree.cc.isPressDown(sel));
idAntiPhase = sel(isAntiPhase & isHighAmp);
idAntiPhaseLickUp = sel(isAntiPhase & isHighAmp & metaArtiFree.cc.isLickUp(sel));
idAntiPhaseLickFlat = sel(isAntiPhase & isHighAmp & ~metaArtiFree.cc.isLickResponsive(sel));
idAntiPhaseLickDown = sel(isAntiPhase & isHighAmp & metaArtiFree.cc.isLickDown(sel));
idAntiPhasePressUp = sel(isAntiPhase & isHighAmp & metaArtiFree.cc.isPressUp(sel));
idAntiPhasePressFlat = sel(isAntiPhase & isHighAmp & ~metaArtiFree.cc.isPressResponsive(sel));
idAntiPhasePressDown = sel(isAntiPhase & isHighAmp & metaArtiFree.cc.isPressDown(sel));
idFirstQuarterPhase = sel(isFirstQuarterPhase & isHighAmp);
idFirstQuarterPhaseLickUp = sel(isFirstQuarterPhase & isHighAmp & metaArtiFree.cc.isLickUp(sel));
idFirstQuarterPhaseLickFlat = sel(isFirstQuarterPhase & isHighAmp & ~metaArtiFree.cc.isLickResponsive(sel));
idFirstQuarterPhaseLickDown = sel(isFirstQuarterPhase & isHighAmp & metaArtiFree.cc.isLickDown(sel));
idFirstQuarterPhasePressUp = sel(isFirstQuarterPhase & isHighAmp & metaArtiFree.cc.isPressUp(sel));
idFirstQuarterPhasePressFlat = sel(isFirstQuarterPhase & isHighAmp & ~metaArtiFree.cc.isPressResponsive(sel));
idFirstQuarterPhasePressDown = sel(isFirstQuarterPhase & isHighAmp & metaArtiFree.cc.isPressDown(sel));
idThirdQuarterPhase = sel(isThirdQuarterPhase & isHighAmp);
idThirdQuarterPhaseLickUp = sel(isThirdQuarterPhase & isHighAmp & metaArtiFree.cc.isLickUp(sel));
idThirdQuarterPhaseLickFlat = sel(isThirdQuarterPhase & isHighAmp & ~metaArtiFree.cc.isLickResponsive(sel));
idThirdQuarterPhaseLickDown = sel(isThirdQuarterPhase & isHighAmp & metaArtiFree.cc.isLickDown(sel));
idThirdQuarterPhasePressUp = sel(isThirdQuarterPhase & isHighAmp & metaArtiFree.cc.isPressUp(sel));
idThirdQuarterPhasePressFlat = sel(isThirdQuarterPhase & isHighAmp & ~metaArtiFree.cc.isPressResponsive(sel));
idThirdQuarterPhasePressDown = sel(isThirdQuarterPhase & isHighAmp & metaArtiFree.cc.isPressDown(sel));
colors = getColor([1, 3, 2, 4], 4, 0.6);%getColor(1:4, 4, 0.6);

euSel = euArtiFree(sel);
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

% 5e. Osci lick phase distribution histogram
ax = nexttile(layout.child(3).tl, 1 + layout.child(3).cw(2), [layout.child(3).h(1), layout.child(3).w(2)]);
hold(ax, 'on')
edges = -0:2*pi/32:2*pi;
histogram(ax, rectifiedPhase(isInPhase), edges, FaceColor=colors(1, :), FaceAlpha=1, EdgeAlpha=0.5);
histogram(ax, rectifiedPhase(isFirstQuarterPhase), edges, FaceColor=colors(2, :), FaceAlpha=1, EdgeAlpha=0.5);
histogram(ax, rectifiedPhase(isAntiPhase), edges, FaceColor=colors(4, :), FaceAlpha=1, EdgeAlpha=0.5);
histogram(ax, rectifiedPhase(isThirdQuarterPhase), edges, FaceColor=colors(3, :), FaceAlpha=1, EdgeAlpha=0.5);
xticks(ax, 0:pi:2*pi)
yticks(ax, [0, 15])
xlim(ax, pi*[0, 2])
xticklabels(ax, {'0', '\pi', '2\pi'});
xlabel('Lick phase')
ylabel('units     ')
fontsize(ax, p.fontSize, 'points')

hLetter = text(ax, 0, 0, 'e', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.25, ax.Position(4) + 0.25, 0];

ID = {idInPhase; idFirstQuarterPhase; idAntiPhase; idThirdQuarterPhase};
IDSplit = { ...
    {idInPhaseLickUp, idInPhaseLickFlat, idInPhaseLickDown}, {idInPhasePressUp, idInPhasePressFlat, idInPhasePressDown}; ...
    {idFirstQuarterPhaseLickUp, idFirstQuarterPhaseLickFlat, idFirstQuarterPhaseLickDown}, {idFirstQuarterPhasePressUp, idFirstQuarterPhasePressFlat, idFirstQuarterPhasePressDown}; ...
    {idAntiPhaseLickUp, idAntiPhaseLickFlat, idAntiPhaseLickDown}, {idAntiPhasePressUp, idAntiPhasePressFlat, idAntiPhasePressDown}; ...
    {idThirdQuarterPhaseLickUp, idThirdQuarterPhaseLickFlat, idThirdQuarterPhaseLickDown}, {idThirdQuarterPhasePressUp, idThirdQuarterPhasePressFlat, idThirdQuarterPhasePressDown}; ...
    };
ETAMOVEBOUT = {etaArtiFree.correctLickBoutNorm, etaArtiFree.correctPressBoutNorm};
TASKS = ["lick", "press"];
TASKTITLE = ["Self-timed lick", "Self-timed reach"];
PHASENAME = ["2\pi", "1/2\pi", "\pi", "3/2\pi"];
ICOLOR = [1, 2, 4, 3];

% 7f. Reach vs. Lick (scatter)
sz = 7;
% sel = true(size(eu));
x = metaArtiFree.meta.lickNorm;
y = metaArtiFree.meta.pressNorm;
subselResp = metaArtiFree.cc.isLick;
subselNone = ~metaArtiFree.cc.isLick;

mdl = fitlm(metaArtiFree.meta.lickNorm(subselResp), metaArtiFree.meta.pressNorm(subselResp));
fprintf('press vs. lick (osci): LM slope p<%g.\n', mdl.Coefficients.pValue(2))


% 7f right, scatter press vs lick META, color by lick entrainment phase: 
ax = nexttile(layout.child(3).tl, 1 + layout.child(3).ch(2)*sum(layout.child(3).w) + layout.child(3).cw(2), [layout.child(3).h(2), layout.child(3).w(2)]);
hold(ax, 'on')
h = gobjects(5, 1);
h(5) = scatter(ax, x(~metaArtiFree.cc.isLick), y(~metaArtiFree.cc.isLick), sz-2, [0.2 0.2 0.2], Marker='o', MarkerEdgeAlpha=0.25, DisplayName='not-entrained');
for i = 1:4
    sel = ID{i};
    h(i) = scatter(ax, x(sel), y(sel), sz, colors(ICOLOR(i), :), 'filled', Marker='o', MarkerFaceAlpha=0.75, MarkerEdgeAlpha=1, DisplayName=PHASENAME(i));
end

plot(ax, [-10, 10], [0, 0], 'k:');
plot(ax, [0, 0], [-10, 10], 'k:');
plot(ax, [-10, 10], [-10, 10], 'k:')

axis(ax, 'equal')
xlim(ax, [-2, 5])
ylim(ax, [-2, 5])

fontsize(ax, p.fontSize, 'points')
fontname(ax, 'Arial')

ylabel(ax, 'Peri-reach activity (a.u.)', FontSize=p.fontSize)
xlabel(ax, 'Peri-lick activity (a.u.)', FontSize=p.fontSize)

hLetter = text(ax, 0, 0, 'f', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.3, ax.Position(4) + 0.3, 0];


% 7g
nBoutsDisp = 6;
W = 40*[0.5 + nBoutsDisp/8; 0.5 + 0.8 + (nBoutsDisp-1)/8];
tl = tiledlayout(layout.child(3).tl, 2, max(W));
tl.Layout.Tile = 1 + layout.child(3).cw(3); 
tl.Layout.TileSpan = [sum(layout.child(3).h), layout.child(3).w(3)];

AX = gobjects(2, 1);
% AX(1) = nexttile(tl, 1 + W(2) - W(1), [1, W(1)]);
AX(1) = nexttile(tl, 1, [1, W(2)]);
AX(2) = nexttile(tl, 1 + max(W), [1, W(2)]);

% 7d (left) and 7g (right)
% AX = gobjects(4, 3);
for iAx = 2
    iEu = ID{iAx};
    for iTask = 1:2
        % First lick
        % ax = nexttile(layout.right.bottom.tl, (iAx-1)*sum(W) + 1 + sum(W(1:iTask)), [1, W(iTask + 1)]); 
        ax = AX(iTask); 
        hold(ax, 'on')
        t = ETAMOVEBOUT{iTask}.t;
        switch TASKS(iTask)
            case "press"
                t(t>0 & t<=2*pi) = t(t>0 & t<=2*pi) / (2*pi) * 0.8;
                t(t>2*pi) = (t(t>2*pi) - 2*pi) ./ (2*pi) / 8 + 0.8;
            case "lick"
                t(t>0) = t(t>0) ./ (2*pi) / 8;
        end

        for iDir = [1, 3]
            iEuDir = IDSplit{iAx, iTask}{iDir};
            X = ETAMOVEBOUT{iTask}.X(iEuDir, :);
            if isempty(X)
                continue
            end
            X = smoothdata(X, 2, 'gaussian', 5);
            plot(ax, t, mean(X, 1, 'omitnan'), Color=colors(ICOLOR(iAx), :), LineWidth=1.5)
            plot(ax, t, X, Color=[0.15, 0.15, 0.15, 0.05], LineWidth=0.5)
        end

        ylim(ax, [-0.5, 2])
        yticks(ax, [0, 1, 2])

        switch TASKS(iTask)
            case "press"
                xline(ax, [0, 0.8,  0.8+(1:(nBoutsDisp-1))/8], LineStyle=':')
                xlim(ax, [-0.5, 0.8+(nBoutsDisp-1)/8])
                xticks(ax, [-1, -0.5, 0, 0.8,  0.8+(1:(nBoutsDisp-1))/8])
                xticklabels(ax, {'-1', '-0.5', 'reach', '2\pi', '', '', '', '', '12\pi'})
            case "lick"    
                xline(ax, (0:nBoutsDisp)/8, LineStyle=':')       
                xlim(ax, [-0.5-0.675, nBoutsDisp/8])
                xticks(ax, [-1, -0.5, 0, (1:nBoutsDisp)/8])
                xticklabels(ax, {'-1', '-0.5', 'lick', '', '', '', '', '', '12\pi'}) 
        end
        ax.XAxis.TickLabelRotation = 0;

        if iAx == 1
            title(ax, TASKTITLE(iTask))
        end

        % ax.XGrid = 'on';
        hold(ax, 'off')
        fontsize(ax, p.fontSize, 'points')
        % text(ax, 0.025, 1, sprintf('inc(n=%i)', nnz(IDSplit{iAx, iTask}{1})), Unit='normalized', HorizontalAlignment='left', VerticalAlignment='top', Interpreter='none', FontSize=p.fontSize-1)
        % text(ax, 0.025, 0.025, sprintf('dec(n=%i)', nnz(IDSplit{iAx, iTask}{3})), Unit='normalized', HorizontalAlignment='left', VerticalAlignment='bottom', Interpreter='none', FontSize=p.fontSize-1)
        yline(ax, 0, '--')
    end

    % % 6d. Any bout
    % ax = nexttile(layout.left.bottom.tl); AX(iAx, 1) = ax;
    % X = etaArtiFree.lickBoutNaiveNorm.X(iEu, :);
    % X = smoothdata(X, 2, 'gaussian', 5);
    % plot(ax, etaArtiFree.lickBoutNaiveNorm.t, mean(X, 1, 'omitnan'), Color=colors(ICOLOR(iAx), :), LineWidth=1.5)
    % hold(ax, 'on')
    % plot(ax, etaArtiFree.lickBoutNaiveNorm.t, X, Color=[0.15, 0.15, 0.15, 1./nnz(iEu)], LineWidth=0.5)
    % xticks(ax, 0:2*pi:8*pi)
    % xticklabels(ax, [{'0'}, arrayfun(@(x) sprintf('%i\\pi', x), 2:2:8, UniformOutput=false)]);
    % xlim(ax, [0, 8*pi])
    % ylim(ax, [-4, 4])
    % yticks(ax, [-3, 3])
    % ax.XGrid = 'on';
    % fontsize(ax, p.fontSize, 'points')
    % text(ax, 0.05, -0.025, sprintf('n=%i', nnz(iEu)), Unit='normalized', HorizontalAlignment='left', VerticalAlignment='bottom', Interpreter='none', FontSize=p.fontSize-1.5)
    % yline(ax, 0, '--')
end
% xlabel(layout.right.bottom.tl, 'Time from bar/spout contact (s) & lick phase', FontSize=p.fontSize);
% ylabel(layout.right.bottom.tl, 'Normalized spike rate (a.u.)', FontSize=p.fontSize)
% xlabel(layout.left.bottom.tl, 'Lick phase', FontSize=p.fontSize)
% ylabel(layout.left.bottom.tl, 'Normalized spike rate (a.u.)', FontSize=p.fontSize)
ylabel(tl, 'Norm spike rate (a.u.)', FontSize=p.fontSize)
xlabel(tl, '   Time (s)           Lick phase', FontSize=p.fontSize)


ax = AX(1);
hLetter = text(ax, 0, 0, 'g', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.25, ax.Position(4) + 0.25, 0];

% Colorbars need to be done at the end to avoid tiledlayout recursion nonsense
hCb = colorbar(axc);
hCb.Label.String = 'Norm spike rate (a.u.)';
hCb.Label.Position(1) = 0;
hCb.Label.VerticalAlignment = 'bottom';
% hCb.Layout.Tile = 'east';

hCb = axc2.Colorbar;
hCb.Label.String = 'Norm spike rate (a.u.)';
hCb.Label.Position(1) = 0;
hCb.Label.VerticalAlignment = 'bottom';
delete(hCb)

copygraphics(fig, ContentType='vector', BackgroundColor='none')



%% Alternate panel f
fig = figure(Units='inches', Position=[1, 1, 2, 2]);
thisTL = tiledlayout(fig, 1, 2, TileSpacing='compact');
% 5b/metaArtiFree.cc. PETH Reach vs. Lick vs. Osci Lick
% Reach, Lick, Osci Lick, sorted amongst themselves
% close all
ax = gobjects(1, 2);
ax(1) = nexttile(thisTL);
ax(2) = nexttile(thisTL);
% ax(3) = nexttile(layout.middle.right.tl);
% ax(4) = nexttile(layout.middle.right.tl);


sel = metaArtiFree.cc.isLick;
% sel = 1:length(eu);

metaArtiFree.cc.isLickUnresponsiveButUp = ~metaArtiFree.cc.isLickResponsive & metaArtiFree.meta.lickNorm > 0;
metaArtiFree.cc.isLickUnresponsiveButDown = ~metaArtiFree.cc.isLickResponsive & metaArtiFree.meta.lickNorm < 0;
metaArtiFree.cc.isPressUnresponsiveButUp = ~metaArtiFree.cc.isPressResponsive & metaArtiFree.meta.pressNorm > 0;
metaArtiFree.cc.isPressUnresponsiveButDown = ~metaArtiFree.cc.isPressResponsive & metaArtiFree.meta.pressNorm < 0;

groupVar = NaN(length(sel), 1);
groupVar(metaArtiFree.cc.isPressUnresponsiveButDown & metaArtiFree.cc.isLickUp) = 0;
groupVar(metaArtiFree.cc.isPressUnresponsiveButDown & metaArtiFree.cc.isLickUnresponsiveButUp) = 0;
groupVar(metaArtiFree.cc.isPressDown & metaArtiFree.cc.isLickUp) = 0;
groupVar(metaArtiFree.cc.isPressDown & metaArtiFree.cc.isLickUnresponsiveButUp) = 0;
groupVar(metaArtiFree.cc.isPressUp & metaArtiFree.cc.isLickUnresponsiveButDown) = 1;
groupVar(metaArtiFree.cc.isPressUp & metaArtiFree.cc.isLickDown) = 1;
groupVar(metaArtiFree.cc.isPressUnresponsiveButUp & metaArtiFree.cc.isLickUnresponsiveButDown) = 1;
groupVar(metaArtiFree.cc.isPressUnresponsiveButUp & metaArtiFree.cc.isLickDown) = 1;

groupVar(metaArtiFree.cc.isPressUnresponsiveButDown & metaArtiFree.cc.isLickUnresponsiveButDown) = 2;
groupVar(metaArtiFree.cc.isPressUnresponsiveButDown & metaArtiFree.cc.isLickDown) = 2;
groupVar(metaArtiFree.cc.isPressDown & metaArtiFree.cc.isLickUnresponsiveButDown) = 2;
groupVar(metaArtiFree.cc.isPressDown & metaArtiFree.cc.isLickDown) = 2;
groupVar(metaArtiFree.cc.isPressUp & metaArtiFree.cc.isLickUp) = 3;
groupVar(metaArtiFree.cc.isPressUp & metaArtiFree.cc.isLickUnresponsiveButUp) = 3;
groupVar(metaArtiFree.cc.isPressUnresponsiveButUp & metaArtiFree.cc.isLickUp) = 3;
groupVar(metaArtiFree.cc.isPressUnresponsiveButUp & metaArtiFree.cc.isLickUnresponsiveButUp) = 3;
groupVar = groupVar(sel);
[~, order] = EphysUnit.plotETA(ax(1), etaArtiFree.pressNorm, sel, sortGroup=groupVar, ...
    clim=[-1.5, 1.5], xlim=[-2.5, 0.5], sortWindow=[-3, 0], signWindow=[-0.3, 0], ...
    sortThreshold=0.25, negativeSortThreshold=0.25, hidecolorbar=true);
EphysUnit.plotETA(ax(2), etaArtiFree.lickNorm, sel, order=order, sortGroup=groupVar, ...
    clim=[-1.5, 1.5], xlim=[-2.5, 0.5], hidecolorbar=true);

for iAx = 1:2
    applyCustomColormap(ax(iAx), [-1.5, 1.5], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);
    % applyCustomColormap(ax(iAx), [-1.5, 1.5], hlim=[0.375, 0, -0.375, -0.375], llim=[0.1, 1, 1, 0.3], hpwr=.3, lpwr=0.25, h0=0.33);
    % applyCustomColormap(ax(iAx), [-1.5, 1.5], hlim=[0.375, 0, 0, -0.375], llim=[0.125, 0.5, 0.5, 0.25], hpwr=.5, lpwr=1, h0=0.33);
end

N = histcounts(groupVar, -0.5:2:3.5);
yline(ax(1), cumsum(N(1:end-1)) + 1, 'k--');
yline(ax(2), cumsum(N(1:end-1)) + 1, 'k--');

xline(ax(1), 0, 'k--')
xline(ax(2), 0, 'k--')
ylim(ax(1:2), [0, nnz(sel)+1])
yt = 0:100:nnz(sel);
yt(1) = 1;
if round(yt(end)./100) == round(nnz(sel)./100)
    yt(end) = nnz(sel);
else
    yt(end + 1) = nnz(sel);
end
yt = unique(yt);
yticks(ax(1), yt)
yticks(ax(2), [])

title(ax(1), 'Reach')
title(ax(2), 'Lick')
ylabel(thisTL, 'Unit', FontSize=p.fontSize)
ylabel(ax, '')
xlabel(thisTL, 'Time to bar/spout contact (s)', FontSize=p.fontSize)
xlabel(ax, '')
fontsize(ax, p.fontSize, 'points')
fontname(ax, 'Arial')
axc = ax(2);

hLetter = text(ax(1), 0, 0, 'f', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax(1).Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.3, ax(1).Position(4) + 0.3, 0];

clear yt

copygraphics(fig, ContentType='vector', BackgroundColor='none')
