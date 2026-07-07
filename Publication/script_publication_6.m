%% Set root
% ROOTPATH = "C:\SERVER";
ROOTPATH = 'E:\DATA';

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
load(fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot", "20260617_metaRasterData.mat"));
load(fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot\LickVsReach_DLC_dta_rta_25_100_200to800ms_units1to1225_0boots_20260613.mat"));
load(fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot\LickVsReach_DLC_miBoot_200to800ms_1225units_10000boots_20260614.mat")); % contains updated `p`
load(fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot\LickVsReach_DLC_sdBoot_1225units_1000boots_20260703.mat")); % Load bootstrapped statistics for fig6 panels (sd version of cleanclusters boot)
load(fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot\LickVsReach_DLC_maxRectified_1225units_1000boots_20260703.mat")); % Load bootstrapped statistics for fig6 panels (maxRectified xta, mi)
% %% Count number of "clean clusters" by unit
p.mi.minNumTrialsPerCluster = 5; % This is just for example units in Fig 6d, so that small clusters are not culled
p.mi.ignoreSpine = true;
count_dip_triggered_average_units_clusters
load(fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot\LickVsReach_DLC_maxRectifiedPVal_1225units_1000boots_20260703.mat")); % pVal of observing mimr under null bootstrap
p.mi.requireCleanClusters = true;
p.mi.minNumTrialsPerCluster = 10;
p.mi.mergeDirections = true; % false: 6+1 clusters, true: 3+1 clusters (e.g., merging "jaw open" "jaw close" using `or`)


%%
clearvars -except cc clusterSize etaArtiFree euArtiFree kinematics metaArtiFree metaRasterData mi miBoot miMaxRectified miObs nUnits p pVal ROOTPATH sdBoot tests xta xtaMaxRectified nUnits
nUnits = length(xta.dip);
if ~exist('p', 'var') || ~isfield(p, 'fontSize')
    p.fontSize = 8;
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
% %%
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


%%
alpha = 0.05;
less = struct(dip=[], rise=[]);
more = struct(dip=[], rise=[]);
for cn = ["clu1", "all"]
    less.dip = pVal.(cn).dip < alpha/2; % dips have less movement than chance
    more.dip = pVal.(cn).dip > 1-alpha/2; % dips have more movement than chance
    less.rise = pVal.(cn).rise < alpha/2; % rises have less movement than chance
    more.rise = pVal.(cn).rise > 1-alpha/2; % rises have more movement than chance

    switch cn
        case "clu1"
            fprintf("In cluster 1 (no move):\n")
        case "all"
            fprintf("In all clusters:\n")
    end
    fprintf("\t%i units (%.1f%%) have less movement during dips than chance.\n", nnz(less.dip), 100*nnz(less.dip)./nUnits);
    fprintf("\t%i units (%.1f%%) have more movement during dips than chance.\n", nnz(more.dip), 100*nnz(more.dip)./nUnits);
    fprintf("\t%i units (%.1f%%) have less movement during rises than chance.\n", nnz(less.rise), 100*nnz(less.rise)./nUnits);
    fprintf("\t%i units (%.1f%%) have more movement during rises than chance.\n", nnz(more.rise), 100*nnz(more.rise)./nUnits);


    fprintf("\t%i units (%.1f%%) have less movement during dips and rises than chance.\n", nnz(less.dip & less.rise), 100*nnz(less.dip & less.rise)./nUnits);
    fprintf("\t%i units (%.1f%%) have more movement during dips and rises than chance.\n", nnz(more.dip & more.rise), 100*nnz(more.dip & more.rise)./nUnits);

    fprintf("\t%i units (%.1f%%) have less movement during dips but more movement during rises than chance.\n", nnz(less.dip & more.rise), 100*nnz(less.dip & more.rise)./nUnits);
    fprintf("\t%i units (%.1f%%) have more movement during dips but less movement during rises than chance.\n", nnz(more.dip & less.rise), 100*nnz(more.dip & less.rise)./nUnits);
end
clear alpha cn

%% Fig 6
% close all

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

% Figure layout
fig = figure(Units='inches', Position=[1, 1, 7.5, 8]);

clear layout
layout.w = [10, 13];
% layout.h = [4, 8, 7];
layout.h = [7, 17, 12];
layout.ch = cumsum([0, layout.h]).*sum(layout.w);
layout.tl = tiledlayout(fig, sum(layout.h), sum(layout.w), TileSpacing='loose', Padding='loose');

% 1st row (left) (examples)
layout.child(1).h = nEgUnits;
layout.child(1).w = W;
layout.child(1).cw = cumsum([0, W]);
layout.child(1).tl = tiledlayout(layout.tl, sum(layout.child(1).h), sum(layout.child(1).w), TileSpacing='tight', Padding='tight');
l = layout.child(1).tl; l.Layout.Tile = 1 + layout.ch(1); l.Layout.TileSpan = [layout.h(1), layout.w(1)];

% 2nd row (left) (heatmap)
layout.child(2).h = [4, 21];
layout.child(2).w = W;
layout.child(2).cw = cumsum([0, W]);
layout.child(2).tlp = tiledlayout(layout.tl, sum(layout.child(2).h), 1, TileSpacing='tight', Padding='tight');
l = layout.child(2).tlp; l.Layout.Tile = 1 + layout.ch(2); l.Layout.TileSpan = [layout.h(2), layout.w(1)];
layout.child(2).tl = tiledlayout(layout.child(2).tlp, 1, sum(layout.child(2).w), TileSpacing='tight', Padding='tight');
l = layout.child(2).tl; l.Layout.Tile = 1 + layout.child(2).h(1); l.Layout.TileSpan = [layout.child(2).h(2), 1];

% 3rd row (osci)
layout.child(3).h = [9, 20];
layout.child(3).w = [3, 3, 3];
layout.child(3).cw = cumsum([0, layout.child(3).w]);
layout.child(3).ch = cumsum([0, layout.child(3).h]);
layout.child(3).tl = tiledlayout(layout.tl, sum(layout.child(3).h), sum(layout.child(3).w), TileSpacing='loose', Padding='compact');
l = layout.child(3).tl; l.Layout.Tile = 1 + layout.ch(3); l.Layout.TileSpan = [layout.h(3), sum(layout.w)];

% upper right (dip triggered average)
layout.child(4).h = [2, 5, 7];
layout.child(4).w = [1];
layout.child(4).ch = cumsum([0, layout.child(4).h]).*sum(layout.child(4).w);
layout.child(4).tl = tiledlayout(layout.tl, sum(layout.child(4).h), sum(layout.child(4).w), TileSpacing='compact', Padding='compact');
l = layout.child(4).tl; l.Layout.Tile = 1 + layout.w(1); l.Layout.TileSpan = [sum(layout.h(1:2)), layout.w(2)];

% upper right, 1st row: Where do dips occure in trial
layout.child(4).child(1).w = [3, 3, 1];
layout.child(4).child(1).cw = cumsum([0, layout.child(4).child(1).w]);
layout.child(4).child(1).h = [1];
% upper right, 2nd row: 2 example untis (dip, rise): feature x cluster
layout.child(4).child(2).w = [1, 1];
layout.child(4).child(2).h = [1];
% upper right, 3rd row: 2 subrows (dip vs. rise), 3 subcolumns (nuwcc, ncc, xtamr): 
layout.child(4).child(3).w = [2, 1, 2];
layout.child(4).child(3).h = [1, 1];
layout.child(4).child(3).cw = cumsum([0, layout.child(4).child(3).w]);

for iChild = 1:length(layout.child(4).child)
    layout.child(4).child(iChild).tl = tiledlayout(layout.child(4).tl, sum(layout.child(4).child(iChild).h), sum(layout.child(4).child(iChild).w), TileSpacing='loose', Padding='compact');
    l = layout.child(4).child(iChild).tl; l.Layout.Tile = 1 + layout.child(4).ch(iChild); l.Layout.TileSpan = [layout.child(4).h(iChild), layout.child(4).w(1)];
end

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
ylabel(layout.child(1).tl, 'spike rate (sp/s)', FontSize=p.fontSize, Color=[0.2, 0.2, 0.8, 1.0])
xlabel(layout.child(1).tl, 'time (s)', FontSize=p.fontSize)
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

% Plot color bar on top (using patch instead of colorbar so we can make it skinny)
axc2 = nexttile(layout.child(2).tlp, 1, [layout.child(2).h(1), 1]);
patch(axc2, [-1.5, 1.5, 1.5, -1.5], [0, 0, 1, 1], [-1.5, 1.5, 1.5, -1.5], EdgeColor='k');
applyCustomColormap(axc2, [-1.5, 1.5], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);  
xlim(axc2, [-1.5, 1.5])
ylim(axc2, [0, 1])
xticks(axc2, [-1.5, 0, 1.5])
yticks(axc2, [])
title(axc2, 'norm spike rate (a.u.)', FontWeight='normal')
fontsize(axc2, 7, 'points')

ax = gobjects(1, length(XLIM));
for i = 1:length(XLIM)
    ax(i) = nexttile(layout.child(2).tl, CW(i)+1, [1, W(i)]);
end
for iAx = 1:length(ETA)
    hidecb = true;
    EphysUnit.plotETA(ax(iAx), ETA{iAx}, selUnits, xlim=XLIM{iAx}, clim=[-1.5, 1.5], order=sortOrder, hidecolorbar=hidecb);
    applyCustomColormap(ax(iAx), [-1.5, 1.5], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);    
    if ~hidecb
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
xlabel(layout.child(2).tl, "time (s)", FontSize=p.fontSize)
ylabel(layout.child(2).tl, "unit", FontSize=p.fontSize)
% clear etaCombined nDims coeff score explained sortOrder fig tl ax iAx ETASORT NAME ZEROLABEL XLIM hidecp

ax = ax(1);
hLetter = text(ax, 0, 0, 'b', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.25, ax.Position(4) + 0.75, 0];

% 6c. Where do dips occur in reach/lick trials
tl = layout.child(4).child(1).tl;
ax = gobjects(1, 2);
trialTypes = ["press", "lick"];
trialTypeDispNames = ["reach", "lick"];
dirs = ["dip", "rise"];
colors = 'br';
fprintf('Fig6c:');
for iAx = 1:2
    trialType = trialTypes(iAx);
    ax(iAx) = nexttile(tl, 1+layout.child(4).child(1).cw(iAx), [1, layout.child(4).child(1).w(iAx)]);
    hold(ax(iAx), 'on')
    edges = -4:0.1:1;
    for iDir = 1:2
        histogram(ax(iAx), [metaRasterData.(trialType).(dirs(iDir)).t], edges, Normalization='pdf', ...
            FaceColor='none', EdgeColor=colors(iDir), EdgeAlpha=0.8, ...
            DisplayName=dirs(iDir), DisplayStyle='stairs', LineWidth=1);
        fprintf(" %s: %i %ss;", trialTypeDispNames(iAx), nnz(isin([metaRasterData.(trialType).(dirs(iDir)).t], [-4, 1])), dirs(iDir));
    end
    title(ax(iAx), trialTypeDispNames(iAx))
end
fprintf('\n')
% Find total no. reach vs. lick trials
expName = arrayfun(@(rd) strsplit(string(rd.name), '_'), metaRasterData.press.spike, UniformOutput=false);
expName = cellfun(@(c) strjoin(c(1:2), "_"), expName);
[~, ia, ~] = unique(expName);
nTrialsPress = sum(arrayfun(@length, [metaRasterData.press.spike(ia).duration]));
nTrialsLick = sum(arrayfun(@length, [metaRasterData.lick.spike(ia).duration]));
clear expName ia
fprintf('In total: %i units (%i animals, %i sessions, %i press trials, %i lick trials): %i dips, %i rises\n', length(xta.dip), 6, length(unique([xta.dip.iExp])), nTrialsPress, nTrialsLick, length([xta.dip.t0]), length([xta.rise.t0]))

yticks(ax, [])
xlim(ax, [-4, 1])
xlabel(ax(1), "time to bar contact (s)")
xlabel(ax(2), "time to spout contact (s)")
ylabel(tl, 'prob')
fontsize(tl, p.fontSize, 'points')
lgd = legend(ax(2), FontSize=7, Location='eastoutside', IconColumnWidth=9);
lgd.ItemTokenSize = [9, 9];
lgd.Layout.Tile = layout.child(4).child(1).cw(4);

ax = ax(1);
hLetter = text(ax, 0, 0, 'c', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.1, ax.Position(4) + 0.15, 0];

% 6d. Two example units, dip/rise triggered average
egUnitIndices = [354, 374];
egUnitDirs = ["dip", "rise"];
% egUnitNames = ["daisy29_20251024_Channel56_Unit2", "daisy29_20251025_Channel346_Unit1"];
features = ["spikerate"; "Jaw"; "HandL"; "HandR"; "Spine"];
% featureDispName = ["spike"; "jaw"; "handL"; "handR"; "spine"];
clusterDispNameShort = ["~", "jo", "jc", "lhf", "lhb", "rhf", "rhb"]; % cc features are in semantic order already
% clusterDispName = ["no.move", "jaw.open", "jaw.close", "l.hand.fwd", "r.hand.back", "r.hand.fwd", "r.hand.back"]; % cc features are in semantic order already
clusterDispName = string(1:7);
clusterDispNameLong = ["no move", "jaw open", "jaw close", "left hand forward", "left hand back", "right hand forward", "right hand back"]; % cc features are in semantic order already
clusterDispNameMerged = ["jaw", "l.hand", "r.hand"];
assert(isequal(cc.clusterNames, ["no move", "lick start", "lick stop", "left hand reach", "left hand retract", "right hand reach", "right hand retract"]))
featureDispNameShort = ["spk"; "jaw"; "lh"; "rh"; "spn"];
featureDispName = ["spike"; "jaw"; "l.hand"; "r.hand"; "spine"];
featureUnits = ["(a.u.)"; "(DV a.u.)"; "(AP a.u.)"; "(AP a.u.)"; "(DV a.u.)"];
featureSign = [1; -1; 1; 1; -1];
featureAxisDir = ["normal"; "reverse"; "normal"; "normal"; "normal"];
% ylims = {5; 5; 5; 5; 5};
% ylims = cellfun(@(y) y*[-2/3, 1], ylims, UniformOutput=false);
ylims{featureAxisDir=="reverse"} = [-1, 2/3]*5;
ylims = {4; 4; 4; 4; 4};
ylims = cellfun(@(y) y*[-1, 1], ylims, UniformOutput=false);
ylims{featureAxisDir=="reverse"} = [-1, 1]*4;
xt = [0, 1];

nClusters = length(p.mi.semanticClusterOrder);
nFeatures = length(features);

tlp = layout.child(4).child(2).tl;
tl = gobjects(1, 2);
hXLabel = gobjects(1, 2);
hYLabel = gobjects(1, 2);
fprintf("Fig6d:")
for iTl = 1:2
    iUnit = egUnitIndices(iTl);
    fprintf(" Example unit %i: ", iTl)
    dir = egUnitDirs(iTl);
    tl(iTl) = tiledlayout(tlp, nFeatures, nClusters+1, TileSpacing='compact', Padding='tight', TileIndexing='columnmajor', DefaultAxesFontSize=8);
    tl(iTl).Layout.Tile = iTl;

    ax = gobjects(nFeatures, nClusters);
    for iClu = 1:nClusters
        for iFeat = 1:nFeatures
            ax(iFeat, iClu) = nexttile(tl(iTl));
            ax(iFeat, iClu).Box = 'off';
            ax(iFeat, iClu).XAxis.Visible = 'off';
            ax(iFeat, iClu).YAxis.Visible = 'off';
            ax(iFeat, iClu).YAxis.Label.Visible = 'on';
            ax(iFeat, iClu).Color = 'none';
            if iTl == 1 && iClu == 1 && iFeat == 1
                axLetter = ax(iFeat, iClu);
            end
        end
    end

    idx = mi.(dir)(iUnit).idx;
    
    for iClu = 1:length(p.mi.semanticClusterOrder)
        k = p.mi.semanticClusterOrder(iClu); % raw cluster index
        nTrials = clusterSize(iUnit).(dir).n(k);
        if isnan(nTrials)% || nTrials < p.mi.minNumTrialsPerCluster
            set(ax(:, iClu), Visible=false);
            continue
        end

        set(ax(:, iClu), Visible=true);

        for iFeat = 1:length(features) % k: semantic cluster index
            hold(ax(iFeat, iClu), 'on')
            fn = features(iFeat);
            sel = idx==k;
            % c = getColor(iFeat, length(features), 0.7);
            c = 'k';
            X = xta.(dir)(iUnit).(fn).X(sel, :);
            t = xta.(dir)(iUnit).(fn).t;
            mu = featureSign(iFeat)*mean(X, 1, 'omitnan');
            
            % Check significance of cluster movement index against bootstrap
            if iFeat > 1
                hVal = cc.data(iUnit).(dir).h(iClu, iFeat-1);
                nStars = nnz(hVal);
            else
                fprintf("%s=%i;", clusterDispName(iClu), nnz(sel))
                switch dir
                    case "dip"
                        text(ax(iFeat, iClu), 0.3, 1, string(nnz(sel)), FontSize=5, ...
                            HorizontalAlignment='center', VerticalAlignment='bottom')                        
                    case "rise"
                        text(ax(iFeat, iClu), 0.3, ylims{iFeat}(2)/2, string(nnz(sel)), FontSize=5, ...
                            HorizontalAlignment='center', VerticalAlignment='bottom')
                end
            end

            plot(ax(iFeat, iClu), t, mu, Color=c, LineStyle='-', LineWidth=1);
            if iFeat > 1 && nStars >= 1
                selT = isin(t, p.mi.windowPre) | isin(t, p.mi.windowPost);
                plot(ax(iFeat, iClu), t(selT), mu(selT), Color='red', LineStyle='-', LineWidth=1);
                clear selT
            end
            clear nStars
            plot(ax(iFeat, iClu), [0, 0], ylims{iFeat}, '-', Color=[0.15, 0.15, 0.15, 0.25])
            plot(ax(iFeat, iClu), [-0.5, 0.5], [0, 0], '-', Color=[0.15, 0.15, 0.15, 0.25])
            if iClu == 1
                ylabel(ax(iFeat, iClu), sprintf("%s", featureDispName(iFeat)), Rotation=30, FontSize=p.fontSize)
            end
            if iFeat == 1
                title(ax(iFeat, iClu), strsplit(clusterDispName(iClu), "\\n"), Interpreter='none', FontWeight='normal', Color='black')
            end
            % title(tl(iTl), 'clusters', FontWeight='normal', Color='black', FontSize=p.fontSize)
            ylim(ax(iFeat, iClu), ylims{iFeat})
            hold(ax(iFeat, iClu), 'off')
            ax(iFeat, iClu).YAxis.Direction = featureAxisDir(iFeat);
        end
    end
    xticks(ax, [])
    yticks(ax, [])
    xlim(ax, [-0.5, 0.5])

    % Add scale bar
    ax = nexttile(tl(iTl), (nClusters+1)*nFeatures, [1, 1]);
    ylim(ax, ylims{1})
    xlim(ax, [-0.5, 0.5])
    xticks(ax, [])
    yticks(ax, [])
    ax.YAxisLocation = 'right';
    ax.Color = 'none';
    xlabel(ax, "1 s", FontSize=8)
    ylabel(ax, sprintf("%i sd", diff(ylims{1})), FontSize=8)

    % hXLabel(iTl) = xlabel(tl(iTl), 'cluster', FontSize=8);
    % hYLabel(iTl) = ylabel(tl(iTl), 'bodypart', FontSize=8);
end

fprintf('\n')
% fontsize(tlp, p.fontSize, 'points')

hLetter = text(axLetter, 0, 0, 'd', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
axLetter.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.1, axLetter.Position(4) + 0.3, 0];

% 6e(l), 6e(r), 6f shares a tiled layout, we'll make the axes here:
tl = layout.child(4).child(3).tl;
tl.TileSpacing = 'tight';
tl.TileIndexing = 'rowmajor';
for iRow = 1:2
    for iCol = 1:3
        ax(iRow, iCol) = nexttile(tl, (iRow-1)*layout.child(4).child(3).cw(end) + 1 + layout.child(4).child(3).cw(iCol), [layout.child(4).child(3).h(iRow), layout.child(4).child(3).w(iCol)]);
    end
end

% 6e (left). Histogram count of units per movement type (xticks are cluster names)
nUnits = length(xta.dip);
dirs = ["dip", "rise"];
colors = [0, 0, 1; 1, 0, 0];
iCol = 1;
lgd = gobjects(1, 2);
yl = [0, 1];
for iRow = 1:2
    hold(ax(iRow, iCol), 'on')
    dir = dirs(iRow);
    if p.mi.mergeDirections
        x = 1:3;
    else
        x = 1:length(clusterDispName);
    end
    % nUnits with clean cluster
    if p.mi.requireCleanClusters
        nuwcc = arrayfun(@(ccData) ccData.(dir).n(:, 1) .* uint16(ccData.(dir).clean), cc.data, UniformOutput=false);
    else
        nuwcc = arrayfun(@(ccData) ccData.(dir).n(:, 1), cc.data, UniformOutput=false);
    end
    nuwcc = cat(2, nuwcc{:})';
    nuwcc = nuwcc>=p.mi.minNumTrialsPerCluster;
    if p.mi.mergeDirections
        nuwcc = cat(2, any(nuwcc(:, 2:3), 2), any(nuwcc(:, 4:5), 2), any(nuwcc(:, 6:7), 2));
    end
    y = sum(nuwcc, 1)./nUnits;

    % nUnits with clean cluster, bootstrapped
    if p.mi.requireCleanClusters
        nuwccBoot = arrayfun(@(ccData) ccData.(dir).n .* uint16(ccData.(dir).clean), sdBoot, UniformOutput=false);
    else
        nuwccBoot = arrayfun(@(ccData) ccData.(dir).n, sdBoot, UniformOutput=false);
    end
    nuwccBoot = cat(3, nuwccBoot{:});
    nuwccBoot = permute(nuwccBoot, [3, 2, 1]); % nUnits x nClusters x nBoot
    nuwccBoot = nuwccBoot(:, p.mi.semanticClusterOrder, :); % nUnits x nClusters x nBoot, clusters in semantic order
    nuwccBoot = nuwccBoot>=p.mi.minNumTrialsPerCluster;
    if p.mi.mergeDirections
        nuwccBoot = cat(2, any(nuwccBoot(:, 2:3, :), 2), any(nuwccBoot(:, 4:5, :), 2), any(nuwccBoot(:, 6:7, :), 2));
    end
    yBoot = squeeze(sum(nuwccBoot, 1))./nUnits; % nClusters x nBoot
    muBoot = mean(yBoot, 2, 'omitnan')';
    ciBoot = quantile(yBoot, [0.05/7/2, 1-0.05/7/2], 2)';

    if p.mi.mergeDirections
        % edges = 0:2:nUnits+1;
        % for iClu = 1:3
        %     thisNBoot = histcounts(yBoot(iClu, :), edges);
        %     xVerticesBoot = reshape([thisNBoot; thisNBoot]./max(thisNBoot), [], 1);
        %     yVerticesBoot = reshape([edges(1:end-1); edges(2:end)], [], 1);
        %     patch(ax(iRow, iCol), x(iClu)-xVerticesBoot*0.5, yVerticesBoot, 'k', FaceColor=[0.15, 0.15, 0.15], EdgeColor='none', FaceAlpha=0.5)
        %     patch(ax(iRow, iCol), x(iClu)+xVerticesBoot*0.5, yVerticesBoot, 'k', FaceColor=[0.15, 0.15, 0.15], EdgeColor='none', FaceAlpha=0.5)
        % end
        hBar = bar(ax(iRow, iCol), x, [y; muBoot], FaceAlpha=0.5, EdgeAlpha=0.8, BarWidth=1);
        set(hBar(1), FaceColor=colors(iRow, :), EdgeColor=colors(iRow, :));
        set(hBar(2), FaceColor='none', EdgeColor=[0.15, 0.15, 0.15]);
        hError = errorbar(ax(iRow, iCol), x+1/7, muBoot, muBoot-ciBoot(1, :), ciBoot(2, :)-muBoot, ...
            LineStyle='none', Color=[0.15, 0.15, 0.15, 0.8], CapSize=5, Clipping=false);
    else
        % edges = 0:2:nUnits+1;
        % for iClu = 2:nClusters
        %     thisNBoot = histcounts(yBoot(iClu, :), edges);
        %     xVerticesBoot = reshape([thisNBoot; thisNBoot]./max(thisNBoot), [], 1);
        %     yVerticesBoot = reshape([edges(1:end-1); edges(2:end)], [], 1);
        %     patch(ax(iRow, iCol), x(iClu)-xVerticesBoot*0.5, yVerticesBoot, 'k', FaceColor=[0.15, 0.15, 0.15], EdgeColor='none', FaceAlpha=0.5)
        %     patch(ax(iRow, iCol), x(iClu)+xVerticesBoot*0.5, yVerticesBoot, 'k', FaceColor=[0.15, 0.15, 0.15], EdgeColor='none', FaceAlpha=0.5)
        % end
        hBar = bar(ax(iRow, iCol), x(2:end), [y(2:end); muBoot(2:end)], FaceAlpha=0.5, EdgeAlpha=0.8, BarWidth=1);
        set(hBar(1), FaceColor=colors(iRow, :), EdgeColor=colors(iRow, :));
        set(hBar(2), FaceColor='none', EdgeColor=[0.15, 0.15, 0.15]);
        hError = errorbar(ax(iRow, iCol), x(2:end)+1/7, muBoot, muBoot(2:end)-ciBoot(1, 2:end), ciBoot(2, 2:end)-muBoot(2:end), ...
            LineStyle='none', Color=[0.15, 0.15, 0.15, 0.8], CapSize=5, Clipping=false);
    end
    if iRow == 1
        hBar(3) = bar(ax(iRow, iCol), NaN, NaN, FaceAlpha=0.5, EdgeAlpha=0.8, BarWidth=1, FaceColor=colors(2, :), EdgeColor=colors(2, :));
        h = [hBar([1, 3, 2]), hError];
        h(1).DisplayName = 'dip';
        h(2).DisplayName = 'rise';
        h(3).DisplayName = 'boot';
        h(4).DisplayName = '95%CI';
    end

    if p.mi.mergeDirections
        % xticks(ax(1, iCol), [])
        % xticks(ax(2, iCol), 1:length(clusterDispNameMerged))
        % xticklabels(ax(2, iCol), ["", "", ""])
        % xtickangle(ax(2, iCol), 0)
        xticks(ax(:, iCol), [])
        for i = 1:length(clusterDispNameMerged)
            text(ax(2, iCol), x(i), -0, clusterDispNameMerged(i), Color=[0.15, 0.15, 0.15], ...
                VerticalAlignment='top', HorizontalAlignment='right', Rotation=30, ...
                FontWeight='normal', FontName='Arial')
        end
        xlim(ax(iRow, iCol), [0.3, length(clusterDispNameMerged)+0.7])
    else
        xticks(ax(1, iCol), [])
        xticks(ax(2, iCol), 1:length(clusterDispName))
        xticklabels(ax(2, iCol), clusterDispName)
        xtickangle(ax(2, iCol), 0)
        xlim(ax(iRow, iCol), [1.3, length(clusterDispName)+0.7])
    end

    ylim(ax(iRow, iCol), yl)
    yticks(ax(iRow, iCol), [0, 0.5, 1])
    ylabel(ax(iRow, iCol), "fraction of units")
    hold(ax(iRow, iCol), 'off')

    lgd = legend(ax(1, iCol), h, Location='north', Orientation='horizontal', FontSize=p.fontSize-1, IconColumnWidth=7, NumColumns=2);
    lgd.ItemTokenSize = [7, 7];
    lgd.Position(1:2) = lgd.Position(1:2) + [0.05, 0.16];
end
% xlabel(ax(2, iCol), "bodyparts")

% 6e (right). Histogram count of no. movement types for each unit (xticks are numCleanClusters 0-6)
iCol = 2;
if p.mi.mergeDirections
    edges = -0.5:1:3.5;
    centers = 0:3;
    maxC = 0.2;
else
    edges = -0.5:1:6.5;
    centers = 0:6;
    maxC = 0.1;
end
ncc = struct(dip=[], rise=[]);
for dir = dirs
    if p.mi.requireCleanClusters
        ncc.(dir) = arrayfun(@(ccData) ccData.(dir).clean & ccData.(dir).n(:, 1)>=p.mi.minNumTrialsPerCluster, cc.data, UniformOutput=false);
    else
        ncc.(dir) = arrayfun(@(ccData) ccData.(dir).n(:, 1)>=p.mi.minNumTrialsPerCluster, cc.data, UniformOutput=false);
    end
    ncc.(dir) = cat(2, ncc.(dir){:})';
    if p.mi.mergeDirections
        ncc.(dir) = [any(ncc.(dir)(:, 2:3), 2), any(ncc.(dir)(:, 4:5), 2), any(ncc.(dir)(:, 6:7), 2)];
        ncc.(dir) = sum(ncc.(dir), 2);
    else
        ncc.(dir) = sum(ncc.(dir)(:, 2:end), 2);
    end
end
nccBoot = struct(dip=[], rise=[]);
for dir = dirs
    if p.mi.requireCleanClusters
        nccBoot.(dir) = arrayfun(@(ccData) ccData.(dir).n>=p.mi.minNumTrialsPerCluster & ccData.(dir).clean, sdBoot, UniformOutput=false);
    else
        nccBoot.(dir) = arrayfun(@(ccData) ccData.(dir).n>=p.mi.minNumTrialsPerCluster, sdBoot, UniformOutput=false);
    end
    nccBoot.(dir) = permute(cat(3, nccBoot.(dir){:}), [3, 2, 1]); % units x clusters x boots

    if p.mi.mergeDirections
        nccBoot.(dir) = nccBoot.(dir)(:, p.mi.semanticClusterOrder, :);
        nccBoot.(dir) = cat(2, any(nccBoot.(dir)(:, 2:3, :), 2), any(nccBoot.(dir)(:, 4:5, :), 2), any(nccBoot.(dir)(:, 6:7, :), 2));
        nccBoot.(dir) = squeeze(sum(nccBoot.(dir), 2));
    else
        nccBoot.(dir) = squeeze(sum(nccBoot.(dir)(:, 2:end, :), 2));
    end
end
for iRow = 1:2
    hold(ax(iRow, iCol), 'on')
    dir = dirs(iRow);
    n = histcounts(ncc.(dir), edges)./nUnits;
    nBoot = arrayfun(@(iBoot) histcounts(nccBoot.(dir)(:, iBoot), edges), 1:p.sd.nBoot, UniformOutput=false);
    nBoot = cat(1, nBoot{:})./nUnits;
    ciBoot = quantile(nBoot, [0.05/2, 1-0.05/2], 1);
    nBoot = mean(nBoot, 1);
    histogram(ax(iRow, iCol), BinCounts=n, BinEdges=edges, FaceColor=colors(iRow, :), FaceAlpha=0.5, EdgeColor=colors(iRow, :));
    histogram(ax(iRow, iCol), BinCounts=nBoot, BinEdges=edges, EdgeColor=[0.15, 0.15, 0.15], EdgeAlpha=1, DisplayStyle='stairs');
    errorbar(ax(iRow, iCol), centers, nBoot, nBoot-ciBoot(1, :), ciBoot(2, :)-nBoot, LineStyle='none', Color=[0.15, 0.15, 0.15, 0.8], CapSize=5, Clipping=false)

    % ylabel(ax(iRow, iCol), "no. units")
    hold(ax(iRow, iCol), 'off')
end
xticks(ax(1, iCol), [])
xticks(ax(2, iCol), centers)
xtickangle(ax(2, iCol), 0)
if p.mi.mergeDirections
    xlabel(ax(2, iCol), "no. bodyparts")
else
    xlabel(ax(2, iCol), "no. clusters")
end
xlim(ax(:, iCol), [edges(1), edges(end)])
ylim(ax(:, iCol), yl)
yticks(ax(:, iCol), [0, 0.5, 1])

% 6f: xtamr (maxRectifiedXTA)
alpha = 0.05;
aggr(nUnits) = struct(dip=[], rise=[]);
for iUnit = 1:nUnits
    for dir = ["dip", "rise"]
        aggr(iUnit).(dir).xObs = mean(xtaMaxRectified(iUnit).(dir).clu1.X, 1, 'omitnan');
        aggr(iUnit).(dir).tObs = xtaMaxRectified(iUnit).(dir).clu1.t;

        if isfield(xtaMaxRectified(iUnit).(dir).clu1, 'XBoot')
            XBoot = xtaMaxRectified(iUnit).(dir).clu1.XBoot';
            aggr(iUnit).(dir).tBoot = xtaMaxRectified(iUnit).(dir).clu1.tBoot;
            aggr(iUnit).(dir).muBoot = mean(XBoot, 1, 'omitnan');
            aggr(iUnit).(dir).ciBoot = quantile(XBoot, [alpha/2, 1-alpha/2], 1);
        else
            aggr(iUnit).(dir).tBoot = [];
            aggr(iUnit).(dir).muBoot = NaN(size(aggr(1).dip.muBoot), 'single');
            aggr(iUnit).(dir).ciBoot = NaN(size(aggr(1).dip.ciBoot), 'single');
        end
    end
end
clear iUnit dir XBoot
iCol = 3;
h = gobjects(4, 1);
xl = [-0.1, 0.3];
yl = [0, 0.4];
for iRow = 1:2
    hold(ax(iRow, iCol), 'on')
    dir = dirs(iRow);
    tObs = aggr(1).(dir).tObs;
    tBoot = aggr(1).(dir).tBoot;
    selLo =  pVal.clu1.(dir) < alpha/2;
    selHi =  pVal.clu1.(dir) > 1 - alpha/2;

    xObs = arrayfun(@(aggr) aggr.(dir).xObs, aggr, UniformOutput=false);
    xObs = cat(1, xObs{:}); % untis x timestamps
    xObsLo = mean(xObs(selLo, :), 1, 'omitnan');
    xObsHi = mean(xObs(selHi, :), 1, 'omitnan');
    xObsNull = mean(xObs(~selLo&~selHi, :), 1, 'omitnan');
    xObs = mean(xObs, 1, 'omitnan');

    ciBoot = arrayfun(@(aggr) aggr.(dir).ciBoot, aggr, UniformOutput=false);
    ciBoot = cat(3, ciBoot{:}); % quantile x timestamps x units (2x60x1225)
    ciBootLo = mean(ciBoot(:, :, selLo), 3, 'omitnan');
    ciBootHi = mean(ciBoot(:, :, selHi), 3, 'omitnan');
    ciBootHiNull = mean(ciBoot(:, :, ~selLo&~selHi), 3, 'omitnan');
    ciBoot = mean(ciBoot, 3, 'omitnan');

    h(1) = patch(ax(iRow, iCol), [tBoot, flip(tBoot)], [ciBoot(1, :), flip(ciBoot(2, :))], [0.15, 0.15, 0.15], FaceAlpha=0.1, EdgeAlpha=0.5, DisplayName=sprintf('%i%% CI', round(100*(1-alpha))));
    h(2) = plot(ax(iRow, iCol), tObs, xObs, 'k-', DisplayName="all units", LineWidth=1.5);
    h(3) = plot(ax(iRow, iCol), tObs, xObsHi, 'r:', DisplayName="move+", LineWidth=1);
    h(4) = plot(ax(iRow, iCol), tObs, xObsLo, 'b:', DisplayName="move-", LineWidth=1);
    % h(5) = plot(ax(iRow, iCol), tObs, xObsNull, 'k--', DisplayName="move~", LineWidth=1.5);
    xline(ax(iRow, iCol), 0, 'k--')
    yline(ax(iRow, iCol), 0, 'k--')
    if iRow == 1
        lgd = legend(ax(iRow, iCol), h, Location="north", FontSize=p.fontSize-1, IconColumnWidth=7, NumColumns=2);
        lgd.ItemTokenSize = [7, 7];
    end


    fprintf("xtamr (%s): ", dir);
    fprintf("\tmove+: %i units (%.1f%%);", nnz(selHi), nnz(selHi)/nUnits*100)
    fprintf("\tmove-: %i units (%.1f%%);", nnz(selLo), nnz(selLo)/nUnits*100)
    fprintf("\tmove~: %i units (%.1f%%);", nnz(~selLo&~selHi), nnz(~selLo&~selHi)/nUnits*100)
    fprintf("\n")

    ylabel(ax(iRow, iCol), "movement")

    % ax(iRow, iCol).XLimMode = 'manual';
    % ax(iRow, iCol).YLimMode = 'manual';
    xlim(ax(iRow, iCol), xl)
    ylim(ax(iRow, iCol), yl)
    hold(ax(iRow, iCol), 'off')
end
xticks(ax(1, iCol), [])
xlabel(ax(2, iCol), "time to dip/rise onset (s)")
lgd.Position(2) = lgd.Position(2) + 0.16;
% title(ax(1, :), 'dip')
% title(ax(2, :), 'rise')


% 6e/f - common labels/fonts
fontsize(tl, 8, 'points')

axLetter = ax(1, 1);
hLetter = text(axLetter, 0, 0, 'e', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
axLetter.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.15, axLetter.Position(4) + 0.35, 0];

axLetter = ax(1, 3);
hLetter = text(axLetter, 0, 0, 'f', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
axLetter.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.1, axLetter.Position(4) + 0.35, 0];

% 6h. Correct lick bout lick histogram
nBoutsDisp = 6;
ax = nexttile(layout.child(3).tl, [layout.child(3).h(1), layout.child(3).w(1)]);
counts = lickHist.correctLickOsci.count;
smoothedCounts = smoothdata(counts, 'gaussian', 25);
[pks, locs] = findpeaks(smoothedCounts, lickHist.correctLickOsci.t, MinPeakProminence=0.5);
histogram(ax, BinEdges=lickHist.correctLickOsci.edges, BinCounts=counts, Normalization='probability', EdgeColor='none', FaceColor='black', FaceAlpha=0.9)
hold(ax, 'on')
xline(ax, locs, LineStyle=':', Color=[0, 0.5, 0.5, 1], LineWidth=2.5)
ylabel(ax, ' lick\newlineprob')
xlabel(ax, 'time from first lick (s)')
xticks(ax, -1:0.5:2)
xlim(ax, [0, 1/8*nBoutsDisp])
yticks(ax, [])
ax.Box = 'off';
fontsize(ax, p.fontSize, 'points')

hLetter = text(ax, 0, 0, 'g', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.35, ax.Position(4) + 0.1, 0];


% 6h ETA Heatmap osci lick
tlp = tiledlayout(layout.child(3).tl, 1, 6+1, TileSpacing='tight');
tlp.Layout.Tile = 1 + layout.child(3).ch(2)*sum(layout.child(3).w); tlp.Layout.TileSpan = [layout.child(3).h(2), layout.child(3).w(1)];

ax = nexttile(tlp, 1, [1, 6]);

etaArtiFree.lickBoutNaiveNorm = etaArtiFree.lickBoutNaive;
etaArtiFree.lickBoutNaiveNorm.X = normalize(etaArtiFree.lickBoutNaive.X, 2, 'zscore', 'robust');

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
title(ax, '')
ylabel(ax, 'unit')
xlabel(ax, 'lick phase')
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

hLetter = text(ax, 0, 0, 'h', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
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

% 6i. Osci lick phase distribution histogram
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
xlabel('lick phase')
ylabel('units     ')
fontsize(ax, p.fontSize, 'points')

hLetter = text(ax, 0, 0, 'i', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
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

% 6j. Reach vs. Lick (scatter)
sz = 7;
% sel = true(size(eu));
x = metaArtiFree.meta.lickNorm;
y = metaArtiFree.meta.pressNorm;
subselResp = metaArtiFree.cc.isLick;
subselNone = ~metaArtiFree.cc.isLick;

mdl = fitlm(metaArtiFree.meta.lickNorm(subselResp), metaArtiFree.meta.pressNorm(subselResp));
fprintf('Fig 6j: press vs. lick (osci): LM slope p<%g.\n', mdl.Coefficients.pValue(2))


% right, scatter press vs lick META, color by lick entrainment phase: 
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

ylabel(ax, 'peri-reach activity (a.u.)', FontSize=p.fontSize)
xlabel(ax, 'peri-lick activity (a.u.)', FontSize=p.fontSize)

hLetter = text(ax, 0, 0, 'j', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.3, ax.Position(4) + 0.3, 0];


% 6k
nBoutsDisp = 6;
W = 40*[0.5 + nBoutsDisp/8; 0.5 + 0.8 + (nBoutsDisp-1)/8];
tl = tiledlayout(layout.child(3).tl, 2, max(W));
tl.Layout.Tile = 1 + layout.child(3).cw(3); 
tl.Layout.TileSpan = [sum(layout.child(3).h), layout.child(3).w(3)];

AX = gobjects(2, 1);
% AX(1) = nexttile(tl, 1 + W(2) - W(1), [1, W(1)]);
AX(1) = nexttile(tl, 1, [1, W(2)]);
AX(2) = nexttile(tl, 1 + max(W), [1, W(2)]);

% 6d (left) and 7g (right)
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
end
ylabel(tl, 'norm spike rate (a.u.)', FontSize=p.fontSize)
xlabel(tl, '   time (s)           lick phase', FontSize=p.fontSize)


ax = AX(1);
hLetter = text(ax, 0, 0, 'k', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.25, ax.Position(4) + 0.25, 0];



% Colorbars need to be done at the end to avoid tiledlayout recursion nonsense
hCb = colorbar(axc);
hCb.Label.String = 'norm spike rate (a.u.)';
hCb.Label.Position(1) = 0;
hCb.Label.VerticalAlignment = 'bottom';

fontname(fig, 'Arial')

% drawnow()


% Adjust locations of some tiledlayout x/y labels at the end
% for iTl = 1:2
%     hXLabel(iTl).Text.Position = hXLabel(iTl).Text.Position + [-0.025, 0.4, 0];
% end
% for iTl = 1:2
%     hYLabel(iTl).Text.Position = hYLabel(iTl).Text.Position + [0.02, -0.025, 0];
% end


% exportgraphics(fig, 'test.emf', ContentType='vector', BackgroundColor='none', Padding='Figure')
copygraphics(fig, ContentType='vector', BackgroundColor='none')


%% S6

% from 6e
nUnits = length(xta.dip);
tl = layout.child(4).child(3).tl;
tl.TileSpacing = 'tight';
ax = gobjects(2, 2);
dirs = ["dip", "rise"];
colors = [0, 0, 1; 1, 0, 0];
yl = [0, 100];
iRow = 1;
lgd = gobjects(1, 2);
hPatch = gobjects(2, 2);
for iCol = 1:2
    dir = dirs(iCol);
    ax(iRow, iCol) = nexttile(tl);
    hold(ax(iRow, iCol), 'on')

    x = 1:length(clusterDispName);

    if p.mi.requireCleanClusters
        cs = arrayfun(@(ccData) ccData.(dir).n(:, 1) .* uint16(ccData.(dir).clean), cc.data, UniformOutput=false);
    else
        cs = arrayfun(@(ccData) ccData.(dir).n(:, 1), cc.data, UniformOutput=false);
    end
    cs = double(cat(2, cs{:})');
    cs(cs<p.mi.minNumTrialsPerCluster) = NaN;
    % mu = mean(cs, 1, 'omitnan');

    % Bootstrapped cluster size (null)
    if p.mi.requireCleanClusters
        csBoot = arrayfun(@(d) d.(dir).n .* uint16(d.(dir).clean), sdBoot, UniformOutput=false);
    else
        csBoot = arrayfun(@(d) d.(dir).n, sdBoot, UniformOutput=false);
    end
    csBoot = double(cat(3, csBoot{:}));
    csBoot = permute(csBoot, [1, 3, 2]); % nBoot x nUnits x nClusters
    csBoot(csBoot<p.mi.minNumTrialsPerCluster) = NaN;
    csBoot = csBoot(:, :, p.mi.semanticClusterOrder); % put clusters in semantic order
    % muBoot = reshape(mean(csBoot, [1, 2], 'omitnan'), 1, nClusters);
    % ciBoot = quantile(squeeze(mean(csBoot, 2, 'omitnan')), [0.05/7/2, 1-0.05/7/2], 1);

    edges = 0:2:1226;
    for iClu = 2:nClusters
        thisN = histcounts(cs(:, iClu), edges);
        thisNBoot = histcounts(reshape(csBoot(:, :, iClu), [], 1), edges);
        xVertices = reshape([thisN; thisN]./max(thisN), [], 1);
        yVertices = reshape([edges(1:end-1); edges(2:end)], [], 1);
        xVerticesBoot = reshape([thisNBoot; thisNBoot]./max(thisNBoot), [], 1);
        yVerticesBoot = reshape([edges(1:end-1); edges(2:end)], [], 1);
        xVertices = x(iClu) + xVertices*0.5;
        xVerticesBoot = x(iClu) - xVerticesBoot*0.5;
        hPatch(1, iCol) = patch(ax(iRow, iCol), xVertices, yVertices, 'k', FaceColor=colors(iCol, :), EdgeColor='none', FaceAlpha=0.5, DisplayName='obs');
        hPatch(2, iCol) = patch(ax(iRow, iCol), xVerticesBoot, yVerticesBoot, 'k', FaceColor=[0.15, 0.15, 0.15], EdgeColor='none', FaceAlpha=0.5, DisplayName='boot');
    end

    % % The first cluster is too tall to draw, we do with clipping and write n+-sd on top.
    % violinplot(ax(iRow, iCol), x(1), cs(:, 1), DensityDirection='positive', FaceColor=colors(iCol, :), EdgeColor=colors(iCol, :), FaceAlpha=0.5, Clipping='on')
    % violinplot(ax(iRow, iCol), x(1), reshape(csBoot(:, :, 1), [], 1), DensityDirection='negative', FaceColor='none', EdgeColor=[0.15, 0.15, 0.15], Clipping='on')
    % % Draw clusters 2:7 without clipping so errorbars display correctly
    % violinplot(ax(iRow, iCol), x(2:end), cs(:, 2:end), DensityDirection='positive', FaceColor=colors(iCol, :), EdgeColor=colors(iCol, :), FaceAlpha=0.5, Clipping='on')
    % violinplot(ax(iRow, iCol), x(2:end), reshape(csBoot(:, :, 2:end), [], nClusters-1), DensityDirection='negative', FaceColor='none', EdgeColor=[0.15, 0.15, 0.15], Clipping='on')

    ylabel(ax(iRow, iCol), sprintf("%ss", dir))
    xlim(ax(iRow, iCol), [1.3, length(clusterDispName)+0.7])
    ylim(ax(iRow, iCol), yl)
    hold(ax(iRow, iCol), 'off')
    lgd(iCol) = legend(ax(iRow, iCol), hPatch(:, iCol), Location='northeast', Orientation='horizontal', IconColumnWidth=7, FontSize=7);
    lgd(iCol).ItemTokenSize = [7, 7];
    lgd(iCol).Position(2) = 0.72;
end
lgd(1).Position(1) = 0.26;
lgd(2).Position(1) = 0.72;

%from 6f
iAx = 3;
[nBoot, xEdges, yEdges] = histcounts2(mean(nccBoot.dip, 2, 'omitnan'), mean(nccBoot.rise, 2, 'omitnan'), edges, edges);
nBoot = nBoot./nUnits;
histogram2(ax(iAx), XBinEdges=xEdges, YBinEdges=yEdges, BinCounts=nBoot, ...
    DisplayStyle='tile', ShowEmptyBins=true);
applyCustomColormap(ax(iAx), [0, maxC], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.5, h0=0.33);
title(ax(iAx), 'bootstrap')

iAx = 4;
[n, xEdges, yEdges] = histcounts2(ncc.dip, ncc.rise, edges, edges);
n = n./nUnits;
histogram2(ax(iAx), XBinEdges=xEdges, YBinEdges=yEdges, BinCounts=n, ...
    DisplayStyle='tile', ShowEmptyBins=true);
applyCustomColormap(ax(iAx), [0, maxC], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.5, h0=0.33);
title(ax(iAx), 'observed')

iAx = 5;
[n, xEdges, yEdges] = histcounts2(ncc.dip, ncc.rise, edges, edges);
n = n./nUnits;
n = log((n + 1/nUnits)./(nBoot + 1/nUnits));
imagesc(ax(iAx), n);
applyCustomColormap(ax(iAx), [-2, 2], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.5, h0=0.33);
title(ax(iAx), 'log(observed/bootstrap)')

axis(ax(3:5), 'equal')
xlabel(ax(3:5), 'dip')
ylabel(ax(3:5), 'rise')
xticks(ax(3:5), 0:max(xEdges))
xtickangle(ax(3:5), 0)
if p.mi.mergeDirections
    yticks(ax(3:5), 0:max(xEdges))
else
    yticks(ax(3:5), 0:2:6)
end

% Colorbar right
iAx = 6;
axc3 = ax(iAx);
axc3.Layout.Tile = layout.child(4).child(3).cw(5)-1;
axc3.Layout.TileSpan = [1, 4];
patch(axc3, [0, 1, 1, 0], [0, 0, maxC, maxC], [0, 0, maxC, maxC], EdgeColor='k');
applyCustomColormap(axc3, [0, maxC], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.5, h0=0.33);
xlim(axc3, [0, 1])
ylim(axc3, [0, maxC])
xticks(axc3, [])
yticks(axc3, [0, maxC])
yticklabels(axc3, ["0", sprintf("%g", maxC*100)]);
text(axc3, 0.5, maxC/2, '% units', Rotation=90, FontSize=p.fontSize-1, VerticalAlignment='middle', HorizontalAlignment='center')
axc3.YAxisLocation = 'right';
