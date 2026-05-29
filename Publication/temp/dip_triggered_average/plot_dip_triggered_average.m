%% Set root
% ROOTPATH = "E:\DATA";
ROOTPATH = 'C:\SERVER';

%% Clear temp vars
clearvars -except xta p kinematics ROOTPATH

%% Load data
% load(fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot\LickVsReach_DLC_dta_rta_25_25_200to800ms_units1to1443_100boots.mat"));
load(fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot\LickVsReach_DLC_dta_rta_25_100_200to800ms_units1to1443_100boots.mat"));
load(fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot\LickVsReach_DLC_miBoot_1443units_1000boots.mat"));

%% Combine movement indices (mi) across dips from all units, then cluster them
features = ["spikerate", "Jaw", "Tongue", "HandL", "HandR", "Spine"];
featureDispName = ["spike rate", "jaw", "tongue", "left paw", "right paw", "spine"];
featureUnits = ["spike rate (a.u.)", "DV pos (a.u.)", "protrusion prob", "AP pos (a.u.)", "AP pos (a.u.)", "DV pos (a.u.)"];
featureSign = [1, -1, 1, 1, 1, -1];
featureAxisDir = ["normal", "reverse", "normal", "normal", "normal", "normal"];
ylims = {5, 3, 1, 3, 3, 3};
ylims = cellfun(@(y) y*[-2/3, 1], ylims, UniformOutput=false);
ylims{featureAxisDir=="reverse"} = [-1, 2/3]*3;
p.mi.features = ["Jaw", "Tongue", "HandL", "HandR", "Spine"];
p.mi.windowPre = [-1, -0.3];
p.mi.windowPost = [0, 0.6];
p.mi.nClusters = 7;
p.mi.clusterMethod = "kmeans"; % "gaussian", "kmeans"
p.mi.clusterDimensions = 4;
p.mi.clusterSeed = 42; % 2 is also good
p.mi.semanticClusterOrder = [1, 2, 3, 4, 7, 5, 6];
% p.mi.semanticClusterOrder = 1:7;
p.mi.semanticClusterLabels = ["no move", "lick start", "lick stop", "left hand retract", "left hand reach", "right hand retract", "right hand reach"];
p.mi.semanticClusterSign = [0, 1, -1, -1, 1, -1, 1]; % 0: two-tailed, 1: right, 2: left
v = matlabRelease();
if v.Date >= datetime(2026, 5, 5) % umap was added in 2026a
    p.mi.dimensionReductionMethod = "manual+umap"; % "pca", "tsne", "umap", "manual", "manual+umap"... manual: avg(4 limbs) vs. avg(tongue/jaw) vs. spine
else
    p.mi.dimensionReductionMethod = "manual"; % "pca", "tsne", "umap", "manual", "manual+umap"... manual: avg(4 limbs) vs. avg(tongue/jaw) vs. spine
end
clear v
p.mi.displayDimensions = 2;
switch p.mi.dimensionReductionMethod
    case "pca"
        p.mi.displayMethodName = "Movement Profile (PCA)";
    case "tsne"
        p.mi.displayMethodName = "Movement Profile (t-SNE)";
    case "umap"
        p.mi.displayMethodName = "Movement Profile (UMAP)";
    case "manual"
        p.mi.displayMethodName = "Movement Profile (Manual)";
        assert(p.mi.clusterDimensions == 4)
    case "manual+umap"
        p.mi.displayMethodName = "Movement Profile (Manual + UMAP)";
        assert(p.mi.clusterDimensions == 4)
    otherwise
        error("unknown method %s", p.mi.dimensionReductionMethod)
end
tTic = tic();
clear mi
for dir = ["dip", "rise"]
    mi.(dir)(length(xta.(dir))) = struct(HandR=[], HandL=[], Spine=[], Jaw=[], Tongue=[]);
end
fprintf("Calculating movement index...")
for iUnit = 1:length(xta.dip)
    for dir = ["dip", "rise"]
        for fn = p.mi.features
            X = xta.(dir)(iUnit).(fn).X; % trials x timestamps
            t = xta.(dir)(iUnit).(fn).t;
            XPre = mean(X(:, isin(t, p.mi.windowPre)), 2, 'omitnan');
            XPost = mean(X(:, isin(t, p.mi.windowPost)), 2, 'omitnan');
            mi.(dir)(iUnit).(fn) = XPost - XPre;
        end
    end
end
fprintf("(%.1fs)\n", toc(tTic));

% Cluster movements
fprintf("Concatenating...")
iFeat = 0;
nTrials.dip = cellfun(@length, {mi.dip.HandR});
nTrials.rise = cellfun(@length, {mi.rise.HandR});
X = NaN(sum(nTrials.dip) + sum(nTrials.rise), length(p.mi.features));
for fn = p.mi.features
    iFeat = iFeat + 1;
    X(:, iFeat) = vertcat(vertcat(mi.dip.(fn)), vertcat(mi.rise.(fn)));
end
if ismember(p.mi.dimensionReductionMethod, ["manual", "manual+umap"])
    Y = NaN(size(X, 1), 4);
    Y(:, 1) = mean(X(:, ismember(p.mi.features, ["HandL"])), 2, 'omitnan');
    Y(:, 2) = mean(X(:, ismember(p.mi.features, ["Tongue", "Jaw"])), 2, 'omitnan');
    Y(:, 3) = mean(X(:, ismember(p.mi.features, ["HandR"])), 2, 'omitnan');
    Y(:, 4) = mean(X(:, ismember(p.mi.features, ["Spine"])), 2, 'omitnan');
    X = Y;
    clear Y
end
X(isnan(X)) = 0;
fprintf("(%.1fs)\n", toc(tTic));
fprintf("%s...", p.mi.displayMethodName)
switch p.mi.dimensionReductionMethod
    case "pca"
        [~, pcaScoreMerge, ~, ~] = pca(X);
        dispScoreMerge = pcaScoreMerge(:, 1:p.mi.displayDimensions);
    case "tsne"
        [~, pcaScoreMerge, ~, ~] = pca(X);
        dispScoreMerge = tsne(X, NumDimensions=p.mi.displayDimensions);
    case "umap"
        [~, pcaScoreMerge, ~, ~] = pca(X);
        dispScoreMerge = umap(X, NumDimensions=p.mi.displayDimensions);
    case "manual"
        pcaScoreMerge = X;
        dispScoreMerge = X;
    case "manual+umap"
        pcaScoreMerge = X;
        dispScoreMerge = umap(X, NumDimensions=p.mi.displayDimensions);
    otherwise
        error("unknown method %s", p.mi.dimensionReductionMethod)
end
fprintf("(%.1fs)\n", toc(tTic));

fprintf("Clustering ")
switch p.mi.clusterMethod
    case "kmeans"
        fprintf("(kmeans)...")
        % Fix the seed so we don't have to manually assign semantic cluster labels repeatedly
        rng(p.mi.clusterSeed);
        idxMerge = kmeans(pcaScoreMerge(:, 1:p.mi.clusterDimensions), p.mi.nClusters);
    case "gaussian"
        fprintf("(Gaussian mixture)...")
        rng(p.mi.clusterSeed);
        gm = fitgmdist(pcaScoreMerge(:, 1:p.mi.clusterDimensions), p.mi.nClusters);
        idxMerge = cluster(gm, pcaScoreMerge(:, 1:p.mi.clusterDimensions));
        clear gm
    otherwise
        error("unknown method %s", p.mi.clusterMethod)
end
% Cluster indices should reflect cluster size
n = arrayfun(@(k) nnz(idxMerge==k), 1:max(idxMerge));
[~, clusterOrder] = sort(n, 'descend');
idxMerge = changem(idxMerge, 1:max(idxMerge), clusterOrder);
clear n clusterOrder

fprintf("(%.1fs)\n", toc(tTic));

% Reassign cluster indices to units
clear idx pcaScore dispScore
idx.dip = idxMerge(1:sum(nTrials.dip));
idx.rise = idxMerge(sum(nTrials.dip)+1:end);
pcaScore.dip = pcaScoreMerge(1:sum(nTrials.dip), :);
pcaScore.rise = pcaScoreMerge(sum(nTrials.dip)+1:end, :);
dispScore.dip = dispScoreMerge(1:sum(nTrials.dip), :);
dispScore.rise = dispScoreMerge(sum(nTrials.dip)+1:end, :);

i0.dip = 0;
i0.rise = 0;
for iUnit = 1:length(xta.dip)
    for dir = ["dip", "rise"]
        n = nTrials.(dir)(iUnit);
        mi.(dir)(iUnit).idx = idx.(dir)(i0.(dir)+1 : i0.(dir)+n);
        mi.(dir)(iUnit).pcaScore = pcaScore.(dir)(i0.(dir)+1 : i0.(dir)+n, :);
        mi.(dir)(iUnit).dispScore = dispScore.(dir)(i0.(dir)+1 : i0.(dir)+n, :);
        i0.(dir) = i0.(dir) + n;
    end
end
assert(i0.dip == sum(nTrials.dip))
assert(i0.rise == sum(nTrials.rise))



% Scatter plot of all movement profiles
close all
layout.h = [4, 1, 1];
fig = figure(Units='normalized', Position=[0.05, 0.05, 0.9, 0.9]);
tlp = tiledlayout(fig, sum(layout.h), 1);
ax = nexttile(tlp, 1, [layout.h(1), 1]);
hold(ax, 'on')
h = gobjects(2, p.mi.nClusters);
iDir = 0;
for dir = ["dip", "rise"]
    iDir = iDir + 1;
    for k = 1:p.mi.nClusters
        sel = idx.(dir)==p.mi.semanticClusterOrder(k);
        switch dir
            case "dip"
                faceColor = "none";
                style = 'o';
            case "rise"
                faceColor = [getColor(k, p.mi.nClusters, 0.7)];
                style = 'o';
        end

        if p.mi.displayDimensions==2
            h(iDir, k) = scatter(ax, dispScore.(dir)(sel, 1), dispScore.(dir)(sel, 2), 1, Marker=style, MarkerEdgeColor=[getColor(k, p.mi.nClusters, 0.7)], MarkerFaceColor=faceColor, DisplayName=sprintf('Clu%i (n=%i %ss)', k, nnz(sel), dir));
        else
            h(iDir, k) = scatter3(ax, dispScore.(dir)(sel, 1), dispScore.(dir)(sel, 2), dispScore.(dir)(sel, 3), 1, Marker=style, MarkerEdgeColor=[getColor(k, p.mi.nClusters, 0.7)], MarkerFaceColor=faceColor, DisplayName=sprintf('Clu%i (n=%i %ss)', k, nnz(sel), dir));
        end
        clear faceColor style
    end
end
% axis(ax, 'equal')
xl = quantile(dispScoreMerge(:, 1), [0, 1]);
yl = quantile(dispScoreMerge(:, 2), [0, 1]);
xlim(ax, xl)
ylim(ax, yl)
plot3(ax, xl, [0, 0], [0, 0], 'k-')
plot3(ax, [0, 0], yl, [0, 0], 'k-')
if p.mi.displayDimensions == 3
    zl = quantile(dispScoreMerge(:, 3), [0.01, 0.99]);
    plot3(ax, [0, 0], [0, 0], zl, 'k-')
    zlim(ax, zl)
end
switch p.mi.dimensionReductionMethod
    case "manual"
        xlabel(ax, 'HandL')
        ylabel(ax, 'Tongue/Jaw')
        if p.mi.displayDimensions == 3
            zlabel(ax, 'HandR')
        end
    otherwise
        xlabel(ax, 'PC1')
        ylabel(ax, 'PC2')
        if p.mi.displayDimensions == 3
            zlabel(ax, 'PC3')
        end
end
hold(ax, 'off')
legend(h, Location='northeast')
title(ax, p.mi.displayMethodName)

% Average movement profiles by cluster
clear mp0 mp
mp0(length(xta.dip), 1) = struct(spikerate=[], HandR=[], HandL=[], Spine=[], Jaw=[], Tongue=[]);
mp = struct(dip=mp0, rise=mp0);
for iUnit = 1:length(xta.dip)
    for dir = ["dip", "rise"]
        for k = 1:p.mi.nClusters
            sel = mi.(dir)(iUnit).idx==k;
            for fn = features
                mp.(dir)(iUnit, k).(fn) = xta.(dir)(iUnit).(fn).X(sel, :);
            end
        end
    end
end
clear mp0 iUnit dir k sel fn

mpMean = struct(dip=[], rise=[]);
for dir = ["dip", "rise"]
    for fn = features
        XCell = arrayfun(@(mp) mp.(fn), mp.(dir), UniformOutput=false);
        t = xta.(dir)(1).(fn).t;
        X = NaN(p.mi.nClusters, length(t));
        for k = 1:p.mi.nClusters
            X(k, :) = mean(cat(1, XCell{:, k}), 1, 'omitnan');
        end
        mpMean.(dir).(fn).X = X;
        mpMean.(dir).(fn).t = t;
        mpMean.(dir).(fn).N = sum(cellfun(@(x) size(x, 1), XCell), 1);
    end
end
clear dir fn XCell X k t

% Plot grand average movement trajectories by cluster
lineStyles = ["-", "-", "--", "--", "-", "--", "-", "--"];
tl = gobjects(2, 1);
for iDir = 1:2
    tl(iDir) = tiledlayout(tlp, 1, length(features));
    tl(iDir).Layout.Tile = sum(layout.h(1:iDir)) + 1;
    tl(iDir).Layout.TileSpan = [layout.h(iDir+1), 1];
end
iClu = 0;
h = gobjects(2, length(features), p.mi.nClusters);
for iFeat = 1:length(features)
    fn = features(iFeat);
    iClu = iClu + 1;
    iDir = 0;
    for dir = ["dip", "rise"]
        iDir = iDir + 1;
        ax = nexttile(tl(iDir));
        hold(ax, 'on')
        for k = 1:p.mi.nClusters
            k0 = p.mi.semanticClusterOrder(k);
            h(iDir, iClu, k) = plot(ax, mpMean.(dir).(fn).t, featureSign(iFeat)*mpMean.(dir).(fn).X(k0, :), Color=[getColor(k, p.mi.nClusters, 0.7)], LineStyle=lineStyles(k), LineWidth=1.5, DisplayName=sprintf('Clu%i (%s, n=%i)', k, p.mi.semanticClusterLabels(k), mpMean.(dir).(fn).N(k0)));
        end
        hold(ax, 'off')
        title(ax, featureDispName(iFeat));
        xlabel(ax, 'time (ms)')
        ylabel(ax, featureUnits(iClu))
        yline(ax, 0, Color=[0.15, 0.15, 0.15, 0.5], LineStyle=':')
        xline(ax, 0, Color=[0.15, 0.15, 0.15, 0.5], LineStyle=':')
        ylim(ax, ylims{iClu})
        ax.YAxis.Direction = featureAxisDir(iFeat);
    end
end
legend(h(1, length(features), :), Location='eastoutside')
legend(h(2, length(features), :), Location='eastoutside')

fontsize(fig, 9, 'points')


p.mi.minNumTrialsPerCluster = 5;
clear clusterSize
clusterSize(length(xta.dip)) = struct(dip=[], rise=[]);
for iUnit = 1:length(xta.dip)
    for dir = ["dip", "rise"]
        idx = mi.(dir)(iUnit).idx;
        clusterSize(iUnit).(dir).n = zeros(1, p.mi.nClusters);
        clusterSize(iUnit).(dir).nRaw = zeros(1, p.mi.nClusters);
        clusterSize(iUnit).(dir).nTotal = length(idx);
        clusterSize(iUnit).(dir).threshold = p.mi.minNumTrialsPerCluster;
        for k = 1:p.mi.nClusters
            n = nnz(idx==k);
            clusterSize(iUnit).(dir).nRaw(k) = n;
            if n < clusterSize(iUnit).(dir).threshold
                clusterSize(iUnit).(dir).n(k) = NaN;
            else
                clusterSize(iUnit).(dir).n(k) = n;
            end
        end
    end
end
clear iUnit dir idx k n


clear tTic iUnit dir fn X t XPre XPost
clear iFeat fn pcaScoreMerge pcaExplained nClusters ax k sel
clear i0 iUnit dir n h
clear fig ax tl tlp layout iClu fn k faceColor dir iDir xl yl zl
clear dispScore dispScoreMerge featureAxisDir featureDispName features featureSign featureUnits idxMerge k0 nTrials pcaScore ylims
%% Boot
% boot_dta_clustered_movement_index;

%% Count number of "clean clusters" by unit
% A clean cluster is a cluster of dips/rises where one bodypart moved but
% nothing else (e.g. for right-hand-reach cluster, tongue/left-hand/spine
% must be stationary)

clear cc
alpha = 0.01;
semanticLabels = ["no move", "lick start", "lick stop", "left hand retract", "left hand reach", "right hand retract", "right hand reach"];
assert(isequal(semanticLabels, p.mi.semanticClusterLabels), "these must match or we could not use semanticClusterOrder below.")
mustMove = {"", "Jaw", "Jaw", "HandL", "HandL", "HandR", "HandR"};
mustNotMove = {["Jaw", "HandL", "HandR", "Spine"], ["HandL", "HandR", "Spine"], ["HandL", "HandR", "Spine"], ["Jaw", "HandR", "Spine"], ["Jaw", "HandR", "Spine"], ["Jaw", "HandL", "Spine"], ["Jaw", "HandL", "Spine"]};
% k: cluster displayOrder
% k0: real cluster id
nUnits = length(xta.dip);
for i = 1:length(semanticLabels)
    cc(i) = struct(idx=p.mi.semanticClusterOrder(i), label=semanticLabels(i), expectedSign=p.mi.semanticClusterSign(i), mustMove=mustMove{i}, mustNotMove=mustNotMove{i}, n=struct(dip=NaN(nUnits, 1), rise=NaN(nUnits, 1)));
end
clear semanticLabels mustMove mustNotMove i
fprintf('start\n')
% For each unit, count number of dips/rises that belong to each clean cluster
isMissingData = struct(dip=false(nUnits, length(cc)), rise=false(nUnits, length(cc)));
for i = 1:length(cc)
    k = cc(i).idx;
    iMustMove = find(ismember(p.mi.boot.features, cc(i).mustMove)); % this is one index, sometimes empty
    iMustNotMove = find(ismember(p.mi.boot.features, cc(i).mustNotMove)); % this is often an array
    assert(length(iMustMove) <= 1)
    assert(length(iMustNotMove) >= 1)
    nAborts = 0;
    for iUnit = 1:nUnits
        abortThisDir = false;
        for dir = ["dip", "rise"]
            if abortThisDir
                nAborts = nAborts + 1;
                break
            end
            n = clusterSize(iUnit).(dir).n(i);
            if isnan(n) || n == 0
                continue
            end
            % For mustMove: do one-tailed test depending on expectedSign.
            if isempty(iMustMove)
                mustMoveSatisfied = true;
            else
                iFeat = iMustMove;
                xObs = miObs(iUnit).(dir)(k, iFeat);
                xBoot = miBoot(iUnit).(dir)(:, k, iFeat);
                if sum(isnan(xBoot))/length(xBoot)>0.5 || isnan(xObs)
                    warning('\tcc(%i)-"%s", mustMove: Unit %i had all NaNs for feature %i-"%s"', i, cc(i).label, iUnit, iFeat, p.mi.boot.features(iFeat))
                    abortThisDir = true;
                    cc(i).n.(dir)(iUnit) = NaN;
                    isMissingData.(dir)(iUnit, i) = true;
                    continue
                end
                switch cc(i).expectedSign
                    case 1 % "right-sided test, expect increase"
                        pVal = sum(xBoot > xObs) / p.mi.boot.nBoot;
                    case -1 % "left-sided test, expect decrease"
                        pVal = sum(xBoot < xObs) / p.mi.boot.nBoot;
                    case 0
                        error("Never should have come here!")
                end
                mustMoveSatisfied = pVal < alpha;
            end

            % For mustNotMove: do two-tailed test (any is bad)
            mustNotMoveSatisfied = true;
            for iFeat = iMustNotMove(:)'
                xObs = miObs(iUnit).(dir)(k, iFeat);
                xBoot = miBoot(iUnit).(dir)(:, k, iFeat);
                if sum(isnan(xBoot))/length(xBoot)>0.5 || isnan(xObs)
                    warning('\tcc(%i)-"%s", mustNotMove: Unit %i had all NaNs for feature %i-"%s"', i, cc(i).label, iUnit, iFeat, p.mi.boot.features(iFeat))
                    abortThisDir = true;
                    cc(i).n.(dir)(iUnit) = NaN;
                    isMissingData.(dir)(iUnit, i) = true;
                    break
                end
                switch cc(i).expectedSign
                    case 1 % "right-sided test, expect increase"
                        pVal = sum(xBoot > xObs) / p.mi.boot.nBoot;
                        moveDetected = pVal < alpha;
                    case -1 % "left-sided test, expect decrease"
                        pVal = sum(xBoot < xObs) / p.mi.boot.nBoot;
                        moveDetected = pVal < alpha;
                    case 0
                        pValRight = sum(xBoot > xObs) / p.mi.boot.nBoot;
                        pValLeft = sum(xBoot < xObs) / p.mi.boot.nBoot;
                        moveDetected = pValRight < alpha/2 || pValLeft < alpha/2;
                end
                if moveDetected
                    mustNotMoveSatisfied = false;
                    break
                end
            end

            if mustMoveSatisfied && mustNotMoveSatisfied
                cc(i).n.(dir)(iUnit) = n;
            else
                cc(i).n.(dir)(iUnit) = 0;
            end
        end
    end

    if nAborts > 0
        warning("Cluster %i has %i unitsxdirs aborted", i, nAborts)
    end
end

% Units with missing data for certain clusters can be skipped over, optionally
% Because we do not know what happens.
isGoodUnit = struct(dip=[], rise=[], both=[]);
for dir = ["dip", "rise"]
    isGoodUnit.(dir) = all(~isMissingData.(dir), 2);
end
isGoodUnit.both = isGoodUnit.dip & isGoodUnit.rise;
isMissingData.both = isMissingData.dip | isMissingData.rise; % missing either rise or dip

close all
fig = figure;
tl = tiledlayout(fig, 1, 3);
for dir = ["dip", "rise", "both"]
    ax = nexttile(tl); 
    [~, I] = sort(isGoodUnit.(dir));
    imagesc(ax, ~isMissingData.(dir)(I, :))
    colormap(ax, 'gray')
    xticks(ax, 1:7)
    xticklabels(ax, [cc.label])
    title(ax, dir)
end
title(tl, "black bars are missing data, unit x cluster")
clear fig tl dir ax I

clear i k iMustMove iMustNotMove iUnit dir n mustMoveSatisfied mustNotMoveSatisfied iFeat xObs xBoot pVal pValLeft pValRight moveDetected
clear nBoot nAborts abort abortThisDir

% Count units
% To say: this unit moved one body part and nothing else, for that
% cluster, we must not miss any data for any bodypart.

% But: a unit need not contain all syllables! We're probably conservative
% about the number of units with >=2 clean clusters.

% For each unit, we ask:
% 1) how many clean clusters do you contain (except cluster 1)?
% 2) bodypart specificity: does it have clean clusters belonging to more than one of the following categories?
%   - ["lick start", "lick stop"]
%   - ["left hand reach", "left hand retract"]
%   - ["right hand reach", "right hand retract"]
% 3) start vs. stop specificity: does it have clean clusters belonging to more than one of the following categories?
%   - ["lick start", "left hand reach", "right hand reach"]
%   - ["lick stop", "left hand retract", "right hand retract"]
% 4) locomotion vs. consumption specificity: does it have clean clusters belonging to more than one of the following categories?
%   - ["lick start", "left hand retract", "right hand retract"]
%   - ["lick stop", "left hand reach", "right hand reach"]

testNames = ["clean clusters", "bodypart+movement specificity", "lick vs. reach vs. retract", "bodypart specificity", "start vs. stop specificity", "locomotion vs. consumption specificity"];
testGroups = { ... 
    {"no move", "lick start", "lick stop", "left hand reach", "left hand retract", "right hand reach", "right hand retract"}, ...
    {"lick start", "lick stop", "left hand reach", "left hand retract", "right hand reach", "right hand retract"}, ...
    {["lick start", "lick stop"], "left hand reach", "left hand retract", "right hand reach", "right hand retract"}, ...
    {["lick start", "lick stop"], ["left hand reach", "left hand retract"], ["right hand reach", "right hand retract"]}, ...
    {["lick start", "left hand reach", "right hand reach"], ["lick stop", "left hand retract", "right hand retract"]}, ...
    {["lick start", "left hand retract", "right hand retract"], ["lick stop", "left hand reach", "right hand reach"]}, ...
    };
clear tests
tests(length(testGroups)) = struct(labels=[]);
for iTest = 1:length(tests)
    tests(iTest).name = testNames(iTest);
    tests(iTest).labels = testGroups{iTest};
    for dir = ["dip", "rise"]
        for iGrp = 1:length(tests(iTest).labels)
            selCC = ismember([cc.label], tests(iTest).labels{iGrp});
            hasData = any(~isMissingData.(dir)(:, selCC), 2);
            found = arrayfun(@(c) c.n.(dir)>=p.mi.minNumTrialsPerCluster, cc(selCC), UniformOutput=false); % if cc.n.(dir) is nan (missing video data for unit-cluster), will return false, so that works out
            found = any(cat(2, found{:}), 2); % nUnits x nCategories -> nUnitsx1, any categories
            tests(iTest).(dir).found(:, iGrp) = found;
            tests(iTest).(dir).hasData(:, iGrp) = hasData;
        end
        tests(iTest).(dir).nCleanClustersFound = sum(tests(iTest).(dir).found, 2);
        tests(iTest).(dir).nCleanClustersHasData = sum(tests(iTest).(dir).hasData, 2);
        tests(iTest).(dir).prcCleanClustersFound = tests(iTest).(dir).nCleanClustersFound ./ tests(iTest).(dir).nCleanClustersHasData;
        for n = 1:length(tests(iTest).labels)
            tests(iTest).(dir).nUnitsWithNCleanClusters(n) = sum(tests(iTest).(dir).nCleanClustersFound==n);
            tests(iTest).(dir).nUnitsWithNNonMissingClusters(n) = sum(tests(iTest).(dir).nCleanClustersHasData==n);
            tests(iTest).(dir).prcUnitsWithNCleanClusters(n) = tests(iTest).(dir).nUnitsWithNCleanClusters(n) ./ nUnits;
            tests(iTest).(dir).prcUnitsWithNCleanClustersNonMissing(n) = tests(iTest).(dir).nUnitsWithNCleanClusters(n) ./ tests(iTest).(dir).nUnitsWithNNonMissingClusters(n);
        end
        for n = 1:length(tests(iTest).labels)
            tests(iTest).(dir).nUnitsWithNPlusCleanClusters(n) = sum(tests(iTest).(dir).nCleanClustersFound>=n); % Of all units, how many has at least 2 clean clusters?
            tests(iTest).(dir).nUnitsWithNPlusNonMissingClusters(n) = sum(tests(iTest).(dir).nCleanClustersHasData>=n); % Of all units, how many has at least 2 clusters not missing data?
            tests(iTest).(dir).prcUnitsWithNPlusCleanClusters(n) = tests(iTest).(dir).nUnitsWithNPlusCleanClusters(n) ./ nUnits;
            tests(iTest).(dir).prcUnitsWithNPlusCleanClustersNonMissing(n) = tests(iTest).(dir).nUnitsWithNPlusCleanClusters(n) ./ tests(iTest).(dir).nUnitsWithNPlusNonMissingClusters(n);
        end
    end
end
clear iTest dir iGrp testNames testGroups selCC hasData found

clc
for iTest = 1:length(tests)
    fprintf("Test %i (%s):\n%s:\n", ...
        iTest, tests(iTest).name, ...
        strjoin(cellfun(@(labels) sprintf("[%s]", strjoin(labels, ", ")), tests(iTest).labels), " vs. ") ...
        )
    for n = 1
        for dir = ["dip", "rise"]
            fprintf("\t%s\t%2i%% (%4i/%i total units)\t%2i%% (%4i/%i valid units)\thas exactly\t%i clean cluster(s).\n", ...
                dir, ...
                round(100*tests(iTest).(dir).prcUnitsWithNPlusCleanClusters(n)), ...
                tests(iTest).(dir).nUnitsWithNPlusCleanClusters(n), nUnits, ...
                round(100*tests(iTest).(dir).prcUnitsWithNPlusCleanClustersNonMissing(n)), ...
                tests(iTest).(dir).nUnitsWithNPlusCleanClusters(n), tests(iTest).(dir).nUnitsWithNPlusNonMissingClusters(n), ...
                n);
        end
    end
    for n = 2:length(tests(iTest).labels)
        for dir = ["dip", "rise"]
            fprintf("\t%s\t%2i%% (%4i/%i total units)\t%2i%% (%4i/%i valid units)\thas\t%i+ clean clusters.\n", ...
                dir, ...
                round(100*tests(iTest).(dir).prcUnitsWithNPlusCleanClusters(n)), ...
                tests(iTest).(dir).nUnitsWithNPlusCleanClusters(n), nUnits, ...
                round(100*tests(iTest).(dir).prcUnitsWithNPlusCleanClustersNonMissing(n)), ...
                tests(iTest).(dir).nUnitsWithNPlusCleanClusters(n), tests(iTest).(dir).nUnitsWithNPlusNonMissingClusters(n), ...
                n);
        end
    end
end
clear iTest n dir

%% Plot individual units (each axis is a bodypart, traces for clusters)
close all

features = ["spikerate", "Jaw", "HandL", "HandR", "Spine"];
featureDispName = ["spike rate", "jaw", "left paw", "right paw", "spine"];
featureUnits = ["spike rate (a.u.)", "DV pos (a.u.)", "AP pos (a.u.)", "AP pos (a.u.)", "DV pos (a.u.)"];
featureSign = [1, -1, 1, 1, -1];
featureAxisDir = ["normal", "reverse", "normal", "normal", "normal"];
% p.mi.semanticClusterLabels = ["no move", "lick start", "lick stop", "left hand retract", "left hand reach", "right hand retract", "right hand reach"];
% semanticClusterSign = [0, 1, -1, -1, 1, -1, 1]; % 0: two-tailed, 1: right, 2: left
semanticClusterSign = [0, 0, 0, 0, 0, 0, 0]; % 0: two-tailed, 1: right, 2: left
ylims = {5, 3, 3, 3, 3};
ylims = cellfun(@(y) y*[-2/3, 1], ylims, UniformOutput=false);
ylims{featureAxisDir=="reverse"} = [-1, 2/3]*3;
lineStyles = ["-", "-", "--", "--", "-", "--", "-", "--"];

nAx = length(features);
dirs = ["dip", "rise"];
fig = figure(Units='inches', InnerPosition=[2, 2, 2*(2+length(features)), 7]);
tlp = tiledlayout(fig, 2, 1, TileSpacing='compact', Padding='compact');
tl = gobjects(2, 1);
tl(1) = tiledlayout(tlp, 1, nAx, TileSpacing='compact', Padding='compact');
tl(2) = tiledlayout(tlp, 1, nAx, TileSpacing='compact', Padding='compact');
tl(1).Layout.Tile = 1;
tl(2).Layout.Tile = 2;

ax = gobjects(2, nAx);
for iDir = 1:2
    for iClu = 1:nAx
        ax(iDir, iClu) = nexttile(tl(iDir));
    end
end

iTest = 3;
for nCleanClusters = [5, 4, 3, 2]
    exportPath = fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot\Figures", sprintf("LickVsReach_DLC_dta_rta_%i_%i_%ito%ims_std", 100*p.xta.dip.thresholdQuantile, 100*p.xta.dip.thresholdSubQuantile, 100*p.xta.dip.samples(1), 100*p.xta.rise.samples(2)), sprintf("Test%i - %s", iTest, tests(iTest).name), sprintf('%i clean clusters', nCleanClusters));
    if ~exist(exportPath, 'dir')
        mkdir(exportPath)
    end
    selUnits = reshape(find(tests(iTest).dip.nCleanClustersFound == nCleanClusters), 1, []);
    
    for iUnit = selUnits
        for iDir = 1:2
            for iClu = 1:nAx
                cla(ax(iDir, iClu))
            end
        end
        for iDir = 1:2
            dir = dirs(iDir);
            % Check existence
            if isempty(xta.(dir)(iUnit).t0) && xta.(dir)(iUnit).iExp > 0
                for iClu = 1:length(features)
                    cla(ax(iDir, iClu))
                    ax(iDir, iClu).Visible = false;
                end
                continue
            end
    
            idx = mi.(dir)(iUnit).idx;
            
            for iClu = 1:length(features)
                fn = features(iClu);
                if isempty(xta.(dir)(iUnit).(fn))
                    ax(iDir, iClu).Visible = false;
                    continue
                end
    
                % if p.nBoot > 0 && ismember(fn, p.std.features) && isfield(xta.(dir)(iUnit).(fn), 'stats')
                %     prcSTD = quantile(xta.(dir)(iUnit).(fn).stats.stdBoot, [0.95, 0.99, 0.999]);
                %     nStarsSTD = sum(xta.(dir)(iUnit).(fn).stats.std > prcSTD);
                % else
                %     nStarsSTD = 0;
                % end
    
                ax(iDir, iClu).Visible = true;
    
                hold(ax(iDir, iClu), 'on')
                t = 1e3*xta.(dir)(iUnit).(fn).t;
                plot(ax(iDir, iClu), t, featureSign(iClu)*mean(xta.(dir)(iUnit).(fn).X, 1, 'omitnan'), Color=[0.15, 0.15, 0.15, 1], LineWidth=1.5, LineStyle=':');
                h = gobjects(p.mi.nClusters, 1);
                isClusterEmpty = false(1, p.mi.nClusters);
                for k = 1:p.mi.nClusters % k: semantic cluster index
                    k0 = p.mi.semanticClusterOrder(k); % Raw cluster index
                    sel = idx==k0;
                    if nnz(sel) < p.mi.minNumTrialsPerCluster
                        isClusterEmpty(k) = true;
                        continue
                    end
                    c = getColor(k, p.mi.nClusters, 0.7);
                    X = xta.(dir)(iUnit).(fn).X(sel, :);
                    mu = featureSign(iClu)*mean(X, 1, 'omitnan');
                    
                    % Check significance of cluster movement index against boostrap
                    if iClu > 1 % Skip spikerate
                        iFeat = iClu - 1;
                        xObs = miObs(iUnit).(dir)(k0, iFeat); % nClusters x nFeatures
                        xBoot = miBoot(iUnit).(dir)(:, k0, iFeat); % nBoot x nClusters x nFeatures
                        switch semanticClusterSign(k)
                            case 1
                                pVal = nnz(xBoot > xObs) ./ length(xBoot);
                                nStars = sum(xObs > quantile(xBoot, 1 - [0.05, 0.01, 0.001]));
                            case -1
                                pVal = nnz(xBoot < xObs) ./ length(xBoot);
                                nStars = sum(xObs < quantile(xBoot, [0.05, 0.01, 0.001]));
                            case 0
                                pVal = nnz(xBoot > xObs) ./ length(xBoot);
                                if pVal < 0.5 % obs on right tail
                                    nStars = sum(xObs > quantile(xBoot, 1 - 0.5*[0.05, 0.01, 0.001]));
                                else % obs on left tail
                                    pVal = 1 - pVal;
                                    nStars = sum(xObs < quantile(xBoot, 0.5*[0.05, 0.01, 0.001]));
                                end
                                pVal = pVal * 2;
                        end
                        dispName = sprintf("%s%s (n=%i, p<%g)", repmat('*', [1, nStars]), p.mi.semanticClusterLabels(k), nnz(sel), pVal);
                        clear iFeat xObs xBoot pVal nStars
                    else
                        dispName = sprintf("%s (n=%i)", p.mi.semanticClusterLabels(k), nnz(sel));
                    end
    
                    h(k) = plot(ax(iDir, iClu), t, mu, Color=c, LineStyle=lineStyles(k), LineWidth=1.5, DisplayName=dispName);
                    clear dispName
                end
                if p.nBoot > 0 && isfield(xta.(dir)(iUnit).(fn), 'XBoot')
                    prc = quantile(xta.(dir)(iUnit).(fn).XBoot, [p.bootAlpha/2, 1-p.bootAlpha/2], 1);
                    patch(ax(iDir, iClu), [t, flip(t)], featureSign(iClu)*[prc(1, :), flip(prc(2, :))], [0.15, 0.15, 0.15], FaceAlpha=0.05, EdgeColor=[0.15, 0.15, 0.15], EdgeAlpha=0.5);
                end
                xline(ax(iDir, iClu), 1e3*p.std.window, 'k--', Alpha=0.1)
                xline(ax(iDir, iClu), 0, 'k-', Alpha=0.1)
                xticks(ax(iDir, iClu), [-300, 0, 600])
                xtickangle(ax(iDir, iClu), 0)
    
                xlabel(ax(iDir, iClu), 'time (ms)')
                ylabel(ax(iDir, iClu), featureUnits(iClu))
    
                % fnDisp = sprintf("%s %s", featureDispName(iAx), repmat('*', [1, nStarsSTD]));
                fnDisp = sprintf("%s", featureDispName(iClu));
                title(ax(iDir, iClu), fnDisp, Interpreter='none')
                ylim(ax(iDir, iClu), ylims{iClu})
                hold(ax(iDir, iClu), 'off')
                ax(iDir, iClu).YAxis.Direction = featureAxisDir(iClu);
                legend(ax(iDir, iClu), h(~isClusterEmpty), Location='southoutside', AutoUpdate=false);
            end
            % lgd = legend(ax(iDir, iAx), h(~isClusterEmpty), AutoUpdate=false);
            % lgd.Layout.Tile = 'east';
        end
        for iDir = 1:2
            dir = dirs(iDir);
            if ~isempty(xta.(dir)(iUnit).HandR)
                if xta.(dir)(iUnit).iExp > 0
                    title(tl(iDir), sprintf("Unit %i (n=%i %ss)", iUnit, length(xta.(dir)(iUnit).t0), dir), FontWeight='bold')
                else
                    title(tl(iDir), sprintf("%i units (n=%i %ss)", length(xta)-1, length(xta.(dir)(iUnit).t0), dir), FontWeight='bold')
                end
            else
                xlabel(tl(iDir), '')
                title(tl(iDir), '')
            end
        end
        xlim(ax(:, 1:end-2), 1e3*[-0.5, 1])
        fontsize(fig, 9, 'points')
        print(fig, fullfile(exportPath, sprintf("dta_rta_unit_%03i", iUnit)), '-dpng', '-r0')
    end
end
clear exportPath features nAx dirs dimensionReduction fig tlp tl ax iDir iClu iUnit dir
clear idx nClusters mr i fn t selT xx pcScore explained score eva k sel
clear mdm i fn t selT stdObs hash I idxSorted sepHash prcSTD nStarsSTD t X k c mu err prc
clear iDir dir fnDisp
clear isClusterEmpty t h k k0 sel c X mu prc lgd ylims sn t0 tt s lineStyles



% %% Plot individual units (each axis is a cluster, traces for bodyparts)
% close all
% 
% features = ["Jaw", "HandL", "HandR", "Spine"];
% featureDispName = ["jaw", "left paw", "right paw", "spine"];
% featureUnits = ["DV pos (a.u.)", "AP pos (a.u.)", "AP pos (a.u.)", "DV pos (a.u.)"];
% featureSign = [1, 1, 1, -1];
% featureAxisDir = ["reverse", "normal", "normal", "normal"];
% % p.mi.semanticClusterLabels = ["no move", "lick start", "lick stop", "left hand retract", "left hand reach", "right hand retract", "right hand reach"];
% % semanticClusterSign = [0, 1, -1, -1, 1, -1, 1]; % 0: two-tailed, 1: right, 2: left
% semanticClusterSign = [0, 0, 0, 0, 0, 0, 0]; % 0: two-tailed, 1: right, 2: left
% ylims = {5, 3};
% ylims = cellfun(@(y) y*[-2/3, 1], ylims, UniformOutput=false);
% lineStyles = ["-", "-", "-", "-"];
% % lineStyles = ["-", "-", "--", "--", "-", "--", "-", "--"];
% 
% nAx = length(p.mi.semanticClusterOrder); % spikerate + nClusters
% dirs = ["dip", "rise"];
% fig = figure(Units='inches', InnerPosition=[2, 2, 2*(nAx+1), 7]);
% tlp = tiledlayout(fig, 1, nAx+1, TileSpacing='compact', Padding='compact', TileIndexing='rowmajor');
% tll = tiledlayout(tlp, 2, 1, TileSpacing='compact', Padding='compact');
% tll.Layout.Tile = 1;
% tll.Layout.TileSpan = [1, 1];
% tlr = tiledlayout(tlp, 2, 1, TileSpacing='compact', Padding='compact');
% tlr.Layout.Tile = 2;
% tlr.Layout.TileSpan = [1, nAx];
% tl = gobjects(2, 1);
% tl(1) = tiledlayout(tlr, 1, nAx, TileSpacing='compact', Padding='compact');
% tl(2) = tiledlayout(tlr, 1, nAx, TileSpacing='compact', Padding='compact');
% tl(1).Layout.Tile = 1;
% tl(2).Layout.Tile = 2;
% xlabel(tlr, 'Clusters', FontSize=9, FontWeight='bold')
% 
% ax = gobjects(2, nAx);
% for iDir = 1:2
%     for iAx = 1:nAx
%         ax(iDir, iAx) = nexttile(tl(iDir));
%     end
% end
% 
% iTest = 3;
% % for nCleanClusters = [5, 4, 3, 2]
% for nCleanClusters = 1
%         exportPath = fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot\Figures", sprintf("LickVsReach_DLC_dta_rta_%i_%i_%ito%ims_std", 100*p.xta.dip.thresholdQuantile, 100*p.xta.dip.thresholdSubQuantile, 100*p.xta.dip.samples(1), 100*p.xta.rise.samples(2)), "By Cluster", sprintf("Test%i - %s", iTest, tests(iTest).name), sprintf('%i clean clusters', nCleanClusters));
%     if ~exist(exportPath, 'dir')
%         mkdir(exportPath)
%     end
%     selUnits = reshape(find(tests(iTest).dip.nCleanClustersFound == nCleanClusters), 1, []);
% 
%     for iUnit = selUnits
%         for iDir = 1:2
%             for iAx = 1:nAx
%                 cla(ax(iDir, iAx))
%             end
%         end
%         for iDir = 1:2
%             dir = dirs(iDir);
%             % Check existence of dips/rises
%             if isempty(xta.(dir)(iUnit).t0) && xta.(dir)(iUnit).iExp > 0
%                 for iAx = 1:length(features)
%                     cla(ax(iDir, iAx))
%                     ax(iDir, iAx).Visible = false;
%                 end
%                 continue
%             end
% 
%             idx = mi.(dir)(iUnit).idx;
% 
%             for iAx = 1:length(p.mi.semanticClusterOrder)
%                 iClu = p.mi.semanticClusterOrder(iAx);
%                 nTrials = clusterSize(iUnit).(dir).n(iClu);
%                 if isnan(nTrials) || nTrials < p.mi.minNumTrialsPerCluster;
%                     ax(iDir, iAx).Visible = false;
%                     continue
%                 end
% 
%                 ax(iDir, iAx).Visible = true;
% 
%                 hold(ax(iDir, iAx), 'on')
%                 t = 1e3*xta.(dir)(iUnit).HandR.t;
%                 h = gobjects(length(features), 1);
%                 for iFeat = 1:length(features) % k: semantic cluster index
%                     fn = features(iFeat);
%                     sel = idx==iClu;
%                     c = getColor(iFeat, length(features), 0.7);
%                     X = xta.(dir)(iUnit).(fn).X(sel, :);
%                     mu = featureSign(iFeat)*mean(X, 1, 'omitnan');
% 
%                     % Check significance of cluster movement index against bootstrap
%                     xObs = miObs(iUnit).(dir)(iClu, iFeat); % nClusters x nFeatures
%                     xBoot = miBoot(iUnit).(dir)(:, iClu, iFeat); % nBoot x nClusters x nFeatures
%                     switch semanticClusterSign(iAx)
%                         case 1
%                             pVal = nnz(xBoot > xObs) / length(xBoot);
%                             nStars = sum(xObs > quantile(xBoot, 1 - [0.05, 0.01, 0.001]));
%                         case -1
%                             pVal = nnz(xBoot < xObs) / length(xBoot);
%                             nStars = sum(xObs < quantile(xBoot, [0.05, 0.01, 0.001]));
%                         case 0
%                             pVal = nnz(xBoot > xObs) / length(xBoot);
%                             if pVal < 0.5 % obs on right tail
%                                 nStars = sum(xObs > quantile(xBoot, 1 - 0.5*[0.05, 0.01, 0.001]));
%                             else % obs on left tail
%                                 pVal = 1 - pVal;
%                                 nStars = sum(xObs < quantile(xBoot, 0.5*[0.05, 0.01, 0.001]));
%                             end
%                             pVal = pVal * 2;
%                     end
%                     dispName = sprintf("%s%s (p<%g)", repmat('*', [1, nStars]), fn, pVal);
%                     clear xObs xBoot pVal nStars
% 
%                     h(iFeat) = plot(ax(iDir, iAx), t, mu, Color=c, LineStyle=lineStyles(iFeat), LineWidth=1.5, DisplayName=dispName);
%                     clear dispName
%                 end
%                 xline(ax(iDir, iAx), 1e3*p.std.window, 'k--', Alpha=0.1)
%                 xline(ax(iDir, iAx), 0, 'k-', Alpha=0.1)
%                 xticks(ax(iDir, iAx), [-300, 0, 600])
%                 xtickangle(ax(iDir, iAx), 0)
%                 yline(ax(iDir, iAx), 0, 'k-', Alpha=0.1)
% 
%                 xlabel(ax(iDir, iAx), 'time (ms)')
%                 ylabel(ax(iDir, iAx), 'position (a.u.)')
% 
%                 title(ax(iDir, iAx), p.mi.semanticClusterLabels(iAx), Interpreter='none')
%                 ylim(ax(iDir, iAx), ylims{2})
%                 hold(ax(iDir, iAx), 'off')
%                 % ax(iDir, iAx).YAxis.Direction = featureAxisDir(iAx);
%                 legend(ax(iDir, iAx), h, Location='southoutside', AutoUpdate=false);
%             end
%             % lgd = legend(ax(iDir, iAx), h, AutoUpdate=false);
%             % lgd.Layout.Tile = 'east';
%         end
%         for iDir = 1:2
%             dir = dirs(iDir);
%             if ~isempty(xta.(dir)(iUnit).HandR)
%                 if xta.(dir)(iUnit).iExp > 0
%                     title(tl(iDir), sprintf("Unit %i (n=%i %ss)", iUnit, length(xta.(dir)(iUnit).t0), dir), FontWeight='bold')
%                 else
%                     title(tl(iDir), sprintf("%i units (n=%i %ss)", length(xta)-1, length(xta.(dir)(iUnit).t0), dir), FontWeight='bold')
%                 end
%             else
%                 xlabel(tl(iDir), '')
%                 title(tl(iDir), '')
%             end
%         end
%         xlim(ax(:, 1:end-2), 1e3*[-0.5, 1])
%         fontsize(fig, 9, 'points')
%         print(fig, fullfile(exportPath, sprintf("dta_rta_unit_%03i", iUnit)), '-dpng', '-r0')
%     end
% end
% clear exportPath features nAx dirs dimensionReduction fig tlp tll tlr tl ax iDir iAx iUnit dir
% clear idx nClusters mr i fn t selT xx pcScore explained score eva k sel
% clear mdm i fn t selT stdObs hash I idxSorted sepHash prcSTD nStarsSTD t X k c mu err prc
% clear iDir dir fnDisp
% clear isClusterEmpty t h k iClu sel c X mu prc lgd ylims sn t0 tt s lineStyles




%% Plot individual units just-dip matrix (rows are features, columns are clusters)
close all

features = ["spikerate"; "Jaw"; "HandL"; "HandR"; "Spine"];
featureDispName = ["spike rate"; "jaw"; "left hand"; "right hand"; "spine"];
featureUnits = ["(a.u.)"; "(DV a.u.)"; "(AP a.u.)"; "(AP a.u.)"; "(DV a.u.)"];
featureSign = [1; -1; 1; 1; -1];
featureAxisDir = ["normal"; "normal"; "normal"; "normal"; "normal"];
% p.mi.semanticClusterLabels = ["no move", "lick start", "lick stop", "left hand retract", "left hand reach", "right hand retract", "right hand reach"];
% semanticClusterSign = [0, 1, -1, -1, 1, -1, 1]; % 0: two-tailed, 1: right, 2: left
semanticClusterSign = [0, 1, -1, -1, 1, -1, 0]; % 0: two-tailed, 1: right, 2: left
clusterDispOrder = [1, 2, 3, 5, 4, 7, 6];
ylims = {4; 4; 4; 4; 4};
ylims = cellfun(@(y) y*[-2/3, 1], ylims, UniformOutput=false);
% ylims{featureAxisDir=="reverse"} = [-1, 2/3]*3;
xt = 1e3*unique([p.mi.windowPost]);

nClusters = length(p.mi.semanticClusterOrder);
nFeatures = length(features);

fig = figure(Units='inches', InnerPosition=[1, 1, 1.1*nClusters, 1*length(features)]);
tl = tiledlayout(fig, nFeatures, nClusters, TileSpacing='tight', Padding='compact', TileIndexing='columnmajor');

ax = gobjects(nFeatures, nClusters);
for iClu = 1:nClusters
    for iFeat = 1:nFeatures
        ax(iFeat, iClu) = nexttile(tl);
        ax(iFeat, iClu).Box = 'on';
    end
end

iTest = 3;
for nCleanClusters = [5, 4, 3, 2, 1]
    exportPath = fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot\Figures", sprintf("LickVsReach_DLC_dta_rta_%i_%i_%ito%ims_std", 100*p.xta.dip.thresholdQuantile, 100*p.xta.dip.thresholdSubQuantile, 100*p.xta.dip.samples(1), 100*p.xta.rise.samples(2)), "By Feature x Cluster", sprintf("Test%i - %s", iTest, tests(iTest).name), sprintf('%i clean clusters', nCleanClusters));
    if ~exist(exportPath, 'dir')
        mkdir(exportPath)
    end
    selUnits = reshape(find(tests(iTest).dip.nCleanClustersFound == nCleanClusters), 1, []);
    
    for iUnit = selUnits
        for dir = ["dip", "rise"]
            for iClu = 1:nClusters
                for iFeat = 1:nFeatures
                    cla(ax(iFeat, iClu))
                end
            end
    
            % Check existence of dips/rises
            if isempty(xta.(dir)(iUnit).t0)
                continue
            end
    
            idx = mi.(dir)(iUnit).idx;
            
            for iClu = 1:length(p.mi.semanticClusterOrder)
                k = clusterDispOrder(iClu);
                k0 = p.mi.semanticClusterOrder(k); % raw cluster index
                nTrials = clusterSize(iUnit).(dir).n(k0);
                if isnan(nTrials) || nTrials < p.mi.minNumTrialsPerCluster
                    set(ax(:, iClu), Visible=false);
                    continue
                end
    
                set(ax(:, iClu), Visible=true);
    
                for iFeat = 1:length(features) % k: semantic cluster index
                    hold(ax(iFeat, iClu), 'on')
                    fn = features(iFeat);
                    sel = idx==k0;
                    % c = getColor(iFeat, length(features), 0.7);
                    c = 'k';
                    X = xta.(dir)(iUnit).(fn).X(sel, :);
                    t = 1e3*xta.(dir)(iUnit).(fn).t;
                    mu = featureSign(iFeat)*mean(X, 1, 'omitnan');
                    
                    % Check significance of cluster movement index against bootstrap
                    if iFeat > 1
                        xObs = miObs(iUnit).(dir)(k0, iFeat-1); % nClusters x nFeatures
                        xBoot = miBoot(iUnit).(dir)(:, k0, iFeat-1); % nBoot x nClusters x nFeatures
                        switch semanticClusterSign(k)
                            case 1
                                pVal = nnz(xBoot > xObs) / length(xBoot);
                                nStars = sum(xObs > quantile(xBoot, 1 - [0.05, 0.01, 0.001]));
                            case -1
                                pVal = nnz(xBoot < xObs) / length(xBoot);
                                nStars = sum(xObs < quantile(xBoot, [0.05, 0.01, 0.001]));
                            case 0
                                pVal = nnz(xBoot > xObs) / length(xBoot);
                                if pVal < 0.5 % obs on right tail
                                    nStars = sum(xObs > quantile(xBoot, 1 - 0.5*[0.05, 0.01, 0.001]));
                                else % obs on left tail
                                    pVal = 1 - pVal;
                                    nStars = sum(xObs < quantile(xBoot, 0.5*[0.05, 0.01, 0.001]));
                                end
                                pVal = pVal * 2;
                        end
                        text(ax(iFeat, iClu), 1e3*mean(p.mi.windowPost), ylims{iFeat}(2), repmat('*', [1, nStars]), ...
                            HorizontalAlignment='center', VerticalAlignment='top', FontSize=12, FontWeight='bold');
                        clear xObs xBoot pVal nStars
                    end
    
                    plot(ax(iFeat, iClu), t, mu, Color=c, LineStyle='-', LineWidth=1.5);
                    xline(ax(iFeat, iClu), xt, 'k--', Alpha=0.1)
                    xline(ax(iFeat, iClu), 0, 'k-', Alpha=0.1)
                    yline(ax(iFeat, iClu), 0, 'k-', Alpha=0.1)
                    if iFeat == nFeatures
                        xticks(ax(iFeat, iClu), xt)
                        xtickangle(ax(iFeat, iClu), 0)
                    else
                        xticks(ax(iFeat, iClu), [])
                    end
                    if iClu == 1
                        yticks(ax(iFeat, iClu), [0, ylims{iFeat}(2)])
                    else
                        yticks(ax(iFeat, iClu), [])
                    end
    
                    % if iFeat == nFeatures
                    %     xlabel(ax(iFeat, iClu), 'time (ms)')
                    % end
                    if iClu == 1
                        ylabel(ax(iFeat, iClu), sprintf("%s\n%s", featureDispName(iFeat), featureUnits(iFeat)))
                    end
                    if iFeat == 1
                        title(ax(iFeat, iClu), sprintf("%s\n(n=%i %ss)", p.mi.semanticClusterLabels(k), nTrials, dir), Interpreter='none', FontWeight='normal')
                    end
                    ylim(ax(iFeat, iClu), ylims{iFeat})
                    hold(ax(iFeat, iClu), 'off')
                    ax(iFeat, iAx).YAxis.Direction = featureAxisDir(iFeat);
                end
            end
            xlim(ax, 1e3*[-0.5, 1])
            xlabel(tl, "time since spike rate change (ms)")
            title(tl, sprintf("Unit %i (n=%i %ss)", iUnit, length(xta.(dir)(iUnit).t0), dir), FontWeight='bold')
            fontsize(fig, 9, 'points')
            print(fig, fullfile(exportPath, sprintf("featxclus_unit_%03i_%s", iUnit, dir)), '-dpng', '-r0')
        end
    end
end

clear features featureDispName featureUnits featureSign featureAxisDir semanticClusterSign clusterDispOrder ylims xt nClusters nFeatures fig tl ax iClu iFeat iTest nCleanClusters exportPath selUnits iUnit dir idx k k0 nTrials fn sel c X t mu xObs xBoot pVal nStars 