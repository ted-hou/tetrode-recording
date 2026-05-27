%% Set root
% ROOTPATH = "E:\DATA";
ROOTPATH = 'C:\SERVER';

%% Clear temp vars
clearvars -except xta p kinematics ROOTPATH

%% Load data
% load(fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot\LickVsReach_DLC_dta_rta_25_25_200to800ms_units1to1443_100boots.mat"));
load(fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot\LickVsReach_DLC_dta_rta_25_100_200to800ms_units1to1443_100boots.mat"));

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
p.mi.dimensionReductionMethod = "manual+umap"; % "pca", "tsne", "umap", "manual", "manual+umap"... manual: avg(4 limbs) vs. avg(tongue/jaw) vs. spine
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

% Plot grand average movement trajectories by cluster
lineStyles = ["-", "-", "--", "--", "-", "--", "-", "--"];
tl = gobjects(2, 1);
for iDir = 1:2
    tl(iDir) = tiledlayout(tlp, 1, length(features));
    tl(iDir).Layout.Tile = sum(layout.h(1:iDir)) + 1;
    tl(iDir).Layout.TileSpan = [layout.h(iDir+1), 1];
end
iAx = 0;
h = gobjects(2, length(features), p.mi.nClusters);
for iFeat = 1:length(features)
    fn = features(iFeat);
    iAx = iAx + 1;
    iDir = 0;
    for dir = ["dip", "rise"]
        iDir = iDir + 1;
        ax = nexttile(tl(iDir));
        hold(ax, 'on')
        for k = 1:p.mi.nClusters
            k0 = p.mi.semanticClusterOrder(k);
            h(iDir, iAx, k) = plot(ax, mpMean.(dir).(fn).t, featureSign(iFeat)*mpMean.(dir).(fn).X(k0, :), Color=[getColor(k, p.mi.nClusters, 0.7)], LineStyle=lineStyles(k), LineWidth=1.5, DisplayName=sprintf('Clu%i (%s, n=%i)', k, p.mi.semanticClusterLabels(k), mpMean.(dir).(fn).N(k0)));
        end
        hold(ax, 'off')
        title(ax, featureDispName(iFeat));
        xlabel(ax, 'time (ms)')
        ylabel(ax, featureUnits(iAx))
        yline(ax, 0, Color=[0.15, 0.15, 0.15, 0.5], LineStyle=':')
        xline(ax, 0, Color=[0.15, 0.15, 0.15, 0.5], LineStyle=':')
        ylim(ax, ylims{iAx})
        ax.YAxis.Direction = featureAxisDir(iFeat);
    end
end
legend(h(1, length(features), :), Location='eastoutside')
legend(h(2, length(features), :), Location='eastoutside')

fontsize(fig, 9, 'points')

clear tTic iUnit dir fn X t XPre XPost
clear iFeat fn pcaScoreMerge pcaExplained nClusters ax k sel
clear i0 iUnit dir n h
clear fig ax tl tlp layout iAx fn k faceColor dir iDir xl yl zl

%% Count number of movement profiles by unit
p.mi.minNumTrialsPerCluster = 5;
% p.mi.minNumTrialsPerClusterQuantile = 0.05;

clear clusterSize
clusterSize(length(xta.dip)) = struct(dip=[], rise=[]);
for iUnit = 1:length(xta.dip)
    for dir = ["dip", "rise"]
        idx = mi.(dir)(iUnit).idx;
        clusterSize(iUnit).(dir).n = zeros(1, p.mi.nClusters);
        clusterSize(iUnit).(dir).nRaw = zeros(1, p.mi.nClusters);
        clusterSize(iUnit).(dir).nTotal = length(idx);
        % clusterSize(iUnit).(dir).threshold = max(p.mi.minNumTrialsPerCluster, length(idx)*p.mi.minNumTrialsPerClusterQuantile);
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

% p.mi.semanticClusterOrder = [1, 2, 3, 4, 7, 5, 6];
% p.mi.semanticClusterLabels = ["no move", "lick start", "lick stop", "left hand retract", "left hand reach", "right hand retract", "right hand reach"];

clear cns cn
for k = 1:p.mi.nClusters
    p.mi.semanticClusterOrderReversed(k) = find(p.mi.semanticClusterOrder == k);
end
cns.none = ["no move"];
cns.all = ["no move", "lick start", "lick stop", "left hand retract", "left hand reach", "right hand retract", "right hand reach"];
cns.any = ["lick start", "lick stop", "left hand retract", "left hand reach", "right hand retract", "right hand reach"];
cns.hand = ["left hand retract", "left hand reach", "right hand retract", "right hand reach"];
cns.lhand = ["left hand retract", "left hand reach"];
cns.rhand = ["right hand retract", "right hand reach"];
cns.reach = ["left hand reach", "right hand reach"];
cns.retract = ["left hand retract", "right hand retract"];
cns.lick = ["lick start", "lick stop"];

for fn = ["none", "all", "any", "hand", "lick", "lhand", "rhand", "reach", "retract"]
    cn.(fn) = ismember(p.mi.semanticClusterLabels, cns.(fn));
    cn.(fn) = p.mi.semanticClusterOrder(cn.(fn));
end

clear semanticClusterSize
for dir = ["dip", "rise"]
    for selType = ["none", "all", "any", "hand", "lick", "lhand", "rhand", "reach", "retract"]
        n = arrayfun(@(cx) cx.(dir).n(cn.(selType)), clusterSize, UniformOutput=false);
        n = cat(1, n{:});
        semanticClusterSize.(dir).(selType) = array2table(n, VariableNames=cns.(selType));
    end
end

clear nExistingProfiles
for dir = ["dip", "rise"]
    for selType = ["none", "all", "any", "hand", "lick", "lhand", "rhand", "reach", "retract"]
        nExistingProfiles.(dir).(selType) = sum(~isnan(table2array(semanticClusterSize.(dir).(selType))), 2);
    end
end
clear dir selType cn cns k

close all
fig = figure();
tl = tiledlayout(fig, 2, 7);
iDir = 0;
for dir = ["dip", "rise"]
    iDir = iDir + 1;
    iType = 0;
    for selType = ["any", "lick", "hand", "lhand", "rhand", "reach", "retract"]
        iType = iType + 1;
        ax = nexttile(tl);
        histogram(ax, nExistingProfiles.(dir).(selType))
        title(ax, sprintf("%s %s", dir, selType))
        xlabel(ax, "no. movement types")
        ylabel(ax, "no. SNr units")
    end
end
clear fig tl iDir iType dir selType ax

nTotal = length(xta.dip);
for dir = ["dip", "rise"]
    fprintf("%i SNr units:\n", nTotal)

    n = nnz(nExistingProfiles.(dir).none > 0);
    fprintf("\t%i (%.1f%%) units had at least %i %ss each where no consistent movements were observed:\n", n, 100*n/nTotal, p.mi.minNumTrialsPerCluster, dir);

    n = nnz(nExistingProfiles.(dir).any >= 2);
    fprintf("\t%i (%.1f%%) units had at least %i %ss each where 2 or more non-overlapping movements were observed:\n", n, 100*n/nTotal, p.mi.minNumTrialsPerCluster, dir);

    n = nnz(nExistingProfiles.(dir).hand >= 2);
    fprintf("\t%i (%.1f%%) units had at least %i %ss each where 2 or more non-overlapping hand movements were observed:\n", n, 100*n/nTotal, p.mi.minNumTrialsPerCluster, dir);

    n = nnz(nExistingProfiles.(dir).lick >= 2);
    fprintf("\t%i (%.1f%%) units had at least %i %ss each where 2 or more non-overlapping lick movements were observed:\n", n, 100*n/nTotal, p.mi.minNumTrialsPerCluster, dir);
end
clear nTotal dir n

%% Boot
% boot_dta_clustered_movement_index;
load(fullfile("C:\SERVER\LickVsReach_DTA_RTA_boot", sprintf("LickVsReach_DLC_miBoot_%iunits_%iboots.mat", 1443, 1000)))

%% Plot individual units
close all
exportPath = fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot\Figures", sprintf("LickVsReach_DLC_dta_rta_%i_%i_%ito%ims_std", 100*p.xta.dip.thresholdQuantile, 100*p.xta.dip.thresholdSubQuantile, 100*p.xta.dip.samples(1), 100*p.xta.rise.samples(2)));
if ~exist(exportPath, 'dir')
    mkdir(exportPath)
end
selUnits = 1:length(xta.dip);
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
    for iAx = 1:nAx
        ax(iDir, iAx) = nexttile(tl(iDir));
    end
end
for iUnit = selUnits
    for iDir = 1:2
        for iAx = 1:nAx
            cla(ax(iDir, iAx))
        end
    end
    for iDir = 1:2
        dir = dirs(iDir);
        % Check existence
        if isempty(xta.(dir)(iUnit).t0) && xta.(dir)(iUnit).iExp > 0
            for iAx = 1:length(features)
                cla(ax(iDir, iAx))
                ax(iDir, iAx).Visible = false;
            end
            continue
        end

        idx = mi.(dir)(iUnit).idx;
        
        for iAx = 1:length(features)
            fn = features(iAx);
            if isempty(xta.(dir)(iUnit).(fn))
                ax(iDir, iAx).Visible = false;
                continue
            end

            % if p.nBoot > 0 && ismember(fn, p.std.features) && isfield(xta.(dir)(iUnit).(fn), 'stats')
            %     prcSTD = quantile(xta.(dir)(iUnit).(fn).stats.stdBoot, [0.95, 0.99, 0.999]);
            %     nStarsSTD = sum(xta.(dir)(iUnit).(fn).stats.std > prcSTD);
            % else
            %     nStarsSTD = 0;
            % end

            ax(iDir, iAx).Visible = true;

            hold(ax(iDir, iAx), 'on')
            t = 1e3*xta.(dir)(iUnit).(fn).t;
            plot(ax(iDir, iAx), t, featureSign(iAx)*mean(xta.(dir)(iUnit).(fn).X, 1, 'omitnan'), Color=[0.15, 0.15, 0.15, 1], LineWidth=1.5, LineStyle=':');
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
                mu = featureSign(iAx)*mean(X, 1, 'omitnan');
                
                % Check significance of cluster movement index against boostrap
                if iAx > 1 % Skip spikerate
                    iFeat = iAx - 1;
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

                h(k) = plot(ax(iDir, iAx), t, mu, Color=c, LineStyle=lineStyles(k), LineWidth=1.5, DisplayName=dispName);
                clear dispName
            end
            if p.nBoot > 0 && isfield(xta.(dir)(iUnit).(fn), 'XBoot')
                prc = quantile(xta.(dir)(iUnit).(fn).XBoot, [p.bootAlpha/2, 1-p.bootAlpha/2], 1);
                patch(ax(iDir, iAx), [t, flip(t)], featureSign(iAx)*[prc(1, :), flip(prc(2, :))], [0.15, 0.15, 0.15], FaceAlpha=0.05, EdgeColor=[0.15, 0.15, 0.15], EdgeAlpha=0.5);
            end
            xline(ax(iDir, iAx), 1e3*p.std.window, 'k--', Alpha=0.1)
            xline(ax(iDir, iAx), 0, 'k-', Alpha=0.1)
            xticks(ax(iDir, iAx), [-300, 0, 600])
            xtickangle(ax(iDir, iAx), 0)

            xlabel(ax(iDir, iAx), 'time (ms)')
            ylabel(ax(iDir, iAx), featureUnits(iAx))

            % fnDisp = sprintf("%s %s", featureDispName(iAx), repmat('*', [1, nStarsSTD]));
            fnDisp = sprintf("%s", featureDispName(iAx));
            title(ax(iDir, iAx), fnDisp, Interpreter='none')
            ylim(ax(iDir, iAx), ylims{iAx})
            hold(ax(iDir, iAx), 'off')
            ax(iDir, iAx).YAxis.Direction = featureAxisDir(iAx);
            legend(ax(iDir, iAx), h(~isClusterEmpty), Location='southoutside', AutoUpdate=false);
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

clear exportPath features nAx dirs dimensionReduction fig tlp tl ax iDir iAx iUnit dir
clear idx nClusters mr i fn t selT xx pcScore explained score eva k sel
clear mdm i fn t selT stdObs hash I idxSorted sepHash prcSTD nStarsSTD t X k c mu err prc
clear iDir dir fnDisp
clear isClusterEmpty t h k k0 sel c X mu prc lgd ylims sn t0 tt s lineStyles

