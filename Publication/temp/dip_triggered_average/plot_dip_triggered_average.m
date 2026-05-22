%% Set root
% ROOTPATH = "E:\DATA";
ROOTPATH = 'C:\SERVER';

%% Clear temp vars
clearvars -except xta p kinematics ROOTPATH

%% Load data
load(fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot\LickVsReach_DLC_dta_rta_25_25_200to800ms_units1to1443_100boots.mat"));
load(fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot\LickVsReach_DLC_dta_rta_25_100_200to800ms_units1to1443_0boots.mat"));

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
%% Plot individual units
close all
exportPath = fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot\Figures", sprintf("LickVsReach_DLC_dta_rta_%i_%i_%ito%ims_std", 100*p.xta.dip.thresholdQuantile, 100*p.xta.dip.thresholdSubQuantile, 100*p.xta.dip.samples(1), 100*p.xta.rise.samples(2)));
if ~exist(exportPath, 'dir')
    mkdir(exportPath)
end
statFeatures = ["Jaw", "Tongue", "HandL", "HandR", "Spine"];
nAx = length(features);
dirs = ["dip", "rise"];
fig = figure(Units='inches', InnerPosition=[2, 2, 1.5*(2+length(features)), 4.5]);
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
for iUnit = 1:length(xta.dip)
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

        % Movement diversity scatter plot
        [lia, statFeatureOrder] = ismember(statFeatures, p.std.features);
        assert(all(lia), 'Some members of statFeatureOrder are not found.')

        idx = mi.(dir)(iUnit).idx;
        % if p.nBoot > 0 && xta.(dir)(iUnit).iExp > 0
            % iAx = nAx - 1;
            % hold(ax(iDir, iAx), 'on')
            % score = mi.(dir)(iUnit).dispScore;
            % for k = 1:p.mi.nClusters
            %     sel = idx==k;
            %     if nnz(sel) == 0
            %         continue
            %     end
            %     scatter(ax(iDir, iAx), score(sel, 1), score(sel, 2), 10, getColor(k, p.mi.nClusters, 0.7))
            % end
            % xlabel(ax(iDir, iAx), sprintf("PC%i", 1))
            % ylabel(ax(iDir, iAx), sprintf("PC%i", 2))
            % switch p.mi.dimensionReductionMethod   
            %     case "pca"
            %         title(ax(iDir, iAx), "PCA")
            %     case "tsne"
            %         title(ax(iDir, iAx), "t-SNE")
            %     case "umap"
            %         title(ax(iDir, iAx), "UMAP")
            % end
            % xticks(ax(iDir, iAx), [])
            % yticks(ax(iDir, iAx), [])
            % hold(ax(iDir, iAx), 'off')
        % end

        % % Movement diversity matrix
        % if p.nBoot > 0 && xta.(dir)(iUnit).iExp > 0
        %     iAx = nAx;
        %     mdm = NaN(length(xta.(dir)(iUnit).t0), length(p.std.features));
        %     for i = 1:length(p.std.features)
        %         fn = p.std.features(statFeatureOrder(i));
        %         if isempty(xta.(dir)(iUnit).(fn))
        %             continue
        %         end
        %         t = xta.(dir)(iUnit).(fn).t;
        %         selT = t >= p.std.window(1) & t <= p.std.window(2);
        %         stdObs = std(xta.(dir)(iUnit).(fn).X(:, selT), 0, 2, 'omitnan');
        %         mdm(:, i) = arrayfun(@(data) nnz(xta.(dir)(iUnit).(fn).stats.stdBoot < data) ./ length(xta.(dir)(iUnit).(fn).stats.stdBoot), stdObs, UniformOutput=true);
        %     end
        %     mdm(isnan(mdm)) = 0;
        %     hash = sum((mdm > 0.95) .* 2.^(size(mdm, 2)-1:-1:0), 2);
        %     hash = hash + (idx-1) .* 2.^(size(mdm, 2));
        %     [~, I] = sort(hash, 'ascend');
        %     idxSorted = idx(I);
        %     sepHash = arrayfun(@(idx) find(idxSorted==idx, 1, 'last'), 1:max(idx)-1);
        %     imagesc(ax(iDir, iAx), mdm(I, :))
        %     if ~isempty(sepHash)
        %         yline(ax(iDir, iAx), 0.5+sepHash, 'k--')
        %         yticks(ax(iDir, iAx), 0.5+unique([1, sepHash, length(idx)]))
        %         yticklabels(ax(iDir, iAx), string(unique([1, sepHash, length(idx)])))
        %     end
        %     applyCustomColormap(ax(iDir, iAx), [-1, 1], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.05, h0=0.33);
        %     xticks(ax(iDir, iAx), 1:length(p.std.features))
        %     xticklabels(ax(iDir, iAx), p.std.features(statFeatureOrder))
        %     colorbar(ax(iDir, iAx), Orientation='horizontal', Location='southoutside')
        %     ax(iDir, iAx).XAxisLocation = 'top';
        %     ylabel(ax(iDir, iAx), 'Trial')
        %     clear mdm i fn t selT stdObs pObs hash I
        % end
        
        for iAx = 1:length(features)
            fn = features(iAx);
            if isempty(xta.(dir)(iUnit).(fn))
                ax(iDir, iAx).Visible = false;
                continue
            end

            if p.nBoot > 0 && ismember(fn, p.std.features) && isfield(xta.(dir)(iUnit).(fn), 'stats')
                prcSTD = quantile(xta.(dir)(iUnit).(fn).stats.stdBoot, [0.95, 0.99, 0.999]);
                nStarsSTD = sum(xta.(dir)(iUnit).(fn).stats.std > prcSTD);
            else
                nStarsSTD = 0;
            end

            ax(iDir, iAx).Visible = true;

            hold(ax(iDir, iAx), 'on')
            % hold(ax(2+1-iDir, iAx), 'on')
            t = 1e3*xta.(dir)(iUnit).(fn).t;
            plot(ax(iDir, iAx), t, featureSign(iAx)*mean(xta.(dir)(iUnit).(fn).X, 1, 'omitnan'), Color=[0.15, 0.15, 0.15, 1], LineWidth=1.5, LineStyle=':');
            % plot(ax(2+1-iDir, iAx), t, X, Color=[0.15, 0.15, 0.15, 0.1], LineWidth=1.5, LineStyle=':');
            h = gobjects(p.mi.nClusters, 1);
            isClusterEmpty = false(1, p.mi.nClusters);
            for k = 1:p.mi.nClusters
                sel = idx==k;
                if nnz(sel) < p.mi.minNumTrialsPerCluster || nnz(sel) < length(idx)*p.mi.minNumTrialsPerClusterQuantile
                    isClusterEmpty(k) = true;
                    continue
                end
                c = getColor(k, p.mi.nClusters, 0.7);
                X = xta.(dir)(iUnit).(fn).X(sel, :);
                mu = featureSign(iAx)*mean(X, 1, 'omitnan');
                % err = std(X, 0, 1, 'omitnan')./sqrt(size(xta.(dir)(iUnit).(fn).X, 1));
    
                h(k) = plot(ax(iDir, iAx), t, mu, Color=c, LineStyle=lineStyles(k), LineWidth=1.5, DisplayName=sprintf("%s (n=%i)", p.mi.semanticClusterLabels(k), nnz(sel)));
            end
            if p.nBoot > 0 && isfield(xta.(dir)(iUnit).(fn), 'XBoot')
                prc = quantile(xta.(dir)(iUnit).(fn).XBoot, [p.bootAlpha/2, 1-p.bootAlpha/2], 1);
                patch(ax(iDir, iAx), [t, flip(t)], featureSign(iAx)*[prc(1, :), flip(prc(2, :))], [0.15, 0.15, 0.15], FaceAlpha=0.05, EdgeColor=[0.15, 0.15, 0.15], EdgeAlpha=0.5);
            end
            % xline(ax(iRow, iAx), 1e3*p.xta.meanWindow, 'k--', Alpha=0.1)
            xline(ax(iDir, iAx), 1e3*p.std.window, 'k--', Alpha=0.1)
            xline(ax(iDir, iAx), 0, 'k-', Alpha=0.1)
            % xticks(ax(iRow, iAx), 1e3*p.xta.meanWindow)
            xticks(ax(iDir, iAx), [-300, 0, 600])
            xtickangle(ax(iDir, iAx), 0)

            xlabel(ax(iDir, iAx), 'time (ms)')
            ylabel(ax(iDir, iAx), featureUnits(iAx))

            fnDisp = sprintf("%s %s", featureDispName(iAx), repmat('*', [1, nStarsSTD]));
            title(ax(iDir, iAx), fnDisp, Interpreter='none')
            ylim(ax(iDir, iAx), ylims{iAx})
            hold(ax(iDir, iAx), 'off')
            ax(iDir, iAx).YAxis.Direction = featureAxisDir(iAx);
            % hold(ax(2+1-iDir, iAx), 'off')
        end
        lgd = legend(ax(iDir, iAx), h(~isClusterEmpty), AutoUpdate=false);
        lgd.Layout.Tile = 'east';

        % Correlegram
        % 
        % iAx = iAx + 1;
        % r = NaN(length(p.std.features));
        % for i = 1:length(p.std.features)
        %     fni = p.std.features(statFeatureOrder(i));
        %     if isempty(xta.(dir)(iUnit).(fni))
        %         continue
        %     end
        %     for j = 1:length(p.std.features)
        %         fnj = p.std.features(statFeatureOrder(j));
        %         if isempty(xta.(dir)(iUnit).(fnj))
        %             continue
        %         end
        %         selT = xta.(dir)(iUnit).(fni).t >= p.std.window(1) & xta.(dir)(iUnit).(fni).t <= p.std.window(2);
        %         r(i, j) = corr(std(xta.(dir)(iUnit).(fni).X(:, selT), 0, 2, 'omitnan'), std(xta.(dir)(iUnit).(fnj).X(:, selT), 0, 2, 'omitnan'), Rows='complete');
        %     end
        % end
        % r(isnan(r)) = 0;
        % imagesc(ax(iDir, iAx), r);
        % xticks(ax(iDir, iAx), 1:length(p.std.features))
        % yticks(ax(iDir, iAx), 1:length(p.std.features))
        % xticklabels(ax(iDir, iAx), p.std.features(statFeatureOrder))
        % yticklabels(ax(iDir, iAx), p.std.features(statFeatureOrder))
        % xtickangle(ax(iDir, iAx), 90)
        % ax(iDir, iAx).XAxisLocation = 'top';
        % axis(ax(iDir, iAx), 'image')
        % ax(iDir, iAx).XAxis.Direction = 'normal';
        % clim(ax(iDir, iAx), [0, 1])
        % colormap(ax(iDir, iAx), 'gray')
        % colorbar(ax(iDir, iAx), 'eastoutside')
        % applyCustomColormap(ax(iDir, iAx), [-1, 1], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);
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

clear exportPath statFeatures nAx dirs dimensionReduction fig tlp tl ax iDir iAx iUnit dir
clear lia statFeatureOrder idx nClusters mr i fn t selT xx pcScore explained score eva k sel
clear mdm i fn t selT stdObs hash I idxSorted sepHash prcSTD nStarsSTD t X k c mu err prc
clear iDir dir fnDisp