%% Set root
% ROOTPATH = "E:\DATA";
ROOTPATH = 'C:\SERVER';

%% Clear temp vars
clearvars -except xta p kinematics ROOTPATH eu exp expIndices

%% Load data
% Load dip/rise triggered averages, the bootstraps contained within are kind of useless.
% load(fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot\LickVsReach_DLC_dta_rta_25_100_200to800ms_units1to1443_100boots.mat"));
load(fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot\LickVsReach_DLC_dta_rta_25_100_200to800ms_units1to1225_0boots_20260613.mat"));

% Load bootstrapped per-cluster averages. Can skip next step unless you
% want to recluster/rebootstrap
load(fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot\LickVsReach_DLC_miBoot_200to800ms_1225units_10000boots_20260614.mat")); % contains updated `p`

%% Combine movement indices (mi) across dips from all units, then cluster them. 
% Do this before bootstrapping per-cluster averages.
OVERWRITE_ALL = false;
requiredFields = ["features", "windowPre", "windowPost", "nClusters", "clusterMethod", "clusterDimensions", "clusterSeed", "semanticClusterOrder", "semanticClusterLabels", "semanticClusterSign", "dimensionReductionMethod"];
if ~isfield(p, 'mi') || any(~isfield(p.mi, requiredFields)) || OVERWRITE_ALL
    if isfield(p, 'mi') && ~OVERWRITE_ALL
        msg = sprintf("The following required p.mi fields are missing:\n%s\n\nOverwrite p.mi?", strjoin(requiredFields(~isfield(p.mi, requiredFields)), ", "));
        choice = questdlg(msg, "Missing p.mi fields", "Overwrite", "Cancel", "Cancel");
    else
        choice = 'Overwrite';
    end
    if strcmp(choice, "Overwrite")
        p.mi.features = ["Jaw", "Tongue", "HandL", "HandR", "Spine"];
        p.mi.windowPre = [-0.3, 0];
        p.mi.windowPost = [0, 0.3];
        p.mi.nClusters = 7;
        p.mi.clusterMethod = "kmeans"; % "gaussian", "kmeans"
        p.mi.clusterDimensions = 4;
        p.mi.clusterSeed = 42; % [-0.3, -0] vs [0, 0.3]: 43,47,42(hands are correlated a bit); [-0.6, -0] vs [0, 0.6]: 42; [-1, -0.3] vs [0, 0.6]: 42, 2
        p.mi.semanticClusterOrder = [1, 3, 2, 6, 7, 5, 4];
        % p.mi.semanticClusterOrder = 1:7;
        p.mi.semanticClusterLabels = ["no move", "lick start", "lick stop", "left hand reach", "left hand retract", "right hand reach", "right hand retract"];
        p.mi.semanticClusterSign = [0, 1, -1, 1, -1, 1, -1]; % 0: two-tailed, 1: right, 2: left

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
    end
    clear msg choice
else
    fprintf('All parameter fields ["%s"] already exist, we will use existing parameters.\n', strjoin(requiredFields, '", "'));
end
tTic = tic();

if exist('mi', 'var') && ~OVERWRITE_ALL
    choice = questdlg('Variable "mi" already exists in workspace. Recalculate and overwrite?', ...
        'Overwrite mi?', 'Yes', 'No', 'No');
else
    choice = 'Yes';
end
if strcmpi(choice, 'Yes')
    clear mi
    for dir = ["dip", "rise"]
        mi.(dir)(length(xta.(dir))) = struct(HandR=[], HandL=[], Spine=[], Jaw=[], Tongue=[]);
    end
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
else
    fprintf('Using existing variable "mi" which contains calculated movement indices, dimensionality reduction results, and clustering results.\n');
end

if exist('clusterSize', 'var') && ~OVERWRITE_ALL
    choice = questdlg('Recalculate existing variable "clusterSize"?', 'Recalculate clusterSize', 'Yes', 'No', 'No');
else
    choice = 'Yes';
end
if strcmpi(choice, 'Yes')
    if ~isfield(p.mi, 'minNumTrialsPerCluster')
        warning('Field p.mi.minNumTrialsPerCluster not found, setting it to 5.')
        p.mi.minNumTrialsPerCluster = 5;
    end
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
else
    fprintf('Using existing variable "clusterSize".\n')
end

clear requiredFields msg choice


% Scatter plot of all movement profiles
features = ["spikerate", "Jaw", "Tongue", "HandL", "HandR", "Spine"];
featureDispName = ["spike rate", "jaw", "tongue", "left paw", "right paw", "spine"];
featureUnits = ["spike rate (a.u.)", "DV pos (a.u.)", "protrusion prob", "AP pos (a.u.)", "AP pos (a.u.)", "DV pos (a.u.)"];
featureSign = [1, -1, 1, 1, 1, -1];
featureAxisDir = ["normal", "reverse", "normal", "normal", "normal", "normal"];
lineStyles = ["-", "-", "--", "-", "--", "-", "--"];
ylims = {5, 3, 1, 3, 3, 3};
ylims = cellfun(@(y) y*[-2/3, 1], ylims, UniformOutput=false);
ylims{featureAxisDir=="reverse"} = [-1, 2/3]*3;
% close all
layout.h = [4, 1, 1];
fig = figure(Units='normalized', Position=[0.05, 0.05, 0.9, 0.9]);
tlp = tiledlayout(fig, sum(layout.h), 1);
ax = nexttile(tlp, 1, [layout.h(1), 1]);
hold(ax, 'on')
h = gobjects(2, p.mi.nClusters);
iDir = 0;

idx = struct(dip=[], rise=[]);
dispScore = struct(dip=[], rise=[]);
for dir = ["dip", "rise"]
    idx.(dir) = vertcat(mi.(dir).idx);
end
for dir = ["dip", "rise"]
    dispScore.(dir) = vertcat(mi.(dir).dispScore);
end
dispScoreMerge = vertcat(dispScore.dip, dispScore.rise);
clear dir

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


clear tTic iUnit dir fn X t XPre XPost
clear iFeat fn pcaScoreMerge pcaExplained nClusters ax k sel
clear i0 iUnit dir n h
clear fig ax tl tlp layout iClu fn k faceColor dir iDir xl yl zl
clear dispScore dispScoreMerge featureAxisDir featureDispName features featureSign featureUnits idxMerge k0 nTrials pcaScore ylims
clear mpMean mp
%% Recalculate bootstrap (takes about 24hrs)
boot_dta_clustered_movement_index;

%% Count number of "clean clusters" by unit
% A clean cluster is a cluster of dips/rises where one bodypart moved but
% nothing else (e.g. for right-hand-reach cluster, tongue/left-hand/spine
% must be stationary)

% For each unit, calculate a pValue matrix (nClusters x nFeatures, clusters are in semantic order)
% For each cluster (row):
%   - do one-sided tests for the feature of interest in `mustMove`
%   - do two-sided tests for all features in `mustNotMove`
% pValue should be NaN for empty clusters, or failed bootstraps (sometimes XBoot or XObs can have too many NaNs if the bodypart is not always tracked by DeepLabCut)

nUnits = length(xta.dip);
nClusters = length(p.mi.semanticClusterLabels);
nFeatures = length(p.mi.boot.features);

p.mi.boot.nBonferroni = nClusters*nFeatures;
assert(isequal(p.mi.semanticClusterLabels, ["no move", "lick start", "lick stop", "left hand reach", "left hand retract", "right hand reach", "right hand retract"]))
p.mi.mustMove = {"", "Jaw", "Jaw", "HandL", "HandL", "HandR", "HandR"};
p.mi.mustNotMove = {["Jaw", "HandL", "HandR", "Spine"], ["HandL", "HandR", "Spine"], ["HandL", "HandR", "Spine"], ["Jaw", "HandR", "Spine"], ["Jaw", "HandR", "Spine"], ["Jaw", "HandL", "Spine"], ["Jaw", "HandL", "Spine"]};

clear cc ccData ccParams
ccData(nUnits) = struct(dip=struct(p=[], h=[], n=[]), rise=struct(p=[], h=[], n=[]));
for iUnit = 1:nUnits
    for dir = ["dip", "rise"]
        ccData(iUnit).(dir).p = NaN(nClusters, nFeatures, 'single');
        ccData(iUnit).(dir).h = false(nClusters, nFeatures);
        ccData(iUnit).(dir).n = zeros(nClusters, nFeatures, 'uint16');
        ccData(iUnit).(dir).clean = false(nClusters, 1);
    end
end
ccParams(nClusters) = struct(idx=[], label=[], expectedSign=[], mustMove=[], mustNotMove=[]);
% i: cluster displayOrder
% k: real cluster id
for iClu = 1:nClusters
    ccParams(iClu) = struct(idx=p.mi.semanticClusterOrder(iClu), label=p.mi.semanticClusterLabels(iClu), expectedSign=p.mi.semanticClusterSign(iClu), mustMove=p.mi.mustMove{iClu}, mustNotMove=p.mi.mustNotMove{iClu});
end
cc = struct(data=ccData, params=ccParams, featureNames=p.mi.boot.features, clusterNames=p.mi.semanticClusterLabels);
clear iClu iUnit dir pMatrix hMatrix nMatrix ccData ccParams

fid = fopen(fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot", sprintf("errorLog_%s.txt", datetime("now", Format="uuuuMMdd_HHmmss"))),'a');
isErrorLogEmpty = true;
for iUnit = 1:nUnits
    for dir = ["dip", "rise"]
        for iClu = 1:nClusters
            k = cc.params(iClu).idx; % True cluster ID
            if clusterSize(iUnit).(dir).nRaw(k) < p.mi.minNumTrialsPerCluster
                continue
            end
           
            % iMustMove = find(ismember(p.mi.boot.features, cc(iClu).mustMove)); % this is one index, sometimes empty
            % assert(length(iMustMove) <= 1)
            % Must not move
            % iMustNotMove = find(ismember(p.mi.boot.features, cc(iClu).mustNotMove)); % this is often an array
            % assert(length(iMustNotMove) >= 1)

            for iFeat = 1:nFeatures
                fn = p.mi.boot.features(iFeat);

                xObs = miObs(iUnit).(dir)(k, iFeat); % nClusters x nFeatures
                xBoot = miBoot(iUnit).(dir)(:, k, iFeat); % nBoot x nClusters x nFeatures

                % xObs is NaN, either:
                if isnan(xObs)
                    % a. zero trials/not enough trials for this cluster; or
                    if clusterSize(iUnit).(dir).nRaw(k) < p.mi.minNumTrialsPerCluster
                        break
                    % b. for these trials, this particular feature is missing
                    else
                        % % Log the same warning message to errorLog.txt
                        % fprintf(fid, "xObs is NaN: iUnit=%i, dir=%s, fn=%s, clu=%s, nTrials=%i\n", iUnit, dir, fn, p.mi.semanticClusterLabels(iClu), clusterSize(iUnit).(dir).nRaw(k));
                        % isErrorLogEmpty = false;
                        % % Also emit the original warning to MATLAB console
                        % warning("xObs is NaN: iUnit=%i, dir=%s, fn=%s, clu=%s, nTrials=%i", iUnit, dir, fn, p.mi.semanticClusterLabels(iClu), clusterSize(iUnit).(dir).nRaw(k));
                        continue
                    end
                end
                
                % xBoot should not be all NaNs, that's going in the book!
                if all(isnan(xBoot))     
                    % Log the same warning message to errorLog.txt
                    fprintf(fid, "xBoot is all NaN: iUnit=%i, dir=%s, fn=%s, clu=%s, nTrials=%i\n", iUnit, dir, fn, p.mi.semanticClusterLabels(iClu), clusterSize(iUnit).(dir).nRaw(k));
                    isErrorLogEmpty = false;
                    % Also emit the original warning to MATLAB console
                    warning("xBoot is all NaN: iUnit=%i, dir=%s, fn=%s, clu=%s, nTrials=%i", iUnit, dir, fn, p.mi.semanticClusterLabels(iClu), clusterSize(iUnit).(dir).nRaw(k))
                end

                % Feature of interest, one sided test based on expectation
                % of increase/decrease
                if ismember(fn, cc.params(iClu).mustMove)
                    expectedSign = p.mi.semanticClusterSign(iClu);
                % Two sided test for features that should not move
                else
                    expectedSign = 0;
                end

                switch expectedSign
                    case 1 % "right-sided test, expect increase"
                        pVal = sum(xBoot > xObs) / p.mi.boot.nBoot;
                    case -1 % "left-sided test, expect decrease"
                        pVal = sum(xBoot < xObs) / p.mi.boot.nBoot;
                    case 0
                        pValRight = sum(xBoot > xObs) / p.mi.boot.nBoot;
                        pValLeft = sum(xBoot < xObs) / p.mi.boot.nBoot;
                        pVal = min(pValLeft, pValRight) * 2;
                end

                cc.data(iUnit).(dir).p(iClu, iFeat) = pVal;
                cc.data(iUnit).(dir).n(iClu, iFeat) = uint16(clusterSize(iUnit).(dir).n(k));
            end
            clear iFeat fn xObs xBoot expectedSign pValLeft pValRight pVal
        end
        cc.data(iUnit).(dir).h = cc.data(iUnit).(dir).p < p.mi.boot.alpha / p.mi.boot.nBonferroni;
        clean = true(nClusters, 1);
        for iClu = 1:nClusters
            iMustMove = find(ismember(cc.featureNames, cc.params(iClu).mustMove));
            iMustNotMove = find(ismember(cc.featureNames, cc.params(iClu).mustNotMove));
            % Cluster of interest must move
            if ~isempty(iMustMove) && ~cc.data(iUnit).(dir).h(iClu, iMustMove) % This should be rare
                assert(isscalar(iMustMove))
                clean(iClu) = false;
                continue
            end
            % Other clusters must not move
            if any(cc.data(iUnit).(dir).h(iClu, iMustNotMove))
                clean(iClu) = false;
            end
        end
        cc.data(iUnit).(dir).clean = clean;
        clear iClu k clean iMustMove iMustNotMove
    end
    clear dir
end
fpath = fopen(fid);
fclose(fid);
if isErrorLogEmpty
    delete(fpath);
end
clear iUnit fid isErrorLogEmpty fpath




clear nUnits nFeatures nClusters
%% Count units
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

nUnits = length(xta.dip);
testNames = ["all clusters", "bodypart+movement specificity", "lick vs. reach vs. retract", "bodypart specificity", "start vs. stop specificity", "locomotion vs. consumption specificity"];
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
    nGroups = length(testGroups{iTest});
    for dir = ["dip", "rise"]
        tests(iTest).(dir).groupContainsCleanClusters = false(nUnits, nGroups); % nUnits x nGroups
        tests(iTest).(dir).groupContainsClusters = false(nUnits, nGroups); % nUnits x nGroups
        for iGrp = 1:length(tests(iTest).labels)
            selClusters = ismember([cc.params.label], tests(iTest).labels{iGrp});
            
            found = arrayfun(@(d) d.clean(selClusters), vertcat(cc.data.(dir)), UniformOutput=false);
            found = any(cat(2, found{:}), 1); % nClustersInGroup x nUnits -> 1 x nUnits          
            tests(iTest).(dir).groupContainsCleanClusters(:, iGrp) = found';

            hasData = arrayfun(@(d) any(max(d.n(selClusters, :), [], 2) > 0), vertcat(cc.data.(dir)), UniformOutput=true);
            tests(iTest).(dir).groupContainsClusters(:, iGrp) = hasData;
        end
        tests(iTest).(dir).numGroupsClean = sum(tests(iTest).(dir).groupContainsCleanClusters, 2);
        tests(iTest).(dir).numGroupsPresent = sum(tests(iTest).(dir).groupContainsClusters, 2);
        tests(iTest).(dir).prcGroupsClean = tests(iTest).(dir).numGroupsClean ./ tests(iTest).(dir).numGroupsPresent;
        tests(iTest).(dir).valid = tests(iTest).(dir).numGroupsPresent >= 2;
    end
end
clear iTest dir iGrp testNames testGroups selClusters hasData found n

close all
for iTest = 1:length(tests)
    fig = figure(Name=tests(iTest).name, Units='inches', Position=[1, 1, 5, 7]);
    tl = tiledlayout(fig, 4, 2, TileIndexing='columnmajor');
    title(tl, sprintf("%s\n[%s]", tests(iTest).name, strjoin(cellfun(@(labels) strjoin(labels, "/"), tests(iTest).labels), "] vs. [")))
    ax = gobjects(4, 2);
    iCol = 0;
    for dir = ["dip", "rise"]
        iCol = iCol + 1;
        iRow = 0;
        for fn = ["numGroupsClean", "numGroupsPresent", "prcGroupsClean"]
            iRow = iRow + 1;
            ax(iRow, iCol) = nexttile(tl);
            edges = 0:length(tests(iTest).labels)+1;% Each bin includes the leading edge, but does not include the trailing edge, except for the last bin which includes both edges.
            if fn == "prcGroupsClean"
                edges = edges./edges(end);
            end
            n = histcounts(tests(iTest).(dir).(fn), edges);
            % n = histcounts(tests(iTest).(dir).(fn)(tests(iTest).(dir).valid), edges);
            if fn ~= "prcGroupsClean"
                bar(ax(iRow, iCol), edges(1:end-1), n, 1, EdgeColor='black', FaceColor='black', FaceAlpha=0.5)
                xlim(ax(iRow, iCol), [edges(1), edges(end)]-0.5)
            else
                histogram(ax(iRow, iCol), BinEdges=edges, BinCounts=n, EdgeColor='black', FaceColor='black', FaceAlpha=0.5)
                xlim(ax(iRow, iCol), [-0.1, 1.1])
            end
            xlabel(ax(iRow, iCol), fn)
            ylabel(ax(iRow, iCol), "no. units")
            if iRow == 1
                title(ax(iRow, iCol), dir)
            end
        end
        iRow = iRow + 1;
        ax(iRow, iCol) = nexttile(tl);
        [n, xEdges, yEdges] = histcounts2(tests(iTest).(dir).numGroupsPresent, tests(iTest).(dir).numGroupsClean);
        n = n./nUnits;
        histogram2(ax(iRow, iCol), XBinEdges=xEdges, YBinEdges=yEdges, BinCounts=n, ...
            DisplayStyle='tile', ShowEmptyBins=true)
        xlabel(ax(iRow, iCol), "numGroupsPresent")
        ylabel(ax(iRow, iCol), "numGroupsClean")
        applyCustomColormap(ax(iRow, iCol), [0, 0.25], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);
        % colormap(ax(iRow, iCol), 'sky')
        colorbar(ax(iRow, iCol), 'eastoutside')
        axis(ax(iRow, iCol), 'equal')
        xticks(ax(iRow, iCol), xEdges+0.5)
        yticks(ax(iRow, iCol), yEdges+0.5)
        grid(ax(iRow, iCol), 'off');
    end
    fontsize(fig, 9, 'points')
end
clear iTest fig tl ax iCol dir iRow fn centers edges n xEdges yEdges

% Do rises contain more movement types than dips?
fig = figure(Units='inches', Position=[1 1 4 8]);
tl = tiledlayout(fig, 3, 2);
maxC = [0.10, 0.10, 0.15, 0.2, 0.3, 0.3];
for iTest = 1:length(tests)
    ax = nexttile(tl);
    [n, xEdges, yEdges] = histcounts2(tests(iTest).dip.numGroupsClean, tests(iTest).rise.numGroupsClean);
    n = n./nUnits;
    histogram2(ax, XBinEdges=xEdges, YBinEdges=yEdges, BinCounts=n, ...
        DisplayStyle='tile', ShowEmptyBins=true)
    applyCustomColormap(ax, [0, maxC(iTest)], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.5, h0=0.33);
    axis(ax, 'equal')
    title(ax, tests(iTest).name)
    xlabel(ax, 'dip')
    ylabel(ax, 'rise')

    xticks(ax, 0:max(xEdges))
    yticks(ax, 0:max(yEdges))

    cb = colorbar(ax, Location='southoutside');
    % cb.Layout.Tile = 'south';
    cb.Label.String = sprintf("%% of %i units", nUnits);
    cb.Ticks = [0, maxC(iTest)];
    cb.TickLabels = ["0", sprintf("%g%%", maxC(iTest)*100)];

end

fontsize(fig, 9, 'points')
clear fig tl iTest ax n xEdges yEdges cb maxC

%%
clc
for iTest = 1:length(tests)
    fprintf("Test %i (%s):\n%s:\n", ...
        iTest, tests(iTest).name, ...
        strjoin(cellfun(@(labels) sprintf("[%s]", strjoin(labels, ", ")), tests(iTest).labels), " vs. ") ...
        )
    for dir = ["dip", "rise"]
        for n = 1:length(tests(iTest).labels)
            if n > 1
                plural = "s";
            else
                plural = "";
            end
            fprintf("\t%s\t%2i%% (%4i/%i total units)\t%2i%% (%4i/%i valid units)\thas\t%i clean cluster%s.\n", ...
                dir, ...
                round(100*nnz(tests(iTest).(dir).numGroupsClean == n) / nUnits), ...
                nnz(tests(iTest).(dir).numGroupsClean == n), nUnits, ...
                round(100*nnz(tests(iTest).(dir).numGroupsClean == n & tests(iTest).(dir).numGroupsPresent >=n) / nnz(tests(iTest).(dir).numGroupsPresent >=n)), ...
                nnz(tests(iTest).(dir).numGroupsClean == n & tests(iTest).(dir).numGroupsPresent >=n), nnz(tests(iTest).(dir).numGroupsPresent >=n), ...
                n, plural);
        end
    end
    for dir = ["dip", "rise"]
        for n = 2:length(tests(iTest).labels)
            fprintf("\t%s\t%2i%% (%4i/%i total units)\t%2i%% (%4i/%i valid units)\thas\t%i+ clean cluster%s.\n", ...
                dir, ...
                round(100*nnz(tests(iTest).(dir).numGroupsClean >= n) / nUnits), ...
                nnz(tests(iTest).(dir).numGroupsClean >= n), nUnits, ...
                round(100*nnz(tests(iTest).(dir).numGroupsClean >= n & tests(iTest).(dir).numGroupsPresent >=n) / nnz(tests(iTest).(dir).numGroupsPresent >=n)), ...
                nnz(tests(iTest).(dir).numGroupsClean >= n & tests(iTest).(dir).numGroupsPresent >=n), nnz(tests(iTest).(dir).numGroupsPresent >=n), ...
                n, plural);
        end
    end
end
clear iTest n dir


%% Plot individual units just-dip matrix (rows are features, columns are clusters)
close all

features = ["spikerate"; "Jaw"; "HandL"; "HandR"; "Spine"];
% featureDispName = ["spike rate"; "jaw"; "left hand"; "right hand"; "spine"];
% featureUnits = ["(a.u.)"; "(DV a.u.)"; "(AP a.u.)"; "(AP a.u.)"; "(DV a.u.)"];
clusterDispName = ["none", "jaw\nopen", "jaw\nclose", "handL\nreach", "handL\nretract", "handR\nreach", "handR\nretract"];
featureDispName = ["spike"; "jaw"; "handL"; "handR"; "spine"];
featureUnits = ["(a.u.)"; "(DV a.u.)"; "(AP a.u.)"; "(AP a.u.)"; "(DV a.u.)"];
featureSign = [1; -1; 1; 1; -1];
featureAxisDir = ["normal"; "reverse"; "normal"; "normal"; "normal"];
% p.mi.semanticClusterLabels = ["no move", "lick start", "lick stop", "left hand retract", "left hand reach", "right hand retract", "right hand reach"];
% semanticClusterSign = [0, 1, -1, -1, 1, -1, 1]; % 0: two-tailed, 1: right, 2: left
% semanticClusterSign = [0, 1, -1, -1, 1, -1, 0]; % 0: two-tailed, 1: right, 2: left
ylims = {5; 5; 5; 5; 5};
ylims = cellfun(@(y) y*[-2/3, 1], ylims, UniformOutput=false);
ylims{featureAxisDir=="reverse"} = [-1, 2/3]*5;
% xt = unique([p.mi.windowPost, p.mi.windowPre]);
xt = [0, 1];

nClusters = length(p.mi.semanticClusterOrder);
nFeatures = length(features);

fig = figure(Units='inches', InnerPosition=[1, 1, 0.4*nClusters, 0.6*length(features)]);
tl = tiledlayout(fig, nFeatures, nClusters, TileSpacing='tight', Padding='tight', TileIndexing='columnmajor');

ax = gobjects(nFeatures, nClusters);
for iClu = 1:nClusters
    for iFeat = 1:nFeatures
        ax(iFeat, iClu) = nexttile(tl);
        ax(iFeat, iClu).Box = 'on';
    end
end

% for iTest = [3, 2, 4]
for iTest = 2
    % for nCleanGroups = [5, 4, 3, 2, 1]
    for nCleanGroups = 6
        for requirement = ["dip or rise", "both"]
            exportPath = fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot\Figures", sprintf("%s_LickVsReach_DLC_dta_rta_%i_%i_%ito%ims_std", datetime("now", Format="uuuuMMdd"), 100*p.xta.dip.thresholdQuantile, 100*p.xta.dip.thresholdSubQuantile, p.spikeRes*1000*p.xta.dip.samples(1), p.spikeRes*1000*p.xta.dip.samples(2)), "By Feature x Cluster", sprintf("Test%i - %s", iTest, tests(iTest).name), sprintf('%i clean clusters for %s', nCleanGroups, requirement));
            if ~exist(exportPath, 'dir')
                mkdir(exportPath)
            else
                rmdir(exportPath, 's')
                mkdir(exportPath)   
            end
            
            for dir = ["dip", "rise"]
                switch requirement
                    case "dip or rise"
                        selUnits = reshape(find(tests(iTest).(dir).numGroupsClean == nCleanGroups), 1, []);
                    case "both"
                        selUnits = reshape(find(tests(iTest).dip.numGroupsClean == nCleanGroups & tests(iTest).rise.numGroupsClean == nCleanGroups), 1, []);
                    otherwise
                        error("unknown requirement=%s", requirment)
                end
                for iUnit = selUnits
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
                        k = p.mi.semanticClusterOrder(iClu); % raw cluster index
                        nTrials = clusterSize(iUnit).(dir).n(k);
                        if isnan(nTrials) || nTrials < p.mi.minNumTrialsPerCluster
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
                            
                            % Draw n trials
                            if iFeat == 1
                                text(ax(iFeat, iClu), -0.9, 4, sprintf("n=%i", nTrials), HorizontalAlignment='left', ...
                                    VerticalAlignment='top', FontSize=8, FontWeight='normal')
                            % Check significance of cluster movement index against bootstrap
                            else
                                hVal = cc.data(iUnit).(dir).h(iClu, iFeat-1);
                                nStars = nnz(hVal);
                                if featureAxisDir(iFeat) == "reverse"
                                    yPos = ylims{iFeat}(1)+0.5;
                                else
                                    yPos = ylims{iFeat}(2)-0.5;
                                end
                                text(ax(iFeat, iClu), mean(p.mi.windowPost), yPos, repmat('*', [1, nStars]), ...
                                    HorizontalAlignment='center', VerticalAlignment='top', FontSize=12, FontWeight='bold', Color='red');
                            end
            
                            plot(ax(iFeat, iClu), t, mu, Color=c, LineStyle='-', LineWidth=1.5);
                            if iFeat > 1 && nStars >= 1
                                selT = isin(t, p.mi.windowPre) | isin(t, p.mi.windowPost);
                                plot(ax(iFeat, iClu), t(selT), mu(selT), Color='red', LineStyle='-', LineWidth=1.5);
                                clear selT
                            end
                            clear pVal nStars yPos
                            % xline(ax(iFeat, iClu), xt, 'k--', Alpha=0.1)
                            xline(ax(iFeat, iClu), 0, 'k-', Alpha=0.1)
                            yline(ax(iFeat, iClu), 0, 'k-', Alpha=0.1)
                            if iFeat == nFeatures
                                xticks(ax(iFeat, iClu), xt)
                                xtickangle(ax(iFeat, iClu), 0)
                            else
                                xticks(ax(iFeat, iClu), [])
                            end
                            if iClu == 1
                                if featureAxisDir(iFeat) == "reverse"
                                    yticks(ax(iFeat, iClu), [ylims{iFeat}(1), 0])
                                    yticklabels(ax(iFeat, iClu), [-ylims{iFeat}(1), 0])
                                else
                                    yticks(ax(iFeat, iClu), [0, ylims{iFeat}(2)])
                                end
                            else
                                yticks(ax(iFeat, iClu), [])
                            end
            
                            % if iFeat == nFeatures
                            %     xlabel(ax(iFeat, iClu), 'time (ms)')
                            % end
                            if iClu == 1
                                % ylabel(ax(iFeat, iClu), sprintf("%s\n%s", featureDispName(iFeat), featureUnits(iFeat)))
                                ylabel(ax(iFeat, iClu), sprintf("%s", featureDispName(iFeat)))
                            end
                            if iFeat == 1
                                % title(ax(iFeat, iClu), sprintf("%s\n(n=%i %ss)", p.mi.semanticClusterLabels(iClu), nTrials, dir), Interpreter='none', FontWeight='bold', Color='black')
                                title(ax(iFeat, iClu), strsplit(clusterDispName(iClu), "\\n"), Interpreter='none', FontWeight='normal', Color='black')
                            end
                            ylim(ax(iFeat, iClu), ylims{iFeat})
                            hold(ax(iFeat, iClu), 'off')
                            ax(iFeat, iClu).YAxis.Direction = featureAxisDir(iFeat);
                        end
                    end
                    xlim(ax, [-1, 1])
                    xlabel(tl, sprintf("time since %s (s)", dir))
                    title(tl, sprintf("Unit %i (n=%i %ss)", iUnit, length(xta.(dir)(iUnit).t0), dir), FontWeight='bold')
                    fontsize(fig, 8, 'points')
                    print(fig, fullfile(exportPath, sprintf("featxclus_unit_%03i_%s", iUnit, dir)), '-dpng', '-r0')
                end
            end
        end
    end
end

clear features featureDispName featureUnits featureSign featureAxisDir semanticClusterSign clusterDispOrder ylims xt nClusters nFeatures fig tl ax iClu iFeat iTest nCleanGroups exportPath selUnits iUnit dir idx iClu k nTrials fn sel c X t mu xObs xBoot pVal nStars 
