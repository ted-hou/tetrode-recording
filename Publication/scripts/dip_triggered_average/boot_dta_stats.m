% Do bootstrapped statistics for Figure 6 

%% Set root
% ROOTPATH = "E:\DATA";
ROOTPATH = 'C:\SERVER';

%% Clear temp vars
clearvars -except xta p kinematics ROOTPATH eu exp expIndices mi miBoot miObs

%% Load data
% Load dip/rise triggered averages, the bootstraps contained within are kind of useless.
% load(fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot\LickVsReach_DLC_dta_rta_25_100_200to800ms_units1to1443_100boots.mat"));
load(fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot\LickVsReach_DLC_dta_rta_25_100_200to800ms_units1to1225_0boots_20260613.mat"));

% Load bootstrapped per-cluster averages. Can skip next step unless you
% want to recluster/rebootstrap
load(fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot\LickVsReach_DLC_miBoot_200to800ms_1225units_10000boots_20260614.mat")); % contains updated `p`

%% Notes
% Before bootstrapping
% a) validate 3SD as a potential replacement for 10000 boots (where alpha=0.05/nBonferroni)
%    mu,sd from movement indices from all dips+rises
% b) calculate kmeans cluster centroids

% For each bootstrap:
% 1) randomly sample the session instead of picking true dips(rises)
% 2) recalculate dip(rise)-triggered-averages
% 3) calculate movement indices and cluster using existing kmeans centroids
% 4) calculate mean movement indices for each cluster, test for
%    significance using 3SD
% 5) determine for each unit which clusters are present and clean
%    i.e., calculate `cc.data`
% 6) Redo the statistics in Figures 6e-f: from `cc.data`
%    a) `csTrials`: count no. trials (dips/rises) per cluster, avg across units (6e top)
%    b) `csUnits`: count no. units that contain 5+ dip/rises per cluster (6e bottom)
%    c) `histcounts(ncc)`: histcount no. moves per unit, dip, rise, dip vs. rise (6f)

%% a) validate 3SD as a potential replacement for 10000 boots (where alpha=0.05/nBonferroni)
nUnits = length(xta.dip);
tryNSigmas = 1:0.1:3;
results(length(tryNSigmas)) = struct(nMismatchUnits=[], nMismatchFeatures=[], nFalse=[], nFalsePositives=[], nFalseNegatives=[]);

for iTry = 1:length(tryNSigmas)
    nSigmas = tryNSigmas(iTry);
    nMismatchUnits = 0;
    nMismatchFeatures = 0;
    nFalsePositives = 0;
    nFalseNegatives = 0;
    parfor iUnit = 1:nUnits
        isInvalid = logical(arrayfun(@(i) mean(nnz(mi.dip(iUnit).idx==i)<p.mi.minNumTrialsPerCluster), p.mi.semanticClusterOrder));
        hasMismatch = false;
        for ifn = 1:length(p.mi.boot.features)
            fn = p.mi.boot.features(ifn);
            obsSD = arrayfun(@(i) mean(mi.dip(iUnit).(fn)(mi.dip(iUnit).idx==i), 'omitnan'), p.mi.semanticClusterOrder);
            obsBoot = miObs(iUnit).dip(p.mi.semanticClusterOrder, ifn)';
            hLeftSD = arrayfun(@(i) mean(mi.dip(iUnit).(fn)(mi.dip(iUnit).idx==i), 'omitnan'), p.mi.semanticClusterOrder) >= mean(mi.dip(iUnit).(fn), 'omitnan') + nSigmas*std(mi.dip(iUnit).(fn), 0, 'omitnan');
            hLeftBoot = miObs(iUnit).dip(p.mi.semanticClusterOrder, ifn)' >= quantile(miBoot(iUnit).dip(:, p.mi.semanticClusterOrder, ifn), 1-0.05/28, 1);
            hRightSD = arrayfun(@(i) mean(mi.dip(iUnit).(fn)(mi.dip(iUnit).idx==i), 'omitnan'), p.mi.semanticClusterOrder) <= mean(mi.dip(iUnit).(fn), 'omitnan') - nSigmas*std(mi.dip(iUnit).(fn), 0, 'omitnan');
            hRightBoot = miObs(iUnit).dip(p.mi.semanticClusterOrder, ifn)' <= quantile(miBoot(iUnit).dip(:, p.mi.semanticClusterOrder, ifn), 0.05/28, 1);
            obsSD(isInvalid) = NaN;
            hLeftSD(isInvalid) = false;
            hRightSD(isInvalid) = false;
        
            if ~isequal(hLeftSD, hLeftBoot) || ~isequal(hRightSD, hRightBoot)
                % fprintf("\nUnit %i, %s:", iUnit, fn)
                % fprintf("\nobsSD:     "), fprintf("\t%.3f", obsSD)
                % fprintf("\nobsBoot:   "), fprintf("\t%.3f", obsBoot)
                % fprintf("\nhLeftSD:   "), fprintf("\t%i", hLeftSD)
                % fprintf("\nhLeftBoot: "), fprintf("\t%i", hLeftBoot)
                % fprintf("\nhRightSD:  "), fprintf("\t%i", hRightSD)
                % fprintf("\nhRightBoot:"), fprintf("\t%i", hRightBoot)
                % fprintf('\n')
                nMismatchFeatures = nMismatchFeatures + 1;
                nFalsePositives = nFalsePositives + sum(hLeftSD > hLeftBoot) + sum(hRightSD > hRightBoot);
                nFalseNegatives = nFalseNegatives + sum(hLeftSD < hLeftBoot) + sum(hRightSD < hRightBoot);
                hasMismatch = true;
            end
        end
        if hasMismatch
            nMismatchUnits = nMismatchUnits + 1;
        end
    end
    results(iTry) = struct(nMismatchUnits=nMismatchUnits, nMismatchFeatures=nMismatchFeatures, nFalse=nFalsePositives+nFalseNegatives, nFalsePositives=nFalsePositives, nFalseNegatives=nFalseNegatives);
    fprintf('%gSD: %i mismatched units,\t%i mismatched features,\t%i false,\t%i false positives,\t%i false negatives;\n', nSigmas, nMismatchUnits, nMismatchFeatures, nFalsePositives+nFalseNegatives, nFalsePositives, nFalseNegatives)
end


colors = 'yckrb';
fig = figure(Units='inches', Position=[1, 1, 4, 3]);
tl = tiledlayout(fig, 1, 2);
ax = nexttile(tl);
hold(ax, 'on')
nComparisons = nUnits*length(p.mi.boot.clusters)*length(p.mi.boot.features);
i = 0;
for fn = string(fieldnames(results)')
    i = i + 1;
    plot(ax, tryNSigmas, [results.(fn)], '-', Color=colors(i), DisplayName=fn);
end
legend(ax, Location='southoutside')
xlabel(ax, 'threshold (SD)')
ylabel(ax, 'error count')
hold(ax, 'off')
title(ax, sprintf("%i units", nUnits))

ax = nexttile(tl);
hold(ax, 'on')
colors = 'rbk';
fNames = ["nFalsePositives", "nFalseNegatives", "nFalse"];
fDispNames = ["false positive rate", "false negative rate", "error rate"];
for i = 1:3
    plot(ax, tryNSigmas, 100*[results.(fNames(i))]./nComparisons, '-', Color=colors(i), DisplayName=fDispNames(i))
end
hold(ax, 'off')
xlabel(ax, 'threshold (SD)')
ylabel(ax, 'error rate (%)')
ylim(ax, [0, 10])
legend(ax, Location='southoutside')
title(ax, sprintf("%i total comparisons\n(unit x feature x cluster)", nComparisons))

fontsize(fig, 8, 'points')

copygraphics(fig, BackgroundColor='none', ContentType='vector')

% After some thinking I decided to used 3SD as a threshold to approximate
% boostrapping with an alpha of 0.05/28. We will get very few false
% positives (<1%) but more false negatives (8%)
% sometimes bootstrapping can be more sensitive, since in this test, our SD
% is calculated during real dips, so more movement/SD is expected; in
% reality, we'd be generating fake dips from random parts of the session,
% so the SD should be smaller, and we'd get fewer false negatives.
clearvars -except xta p kinematics ROOTPATH eu exp expIndices mi miBoot miObs

%% Notes:
% Question: We know that dip-triggered average kinematics appear flat (i.e.
% mean movement index across dips, mi=0). Is that because it is flat on all
% dips, or because positive and negative movements cancel out?
% In other words, do dips tend to indicate "any movement" better than
% chance?

%% For each dip, rectify the movement index, then take the max among all
% body parts
nUnits = length(xta.dip);
nClusters = length(p.mi.semanticClusterOrder);
nFeatures = length(p.mi.boot.features);
clear miMaxRectified, miMaxRectified(nUnits) = struct(dip=[], rise=[]);
for iUnit = 1:length(xta.dip)
    for dir = ["dip", "rise"]
        X = arrayfun(@(fn) mi.(dir)(iUnit).(fn), p.mi.boot.features, UniformOutput=false);
        X = cat(2, X{:});
        X = max(abs(X), [], 2);
        miMaxRectified(iUnit).(dir).X = X;
        miMaxRectified(iUnit).(dir).idx = mi.(dir)(iUnit).idx;
    end
end
clear iUnit dir X


%% b) calculate kmeans cluster centroids
p.sd.nBoot = 1000;
nUnits = length(xta.dip);
nFeatures = length(p.mi.boot.features);
nClusters = length(p.mi.boot.clusters);
if ~exist('expIndices', 'var')
    expIndices = [xta.dip.iExp];
end

assert(ismember(p.mi.dimensionReductionMethod, ["manual", "manual+umap"]), "We assert this so we don't need to recapitulate the PCA.")

% Calculate existing cluster centroids in manual feature space (2nd dimension is tongue+jaw)
idx = struct(dip=[], rise=[]);
pcaScore = idx;
for dir = ["dip", "rise"]
    idx.(dir) = vertcat(mi.(dir).idx);
    pcaScore.(dir) = vertcat(mi.(dir).pcaScore);
end
idx = vertcat(idx.dip, idx.rise);
pcaScore = vertcat(pcaScore.dip, pcaScore.rise);
centroids = NaN(nClusters, size(pcaScore, 2));
for k = 1:nClusters
    centroids(k, :) = mean(pcaScore(idx==k, :), 1, 'omitnan');
end

clear idx pcaScore k dir


% Calculate mi from shuffled data, and assign to existing clusters
if isempty(gcp('nocreate'))
    pool = parpool('Processes');
else
    pool = gcp('nocreate');
end
% rng(42) % This does nothing in parfor, need to do something kind of streaming
assert(isequal(p.mi.features, ["Jaw", "Tongue", "HandL", "HandR", "Spine"]))
pcFeatures = cellfun(@(f) ismember(p.mi.features, f), {"HandL", ["Tongue", "Jaw"], "HandR", "Spine"}, UniformOutput=false);
miBootFeatures = ismember(p.mi.features, p.mi.boot.features);
clear sdBoot
sdBoot(nUnits) = struct(dip=struct(h=[], n=[]), rise=struct(h=[], n=[]));
parfevalOnAll(pool, @warning, 0, 'off', 'stats:kmeans:FailedToConverge');
%%
tTic = tic();
ll = 0;
hasWarning = false;
for iUnit = 1:nUnits
    if ~hasWarning
        fprintf(repmat('\b', [1, ll]));
    end
    hasWarning = false;
    ll = fprintf("Unit %i/%i, time elapsed: %.1f seconds; estimated remaining: %.1f seconds...\n", iUnit, nUnits, toc(tTic), toc(tTic)/iUnit*(nUnits-iUnit));
    iExp = expIndices(iUnit);
    kine = kinematics(iExp);
    tMax = min(arrayfun(@(fn) kine.(fn).t(end), p.mi.features));
    for dir = ["dip", "rise"]
        nBootTemp = zeros(p.sd.nBoot, nClusters, 'uint16');
        hBootTemp = zeros(p.sd.nBoot, nClusters, nFeatures, 'int8');
        nTrials = length(mi.(dir)(iUnit).idx);
        miMaxRectifiedTemp = zeros(nTrials, p.sd.nBoot, 'single');
        idxTemp = zeros(nTrials, p.sd.nBoot, 'uint16');
        if nTrials < p.mi.minNumTrialsPerCluster || nTrials < nClusters
            sdBoot(iUnit).(dir) = struct(h=hBootTemp, n=nBootTemp);
            warning("Unit %i has too few (%i) %ss, we cannot do kmeans so this unit/dir is skipped.", iUnit, nTrials, dir)
            hasWarning = true;
            continue
        end
        parfor iBoot = 1:p.sd.nBoot
            pcaScore = zeros(nTrials, 4, 'single');
            T0 = -p.mi.windowPre(1) + rand([1, nTrials])*(tMax - p.mi.windowPost(2) + p.mi.windowPre(1));
            miTemp = NaN(nTrials, length(p.mi.features), 'single');
            for iTrial = 1:length(T0)
                for iFeat = 1:length(p.mi.features)
                    fn = p.mi.features(iFeat);
                    [iStartPre, iStopPre] = isin(kine.(fn).t, T0(iTrial) + p.mi.windowPre, true, true);
                    [iStartPost, iStopPost] = isin(kine.(fn).t, T0(iTrial) + p.mi.windowPost, true, true);
                    miTemp(iTrial, iFeat) = mean(kine.(fn).X(iStartPost:iStopPost), 'omitnan') - mean(kine.(fn).X(iStartPre:iStopPre), 'omitnan');
                end
                for iFeat = 1:length(pcFeatures)
                    pcaScore(:, iFeat) = mean(miTemp(:, pcFeatures{iFeat}), 2, 'omitnan');
                end
            end
            miTemp = miTemp(:, miBootFeatures);
            miTemp(isnan(miTemp)) = 0;
            pcaScore(isnan(pcaScore)) = 0;
            idx = kmeans(pcaScore, nClusters, Start=centroids, MaxIter=0); % assign to existing clusters without updating centroids

            % Now we have:
            %   miTemp (nTrials x 4 features, removed tongue)
            %   idx (cluster id)
            % We need:
            %   convert miTemp from magnitude to significance (3SD)
            %   h (nClusters x nFeatures, 7x4) -1, 0, 1 for dec, flat, inc
            %   n (nClusters): numTrialsPerCluster
            %   miMaxRectifiedTemp
            %   idx
            mu = mean(miTemp, 1, 'omitnan');
            sd = std(miTemp, 0, 1, 'omitnan');
            hTemp = zeros(nClusters, nFeatures, 'int8');
            for k = 1:nClusters
                inCluster = idx==k;
                nTrialsInCluster = nnz(inCluster);
                nBootTemp(iBoot, k) = nTrialsInCluster;
                if nTrialsInCluster < p.mi.minNumTrialsPerCluster
                    continue
                end
                miCluster = mean(miTemp(inCluster, :), 1, 'omitnan');
                for ifn = 1:nFeatures
                    if miCluster(ifn) >= mu + 3*sd
                        hTemp(k, ifn) = 1;
                    elseif miCluster(ifn) <= mu - 3*sd
                        hTemp(k, ifn) = -1;
                    end
                end
            end
            hBootTemp(iBoot, :, :) = hTemp;

            % Tack on the miMaxRectified stuff
            miMaxRectifiedTemp(:, iBoot) = max(abs(miTemp), [], 2);
            idxTemp(:, iBoot) = idx;
        end
        sdBoot(iUnit).(dir) = struct(h=hBootTemp, n=nBootTemp, miMaxRectified=miMaxRectifiedTemp, idx=idxTemp);
        ll = ll + fprintf('%i positives.\n', nnz(hBootTemp~=0));
    end
end


parfevalOnAll(pool, @warning, 0, 'on', 'stats:kmeans:FailedToConverge');
clear pool
clear nUnits nFeatures nClusters idx pcaScore centroids k dir
clear pcFeatures miBootFeatures tTic ll iExp kine tMax dir nBootTemp hBootTemp nTrials iBoot pcaScore T0 miTemp iTrial iFeat fn iStartPre iStopPre iStartPost iStopPost idx mu sd hTemp k inCluster nTrialsInCluster miCluster ifn


exportPath = fullfile("C:\SERVER\LickVsReach_DTA_RTA_boot", sprintf("LickVsReach_DLC_sdBoot_%iunits_%iboots_%s.mat", length(xta.dip), p.sd.nBoot, datetime("now", Format="yyyyMMdd")));
save(exportPath, 'sdBoot', 'p', 'miMaxRectified', '-v7.3')
fprintf("Saved to %s\n", exportPath);

warning('sdBoot clusters are in native order, not semantic order')


%% Some plotting heh
for iUnit = 1:nUnits
    for dir = ["dip", "rise"]
        if isfield(sdBoot(iUnit).(dir), 'miMaxRectified')
            miMaxRectified(iUnit).(dir).XBoot = sdBoot(iUnit).(dir).miMaxRectified;
            miMaxRectified(iUnit).(dir).idxBoot = sdBoot(iUnit).(dir).idx;
        else
            miMaxRectified(iUnit).(dir).XBoot = [];
            miMaxRectified(iUnit).(dir).idxBoot = [];
        end
    end
end

%% Do per-unit tests for miMaxRectified
close all
clear pVal
pVal.all = struct(dip=NaN(nUnits, 1), rise=NaN(nUnits, 1));
pVal.clu1 = struct(dip=NaN(nUnits, 1), rise=NaN(nUnits, 1));
for iUnit = 1:nUnits
    for dir = ["dip", "rise"]
        if isempty(miMaxRectified(iUnit).(dir).XBoot)
            continue
        end
        % Cluster 1
        xObs = mean(miMaxRectified(iUnit).(dir).X(miMaxRectified(iUnit).(dir).idx==1), 1, 'omitnan');
        XBoot = arrayfun(@(iBoot) miMaxRectified(iUnit).(dir).XBoot(miMaxRectified(iUnit).(dir).idxBoot(:, iBoot)==1, iBoot), 1:p.sd.nBoot, UniformOutput=false);
        XBoot = cellfun(@(XBoot) mean(XBoot, 1, 'omitnan'), XBoot);
        pVal.clu1.(dir)(iUnit) = nnz(XBoot<=xObs)./length(XBoot);

        % All
        xObs = mean(miMaxRectified(iUnit).(dir).X, 1, 'omitnan');
        XBoot = mean(miMaxRectified(iUnit).(dir).XBoot, 1, 'omitnan');
        pVal.all.(dir)(iUnit) = nnz(XBoot<=xObs)./length(XBoot);   
    end
end

alpha = 0.05;
edges = 0:alpha/2:1;

for src = ["all", "clu1"]

    fig = figure(Units='inches', Position=[3, 3, 5, 3]);
    tl = tiledlayout(fig, 1, 1);
    AX = gobjects(1, 3);
    ax = nexttile(tl); AX(1) = ax;
    n = histcounts2(pVal.(src).dip, pVal.(src).rise, edges, edges);
    histogram2(ax, XBinEdges=edges, YBinEdges=edges, BinCounts=n, ...
        DisplayStyle='tile', ShowEmptyBins=true, EdgeAlpha=0);
    maxC = 100;
    cb = colorbar(ax);
    applyCustomColormap(ax, [0, maxC], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.5, h0=0.33);
    cb.Label.String = 'no. units';
    cb.Label.Position(1) = 0;
    
    xlabel(ax, 'P(x_{dip}\geqx_{boot})')
    ylabel(ax, 'P(x_{rise}\geqx_{boot})')
    xlim(ax, [0, 1])
    ylim(ax, [0, 1])
    axis(ax, 'equal')
    ax.Box = 'off';
    
    ax = nexttile(tl, 'south'); AX(2) = ax;
    histogram(ax, pVal.(src).dip, edges, EdgeColor='none', FaceColor='k', FaceAlpha=1, DisplayName='dip', DisplayStyle='bar')
    xlim(ax, [0, 1])
    ylabel(ax, 'no. units')
    ylim(ax, [0, nUnits])
    yticks(ax, [0, nUnits])
    
    ax = nexttile(tl, 'west'); AX(3) = ax;
    histogram(ax, pVal.(src).rise, edges, EdgeColor='none', FaceColor='k', FaceAlpha=1, DisplayName='rise', DisplayStyle='bar', Orientation='horizontal')
    ylim(ax, [0, 1])
    ax.XDir = 'reverse';
    ax.YAxisLocation = 'right';
    xlabel(ax, 'no. units')
    xlim(ax, [0, nUnits])
    xticks(ax, [0, nUnits])
    xtickangle(ax, 90)

    xticks(AX([1, 2]), [0, 0.5, 1])
    yticks(AX([1, 3]), [0, 0.5, 1])

    switch src
        case "all"
            title(tl, "dips/rises in all clusters", FontWeight='bold')
        case "clu1"
            title(tl, sprintf("dips/rises in cluster 1 (""%s"")", p.mi.semanticClusterLabels(1)), FontWeight='bold')
    end

    fontsize(fig, 8, 'points')
end

%% Plot aggregate distributions
close all
fig = figure(Units='inches', Position=[3, 3, 10, 6]);
tl = tiledlayout(fig, 2, 2, TileIndexing='rowmajor');
ax = gobjects(2, 2);
for i = 1:2
    for j = 1:2
        ax(i, j) = nexttile(tl);
        hold(ax(i, j), 'on')
    end
end

iCol = 0;
for dir = ["dip", "rise"]
    iCol = iCol + 1;
    X = arrayfun(@(data) data.(dir).X, miMaxRectified, UniformOutput=false);
    X = cat(1, X{:});
    XBoot = arrayfun(@(data) data.(dir).XBoot, miMaxRectified, UniformOutput=false);
    XBoot = cat(1, XBoot{:});
    idx = arrayfun(@(data) data.(dir).idx, miMaxRectified, UniformOutput=false);
    idx = cat(1, idx{:});
    idxBoot = arrayfun(@(data) data.(dir).idxBoot, miMaxRectified, UniformOutput=false);
    idxBoot = cat(1, idxBoot{:});

    iRow = 1;
    for iClu = 1:length(p.mi.semanticClusterOrder)
        k = p.mi.semanticClusterOrder(iClu);
        sel = idxBoot==k;
        histogram(ax(iRow, iCol), XBoot(sel), [0:0.1:5, Inf], Normalization='pdf', EdgeColor=getColor(iClu, 7, 0.7), DisplayName=p.mi.semanticClusterLabels(iClu), DisplayStyle='stairs');
    end
    histogram(ax(iRow, iCol), XBoot(:), [0:0.1:5, Inf], Normalization='pdf', EdgeColor='k', DisplayName='all', DisplayStyle='stairs', LineWidth=2)
    title(ax(iRow, iCol), sprintf('Bootstrap (fake %ss)', dir))

    iRow = 2;
    for iClu = 1:length(p.mi.semanticClusterOrder)
        k = p.mi.semanticClusterOrder(iClu);
        sel = idx==k;
        histogram(ax(iRow, iCol), X(sel), [0:0.1:5, Inf], Normalization='pdf', EdgeColor=getColor(iClu, 7, 0.7), DisplayName=p.mi.semanticClusterLabels(iClu), DisplayStyle='stairs');
    end
    histogram(ax(iRow, iCol), X, [0:0.1:5, Inf], Normalization='pdf', EdgeColor='k', DisplayName='all', DisplayStyle='stairs', LineWidth=2)
    title(ax(iRow, iCol), sprintf('Observed (real %ss)', dir))
end

lgd = legend(ax(2, 2));
lgd.Layout.Tile = 'east';
xlabel(tl, 'max rectified movement index')
ylabel(tl, 'pdf')
