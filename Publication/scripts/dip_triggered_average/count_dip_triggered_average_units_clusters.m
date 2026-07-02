
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
% p.mi.mustNotMove = {["Jaw", "HandL", "HandR"], ["HandL", "HandR"], ["HandL", "HandR"], ["Jaw", "HandR"], ["Jaw", "HandR"], ["Jaw", "HandL"], ["Jaw", "HandL"]};

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

%% Do the same for sdBoot, these clusters are in native order, not semantic
nUnits = length(xta.dip);
nClusters = length(p.mi.semanticClusterLabels);
nFeatures = length(p.mi.boot.features);

for iUnit = 1:nUnits
    for dir = ["dip", "rise"]
        n = sdBoot(iUnit).(dir).n;
        n(n < p.mi.minNumTrialsPerCluster) = 0;
        sdBoot(iUnit).(dir).n = n;
        clean = true(p.sd.nBoot, nClusters);
        for iClu = 1:nClusters
            k = cc.params(iClu).idx;
            iMustMove = find(ismember(cc.featureNames, cc.params(iClu).mustMove));
            iMustNotMove = find(ismember(cc.featureNames, cc.params(iClu).mustNotMove));
            expectedSign = cc.params(iClu).expectedSign;
            hTemp = sdBoot(iUnit).(dir).h(:, k, :);
            for iBoot = 1:p.sd.nBoot % THERE'S NO TIME TO OPTIMIZE THIS! % Okay I optimize a bit
                % Cluster of interest must move
                if ~isempty(iMustMove) && hTemp(iBoot, :, iMustMove)~=expectedSign
                    % assert(isscalar(iMustMove))
                    clean(iBoot, k) = false;
                    continue
                end
                % Other clusters must not move
                if any(hTemp(iBoot, :, iMustNotMove))
                    clean(iBoot, k) = false;
                end
            end
            clear hTemp
        end
        sdBoot(iUnit).(dir).clean = clean;
        clear iClu clean iMustMove iMustNotMove n
    end
    clear dir
end
clear iUnit nUnits nFeatures nClusters iBoot k expectedSign

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
clear iTest dir iGrp testNames testGroups selClusters hasData found n nGroups

% %% Plotting - 1 figrue per test
% close all
% for iTest = 1:length(tests)
%     fig = figure(Name=tests(iTest).name, Units='inches', Position=[1, 1, 5, 7]);
%     tl = tiledlayout(fig, 4, 2, TileIndexing='columnmajor');
%     title(tl, sprintf("%s\n[%s]", tests(iTest).name, strjoin(cellfun(@(labels) strjoin(labels, "/"), tests(iTest).labels), "] vs. [")))
%     ax = gobjects(4, 2);
%     iCol = 0;
%     for dir = ["dip", "rise"]
%         iCol = iCol + 1;
%         iRow = 0;
%         for fn = ["numGroupsClean", "numGroupsPresent", "prcGroupsClean"]
%             iRow = iRow + 1;
%             ax(iRow, iCol) = nexttile(tl);
%             edges = 0:length(tests(iTest).labels)+1;% Each bin includes the leading edge, but does not include the trailing edge, except for the last bin which includes both edges.
%             if fn == "prcGroupsClean"
%                 edges = edges./edges(end);
%             end
%             n = histcounts(tests(iTest).(dir).(fn), edges);
%             % n = histcounts(tests(iTest).(dir).(fn)(tests(iTest).(dir).valid), edges);
%             if fn ~= "prcGroupsClean"
%                 bar(ax(iRow, iCol), edges(1:end-1), n, 1, EdgeColor='black', FaceColor='black', FaceAlpha=0.5)
%                 xlim(ax(iRow, iCol), [edges(1), edges(end)]-0.5)
%             else
%                 histogram(ax(iRow, iCol), BinEdges=edges, BinCounts=n, EdgeColor='black', FaceColor='black', FaceAlpha=0.5)
%                 xlim(ax(iRow, iCol), [-0.1, 1.1])
%             end
%             xlabel(ax(iRow, iCol), fn)
%             ylabel(ax(iRow, iCol), "no. units")
%             if iRow == 1
%                 title(ax(iRow, iCol), dir)
%             end
%         end
%         iRow = iRow + 1;
%         ax(iRow, iCol) = nexttile(tl);
%         [n, xEdges, yEdges] = histcounts2(tests(iTest).(dir).numGroupsPresent, tests(iTest).(dir).numGroupsClean);
%         n = n./nUnits;
%         histogram2(ax(iRow, iCol), XBinEdges=xEdges, YBinEdges=yEdges, BinCounts=n, ...
%             DisplayStyle='tile', ShowEmptyBins=true)
%         xlabel(ax(iRow, iCol), "numGroupsPresent")
%         ylabel(ax(iRow, iCol), "numGroupsClean")
%         applyCustomColormap(ax(iRow, iCol), [0, 0.25], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);
%         % colormap(ax(iRow, iCol), 'sky')
%         colorbar(ax(iRow, iCol), 'eastoutside')
%         axis(ax(iRow, iCol), 'equal')
%         xticks(ax(iRow, iCol), xEdges+0.5)
%         yticks(ax(iRow, iCol), yEdges+0.5)
%         grid(ax(iRow, iCol), 'off');
%     end
%     fontsize(fig, 9, 'points')
% end
% clear iTest fig tl ax iCol dir iRow fn centers edges n xEdges yEdges
% 
% % Do rises contain more movement types than dips?
% fig = figure(Units='inches', Position=[1 1 4 8]);
% tl = tiledlayout(fig, 3, 2);
% maxC = [0.10, 0.10, 0.15, 0.2, 0.3, 0.3];
% for iTest = 1:length(tests)
%     ax = nexttile(tl);
%     [n, xEdges, yEdges] = histcounts2(tests(iTest).dip.numGroupsClean, tests(iTest).rise.numGroupsClean);
%     n = n./nUnits;
%     histogram2(ax, XBinEdges=xEdges, YBinEdges=yEdges, BinCounts=n, ...
%         DisplayStyle='tile', ShowEmptyBins=true)
%     applyCustomColormap(ax, [0, maxC(iTest)], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.5, h0=0.33);
%     axis(ax, 'equal')
%     title(ax, tests(iTest).name)
%     xlabel(ax, 'dip')
%     ylabel(ax, 'rise')
% 
%     xticks(ax, 0:max(xEdges))
%     yticks(ax, 0:max(yEdges))
% 
%     cb = colorbar(ax, Location='southoutside');
%     % cb.Layout.Tile = 'south';
%     cb.Label.String = sprintf("%% of %i units", nUnits);
%     cb.Ticks = [0, maxC(iTest)];
%     cb.TickLabels = ["0", sprintf("%g%%", maxC(iTest)*100)];
% 
% end
% 
% fontsize(fig, 9, 'points')
% clear fig tl iTest ax n xEdges yEdges cb maxC

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
clear iTest n dir plural

% %% Save results
% exportPath = fullfile("C:\SERVER\LickVsReach_DTA_RTA_boot", sprintf("LickVsReach_DLC_cc_%iunits_%iboots_%s.mat", length(xta.dip), p.sd.nBoot, datetime("now", Format="yyyyMMdd")));
% save(exportPath, 'p', 'tests', 'cc', 'clusterSize', '-v7.3')
% fprintf("Saved to %s\n", exportPath);
