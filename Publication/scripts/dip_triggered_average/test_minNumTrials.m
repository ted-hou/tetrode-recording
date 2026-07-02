close all
for minNumTrials = 5:20%:20
    p.minNumTrialsPerCluster = minNumTrials;
    fig = figure(Units='inches', Position=[3, 3, 6, 6]);
    tl = tiledlayout(fig, 4, 2, TileSpacing = 'tight');
    ax = gobjects(4, 2);
    for REQUIRE_CLEAN_CLUSTERS = [false, true]
        SKIP_TESTS = true;
        SKIP_COUNTING = true;
        if ~SKIP_COUNTING
            count_dip_triggered_average_units_clusters
        end
    
        nClusters = length(p.mi.semanticClusterOrder);
        nFeatures = length(features);
        nUnits = length(xta.dip);

        dirs = ["dip", "rise"];
        colors = [0, 0, 1; 1, 0, 0];
        yl = [0, 100];
        iRow = 1 + REQUIRE_CLEAN_CLUSTERS*2;
        lgd = gobjects(1, 2);
        hPatch = gobjects(2, 2);
        for iCol = 1:2
            dir = dirs(iCol);
            ax(iRow, iCol) = nexttile(tl);
            hold(ax(iRow, iCol), 'on')
        
            x = 1:length(clusterDispName);

            if REQUIRE_CLEAN_CLUSTERS
                cs = arrayfun(@(ccData) ccData.(dir).n(:, 1) .* uint16(ccData.(dir).clean), cc.data, UniformOutput=false);
            else
                cs = arrayfun(@(ccData) ccData.(dir).n(:, 1), cc.data, UniformOutput=false);
            end
            cs = double(cat(2, cs{:})');
            % cs(cs<p.minNumTrialsPerCluster) = NaN;
            cs(cs<1) = NaN;

        
            % Bootstrapped cluster size (null)
            if REQUIRE_CLEAN_CLUSTERS
                csBoot = arrayfun(@(d) d.(dir).n .* uint16(d.(dir).clean), sdBoot, UniformOutput=false);
            else
                csBoot = arrayfun(@(d) d.(dir).n, sdBoot, UniformOutput=false);
            end
            csBoot = double(cat(3, csBoot{:}));
            csBoot = permute(csBoot, [1, 3, 2]); % nBoot x nUnits x nClusters
            % csBoot(csBoot<p.minNumTrialsPerCluster) = NaN;
            csBoot(csBoot<1) = NaN;
            csBoot = csBoot(:, :, p.mi.semanticClusterOrder); % put clusters in semantic order
            % muBoot = reshape(mean(csBoot, [1, 2], 'omitnan'), 1, nClusters);
            % ciBoot = quantile(squeeze(mean(csBoot, 2, 'omitnan')), [0.05/7/2, 1-0.05/7/2], 1);

    
            edges = 0:2:1000;
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
            % lgd(iCol).Position(2) = 0.72;
        end
        % lgd(1).Position(1) = 0.26;
        % lgd(2).Position(1) = 0.72;
    
    
        iRow = 2 + REQUIRE_CLEAN_CLUSTERS*2;
        yl = [0, nUnits];
        for iCol = 1:2
            dir = dirs(iCol);
            ax(iRow, iCol) = nexttile(tl);
            hold(ax(iRow, iCol), 'on')
            x = 1:length(clusterDispName);
            % nUnits with clean cluster
            if REQUIRE_CLEAN_CLUSTERS
                nuwcc = arrayfun(@(ccData) ccData.(dir).n(:, 1) .* uint16(ccData.(dir).clean), cc.data, UniformOutput=false);
            else
                nuwcc = arrayfun(@(ccData) ccData.(dir).n(:, 1), cc.data, UniformOutput=false);
            end
            nuwcc = cat(2, nuwcc{:})';    
            y = sum(nuwcc>=p.minNumTrialsPerCluster, 1);
        
            % nUnits with clean cluster, bootstrapped
            if REQUIRE_CLEAN_CLUSTERS
                nuwccBoot = arrayfun(@(ccData) ccData.(dir).n .* uint16(ccData.(dir).clean), sdBoot, UniformOutput=false);
            else
                nuwccBoot = arrayfun(@(ccData) ccData.(dir).n, sdBoot, UniformOutput=false);
            end
            nuwccBoot = cat(3, nuwccBoot{:});
            nuwccBoot = permute(nuwccBoot, [3, 2, 1]); % nUnits x nClusters x nBoot
            nuwccBoot = nuwccBoot(:, p.mi.semanticClusterOrder, :); % nUnits x nClusters x nBoot, clusters in semantic order    
            yBoot = squeeze(sum(nuwccBoot>=p.minNumTrialsPerCluster, 1)); % nClusters x nBoot
            % muBoot = mean(yBoot, 2, 'omitnan')';
            % ciBoot = quantile(yBoot, [0.05/7/2, 1-0.05/7/2], 2)';
        
            edges = 0:2:nUnits+1;
            for iClu = 2:nClusters
                thisNBoot = histcounts(yBoot(iClu, :), edges);
                xVerticesBoot = reshape([thisNBoot; thisNBoot]./max(thisNBoot), [], 1);
                yVerticesBoot = reshape([edges(1:end-1); edges(2:end)], [], 1);
                patch(ax(iRow, iCol), x(iClu)-xVerticesBoot*0.5, yVerticesBoot, 'k', FaceColor=[0.15, 0.15, 0.15], EdgeColor='none', FaceAlpha=0.5)
                patch(ax(iRow, iCol), x(iClu)+xVerticesBoot*0.5, yVerticesBoot, 'k', FaceColor=[0.15, 0.15, 0.15], EdgeColor='none', FaceAlpha=0.5)
            end
            bar(ax(iRow, iCol), x(2:end), y(2:end), FaceColor=colors(iCol, :), EdgeColor=colors(iCol, :), FaceAlpha=0.5, EdgeAlpha=0.8, BarWidth=0.33);
        
           
            % hBar = bar(ax(iRow, iCol), x, [y; muBoot], FaceAlpha=0.5, Clipping='on');
            % hBar(1).FaceColor = colors(iCol, :);
            % hBar(1).EdgeColor = colors(iCol, :);
            % hBar(2).FaceColor = 'none';
            % hBar(2).EdgeColor = [0.15, 0.15, 0.15];
            % hBar(2).EdgeAlpha = 0.8;
        
            xticks(ax(iRow, iCol), 1:length(clusterDispName))
            xticklabels(ax(iRow, iCol), clusterDispName)
            xtickangle(ax(iRow, iCol), 0)
            % yticks(ax(iRow, iCol), yl)
            ylabel(ax(iRow, iCol), "units")
            xlim(ax(iRow, iCol), [1.3, length(clusterDispName)+0.7])
            ylim(ax(iRow, iCol), yl)
            hold(ax(iRow, iCol), 'off')
        end
    end
    xticks(ax([1, 3], :), [])
    title(ax(1, 1), 'dips (unclean)')
    title(ax(1, 2), 'rises (unclean)')
    title(ax(3, 1), 'dips (clean)')
    title(ax(3, 2), 'rises (clean)')
    xlabel(tl, 'movement cluster')
    title(tl, sprintf('minNumTrialsPerCluster=%i', p.minNumTrialsPerCluster), FontWeight='bold')
    fontsize(tl, 8, 'points')
    fontsize(lgd, 7, 'points');

    print(fig, sprintf('test_minNumTrials_%i.png', minNumTrials), '-dpng', '-r0')
end
