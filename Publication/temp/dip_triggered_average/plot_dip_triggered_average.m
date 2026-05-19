
%% Load data
load('C:\SERVER\LickVsReach_DTA_RTA_boot\NewData\LickVsReach_DLC_dta_rta_25_25_200to800ms_units1to1443_100boots.mat');

%% Make a metaDTA/metaRTA
for fn = ["spikerate", "HandR", "HandL", "FootR", "FootL", "Spine", "Jaw", "Tongue"]
    X = arrayfun(@(xta) xta.(fn).X, xta.dip(selUnits), UniformOutput=false);
    X = cat(1, X{:});
    xta.dip(length(eu) + 1).(fn) = struct(X=mean(X, 1, 'omitnan'), t=xta.dip(1).(fn).t);
    X = arrayfun(@(xta) xta.(fn).X, xta.rise(selUnits), UniformOutput=false);
    X = cat(1, X{:});
    xta.rise(length(eu) + 1).(fn) = struct(X=mean(X, 1, 'omitnan'), t=xta.rise(1).(fn).t);
    clear X
end
xta.dip(length(eu) + 1).iExp = 0;
xta.rise(length(eu) + 1).iExp = 0;
xta.dip(length(eu) + 1).t0 = [xta.dip(selUnits).t0];
xta.rise(length(eu) + 1).t0 = [xta.rise(selUnits).t0];

%% Combine movement ranges (mr) across dips from all units, then cluster them
statFeatures = ["HandR", "HandL", "FootR", "FootL", "Spine", "Jaw", "Tongue"];


%% Plot dip-triggered average kinematics (STD Version
close all
exportPath = fullfile("C:\SERVER\LickVsReach_DTA_RTA_boot\NewData\Figures", sprintf("LickVsReach_DLC_dta_rta_%i_%i_%ito%ims_std", 100*p.xta.dip.thresholdQuantile, 100*p.xta.dip.thresholdSubQuantile, 100*p.xta.dip.samples(1), 100*p.xta.rise.samples(2)));
if ~exist(exportPath, 'dir')
    mkdir(exportPath)
end
features = ["spikerate", "HandR", "HandL", "FootR", "FootL", "Spine", "Jaw", "Tongue"];
featureUnits = ["spike rate (a.u.)", "AP pos (a.u.)", "AP pos (a.u.)", "AP pos (a.u.)", "AP pos (a.u.)", "DV pos (a.u.)", "DV pos (a.u.)", "prob"];
statFeatures = ["HandR", "HandL", "FootR", "FootL", "Spine", "Jaw", "Tongue"];
nAx = length(features) + 2;
dirs = ["dip", "rise"];
dimensionReduction = "tsne"; % pca, tsne, umap
% yl = {[-2.5, 5], [-0.1, 5.1], [-0.1, 5.1], [-0.1, 5.1], [-0.1, 5.1], [-0.1, 5.1], [-0.1, 5.1], [-0.1, 1.1]};
yl = {[-2.5, 5], [-2, 2], [-2, 2], [-2, 2], [-2, 2], [-2, 2], [-2, 2], [0, 1]};
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
for iUnit = 1:2%length(xta.dip)
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
        
        if p.nBoot > 0 && xta.(dir)(iUnit).iExp > 0
            iAx = nAx - 1;
            if length(xta.(dir)(iUnit).t0) <= 2
                idx = ones(length(xta.(dir)(iUnit).t0), 1);
                nClusters = 1;
            else
                hold(ax(iDir, iAx), 'on')
                mr = NaN(length(xta.(dir)(iUnit).t0), length(p.std.features));% movement range: nTrials x nFeatures
                for i = 1:length(p.std.features)
                    fn = p.std.features(statFeatureOrder(i));
                    if isempty(xta.(dir)(iUnit).(fn))
                        continue
                    end
                    t = xta.(dir)(iUnit).(fn).t;
                    selT = t >= p.std.window(1) & t <= p.std.window(2);
                    xx = xta.(dir)(iUnit).(fn).X(:, selT);
                    mr(:, i) = std(xx, 0, 2, 'omitnan');
                end
                mr(isnan(mr)) = 0;
                [~, pcScore, ~, ~, explained] = pca(mr);
                switch dimensionReduction   
                    case "pca"
                        score = pcScore;
                    case "umap"
                        score = umap(mr, NumDimensions=2);
                    case "tsne"
                        score = tsne(mr);
                end
                eva = evalclusters(pcScore(:, 1:3), 'kmeans', 'CalinskiHarabasz', KList=1:3);
                nClusters = eva.OptimalK;
                clear eva
                idx = kmeans(pcScore(:, 1:3), nClusters);
                for k = 1:nClusters
                    sel = idx==k;
                    scatter(ax(iDir, iAx), score(sel, 1), score(sel, 2), 10, getColor(k, 2))
                end
                switch dimensionReduction   
                    case "pca"
                        xlabel(ax(iDir, iAx), sprintf("PC%i (%.1f%%)", 1, explained(1)))
                        ylabel(ax(iDir, iAx), sprintf("PC%i (%.1f%%)", 2, explained(2)))
                        title(ax(iDir, iAx), "PCA")
                    case "tsne"
                        xlabel(ax(iDir, iAx), sprintf("PC%i", 1))
                        ylabel(ax(iDir, iAx), sprintf("PC%i", 2))
                        title(ax(iDir, iAx), "t-SNE")
                    case "umap"
                        xlabel(ax(iDir, iAx), sprintf("PC%i", 1))
                        ylabel(ax(iDir, iAx), sprintf("PC%i", 2))
                        title(ax(iDir, iAx), "UMAP")
                end
                xticks(ax(iDir, iAx), [])
                yticks(ax(iDir, iAx), [])
            end
        end

        % Movement diversity matrix
        if p.nBoot > 0 && xta.(dir)(iUnit).iExp > 0
            iAx = nAx;
            mdm = NaN(length(xta.(dir)(iUnit).t0), length(p.std.features));
            for i = 1:length(p.std.features)
                fn = p.std.features(statFeatureOrder(i));
                if isempty(xta.(dir)(iUnit).(fn))
                    continue
                end
                t = xta.(dir)(iUnit).(fn).t;
                selT = t >= p.std.window(1) & t <= p.std.window(2);
                stdObs = std(xta.(dir)(iUnit).(fn).X(:, selT), 0, 2, 'omitnan');
                mdm(:, i) = arrayfun(@(data) nnz(xta.(dir)(iUnit).(fn).stats.stdBoot < data) ./ length(xta.(dir)(iUnit).(fn).stats.stdBoot), stdObs, UniformOutput=true);
            end
            mdm(isnan(mdm)) = 0;
            hash = sum((mdm > 0.95) .* 2.^(size(mdm, 2)-1:-1:0), 2);
            hash = hash + (idx-1) .* 2.^(size(mdm, 2));
            [~, I] = sort(hash, 'ascend');
            idxSorted = idx(I);
            sepHash = arrayfun(@(idx) find(idxSorted==idx, 1, 'last'), 1:max(idx)-1);
            imagesc(ax(iDir, iAx), mdm(I, :))
            if ~isempty(sepHash)
                yline(ax(iDir, iAx), 0.5+sepHash, 'k--')
                yticks(ax(iDir, iAx), 0.5+unique([1, sepHash, length(idx)]))
                yticklabels(ax(iDir, iAx), string(unique([1, sepHash, length(idx)])))
            end
            applyCustomColormap(ax(iDir, iAx), [-1, 1], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.05, h0=0.33);
            xticks(ax(iDir, iAx), 1:length(p.std.features))
            xticklabels(ax(iDir, iAx), p.std.features(statFeatureOrder))
            colorbar(ax(iDir, iAx), Orientation='horizontal', Location='southoutside')
            ax(iDir, iAx).XAxisLocation = 'top';
            ylabel(ax(iDir, iAx), 'Trial')
            clear mdm i fn t selT stdObs pObs hash I
        end
        
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
            hold(ax(2+1-iDir, iAx), 'on')
            t = 1e3*xta.(dir)(iUnit).(fn).t;
            X = mean(xta.(dir)(iUnit).(fn).X, 1, 'omitnan');
            plot(ax(iDir, iAx), t, X, Color=[0.15, 0.15, 0.15, 1], LineWidth=1.5, LineStyle=':');
            plot(ax(2+1-iDir, iAx), t, X, Color=[0.15, 0.15, 0.15, 0.1], LineWidth=1.5, LineStyle=':');
            for k = 1:nClusters
                % c = getColor(iAx, length(features), 0.7);
                c = getColor(k, nClusters, 0.7);
                % c = [0.15, 0.15, 0.15];
    
                X = xta.(dir)(iUnit).(fn).X(idx==k, :);
                mu = mean(X, 1, 'omitnan');
                err = std(X, 0, 1, 'omitnan')./sqrt(size(xta.(dir)(iUnit).(fn).X, 1));
    
                plot(ax(iDir, iAx), t, mu, Color=c, LineWidth=1.5);
            end
            if p.nBoot > 0 && isfield(xta.(dir)(iUnit).(fn), 'XBoot')
                prc = quantile(xta.(dir)(iUnit).(fn).XBoot, [p.bootAlpha/2, 1-p.bootAlpha/2], 1);
                patch(ax(iDir, iAx), [t, flip(t)], [prc(1, :), flip(prc(2, :))], [0.15, 0.15, 0.15], FaceAlpha=0.05, EdgeColor=[0.15, 0.15, 0.15], EdgeAlpha=0.5);
            end
            % xline(ax(iRow, iAx), 1e3*p.xta.meanWindow, 'k--', Alpha=0.1)
            xline(ax(iDir, iAx), 1e3*p.std.window, 'k--', Alpha=0.1)
            xline(ax(iDir, iAx), 0, 'k-', Alpha=0.1)
            % xticks(ax(iRow, iAx), 1e3*p.xta.meanWindow)
            xticks(ax(iDir, iAx), [-300, 0, 600])
            xtickangle(ax(iDir, iAx), 0)

            xlabel(ax(iDir, iAx), 'time (ms)')
            ylabel(ax(iDir, iAx), featureUnits(iAx))

            fnDisp = sprintf("%s %s", fn, repmat('*', [1, nStarsSTD]));
            title(ax(iDir, iAx), fnDisp, Interpreter='none')
            ylim(ax(iDir, iAx), yl{iAx})
            hold(ax(iDir, iAx), 'off')
            hold(ax(2+1-iDir, iAx), 'off')
        end

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

clear exportPath features featureUnits statFeatures nAx dirs dimensionReduction yl fig tlp tl ax iDir iAx iUnit dir
clear lia statFeatureOrder idx nClusters mr i fn t selT xx pcScore explained score eva k sel
clear mdm i fn t selT stdObs hash I idxSorted sepHash prcSTD nStarsSTD t X k c mu err prc
clear iDir dir