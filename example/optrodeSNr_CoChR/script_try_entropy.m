script_loadPreprocessedData_optrodeSNr_CoChR_VGATCRE_20260729;

%%
clear entropy
entropy(length(exp)) = struct(ctrl=[], stim=[]);
fig = figure();
tl = tiledlayout(fig, length(exp), 1, TileSpacing='tight', Padding='tight');
for iExp = 1:length(exp)
    ax = nexttile(tl);
    hold(ax, 'on')
    for cond = ["ctrl", "stim"]
        t = stimData(iExp).psth.(cond).t;
        X = stimData(iExp).psth.(cond).X;
        if cond == "ctrl"
            X = X(stimData(iExp).psmh.ctrl.press.isRewarded, :, :);
        end
        edges = -2:0.1:5;
        H = getEntropyForPC1(X, t, edges);
        entropy(iExp).(cond).H = H;
        entropy(iExp).(cond).edges = edges;
        switch cond
            case "ctrl"
                color = 'k';
            case "stim"
                color = 'r';
        end
        % histogram(ax, BinEdges=edges, BinCounts=H, DisplayName=cond, DisplayStyle='stairs', EdgeColor=color);
        x = [edges(1:end-1); edges(2:end)];
        x = x(:);
        y = [H'; H'];
        y = y(:);
        plot(ax, x, y, DisplayName=cond, Color=color);
        % ylim(ax, [4, 6.5])
    end
    xticks(ax, [0, 1, 2, 3])
end
lgd = legend(ax, Orientation='horizontal');
lgd.Layout.Tile = 'north';
xlabel(tl, 'time to opto/decoder onset (s)')
ylabel(tl, 'entropy (bits)')

function H = getEntropyForPC1(X, t, edges)
    assert(size(X, 2)==length(t)) % X: trials x time x units
    nUnits = size(X, 3);
    XMerge = reshape(X, [], nUnits); % Merge trials, time => samples x units
    unitsToRemove = any(isnan(XMerge), 1);
    XMerge(:, unitsToRemove) = [];
    X(:, :, unitsToRemove) = [];
    nUnits = size(X, 3);
    [coeff, score] = pca(XMerge);
    minScore = min(score(:, 1));
    maxScore = max(score(:, 1));
    mu = mean(XMerge, 1);
    % score = (Xmerge - mu) * coeff;
    H = NaN(length(edges)-1, 1);
    for iBin = 1:length(edges)-1
        [a, b] = isin(t, edges(iBin:iBin+1), false, true);
        Xi = reshape(X(:, a:b, :), [], nUnits);
        score = (Xi - mu) * coeff; % apply PCA transform
        N = histcounts(score(:, 1), linspace(minScore, maxScore, 21));
        P = (N+1)./sum(N);
        H(iBin) = -sum(P.*log2(P));
    end
end