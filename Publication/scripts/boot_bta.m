p.bootAlpha = 0.01;
%% 3. Compare ETA (binned by trial length) (euclidean distance, bootstrap). To show there is no baseline/ramp slope/peak differences

[bta.pressUpRaw.X, bta.pressUpRaw.T, bta.pressUpRaw.N, bta.pressUpRaw.S, bta.pressUpRaw.B] = eu(c.isPressUp).getBinnedTrialAverage('rate', p.binnedTrialEdges, 'press', 'window', [-10, 1], 'normalize', false);
[bta.pressDownRaw.X, bta.pressDownRaw.T, bta.pressDownRaw.N, bta.pressDownRaw.S, bta.pressDownRaw.B] = eu(c.isPressDown).getBinnedTrialAverage('rate', p.binnedTrialEdges, 'press', 'window', [-10, 1], 'normalize', false);

[bta.lickUpRaw.X, bta.lickUpRaw.T, bta.lickUpRaw.N, bta.lickUpRaw.S, bta.lickUpRaw.B] = eu(c.isLickUp).getBinnedTrialAverage('rate', p.binnedTrialEdges, 'lick', 'window', [-10, 1], 'normalize', false);
[bta.lickDownRaw.X, bta.lickDownRaw.T, bta.lickDownRaw.N, bta.lickDownRaw.S, bta.lickDownRaw.B] = eu(c.isLickDown).getBinnedTrialAverage('rate', p.binnedTrialEdges, 'lick', 'window', [-10, 1], 'normalize', false);

%
binEdges = 2:2:10;
bootBTA = bootstrapBTA(10000, eu, c.isPressResponsive, alpha=p.bootAlpha, trialType='press', binEdges=binEdges, distWindow=[-2, 0]);
c.isPressBTADifferent = reshape(bootBTA.distH == 1, 1, []);
fprintf(1, '\n\nOf %d press responsive units, %d showed significantly different responses for different length trials (p<0.01).\n', nnz(c.isPressResponsive), nnz(c.isPressBTADifferent))


%% Functions
function boot = bootstrapBTA(nboot, eu, varargin)
    p = inputParser();
    p.addRequired('nboot', @isnumeric)
    p.addRequired('eu', @(x) isa(x, 'EphysUnit'))
    p.addOptional('sel', [], @(x) islogical(x) || isnumeric(x))
    p.addParameter('trialType', 'press', @(x) ismember(x, {'press', 'lick'}))
    p.addParameter('binEdges', [2, 4, 6, 10]);
%     p.addParameter('metric', 'dist', @(x) ismember(x, {'dist', 'point'}))
    p.addParameter('distWindow', [-2, 0], @(x) isnumeric(x) && length(x) == 2)
%     p.addParameter('pointTimestamp', -1, @isnumeric)
    p.addParameter('alpha', 0.05, @isnumeric)

    p.parse(nboot, eu, varargin{:})
    r = p.Results;
    nboot = r.nboot;
    eu = r.eu;
    trialType = r.trialType;
    binEdges = r.binEdges;
    nBins = length(binEdges) - 1;

    if isempty(r.sel)
        euIndices = 1:length(eu);
    elseif islogical(r.sel)
        euIndices = reshape(find(r.sel), 1, []);
    else
        euIndices = reshape(sel, 1, []);
    end

    ii = 0;
    boot.distH = NaN(length(eu), 1);
    boot.distCI = NaN(length(eu), 2);
    boot.distObs = NaN(length(eu), 1);
    for iEu = euIndices
        ii = ii + 1;
        fprintf(1, '%d/%d ', ii, length(euIndices))
        if mod(ii, 10) == 0
            fprintf(1, '\n')
        end
        [xx, tt] = eu(iEu).getTrialAlignedData('count', [-10, 0], trialType, allowedTrialDuration=[0, Inf], alignTo='stop', resolution=0.1, includeInvalid=false);
        dd = eu(iEu).getTrials(trialType).duration;
        assert(length(dd) == size(xx, 1))

        selTrials = dd >= binEdges(1) & dd <= binEdges(end);
        selTime = tt >= r.distWindow(1) & tt <= r.distWindow(2);
%         [~, selPoint] = min(abs(tt - r.pointTimestamp));
        xx = xx(selTrials, selTime);
        dd = dd(selTrials);
        tt = tt(selTime);

        [N, ~, bins] = histcounts(dd, binEdges);
        xxMean = NaN(nBins, nnz(selTime));
        for iBin = 1:length(binEdges) - 1
            xxMean(iBin, :) = mean(xx(bins == iBin, :), 1, 'omitnan');
        end
        
        pairs = nchoosek(1:nBins, 2);
        nPairs = size(pairs, 1);
        dist = NaN(nPairs, 1);
        for iPair = 1:nPairs
            x1 = xxMean(pairs(iPair, 1), :);
            x2 = xxMean(pairs(iPair, 2), :);
            sqrs = (x1 - x2).^2;
            dist(iPair) = mean(sqrs, 'omitnan');
        end
        dist = mean(dist);
        
        bsample = zeros(sum(N), nboot);
        for iboot = 1:nboot
            bsample(:, iboot) = randperm(sum(N));
        end
        
        bootBinEdges = [0, cumsum(N)];
        xxMeanBoot = NaN(nBins, nnz(selTime), nboot);
        for iBin = 1:length(binEdges) - 1
            for iboot = 1:nboot
                selTrialsBoot = bsample(bootBinEdges(iBin)+1:bootBinEdges(iBin+1), iboot);
                xxMeanBoot(iBin, :, iboot) = mean(xx(selTrialsBoot, :), 1, 'omitnan');
            end
        end
        
        distBoot = NaN(nPairs, nboot);
        for iPair = 1:nPairs
            x1 = squeeze(xxMeanBoot(pairs(iPair, 1), :, :));
            x2 = squeeze(xxMeanBoot(pairs(iPair, 2), :, :));
            sqrs = (x1 - x2).^2;
            distBoot(iPair, :) = mean(sqrs, 1, 'omitnan');
        end
        distBoot = mean(distBoot, 1, 'omitnan');
        distCI = prctile(distBoot, [50*r.alpha, 100-50*r.alpha]);
        distH = dist >= distCI(2) || dist <= distCI(1);
        boot.distH(iEu) = distH;
        boot.distCI(iEu, 1:2) = distCI;
        boot.distObs(iEu) = dist;
    end

end
