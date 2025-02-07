p.bootAlpha = 0.01;
%% 3. Compare ETA (binned by trial length) (euclidean distance, bootstrap). To show there is no baseline/ramp slope/peak differences
p.binnedTrialEdges = [1, 2, 4, 6, 10];
p.startBlankWindow = [0, 0.5];
[bta.pressUpRaw.X, bta.pressUpRaw.T, bta.pressUpRaw.N, bta.pressUpRaw.S, bta.pressUpRaw.B] = eu(c.isPressUp).getBinnedTrialAverage('count', p.binnedTrialEdges, 'press', window=[-10, 1], normalize=false, resolution=0.100, startBlankWindow=[0, 0.5]);
[bta.pressDownRaw.X, bta.pressDownRaw.T, bta.pressDownRaw.N, bta.pressDownRaw.S, bta.pressDownRaw.B] = eu(c.isPressDown).getBinnedTrialAverage('count', p.binnedTrialEdges, 'press', window=[-10, 1], normalize=false, resolution=0.100, startBlankWindow=[0, 0.5]);
[bta.lickUpRaw.X, bta.lickUpRaw.T, bta.lickUpRaw.N, bta.lickUpRaw.S, bta.lickUpRaw.B] = eu(c.isLickUp).getBinnedTrialAverage('count', p.binnedTrialEdges, 'lick', window=[-10, 1], normalize=false, resolution=0.100, startBlankWindow=[0, 0.5]);
[bta.lickDownRaw.X, bta.lickDownRaw.T, bta.lickDownRaw.N, bta.lickDownRaw.S, bta.lickDownRaw.B] = eu(c.isLickDown).getBinnedTrialAverage('count', p.binnedTrialEdges, 'lick', window=[-10, 1], normalize=false, resolution=0.100, startBlankWindow=[0, 0.5]);
bta.pressUpRaw.X = bta.pressUpRaw.X./0.1;
bta.pressUpRaw.S = bta.pressUpRaw.S./0.1;
bta.pressDownRaw.X = bta.pressDownRaw.X./0.1;
bta.pressDownRaw.S = bta.pressDownRaw.S./0.1;
bta.lickUpRaw.X = bta.lickUpRaw.X./0.1;
bta.lickUpRaw.S = bta.lickUpRaw.S./0.1;
bta.lickDownRaw.X = bta.lickDownRaw.X./0.1;
bta.lickDownRaw.S = bta.lickDownRaw.S./0.1;


% [btaSmooth.pressUpRaw.X, btaSmooth.pressUpRaw.T, btaSmooth.pressUpRaw.N, btaSmooth.pressUpRaw.S, btaSmooth.pressUpRaw.B] = eu(c.isPressUp).getBinnedTrialAverage('rate', p.binnedTrialEdges, 'press', window=[-4, 0], normalize=false, startBlankWindow=[0, 0.5]);
% [btaSmooth.pressDownRaw.X, btaSmooth.pressDownRaw.T, btaSmooth.pressDownRaw.N, btaSmooth.pressDownRaw.S, btaSmooth.pressDownRaw.B] = eu(c.isPressDown).getBinnedTrialAverage('rate', p.binnedTrialEdges, 'press', window=[-4, 0], normalize=false, startBlankWindow=[0, 0.5]);

%% Bootstrap
p.binnedTrialEdgesFine = 1:1:10;

bootBTA = bootstrapBTAMono(10000, eu, c.isPressResponsive, alpha=p.bootAlpha, trialType='press', binEdges=p.binnedTrialEdgesFine, distWindow=[-2, -0.2], startBlankWindow=[0, 0.5]);
c.isPressBTADifferentUp = reshape(bootBTA.distObs > bootBTA.distCI(:, 2), 1, []);
c.isPressBTADifferentDown = reshape(bootBTA.distObs < bootBTA.distCI(:, 1), 1, []);
c.isPressBTADifferent = c.isPressBTADifferentUp | c.isPressBTADifferentDown;
fprintf(1, '\n\nOf %d press responsive units, %d (%d up, %d down) showed significantly different responses for different length trials (p<0.01).\n', nnz(c.isPressResponsive), nnz(c.isPressBTADifferent), nnz(c.isPressBTADifferentUp), nnz(c.isPressBTADifferentDown))

% Make individual figures for the significantly different units

[btaSig.X, btaSig.T, btaSig.N, btaSig.S, btaSig.B] = eu(c.isPressResponsive & c.isPressBTADifferent).getBinnedTrialAverage('count', p.binnedTrialEdgesFine, 'press', window=[-10, 1], normalize=false, resolution=0.100, startBlankWindow=[0, 0.5]);
[btaNul.X, btaNul.T, btaNul.N, btaNul.S, btaNul.B] = eu(c.isPressResponsive & ~c.isPressBTADifferent).getBinnedTrialAverage('count', p.binnedTrialEdgesFine, 'press', window=[-10, 1], normalize=false, resolution=0.100, startBlankWindow=[0, 0.5]);
btaSig.X = btaSig.X./0.1;
btaSig.S = btaSig.S./0.1;
btaNul.X = btaNul.X./0.1;
btaNul.S = btaNul.S./0.1;

[btaUp.X, btaUp.T, btaUp.N, btaUp.S, btaUp.B] = eu(c.isPressUp & c.isPressBTADifferentUp).getBinnedTrialAverage('count', p.binnedTrialEdgesFine, 'press', window=[-10, 1], normalize=false, resolution=0.100, startBlankWindow=[0, 0.5]);
[btaDown.X, btaDown.T, btaDown.N, btaDown.S, btaDown.B] = eu(c.isPressDown & c.isPressBTADifferentDown).getBinnedTrialAverage('count', p.binnedTrialEdgesFine, 'press', window=[-10, 1], normalize=false, resolution=0.100, startBlankWindow=[0, 0.5]);
btaUp.X = btaUp.X./0.1;
btaUp.S = btaUp.S./0.1;
btaDown.X = btaDown.X./0.1;
btaDown.S = btaDown.S./0.1;

[btaNulUp.X, btaNulUp.T, btaNulUp.N, btaNulUp.S, btaNulUp.B] = eu(c.isPressUp & ~c.isPressBTADifferentUp).getBinnedTrialAverage('count', p.binnedTrialEdgesFine, 'press', window=[-10, 1], normalize=false, resolution=0.100, startBlankWindow=[0, 0.5]);
[btaNulDown.X, btaNulDown.T, btaNulDown.N, btaNulDown.S, btaNulDown.B] = eu(c.isPressDown & ~c.isPressBTADifferentDown).getBinnedTrialAverage('count', p.binnedTrialEdgesFine, 'press', window=[-10, 1], normalize=false, resolution=0.100, startBlankWindow=[0, 0.5]);
btaNulUp.X = btaNulUp.X./0.1;
btaNulUp.S = btaNulUp.S./0.1;
btaNulDown.X = btaNulDown.X./0.1;
btaNulDown.S = btaNulDown.S./0.1;

sel = find(c.isPressResponsive & ~c.isPressBTADifferent);
sel = sel(randi(nnz(sel), [1, nnz(c.isPressResponsive & c.isPressBTADifferent)]));
[btaNulSub.X, btaNulSub.T, btaNulSub.N, btaNulSub.S, btaNulSub.B] = eu(sel).getBinnedTrialAverage('count', p.binnedTrialEdgesFine, 'press', window=[-10, 1], normalize=false, resolution=0.100, startBlankWindow=[0, 0.5]);
btaNulSub.X = btaNulSub.X./0.1;
btaNulSub.S = btaNulSub.S./0.1;
fprintf(1, '\n\nOf %d press responsive units, %d (%d up, %d down) showed significantly different responses for different length trials (p<0.01).\n', nnz(c.isPressResponsive), nnz(c.isPressBTADifferent), nnz(c.isPressBTADifferentUp), nnz(c.isPressBTADifferentDown))

%% Bootstrap (lick)
% p.binnedTrialEdgesFine = 1:1:10;

bootBTALick = bootstrapBTAMono(10000, eu, c.isLickResponsive & c.hasPress, alpha=p.bootAlpha, trialType='lick', binEdges=p.binnedTrialEdgesFine, distWindow=[-2, 0], startBlankWindow=[0, 0.5]);
c.isLickBTADifferentUp = reshape(bootBTALick.distObs > bootBTALick.distCI(:, 2), 1, []);
c.isLickBTADifferentDown = reshape(bootBTALick.distObs < bootBTALick.distCI(:, 1), 1, []);
c.isLickBTADifferent = c.isLickBTADifferentUp | c.isLickBTADifferentDown;

% Make individual figures for the significantly different units

[btaLickSig.X, btaLickSig.T, btaLickSig.N, btaLickSig.S, btaLickSig.B] = eu(c.hasPress & c.isLickResponsive & c.isLickBTADifferent).getBinnedTrialAverage('count', p.binnedTrialEdgesFine, 'lick', window=[-10, 1], normalize=false, resolution=0.100, startBlankWindow=[0, 0.5]);
[btaLickNul.X, btaLickNul.T, btaLickNul.N, btaLickNul.S, btaLickNul.B] = eu(c.hasPress & c.isLickResponsive & ~c.isLickBTADifferent).getBinnedTrialAverage('count', p.binnedTrialEdgesFine, 'lick', window=[-10, 1], normalize=false, resolution=0.100, startBlankWindow=[0, 0.5]);
btaLickSig.X = btaLickSig.X./0.1;
btaLickSig.S = btaLickSig.S./0.1;
btaLickNul.X = btaLickNul.X./0.1;
btaLickNul.S = btaLickNul.S./0.1;

[btaLickUp.X, btaLickUp.T, btaLickUp.N, btaLickUp.S, btaLickUp.B] = eu(c.hasPress & c.isLickUp & c.isLickBTADifferentUp).getBinnedTrialAverage('count', p.binnedTrialEdgesFine, 'lick', window=[-10, 1], normalize=false, resolution=0.100, startBlankWindow=[0, 0.5]);
[btaLickDown.X, btaLickDown.T, btaLickDown.N, btaLickDown.S, btaLickDown.B] = eu(c.hasPress & c.isLickDown & c.isLickBTADifferentDown).getBinnedTrialAverage('count', p.binnedTrialEdgesFine, 'lick', window=[-10, 1], normalize=false, resolution=0.100, startBlankWindow=[0, 0.5]);
btaLickUp.X = btaLickUp.X./0.1;
btaLickUp.S = btaLickUp.S./0.1;
btaLickDown.X = btaLickDown.X./0.1;
btaLickDown.S = btaLickDown.S./0.1;

[btaLickNulUp.X, btaLickNulUp.T, btaLickNulUp.N, btaLickNulUp.S, btaLickNulUp.B] = eu(c.hasPress & c.isLickUp & ~c.isLickBTADifferentUp).getBinnedTrialAverage('count', p.binnedTrialEdgesFine, 'lick', window=[-10, 1], normalize=false, resolution=0.100, startBlankWindow=[0, 0.5]);
[btaLickNulDown.X, btaLickNulDown.T, btaLickNulDown.N, btaLickNulDown.S, btaLickNulDown.B] = eu(c.hasPress & c.isLickDown & ~c.isLickBTADifferentDown).getBinnedTrialAverage('count', p.binnedTrialEdgesFine, 'lick', window=[-10, 1], normalize=false, resolution=0.100, startBlankWindow=[0, 0.5]);
btaLickNulUp.X = btaLickNulUp.X./0.1;
btaLickNulUp.S = btaLickNulUp.S./0.1;
btaLickNulDown.X = btaLickNulDown.X./0.1;
btaLickNulDown.S = btaLickNulDown.S./0.1;

sel = find(c.hasPress & c.isLickResponsive & ~c.isLickBTADifferent);
sel = sel(randi(nnz(sel), [1, nnz(c.hasPress & c.isLickResponsive & c.isLickBTADifferent)]));
[btaNulSub.X, btaNulSub.T, btaNulSub.N, btaNulSub.S, btaNulSub.B] = eu(sel).getBinnedTrialAverage('count', p.binnedTrialEdgesFine, 'lick', window=[-10, 1], normalize=false, resolution=0.100, startBlankWindow=[0, 0.5]);
btaNulSub.X = btaNulSub.X./0.1;
btaNulSub.S = btaNulSub.S./0.1;
fprintf(1, '\n\nOf %d lick responsive units, %d (%d up, %d down) showed significantly different responses for different length trials (p<0.01).\n', nnz(c.hasPress & c.isLickResponsive), nnz(c.hasPress & c.isLickBTADifferent), nnz(c.hasPress & c.isLickBTADifferentUp), nnz(c.hasPress & c.isLickBTADifferentDown))

%%
save('C:\SERVER\boot_bta_20250204.mat', 'p', 'bta', 'bootBTA', 'btaSig', 'btaNul', 'btaSmooth')

%% Functions
function boot = bootstrapBTAMono(nboot, eu, varargin)
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
    p.addParameter('startBlankWindow', [0, 0.5])

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

    rng(42)

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
        [xx, tt] = eu(iEu).getTrialAlignedData('count', [-10, 0], trialType, allowedTrialDuration=[0, Inf], alignTo='stop', resolution=0.1, includeInvalid=false, startBlankWindow=r.startBlankWindow);
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
        
        pairs = [1:nBins-1; 2:nBins]';
        nPairs = size(pairs, 1);
        dist = NaN(nPairs, 1);
        for iPair = 1:nPairs
            x1 = xxMean(pairs(iPair, 1), :);
            x2 = xxMean(pairs(iPair, 2), :);
            dist(iPair) = mean(x2 - x1, 'omitnan');
        end
        dist = mean(dist, 'omitnan');
        
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
            distBoot(iPair, :) = mean(x2 - x1, 1, 'omitnan');
        end
        distBoot = mean(distBoot, 1, 'omitnan');
        distCI = prctile(distBoot, [50*r.alpha, 100-50*r.alpha]);
        distH = dist >= distCI(2) || dist <= distCI(1);
        boot.distH(iEu) = distH;
        boot.distCI(iEu, 1:2) = distCI;
        boot.distObs(iEu) = dist;

        assert(~isnan(distH))
        assert(~isnan(dist))
        assert(all(~isnan(distCI)))
    end

end
