p.baselineWindow = [-4, -2];
nboot = 100000;
alpha = 0.05;

sr = eta.pressRaw.X([find(c.isPressUp), find(c.isPressDown)], :)./0.1;
bsr = eta.pressRaw.X(:, eta.pressRaw.t >= p.baselineWindow(1) & eta.pressRaw.t <= p.baselineWindow(2))./0.1;
bsrUp = bsr(c.isPressUp, :);
bsrDown = bsr(c.isPressDown, :);
bsr = vertcat(bsrUp, bsrDown);
nUp = size(bsrUp, 1); assert(nUp == nnz(c.isPressUp));
nDown = size(bsrDown, 1); assert(nDown == nnz(c.isPressDown));
nTotal = size(bsr, 1);

% Do permutation test
rng(42)
bsample = zeros(nTotal, nboot);
for iboot = 1:nboot
    bsample(:, iboot) = randperm(nTotal);
end

% First test bsrMean (baseline spike rate, averaged across trial and time)
bsrMean = mean(bsr, 2);
upSamples = bsrMean(bsample(1:nUp, :));
downSamples = bsrMean(bsample(nUp+1:end, :));

muDiffObs = mean(bsrMean(1:nUp)) - mean(bsrMean(nUp+1:end));
muDiffBoot = mean(upSamples, 1) - mean(downSamples, 1);
muDiffCI = quantile(muDiffBoot, [alpha*0.5, 1 - alpha*0.5]);
muDiffP = nnz(muDiffBoot<muDiffObs)./nnz(muDiffBoot);
if muDiffP > 0.5
    muDiffP = 1 - muDiffP;
end
assert(nnz(muDiffObs==muDiffBoot)==0)

boot.pressBaseline.muDiff.obs = muDiffObs;
boot.pressBaseline.muDiff.h = muDiffP > alpha;
boot.pressBaseline.muDiff.ci = muDiffCI;
boot.pressBaseline.muDiff.p = muDiffP;
boot.pressBaseline.muDiff.alpha = alpha;
boot.pressBaseline.muDiff.nboot = nboot;

%% Okay this bits tries to estimated shuffled "error" relative to the shuffled mean.
% We can add this error back onto the observed mean for error bar
% visualization. But that's a bit weird so let's just use SEMs maybe?
boot.pressBaseline.muErrUp.ci = NaN(2, size(sr, 2));
boot.pressBaseline.muErrDown.ci = NaN(2, size(sr, 2));
for it = 1:size(sr, 2)
    upSamples = sr(:, it);
    downSamples = sr(:, it);
    upSamples = mean(upSamples(bsample(1:nUp, :)), 1) - mean(mean(upSamples(bsample(1:nUp, :)), 1));
    downSamples = mean(downSamples(bsample(nUp+1:end, :)), 1) - mean(mean(downSamples(bsample(nUp+1:end, :)), 1));
    boot.pressBaseline.muErrUp.ci(:, it) = quantile(upSamples, [alpha*0.5, 1-alpha*0.5]);
    boot.pressBaseline.muErrDown.ci(:, it) = quantile(downSamples, [alpha*0.5, 1-alpha*0.5]);
end