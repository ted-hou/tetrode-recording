%%
clear pSTE
pSTE.kernel = struct(type='exponential', params=struct(lambda1=5, lambda2=10, resolution=0.01, width=1));
pSTE.minTrialLength = 1.5;
pSTE.minNumTrials = 20;
pSTE.cueBlankWindow = 0.5;
pSTE.fitWindow = [-0.5, 0];
pSTE.baselineWindow = [-4, -2];
pSTE.nBoot = 10000;
pSTE.nTrialsPerBin = 4;

close all
fig = figure(Units='inches', Position=[3 3 12 6]);
tl = tiledlayout(fig, 2, 3);
ax = gobjects(2, 3);
for i = 1:2
    for j = 1:3
        ax(i, j) = nexttile(tl);
    end
end

clear obs boot ci95 ci99 pValue

trialType = 'lick';
euIndices = find(c.hasPress & c.hasLick & c.isLickDown);
isUp = false;
for i = 1:length(euIndices)
    for iAx = 1:6
        cla(ax(iAx))
    end

    iEu = euIndices(i);
    
    % Get PETH
    trials = eu(iEu).getTrials(trialType);
    trials = trials(trials.duration() >= pSTE.minTrialLength);
    assert(length(trials) >= pSTE.minNumTrials)
    [X, t] = eu(iEu).getTrialAlignedData('rate', window=[-4, 2], trialType=trialType, resolution=0.01, trials=trials, allowedTrialDuration=[0, Inf], includeInvalid=false, startBlankWindow=[-Inf, pSTE.cueBlankWindow], kernel=pSTE.kernel);
    z = trials.duration();
    
    % Normalize PETH
    inBaselineWindow = t >= pSTE.baselineWindow(1) & t <= pSTE.baselineWindow(2);
    inFitWindow = t >= pSTE.fitWindow(1) & t <= pSTE.fitWindow(2);
    mu = mean(X(:, inBaselineWindow), 'all', 'omitnan');
    sigma = std(X(:, inBaselineWindow), 0, 'all', 'omitnan');
    X = (X - mu) ./ sigma;

    % Have to deal with negative numbers using loglinear fit
    % Here we just take the absolute (which removes the need to reverse the
    % sign for tau for decrease units)
    if ~isUp
        X = -X;
    end
    X(X<=0) = 0.1;

%     X = (X - min(X, [], 2, 'omitnan')) ./ (max(X, [], 2, 'omitnan') - min(X, [], 2, 'omitnan'));
    
    tic;
    [obs, boot, ci95, ci99, pValue] = bootSingleTrialExponential(X(:, inFitWindow), t(:, inFitWindow), z, nBoot=pSTE.nBoot, nTrialsPerBin=pSTE.nTrialsPerBin);
    toc;

    nBins = length(obs.z);
    scatter(ax(1, 1), obs.z, obs.tau, 25, getColor(1:nBins, nBins, 0.8))
    title(ax(1, 1), sprintf('tau\nobs=%g, ci95=[%g, %g]', obs.slope.tau, ci95.tau(1), ci95.tau(2)))

    scatter(ax(2, 1), obs.z, obs.a, 25, getColor(1:nBins, nBins, 0.8))
    title(ax(2, 1), sprintf('a\nobs=%g, ci95=[%g, %g]', obs.slope.a, ci95.a(1), ci95.a(2)))

    hold(ax(:, 2:3), 'on')
    for iBin = 1:nBins
        plot(ax(1, 2), t(inFitWindow)', log(obs.X(iBin, :)'), Color=[getColor(iBin, nBins, 0.8), 0.5], LineWidth=1.5)
        plot(ax(2, 2), t(inFitWindow)', log(obs.a(iBin)') + t(inFitWindow)'./obs.tau(iBin)', Color=[getColor(iBin, nBins, 0.8), 0.5], LineWidth=1.5)
        plot(ax(1, 3), t(inFitWindow)', obs.X(iBin, :)', Color=[getColor(iBin, nBins, 0.8), 0.5], LineWidth=1.5)
        plot(ax(2, 3), t(inFitWindow)', obs.a(iBin)'.*exp(t(inFitWindow)'./obs.tau(iBin)'), Color=[getColor(iBin, nBins, 0.8), 0.5], LineWidth=1.5)
    end
    hold(ax(:, 2:3), 'off')
%     ylim(ax(:, 3), [0, 5])

    title(ax(1, 2), 'Log Scale')
    title(ax(1, 3), 'Real Scale')
end
%
% X = spike rate (nTrials, nTimestampsPerTrial), t = timestamp (1, nTimestampsPerTrial), z = movement time (nTrials, 1)
function [obs, boot, ci95, ci99, pValue] = bootSingleTrialExponential(X, t, z, varargin)
    p = inputParser();
    p.addRequired('X', @isnumeric)
    p.addRequired('t', @isnumeric)
    p.addRequired('z', @isnumeric)
    p.addParameter('nBoot', 100, @(n) n>=3);
    p.addParameter('nTrialsPerBin', 5, @isnumeric);
    p.parse(X, t, z, varargin{:})
    r = p.Results;
    k = r.nTrialsPerBin;
    nBoot = r.nBoot;

    nTrials = size(X, 1);
    nTimestamps = size(X, 2);
    assert(length(t) == nTimestamps);
    assert(length(z) == nTrials);
    t = reshape(t, 1, []);
    z = reshape(z, [], 1);

    % Order trials by movement time
    [~, trialOrderSorted] = sort(z, 'ascend');

    % Bin and average trials by movement time
    nBins = ceil(nTrials/k);
    XBinned = NaN(nBins, nTimestamps);
    zBinned = NaN(nBins, 1);
    for iBin = 1:nBins
        iFirstTrialInBin = 1 + (iBin-1)*k;
        selTrials = trialOrderSorted(iFirstTrialInBin:min(iFirstTrialInBin+k, nTrials));
        XBinned(iBin, :) = mean(X(selTrials, :), 1, 'omitnan');
        zBinned(iBin) = mean(z(selTrials));
    end

    % Fit exponential to each binned PETH
    tau=NaN(nBins, 1);
    a=NaN(nBins, 1);
    for iBin = 1:nBins
        isValid = ~isnan(XBinned(iBin, :));
        tt = t(isValid)';
        xx = XBinned(iBin, isValid)';
        coeff = fastfitlm([ones(size(tt)), tt], log(xx));
        a(iBin) = exp(coeff(1));
        tau(iBin) = 1/coeff(2);
    end

    obs = struct(X=XBinned, z=zBinned, tau=tau, a=a, slope=struct(tau=[], a=[]));

    % Estimate the slope of tau/a/c vs. z
    phi = [ones(size(obs.z)), obs.z]; % design matrix with column of ones and movement times, i.e. phi
    assert(nnz(isnan(obs.z)) == 0)
    assert(nnz(isnan(obs.tau)) == 0)
    assert(nnz(isnan(obs.a)) == 0)
    coeff = fastfitlm(phi, obs.tau);
    obs.slope.tau = coeff(2);
    coeff = fastfitlm(phi, obs.a);
    obs.slope.a = coeff(2);

    % Bootstrap(permute) by shuffling z and re-estimating slopes
    tau=NaN(nBins, nBoot);
    a=NaN(nBins, nBoot);
    boot = struct(tau=NaN(1, nBoot), a=NaN(1, nBoot));
    for iBoot = 1:nBoot
        zShuffled = z(randperm(nTrials));
        [~, trialOrderShuffled] = sort(zShuffled, 'ascend');

        for iBin = 1:nBins
            iFirstTrialInBin = 1 + (iBin-1)*k;
            selTrials = trialOrderShuffled(iFirstTrialInBin:min(iFirstTrialInBin+k, nTrials));
            XBinned(iBin, :) = mean(X(selTrials, :), 1, 'omitnan');
            zBinned(iBin) = mean(zShuffled(selTrials));
        end
        
        for iBin = 1:nBins
            isValid = ~isnan(XBinned(iBin, :));
            tt = t(isValid)';
            xx = XBinned(iBin, isValid)';
            coeff = fastfitlm([ones(size(tt)), tt], log(xx));
            a(iBin, iBoot) = exp(coeff(1));
            tau(iBin, iBoot) = 1/coeff(2);
        end

        phi = [ones(size(zBinned)), zBinned]; % design matrix with column of ones and movement times, i.e. phi

        sel = ~isinf(tau(:, iBoot)) & ~isnan(tau(:, iBoot)) & ~isinf(a(:, iBoot)) & ~isnan(a(:, iBoot));
        coeff = fastfitlm(phi, tau(sel, iBoot));
        boot.tau(iBoot) = coeff(2);

        coeff = fastfitlm(phi, a(sel, iBoot));
        boot.a(iBoot) = coeff(2);
    end

    ci95.tau = quantile(boot.tau, [0.025, 0.975]);
    ci99.tau = quantile(boot.tau, [0.005, 0.995]);
    if obs.slope.tau > median(boot.tau)
        pValue.tau = mean(boot.tau > obs.slope.tau) * 2;
    else
        pValue.tau = mean(boot.tau < obs.slope.tau) * 2;
    end

    ci95.a = quantile(boot.a, [0.025, 0.975]);
    ci99.a = quantile(boot.a, [0.005, 0.995]);
    if obs.slope.a > median(boot.a)
        pValue.a = mean(boot.a > obs.slope.a) * 2;
    else
        pValue.a = mean(boot.a < obs.slope.a) * 2;
    end
end

%%
function coeff = fastfitlm(X, y)
    coeff = (X'*X)\X'*y;
end