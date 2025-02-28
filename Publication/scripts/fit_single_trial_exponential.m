%%
clear pSTE
pSTE.kernel = struct(type='exponential', params=struct(lambda1=5, lambda2=10, resolution=0.01, width=1));
pSTE.minTrialLength = 1.5;
pSTE.minNumTrials = 20;
pSTE.cueBlankWindow = 0.5;
pSTE.fitWindow = [-2, 0];
pSTE.baselineWindow = [-4, -2];
pSTE.nBoot = 10000;
pSTE.nTrialsPerBin = 5;

iEu = find(c.hasPress & c.hasLick & c.isLickUp, 1, 'first');

% Get PETH
trials = eu(iEu).getTrials('lick');
trials = trials(trials.duration() >= pSTE.minTrialLength);
assert(length(trials) >= pSTE.minNumTrials)
[X, t] = eu(iEu).getTrialAlignedData('rate', window=[-4, 2], trialType=trialType, resolution=0.01, trials=trials, allowedTrialDuration=[0, Inf], includeInvalid=false, startBlankWindow=[-Inf, pSTE.cueBlankWindow], kernel=pSTE.kernel);
z = trials.duration();

% Normalize PETH
inBaselineWindow = t >= -4 & t <= -2;
inFitWindow = t >= -2 & t <= 0;
mu = mean(X(:, inBaselineWindow), 'all', 'omitnan');
sigma = std(X(:, inBaselineWindow), 0, 'all', 'omitnan');
X = (X - mu) ./ sigma;

tic;
[obs, boot, ci95, ci99] = bootSingleTrialExponential(X(:, inFitWindow), t(:, inFitWindow), z, nBoot=pSTE.nBoot, nTrialsPerBin=pSTE.nTrialsPerBin, isUp=true);
toc;
%%
% X = spike rate (nTrials, nTimestampsPerTrial), t = timestamp (1, nTimestampsPerTrial), z = movement time (nTrials, 1)
function [obs, boot, ci95, ci99, pValue] = bootSingleTrialExponential(X, t, z, varargin)
    p = inputParser();
    p.addRequired('X', @isnumeric)
    p.addRequired('t', @isnumeric)
    p.addRequired('z', @isnumeric)
    p.addParameter('nBoot', 100, @(n) n>=3);
    p.addParameter('nTrialsPerBin', 5, @isnumeric);
    p.addParameter('isUp', true, @islogical);
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

    if r.isUp
        modelSpec = fittype(@(a, b, c, t) a*exp(b*t)+c, dependent='x', independent='t');
    else
        modelSpec = fittype(@(a, b, c, t) -a*exp(b*t)+c, dependent='x', independent='t');
    end

%     bootSlope.args = struct(nBoot=nBoot, nTrialsPerBin=r.nTrialsPerBin);
%     bootSlope.tau = struct(ci95=[], ci99=[], p=[]);
%     bootSlope.a = struct(ci95=[], ci99=[], p=[]);
%     bootSlope.c = struct(ci95=[], ci99=[], p=[]);

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
    c=NaN(nBins, 1);
    for iBin = 1:nBins
        isValid = ~isnan(XBinned(iBin, :));
        mdl = fit(t(isValid)', XBinned(iBin, isValid)', modelSpec, StartPoint=[1, 10, 0], Lower=[0.1, 0.1, -5], Upper=[10, 100, 5]); % Starting points and bounds for a and c assumes normalized spike rates
        tau(iBin) = 1./mdl.b;
        a(iBin) = mdl.a;
        c(iBin) = mdl.c;
    end

    obs = struct(X=XBinned, z=zBinned, tau=tau, a=a, c=c, slope=struct(tau=[], a=[], c=[]));

    % Estimate the slope of tau/a/c vs. z
    paramNames = ["tau", "a", "c"];
    for iParam = 1:3
        paramName = paramNames(iParam);
        mdl = fitlm(obs.z, obs.(paramName), VarNames=["z", paramName]);
        obs.slope.(paramName) = table2array(mdl.Coefficients('z', 'Estimate'));
    end

    % Bootstrap(permute) by shuffling z and re-estimating slopes
    tau=NaN(nBins, nBoot);
    a=NaN(nBins, nBoot);
    c=NaN(nBins, nBoot);
    boot = struct(tau=NaN(1, nBoot), a=NaN(1, nBoot), c=NaN(1, nBoot));
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
            mdl = fit(t(isValid)', XBinned(iBin, isValid)', modelSpec, StartPoint=[1, 10, 0], Lower=[0.1, 0.1, -5], Upper=[10, 100, 5]); % Starting points and bounds for a and c assumes normalized spike rates
            tau(iBin, iBoot) = 1./mdl.b;
            a(iBin, iBoot) = mdl.a;
            c(iBin, iBoot) = mdl.c;
        end

        mdl = fitlm(zBinned, tau(:, iBoot), VarNames=["z", "tau"]);
        boot.tau(iBoot) = table2array(mdl.Coefficients('z', 'Estimate'));
        mdl = fitlm(zBinned, a(:, iBoot), VarNames=["z", "a"]);
        boot.a(iBoot) = table2array(mdl.Coefficients('z', 'Estimate'));
        mdl = fitlm(zBinned, c(:, iBoot), VarNames=["z", "c"]);
        boot.c(iBoot) = table2array(mdl.Coefficients('z', 'Estimate'));
    end

    for iParam = 1:3
        paramName = paramNames(iParam);
        ci95.(paramName) = quantile(boot.(paramName), [0.025, 0.975]);
        ci99.(paramName) = quantile(boot.(paramName), [0.005, 0.995]);
        if obs.slope.(paramName) > median(boot.(paramName))
            pValue.(paramName) = mean(boot.(paramName) > obs.slope.(paramName));
        else
            pValue.(paramName) = mean(boot.(paramName) < obs.slope.(paramName));
        end
    end
end