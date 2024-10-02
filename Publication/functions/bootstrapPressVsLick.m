function [h, p, ci, obs] = bootstrapPressVsLick(eu, varargin)
    p = inputParser();
    p.addRequired('eu', @(x) length(x) >=1 && isa(x, 'EphysUnit'))
    p.addParameter('nboot', 100000, @isnumeric)
    p.addParameter('baselineWindow', [-4, -2], @(x) isnumeric(x) && length(x) == 2)
    p.addParameter('pressWindow', [-0.3, -0]);
    p.addParameter('lickWindow', [-0.3, 0]);
    p.addParameter('allowedTrialDuration', [2, Inf], @(x) isnumeric(x) && length(x) >= 2 && x(2) >= x(1))
    p.addParameter('withReplacement', false, @islogical)
    p.addParameter('alpha', 0.01, @isnumeric)    
    p.parse(eu, varargin{:});
    r = p.Results;
    eu = r.eu;

    rng(42);

    h = NaN(length(eu), 1);
    p = NaN(length(eu), 1);
    ci = NaN(length(eu), 2);
    obs = NaN(length(eu), 1);

    for iEu = 1:length(eu)
        fprintf(1, '%d/%d ', iEu, length(eu))
        if mod(iEu, 15) == 0
            fprintf(1, '\n')
        end

        pressTrials = eu(iEu).Trials.Press;
        lickTrials = eu(iEu).Trials.Lick;
        pressTrials = pressTrials(pressTrials.duration >= r.allowedTrialDuration(1) & pressTrials.duration <= r.allowedTrialDuration(2));
        lickTrials = lickTrials(lickTrials.duration >= r.allowedTrialDuration(1) & lickTrials.duration <= r.allowedTrialDuration(2));
        nPress = length(pressTrials);
        nLick = length(lickTrials);
        [srPress, ~] = eu(iEu).getTrialAlignedData('count', r.pressWindow, trials=pressTrials, alignTo='stop', resolution=diff(r.pressWindow));
        [blPress, ~] = eu(iEu).getTrialAlignedData('count', r.baselineWindow, trials=pressTrials, alignTo='stop', resolution=diff(r.baselineWindow));
        [srLick, ~] = eu(iEu).getTrialAlignedData('count', r.lickWindow, trials=lickTrials, alignTo='stop', resolution=diff(r.lickWindow));
        [blLick, ~] = eu(iEu).getTrialAlignedData('count', r.baselineWindow, trials=lickTrials, alignTo='stop', resolution=diff(r.baselineWindow));
        sr = [srPress./diff(r.pressWindow) - blPress./diff(r.baselineWindow); srLick./diff(r.lickWindow) - blLick./diff(r.baselineWindow)];

        if r.withReplacement
            [~, bsample] = bootstrp(r.nboot, [], sr);
        else
            bsample = zeros(length(sr), r.nboot);
            for iboot = 1:r.nboot
                bsample(:, iboot) = randperm(length(sr));
            end
        end

        pressSamples = sr(bsample(1:nPress, :));
        lickSamples = sr(bsample(nPress + 1:end, :));

        obs(iEu) = mean(sr(1:nPress)) - mean(sr(nPress + 1:end));
        boot = mean(pressSamples, 1) - mean(lickSamples, 1);

        ci(iEu, :) = prctile(boot, [r.alpha*50, 100 - r.alpha*50]);
        p(iEu) = nnz(boot<=-abs(obs(iEu)) | boot>=abs(obs(iEu))) ./ r.nboot;
        h(iEu) = p(iEu) <= r.alpha;
    end
end