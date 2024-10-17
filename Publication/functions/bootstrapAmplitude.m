function [h, p, ci, obs] = bootstrapAmplitude(eu, varargin)
    p = inputParser();
    p.addRequired('eu', @(x) length(x) >=1 && isa(x, 'EphysUnit'))
    p.addRequired('trialsA', @(x) isa(x, 'Trial') || iscell(x) || ismember(x, {'press', 'lick'}))
    p.addRequired('trialsB', @(x) isa(x, 'Trial') || iscell(x) || ismember(x, {'press', 'lick'}))
    p.addParameter('correctionA', {}, @iscell)
    p.addParameter('correctionB', {}, @iscell)
    p.addParameter('nboot', 100000, @isnumeric)
    p.addParameter('baselineWindow', [-4, -2], @(x) isnumeric(x) && length(x) == 2 || isempty(x))
    p.addParameter('responseWindowA', [-0.3, -0]);
    p.addParameter('responseWindowB', [-0.3, 0]);
    p.addParameter('allowedTrialDuration', [2, Inf], @(x) isnumeric(x) && length(x) >= 2 && x(2) >= x(1))
    p.addParameter('withReplacement', false, @islogical)
    p.addParameter('alpha', 0.01, @isnumeric)    
    p.parse(eu, varargin{:});
    r = p.Results;
    eu = r.eu;
    trialsA = r.trialsA;
    trialsB = r.trialsB;

    rng(42);

    h = NaN(length(eu), 1);
    p = NaN(length(eu), 1);
    ci = NaN(length(eu), 2);
    obs = NaN(length(eu), 1);

    if iscell(trialsA)
        TRIALSA = trialsA;
        useCellTrialsA = true;
    else
        useCellTrialsA = false;
    end
    if iscell(trialsB)
        TRIALSB = trialsB;
        useCellTrialsB = true;
    else
        useCellTrialsB = false;
    end
    if ~isempty(r.correctionA)
        CORRECTIONA = r.correctionA;
        useCellCorrectionA = true;
    else
        useCellCorrectionA = false;
    end
    if ~isempty(r.correctionB)
        CORRECTIONB = r.correctionB;
        useCellCorrectionB = true;
    else
        useCellCorrectionB = false;
    end

    for iEu = 1:length(eu)
        fprintf(1, '%d/%d ', iEu, length(eu))
        if mod(iEu, 15) == 0
            fprintf(1, '\n')
        end

        if useCellTrialsA
            trialsA = TRIALSA{iEu};
        end
        if useCellTrialsB
            trialsB = TRIALSB{iEu};
        end
        if useCellCorrectionA
            correctionA = CORRECTIONA{iEu};
        else
            correctionA = [];
        end
        if useCellCorrectionB
            correctionB = CORRECTIONB{iEu};
        else
            correctionB = [];
        end

        if ischar(trialsA)
            trialsA = eu(iEu).getTrials(trialsA);
        end
        if ischar(trialsB)
            trialsB = eu(iEu).getTrials(trialsB);
        end
        trialsA = trialsA(trialsA.duration >= r.allowedTrialDuration(1) & trialsA.duration <= r.allowedTrialDuration(2));
        trialsB = trialsB(trialsB.duration >= r.allowedTrialDuration(1) & trialsB.duration <= r.allowedTrialDuration(2));
        nTrialsA = length(trialsA);
%         nTrialsB = length(trialsB);
        [responseA, ~] = eu(iEu).getTrialAlignedData('count', r.responseWindowA, trials=trialsA, alignTo='stop', resolution=diff(r.responseWindowA), correction=correctionA, allowedTrialDuration=r.allowedTrialDuration, correctionAdvancedValidation=false);
        [responseB, ~] = eu(iEu).getTrialAlignedData('count', r.responseWindowB, trials=trialsB, alignTo='stop', resolution=diff(r.responseWindowB), correction=correctionB, allowedTrialDuration=r.allowedTrialDuration, correctionAdvancedValidation=false);
        if isempty(r.baselineWindow)
            deltaResponseAll = [responseA./diff(r.responseWindowA); responseB./diff(r.responseWindowB)];
        else
            [baselineA, ~] = eu(iEu).getTrialAlignedData('count', r.baselineWindow, trials=trialsA, alignTo='stop', resolution=diff(r.baselineWindow), correction=correctionA, allowedTrialDuration=r.allowedTrialDuration, correctionAdvancedValidation=false);
            [baselineB, ~] = eu(iEu).getTrialAlignedData('count', r.baselineWindow, trials=trialsB, alignTo='stop', resolution=diff(r.baselineWindow), correction=correctionB, allowedTrialDuration=r.allowedTrialDuration, correctionAdvancedValidation=false);
            deltaResponseAll = [responseA./diff(r.responseWindowA) - baselineA./diff(r.baselineWindow); responseB./diff(r.responseWindowB) - baselineB./diff(r.baselineWindow)];
        end
        % "deltaResponse = response - baseline"

        if r.withReplacement
            [~, bsample] = bootstrp(r.nboot, [], deltaResponseAll);
        else
            bsample = zeros(length(deltaResponseAll), r.nboot);
            for iboot = 1:r.nboot
                bsample(:, iboot) = randperm(length(deltaResponseAll));
            end
        end

        sampleA = deltaResponseAll(bsample(1:nTrialsA, :));
        sampleB = deltaResponseAll(bsample(nTrialsA + 1:end, :));

        % H0: group A mean  - group B mean = 0
        % H1: group A mean >> group B mean || group A mean << group B mean
        obs(iEu) = mean(deltaResponseAll(1:nTrialsA)) - mean(deltaResponseAll(nTrialsA + 1:end));
        boot = mean(sampleA, 1) - mean(sampleB, 1);

        ci(iEu, :) = prctile(boot, [r.alpha*50, 100 - r.alpha*50]);
        p(iEu) = nnz(boot<=-abs(obs(iEu)) | boot>=abs(obs(iEu))) ./ r.nboot;
        h(iEu) = p(iEu) <= r.alpha;
    end
end