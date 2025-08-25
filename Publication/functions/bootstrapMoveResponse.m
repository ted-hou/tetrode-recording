function [h, muDiffCI, muDiffObs] = bootstrapMoveResponse(eu, trialType, varargin)
    p = inputParser();
    p.addRequired('eu', @(x) length(x) >= 1 && isa(x, 'EphysUnit'));
    p.addRequired('trialType');
    p.addParameter('nboot', 100000, @isnumeric)
    p.addParameter('baselineWindow', [-4, -2], @(x) isnumeric(x) && length(x) == 2)
    p.addParameter('responseWindow', [-0.5, -0.2], @(x) isnumeric(x) && length(x) == 2)
    p.addParameter('alignTo', 'stop', @(x) isstruct(x) || ismember(lower(x), {'start', 'stop'}))
    p.addParameter('allowedTrialDuration', [2, Inf], @(x) isstruct(x) || (isnumeric(x) && length(x) >= 2 && x(2) >= x(1)))
    p.addParameter('trialDurationError', 1e-3, @isnumeric) % Used for opto, error allowed when finding identical trial durations.
    p.addParameter('alpha', 0.01, @isnumeric)
    p.addParameter('withReplacement', false, @islogical)
    p.addParameter('oneSided', false, @islogical)
    p.addParameter('correction', {}, @(x) iscell(x) || isstruct(x)) % correction for movement onset time
    p.addParameter('trials', {}, @(x) iscell(x) || isstruct(x))
    p.addParameter('artifacts', [], @(x) isempty(x) || all(isfield(x, {'event', 'length', 'lengthUnit', 'direction'})))
    p.parse(eu, trialType, varargin{:});
    r = p.Results;
    eu = r.eu;

    % Optionally specify two different trial types, one for baseline, one
    % for response
    if isstruct(r.trialType)
        assert(isstruct(r.trialType) && all(isfield(r.trialType, {'baseline', 'response'})), 'Different trial type are specified for baseline and response, please check the following parameters: "alignTo", "allowedTrialDuration", "correction", "trials"')
    else
        r.trialType = struct(baseline=r.trialType, response=r.trialType);
    end
    if isstruct(r.alignTo)
        assert(all(isfield(r.alignTo, {'baseline', 'response'})))
    else
        r.alignTo = struct(baseline=r.alignTo, response=r.alignTo);
    end
    if isstruct(r.allowedTrialDuration)
        assert(all(isfield(r.allowedTrialDuration, {'baseline', 'response'})))
    else
        r.allowedTrialDuration = struct(baseline=r.allowedTrialDuration, response=r.allowedTrialDuration);
    end
    if isstruct(r.correction)
        assert(all(isfield(r.correction, {'baseline', 'response'})))
    else
        r.correction = struct(baseline=r.correction, response=r.correction);
    end
    if isstruct(r.trials)
        assert(all(isfield(r.trials, {'baseline', 'response'})))
    else
        r.trials = struct(baseline=r.trials, response=r.trials);
    end

    for var = ["baseline", "response"]
        if isempty(r.correction) || isempty(r.correction.(var))
            correction.(var) = cell(length(eu), 1);
        else
            correction.(var) = r.correction.(var);
            assert(length(correction.(var)) == length(eu));
        end
        if isempty(r.trials) || isempty(r.trials.(var)) 
            trials.(var) = cell(length(eu), 1);
        else
            trials.(var) = r.trials.(var);
            assert(length(trials.(var)) == length(eu))
        end
    end

    rng(42);

    h = NaN(length(eu), 1);
    muDiffCI = NaN(length(eu), 2);
    muDiffObs = NaN(length(eu), 1);
    lineLength = 0;
    tTicAll = tic();
    for iEu = 1:length(eu)
        tTic = tic();
        [sr, t] = eu(iEu).getTrialAlignedData('count', r.responseWindow, r.trialType.response, alignTo=r.alignTo.response, ...
            allowedTrialDuration=r.allowedTrialDuration.response, trialDurationError=r.trialDurationError, ...
            includeInvalid=false, resolution=0.1, correction=correction.response{iEu}, trials=trials.response{iEu}, ...
            artifacts=r.artifacts);

        [srb, tb] = eu(iEu).getTrialAlignedData('count', r.baselineWindow, r.trialType.baseline, alignTo=r.alignTo.baseline, ...
            allowedTrialDuration=r.allowedTrialDuration.baseline, trialDurationError=r.trialDurationError, ...
            includeInvalid=false, resolution=0.1, correction=correction.baseline{iEu}, trials=trials.baseline{iEu}, ...
            artifacts=r.artifacts);

        if isempty(sr) || isempty(srb)
            warning('Spike rate for %d - %s is empty.', iEu, eu(iEu).getName('_'));
            continue
        end
    
        response = mean(sr(:, t >= r.responseWindow(1) & t <= r.responseWindow(2)), 2, 'omitnan');
        nBins = nnz(t >= r.responseWindow(1) & t <= r.responseWindow(2));
        baselineSampleIndices = find(tb >= r.baselineWindow(1) & tb <= r.baselineWindow(2));
        baselineSampleIndices = baselineSampleIndices((1:nBins) + flip(length(baselineSampleIndices)-nBins:-nBins:0)');
        baseline = NaN(size(srb, 1), size(baselineSampleIndices, 1));
        for i = 1:size(baselineSampleIndices, 1)
            baseline(:, i) = mean(srb(:, baselineSampleIndices(i, :)), 2);
        end
        baseline = baseline(:);
        combined = [baseline; response];
        nBase = length(baseline);
        
        % With replacement
        if r.withReplacement
            [~, bsample] = bootstrp(r.nboot, [], combined);
        else
            bsample = zeros(length(combined), r.nboot);
            for iboot = 1:r.nboot
                bsample(:, iboot) = randperm(length(combined));
            end
        end
        baselineSamples = combined(bsample(1:nBase, :));
        responseSamples = combined(bsample(nBase+1:end, :));
        muDiffObs(iEu) = mean(response, 'omitnan') - mean(baseline, 'omitnan');
        if muDiffObs(iEu) > 0
            direction = 1;
        else
            direction = -1;
        end
        muDiffBoot = mean(responseSamples, 1, 'omitnan') - mean(baselineSamples, 1, 'omitnan');
        if r.oneSided
            if direction == 1
                muDiffCI(iEu, :) = prctile(muDiffBoot, [0, 100 - r.alpha*100]);
            elseif direction == -1
                muDiffCI(iEu, :) = prctile(muDiffBoot, [r.alpha*100, 100]);
            end
        else
            muDiffCI(iEu, :) = prctile(muDiffBoot, [r.alpha*50, 100 - r.alpha*50]);
        end
        if direction == 1
            h(iEu) = muDiffObs(iEu) > muDiffCI(iEu, 2);
        elseif direction == -1
            h(iEu) = -(muDiffObs(iEu) < muDiffCI(iEu, 1));
        end
        fprintf(repmat('\b', 1, lineLength))
        lineLength = fprintf(1, '%d/%d (%gs)', iEu, length(eu), toc(tTic));
    end
    fprintf(repmat('\b', 1, lineLength))
    fprintf('Bootstrapped %i units in %gs.\n', length(eu), toc(tTicAll))
end