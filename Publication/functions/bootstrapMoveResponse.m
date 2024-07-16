function [h, muDiffCI, muDiffObs] = bootstrapMoveResponse(eu, trialType, varargin)
    p = inputParser();
    p.addRequired('eu', @(x) length(x) >= 1 && isa(x, 'EphysUnit'));
    p.addRequired('trialType', @(x) ismember(x, {'press', 'lick'}));
    p.addParameter('nboot', 100000, @isnumeric)
    p.addParameter('baselineWindow', [-4, -2], @(x) isnumeric(x) && length(x) == 2)
    p.addParameter('responseWindow', [-0.5, -0.2], @(x) isnumeric(x) && length(x) == 2)
    p.addParameter('alignTo', 'stop', @(x) ischar(x) && ismember(lower(x), {'start', 'stop'}))
    p.addParameter('allowedTrialDuration', [2, Inf], @(x) isnumeric(x) && length(x) >= 2 && x(2) >= x(1))
    p.addParameter('trialDurationError', 1e-3, @isnumeric) % Used for opto, error allowed when finding identical trial durations.
    p.addParameter('alpha', 0.01, @isnumeric)
    p.addParameter('withReplacement', false, @islogical)
    p.addParameter('oneSided', false, @islogical)
    p.parse(eu, trialType, varargin{:});
    r = p.Results;
    eu = r.eu;

    rng(42);

    dataWindow = [min(r.baselineWindow(1), r.responseWindow(1)), max(r.baselineWindow(2), r.responseWindow(2))];

    h = NaN(length(eu), 1);
    p = h;
    muDiffCI = NaN(length(eu), 2);
    muDiffObs = NaN(length(eu), 1);
    for iEu = 1:length(eu)
        fprintf(1, '%d/%d ', iEu, length(eu))
        if mod(iEu, 15) == 0
            fprintf(1, '\n')
        end

        [sr, t] = eu(iEu).getTrialAlignedData('count', dataWindow, r.trialType, alignTo=r.alignTo, ...
            allowedTrialDuration=r.allowedTrialDuration, trialDurationError=r.trialDurationError, ...
            includeInvalid=false, resolution=0.1);

        if isempty(sr)
            warning('Spike rate for %d - %s is empty.', iEu, eu(iEu).getName('_'));
            continue
        end
    
        response = mean(sr(:, t >= r.responseWindow(1) & t <= r.responseWindow(2)), 2, 'omitnan');
        nBins = nnz(t >= r.responseWindow(1) & t <= r.responseWindow(2));
        baselineSampleIndices = find(t >= r.baselineWindow(1) & t <= r.baselineWindow(2));
        baselineSampleIndices = baselineSampleIndices((1:nBins) + flip(length(baselineSampleIndices)-nBins:-nBins:0)');
        baseline = NaN(size(sr, 1), size(baselineSampleIndices, 1));
        for i = 1:size(baselineSampleIndices, 1)
            baseline(:, i) = mean(sr(:, baselineSampleIndices(i, :)), 2);
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
%         if muDiffObs(iEu) > muDiffCI(iEu, 2)
%             h(iEu) = 1;
%         elseif muDiffObs(iEu) < muDiffCI(iEu, 1)
%             h(iEu) = -1;
%         else
%             h(iEu) = 0;
%         end
    end
end