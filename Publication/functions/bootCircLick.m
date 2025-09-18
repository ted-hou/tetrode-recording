
function [magH, magCI, Z] = bootCircLick(eta, varargin)
    parser = inputParser();
    parser.addRequired('eta', @isstruct);
    parser.addParameter('selUnits', [], @(x) islogical(x) || isnumeric(x))
    parser.addParameter('alpha', 0.01, @isnumeric);
    parser.addParameter('nBoot', 10000, @isnumeric);
    parser.addParameter('replace', true, @islogical) % true for bootstrap, false for wda
    parser.addParameter('seed', 42, @isnumeric)
    parser.addParameter('interpFirstBin', false, @islogical)
    parser.addParameter('replaceNansWithMean', true, @islogical)
    parser.parse(eta, varargin{:});
    r = parser.Results;
    X = r.eta.X;
    t = r.eta.t;
    selUnits = r.selUnits;
    alpha = r.alpha;
    nBoot = r.nBoot;
    replace = r.replace;
    seed = r.seed;
    interpFirstBin = r.interpFirstBin;
    replaceNansWithMean = r.replaceNansWithMean;
    rng(seed);

    magH = false(size(X, 1), 1);
    magCI = zeros(size(X, 1), 2);
    Z = NaN(size(X, 1), 1) + 1i*NaN(size(X, 1), 1);

    nUnits = size(X, 1);
    if isempty(selUnits)
        selUnits = 1:nUnits;
    elseif islogical(selUnits)
        selUnits = reshape(find(selUnits), 1, []);
    end
    nBins = length(t);
    tTic = tic();
    fprintf('Bootstrapping %i units...', length(selUnits))
    lineLength = 0;
    for iUnit = selUnits
        fprintf(repmat('\b', 1, lineLength));
        lineLength = fprintf('%i/%i', iUnit, nUnits);
        if replace
            I = randi(nBins, [nBoot, nBins]);
        else
            I = zeros(nBoot, nBins);
            for iBoot = 1:nBoot
                I(iBoot, :) = randperm(nBins);
            end
        end
        x = X(iUnit, :);
        if interpFirstBin
            x(1) = mean(x([2, end])); % First bin has lick artifact usually
        end
        if replaceNansWithMean
            x(isnan(x)) = mean(x, 'omitnan');
            zObs = mean(x.*exp(t*1i));
            zRand = mean(x(I).*exp(t*1i), 2);
        else
            zObs = mean(x.*exp(t*1i));
            zRand = mean(x(I).*exp(t*1i), 2);
        end
        
        magObs = abs(zObs);
        magRand = abs(zRand);
        magCI(iUnit, :) = quantile(magRand, [0, 1 - alpha]);
        magH(iUnit) = magObs > magCI(iUnit, 2);
        Z(iUnit) = zObs;
    end
    fprintf('Done (%.2fs)\n', toc(tTic));
end