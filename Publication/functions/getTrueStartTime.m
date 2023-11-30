function tst = getTrueStartTime(exp, varargin)
    p = inputParser();
    p.addRequired('exp', @(x) isa(x, 'CompleteExperiment'))
    p.addParameter('features', {'handL', 'handR'}, @iscell);
    p.addParameter('stats', {'xPos', 'yPos'}, @iscell);
    p.addParameter('threshold', 2.5, @isnumeric);
    p.addParameter('trialType', 'press', @ischar);
    p.parse(exp, varargin{:})
    r = p.Results;
    exp = r.exp;


    [velocityKernels, ~, ~] = CompleteExperiment.makeConsineKernels(0, width=0.2, overlap=0.5, direction='both');
    
    tst = cell(1, length(exp));
    for iExp = 1:length(exp)
        trials = exp(iExp).eu(1).getTrials(r.trialType);
        trueStartTime = NaN(length(trials), 1);
        for iTrial = 1:length(trials)
            t = trials(iTrial).Start:1/30:trials(iTrial).Stop;
            F = exp(iExp).getFeatures(timestamps=t, features=r.features, stats=r.stats);
            F = CompleteExperiment.convolveFeatures(F, velocityKernels, kernelNames={'_smooth'}, ...
                features=r.features, ...
                stats=r.stats, ...
                mode='replace', normalize='none');
            F.inTrial = [];
            F.t = [];
            F = normalize(F);
            data = table2array(F);
            trueStartIndex = find(all(abs(data) < r.threshold, 2), 1, 'last');
            if ~isempty(trueStartIndex)
                trueStartTime(iTrial) = t(trueStartIndex) - t(end);
            end
        end
        tst{iExp} = trueStartTime;
    end

end