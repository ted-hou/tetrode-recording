
% clearvars -except exp
MDL = cell(length(expAcute), 1);
SRT = cell(length(expAcute), 1);
FT = cell(length(expAcute), 1);
TT = cell(length(expAcute), 1);
NOTNAN = cell(length(expAcute), 1);
% MASK = cell(length(exp), 1);

for iExp = 1:length(expAcute)
    F = expAcute(iExp).getFeatures(sampleRate=30, trialType={'press'}, stats={'xVel', 'yVel'}, ...
        features={'handL', 'handR', 'footL', 'footR', 'nose', 'spine', 'trialStart', 'pressTrialRamp', 'firstPressRamp'}, ...
        likelihoodThreshold=0.95);
%     [~, mask] = exp(iExp).maskFeaturesByTrial(F, NaN, 'press', [-1, 3], ...
%         features={'handL', 'handR', 'footL', 'footR', 'nose', 'spine'}, ...
%         stats={'xVel', 'yVel'}, replace=true);
    [velocityKernels, ~, velocityDelays] = CompleteExperiment.makeConsineKernels(4, width=0.1, overlap=0.5, direction='both');
    velocityDelays = round(velocityDelays * 1000);
    F = CompleteExperiment.convolveFeatures(F, velocityKernels, kernelNames=velocityDelays, ...
        features={'handL', 'handR', 'footL', 'footR', 'nose', 'spine'}, ...
        stats={'xVel', 'yVel'}, ...
        mode='replace', normalize='maxabs');

    [eventKernels, ~, eventDelays] = CompleteExperiment.makeConsineKernels(6, width=0.4, overlap=0.75, direction='both');
    eventDelays = round(eventDelays * 1000);
    F = CompleteExperiment.convolveFeatures(F, eventKernels, kernelNames=eventDelays, ...
        features={'trialStart'}, ...
        mode='replace', normalize='none');

    % Mask velocity predictors befor first movement


    t = F.t;
    inTrial = F.inTrial;
    F.t = [];
    F.inTrial = [];
    Ft = F(inTrial, :);
    tt = t(inTrial);
%     mask = mask(inTrial);


    Ft.constant = ones(height(Ft), 1);

    % Names of predictors
    % Cue and first move
    % Movement velocities
    % Trial-length ramping signals
    % Short pre-movement ramps
    names = Ft.Properties.VariableNames;
    cuePredictors = names(contains(names, {'trialStart'}))';
    velocityPredictors = names(contains(names, {'Vel'}))';
    rampPredictors = names(contains(names, {'firstPressRamp'}))';
    trialProgressPredictors = names(contains(names, 'TrialRamp'))';


    variantPredictors = { ...
        {'constant'}, ...
        [{'constant'}; cuePredictors], ...
        [{'constant'}; cuePredictors; velocityPredictors], ...
        [{'constant'}; cuePredictors; velocityPredictors; rampPredictors], ...
        [{'constant'}; cuePredictors; velocityPredictors; rampPredictors; trialProgressPredictors]};
    variantNames = {'Constant', '+Cue', '+Velocity', '+Ramp', '+Trial-progress'};
    nVariants = length(variantNames);
    
    % All nan if one nan
    notnan = all(~isnan(Ft(:, variantPredictors{end}).Variables), 2);
    prctNotNan = nnz(notnan) ./ height(Ft); 
    data = Ft.Variables;
    data(~notnan, :) = NaN;
    Ft.Variables = data;

    FT{iExp} = Ft;
    TT{iExp} = tt;

    % 2.3 Fit GLM
    % Grab a unit and start fitting glms!
    fprintf(1, 'Fitting %g modelsets...\n', length(expAcute(iExp).eu)); tTicAll = tic();
    mdl = cell(length(expAcute(iExp).eu), nVariants);
    srt = cell(length(expAcute(iExp).eu), 1);
    warning('off','all')
    for iEu = 1:length(expAcute(iExp).eu)
        fprintf(1, '\tFitting modelset %g of %g...\n', iEu, length(expAcute(iExp).eu)); tTic = tic();
        euAcute = expAcute(iExp).eu(iEu);
        srTrialAligned = [0, euAcute.getSpikeRates('gaussian', 0.1, t)]'; 
        srt{iEu} = srTrialAligned(inTrial);
    
        for iVariant = 1:nVariants
            thisF = Ft;
            thisF.SpikeRate = double(srt{iEu});
%             if strcmp(variantNames{iVariant}, '+PreMoveVel')
%                 thisF(mask, :) = [];
%             end
            mdl{iEu, iVariant} = fitglm(thisF, ResponseVar='SpikeRate', PredictorVars=variantPredictors{iVariant}, Distribution='poisson');
            fprintf(1, '\t\t%s R^2 = %.2f\n', variantNames{iVariant}, mdl{iEu, iVariant}.Rsquared.Ordinary);
        end

        % fprintf(1, '\tDone (%.0f%% not nan) in %.2f sec.\n', prctNotNan*100, toc(tTic));
    end
    warning('on','all')
    warning('query','all')
    fprintf(1, 'Fitted %g units in %.2f seconds.\n', length(expAcute(iExp).eu), toc(tTicAll));
    MDL{iExp} = cellfun(@compact, mdl, UniformOutput=false);
    SRT{iExp} = srt;
    NOTNAN{iExp} = notnan;
%     MASK{iExp} = mask;
end
mdl = cat(1, MDL{:});
srt = cat(1, SRT{:});
expIndices = zeros(size(mdl, 1), 1);
i = 0;
for iExp = 1:length(expAcute)
    expIndices(i + 1:i + length(expAcute(iExp).eu)) = iExp;
    i = i + length(expAcute(iExp).eu);
end
% clearvars -except exp mdl srt variantNames nVariants FT TT NOTNAN expIndices

% %% Estimate loss
% L = {};
% eu = [exp.eu];
% curExpIndex = -1;
% for iEu = 1:length(eu)
%     iExp = expIndices(iEu);
%     if curExpIndex ~= iEu
%         curExpIndex = iEu;
%         
%         F = exp(iExp).getFeatures(sampleRate=30, stats={'xVel', 'yVel'}, ...
%             features={'handL', 'handR', 'footL', 'footR', 'nose', 'spine' 'trialStart', 'firstPress', 'firstLick', 'pressTrialRamp', 'lickTrialRamp', 'firstPressRamp', 'firstLickRamp'}, ...
%             likelihoodThreshold=0.95);
%         [velocityKernels, ~, velocityDelays] = CompleteExperiment.makeConsineKernels(4, width=0.1, overlap=0.5, direction='both');
%         velocityDelays = round(velocityDelays * 1000);
%         F = CompleteExperiment.convolveFeatures(F, velocityKernels, kernelNames=velocityDelays, ...
%             features={'handL', 'handR', 'footL', 'footR', 'nose', 'spine'}, ...
%             stats={'xVel', 'yVel'}, ...
%             mode='replace', normalize='maxabs');
% 
%         [eventKernels, ~, eventDelays] = CompleteExperiment.makeConsineKernels(6, width=0.4, overlap=0.75, direction='both');
%         eventDelays = round(eventDelays * 1000);
%         F = CompleteExperiment.convolveFeatures(F, eventKernels, kernelNames=eventDelays, ...
%             features={'trialStart'}, ...
%             mode='replace', normalize='none');
%     
%         [eventKernels, ~, eventDelays] = CompleteExperiment.makeConsineKernels(6, width=0.4, overlap=0.75, direction='left');
%         eventDelays = round(eventDelays * 1000);
%         F = CompleteExperiment.convolveFeatures(F, eventKernels, kernelNames=eventDelays, ...
%             features={'firstPress', 'firstLick'}, ...
%             mode='replace', normalize='none');
%     end
% 
%     for iVariant = 1:nVariants
%         ll = loss(mdl(iEu, iVariant), )
%     end
% 
% end

%% Estimate model performance (R^2)

R2 = NaN(size(mdl));
for iEu = 1:size(mdl, 1)
    X = FT{expIndices(iEu)};
    y = srt{iEu};

    for iVariant = 1:nVariants
        yHat = predict(mdl{iEu, iVariant}, X, Simultaneous=true);
        sel = ~isnan(yHat) & ~isnan(y);
        R2(iEu, iVariant) = corr(yHat(sel), y(sel)) .^ 2;
    end
end

%% Estimate model criterion
crits = ["AIC", "AICc", "BIC", "CAIC"];
for crit = crits
    modelCriterion.(crit) = cellfun(@(mdl) mdl.ModelCriterion.(crit), mdl);
end

%% Calculate trial-average fitted vs. observed for all units
euAcute = vertcat(expAcute.eu);
pAcute.minSpikeRate = 15;
pAcute.minTrialDuration = 2;
pAcute.minNumTrials = 30;
pAcute.etaNorm = [-4, -2];
pAcute.etaWindow = [-4, 2];
pAcute.metaWindowPress = [-0.3, 0];
pAcute.metaWindowLick = [-0.3, 0];
% pAcute.posRespThreshold = 1;
% pAcute.negRespThreshold = -0.5;

cAcute.hasPress = arrayfun(@(e) nnz(e.getTrials('press').duration() >= pAcute.minTrialDuration) >= pAcute.minNumTrials, euAcute)';

% to bootstrap significantly movement-modulated units

pAcute.bootAlpha = 0.01;
bootAcute.press = struct('h', NaN(length(euAcute), 1), 'muDiffCI', NaN(length(euAcute), 2), 'muDiffObs', NaN(length(euAcute), 1));
[bootAcute.press.h(cAcute.hasPress), bootAcute.press.muDiffCI(cAcute.hasPress, :), bootAcute.press.muDiffObs(cAcute.hasPress)] = bootstrapMoveResponse( ...
    euAcute(cAcute.hasPress), 'press', alpha=pAcute.bootAlpha, withReplacement=false, oneSided=false, ...
    responseWindow=[-0.3, 0]);
fprintf(1, '\nAll done\n')

% Report bootstraped movement response direction
assert(nnz(isnan(bootAcute.press.h(cAcute.hasPress))) == 0)

figure, histogram(bootAcute.press.h)
cAcute.isPressUp = bootAcute.press.h' == 1 & cAcute.hasPress;
cAcute.isPressDown = bootAcute.press.h' == -1 & cAcute.hasPress;
cAcute.isPressResponsive = cAcute.isPressUp | cAcute.isPressDown;

% ax = axes(figure()); hold on;

assert(length(euAcute) > 1)
msrAcute = cell(length(euAcute), 1);
msrHatAcute = cell(length(euAcute), 1);
peakAcute = NaN(length(euAcute), 1);
tPeakAcute = NaN(length(euAcute), 1);
tOnsetAcute = NaN(length(euAcute), 1);
peakHatAcute = NaN(length(euAcute), nVariants);
tPeakHatAcute = NaN(length(euAcute), nVariants);
tOnsetHatAcute = NaN(length(euAcute), nVariants);
for iEu = 1:length(euAcute)
    iExp = expIndices(iEu);
    tt = TT{iExp};
    Ft = FT{iExp};

    srtHat = NaN(height(Ft), size(mdl, 2));
    for iVariant = 2:size(mdl, 2)
        srtHat(:, iVariant) = predict(mdl{iEu, iVariant}, Ft);
    end

    trials = expAcute(iExp).eu(1).getTrials('press');

    maxTrialSampleLength = 0;
    for iTrial = 1:length(trials)
        inTrial = tt >= trials(iTrial).Start & tt <= trials(iTrial).Stop;
        maxTrialSampleLength = max(maxTrialSampleLength, nnz(inTrial));
    end

    srTrialAligned = NaN(length(trials), maxTrialSampleLength);
    srTrialAlignedHat = NaN(length(trials), maxTrialSampleLength, size(mdl, 2));
    for iTrial = 1:length(trials)
        inTrial = tt >= trials(iTrial).Start & tt <= trials(iTrial).Stop;
        trialSampleLength = nnz(inTrial);
        srTrialAligned(iTrial, maxTrialSampleLength - trialSampleLength + 1:end) = srt{iEu}(inTrial);
        for iVariant = 2:size(mdl, 2)
            srTrialAlignedHat(iTrial, maxTrialSampleLength - trialSampleLength + 1:end, iVariant) = srtHat(inTrial, iVariant);
        end
    end

    % Trial aligned mean spike rate
    t = (-maxTrialSampleLength + 1:0)*median(diff(tt));
    selT = t >= -4;
    t = t(selT);
    msrAcute{iEu} = mean(srTrialAligned(:, selT), 1, 'omitnan');
    msrHatAcute{iEu} = squeeze(mean(srTrialAlignedHat(:, selT, :), 1, 'omitnan'));

    % Onset/peak of trial-averaged msr/msrHat
    [peakAcute(iEu), tPeakAcute(iEu), tOnsetAcute(iEu), mu, sd] = getPeakAndOnset(srTrialAligned(:, selT), t);
    assert(~isnan(peakAcute(iEu)))

    for iVariant = 2:size(mdl, 2)
        [peakHatAcute(iEu, iVariant), tPeakHatAcute(iEu, iVariant), tOnsetHatAcute(iEu, iVariant)] = ...
            getPeakAndOnset(srTrialAlignedHat(:, selT, iVariant), t, mu=mu, sd=sd);
    end
% 
%     if ismember(iEu, find(c.isPressResponsive))
%         cla(ax);
%         plot(t, msr{iEu}, 'k:');
%         ylim('auto')
%         yl = ax.YLim;
%         plot([tPeak(iEu), tPeak(iEu)], yl, 'k:', LineWidth=1)
%         plot([tOnset(iEu), tOnset(iEu)], yl, 'k:', LineWidth=2)
%         for iVariant = 2:nVariants-1
%             plot(t, msrHat{iEu}(:, iVariant), Color=getColor(iVariant-1, nVariants-1))
% %             plot([tPeakHat(iEu), tPeakHat(iEu)], yl, Color=getColor(iVariant-1, nVariants-1), LineStyle=':', LineWidth=1)
%             plot([tOnsetHat(iEu, iVariant), tOnsetHat(iEu, iVariant)], yl, Color=getColor(iVariant-1, nVariants-1), LineStyle=':', LineWidth=2)
%         end
%         ylim(yl)
%     end
end
msrObs = cat(3, msrAcute{:});
msrHatAcute = cat(3, msrHatAcute{:});

clear ax;
%%
save('C:\SERVER\acute_glm_20241023.mat', 'R2', 'aiAcute', 'bootAcute', 'cAcute', 'fallCorrect', 'fallIncorrect', ...
    'mdl', 'modelCriterion', 't', 'msrAcute', 'msrHatAcute', 'msrObs', 'pAcute', 'peakHatAcute', 'srTrialAligned', 'srTrialAlignedHat', 'srt', 'srtHat', 'tOnsetAcute', 'tOnsetHatAcute', 'tPeakAcute', 'tPeakHatAcute');


%%
function [peak, tPeak, tOnset, mu, sd] = getPeakAndOnset(X, t, varargin)
    p = inputParser();
    p.addParameter('onsetThresholdPos', 0.5, @isnumeric);
    p.addParameter('onsetThresholdNeg', 0.25, @isnumeric);
    p.addParameter('signWindow', [-0.5, -0.2], @isnumeric);
    p.addParameter('peakWindow', [-0.5, 0], @isnumeric);
    p.addParameter('mu', []);
    p.addParameter('sd', []);
    p.addParameter('baselineWindow', [-4, -2], @isnumeric);
    p.parse(varargin{:})
    r = p.Results;

    assert(length(t) == size(X, 2));

    % Calculate baseline for normalization
    if isempty(r.mu) || isempty(r.sd)
        baseX = X(:, t >= r.baselineWindow(1) & t <= r.baselineWindow(2));
        mu = mean(baseX, 'all', 'omitnan');
        sd = std(baseX, 0, 'all', 'omitnan');
    else
        mu = r.mu;
        sd = r.sd;
    end

    normX = (X - mu) ./ sd; % Normalize spike rate to baseline
    normX = mean(normX, 1, 'omitnan'); % Then average across trial

    signX = sign(mean(normX(t >= r.signWindow(1) & t <= r.signWindow(2)), 'omitnan'));
    inPeakWindow = find(t >= r.peakWindow(1) & t <= r.peakWindow(2));
    [~, peakIndex] = max(normX(inPeakWindow) .* signX, [], 'omitnan');
    peakIndex = inPeakWindow(peakIndex);
    peak = mean(X(:, peakIndex), 'omitnan');
    tPeak = t(peakIndex);

    if signX > 0
        isAbove = normX(1:peakIndex) >= r.onsetThresholdPos;
    else
        isAbove = normX(1:peakIndex) <= -abs(r.onsetThresholdNeg);
    end
    isAbove = flip(isAbove);
    isAboveConseq = isAbove & [0, diff(isAbove)] == 0;
    onsetIndex = find(~isAboveConseq, 1, 'first') - 1;
    onsetIndex = peakIndex - onsetIndex + 1;
    if isempty(onsetIndex)
        onsetIndex = peakIndex;
    end
    onsetIndex = min(onsetIndex, length(normX));

    tOnset = t(onsetIndex);
    try
        assert(~isempty(tOnset))
        assert(~isempty(tPeak))
    catch
        disp(1)
    end
end

function c = getColor(i, n, maxHue)
    if nargin < 3
        if n <= 4
            c = 'rgbm';
            c = c(i);
            return
        end
        maxHue = 0.8;
    end
    c = hsl2rgb([maxHue*(i-1)./(n-1), 1, 0.5]);
end