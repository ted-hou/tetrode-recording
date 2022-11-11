%% 1. Load CompleteExperiment data
% AnimalInfo
animalInfo = { ...
    'daisy1', 'wt', 'F', -3.2, -1.6, 'tetrode'; ...
    'daisy2', 'wt', 'F', -3.2, +1.6, 'tetrode'; ...
    'daisy3', 'DAT-Cre', 'F', -3.2, +1.6, 'tetrode'; ...
    'desmond10', 'wt', 'M', -3.28, -1.8, 'double-bundle'; ... % -0.962 for other bunder
    'desmond11', 'wt', 'M', -3.28, +1.8, 'double-bundle'; ... % +0.962 for other bunder
    'daisy4', 'D1-Cre', 'F', -3.28, -1.6, 'bundle'; ...
    'daisy5', 'D1-Cre', 'F', -3.28, +1.6, 'bundle'; ...
    'desmond12', 'DAT-Cre', 'M', -3.2, -1.4, 'bundle'; ...
    'desmond13', 'DAT-Cre', 'M', -3.2, +1.4, 'bundle'; ...
    'desmond15', 'wt', 'M', -3.40, -1.5, 'bundle'; ...
    'desmond16', 'wt', 'M', -3.40, +1.5, 'bundle'; ...
    'desmond17', 'wt', 'M', -3.40, +1.5, 'bundle'; ...
    'desmond18', 'wt', 'M', -3.40, +1.5, 'bundle'; ...
    'desmond20', 'A2A-Cre', 'M', -3.28, +1.6, 'bundle'; ...
    'daisy7', 'A2A-Cre', 'F', -3.28, +1.6, 'bundle'; ...
    'desmond21', 'D1-Cre;Dlx-Flp;Ai80', 'M', -3.28, -1.6, 'bundle'; ...
    'desmond22', 'D1-Cre;Dlx-Flp;Ai80', 'M', -3.28, -1.6, 'bundle'; ...
    'daisy8', 'D1-Cre;Dlx-Flp;Ai80', 'F', -3.28, +1.6, 'bundle'; ...
    'daisy9', 'D1-Cre;Dlx-Flp;Ai80', 'F', -3.28, +1.3, '4shank-neuronexus'; ... % 1.3 = center of 4 shanks -4.8DV tip 900um? wide
    'daisy10', 'D1-Cre;Dlx-Flp;Ai80', 'F', -3.28, -1.3, '4shank-neuronexus'; ... % 1.3 = center of 4 shanks -4.8DV tip 900um? wide
    'daisy12', 'wt', 'F', -3.28, +1.3, '4shank-acute-wide'; ... % 1.3 = center of 4 shanks -4.4DV tip 990um? wide
    'daisy13', 'wt', 'F', -3.28, -1.3, '4shank-acute-wide'; ... % 1.3 = center of 4 shanks -4.2DV tip 990um? wide
    'desmond23', 'D1-Cre;Dlx-Flp;Ai80', 'M', -3.28, -1.3, '4shank-acute'; ... % 1.3 = center of 4 shanks, 450um wide
    'daisy14', 'D1-Cre;Dlx-Flp;Ai80', 'F', -3.28, +1.3, '4shank-acute'; ... % 1.3 = center of 4 shanks
    'desmond24', 'A2A-Cre', 'M', -3.28, +1.3, '4shank-acute'; ... % 1.3 = center of 4 shanks
    'desmond25', 'A2A-Cre', 'M', -3.28, -1.3, '4shank-acute'; ... % 1.3 = center of 4 shanks
    'daisy15', 'A2A-Cre', 'F', -3.28, +1.3, '4shank-acute'; ... % 1.3 = center of 4 shanks
    'daisy16', 'A2A-Cre', 'F', -3.28, -1.3, '4shank-acute'; ... % 1.3 = center of 4 shanks
    'desmond26', 'D1-Cre;Dlx-Flp;Ai80', 'M', -3.28, +1.3, '4shank-acute'; ... % 1.3 = center of 4 shanks
    'desmond27', 'D1-Cre;Dlx-Flp;Ai80', 'M', -3.28, -1.3, '4shank-acute'; ... % 1.3 = center of 4 shanks
    };

ai(size(animalInfo, 1)) = struct('name', '', 'strain', '', 'sex', '', 'ap', [], 'ml', [], 'probe', '');
for i = 1:size(animalInfo, 1)
    ai(i).name = animalInfo{i, 1};
    ai(i).strain = animalInfo{i, 2};
    ai(i).sex = animalInfo{i, 3};
    ai(i).ap = animalInfo{i, 4};
    ai(i).ml = animalInfo{i, 5};
    ai(i).probe = animalInfo{i, 6};
end
clear animalInfo

% 1. Load acute EU objects (duplicates already removed)
eu = EphysUnit.load('C:\SERVER\Units\acute_2cam'); 

% Remove multiunit detected by ISI test.
p.ISIThreshold = 0.0015;
for iEu = 1:length(eu)
    st = eu(iEu).SpikeTimes;
    isi = [NaN, diff(st)];
    st(isi == 0) = [];
    isi = [NaN, diff(st)];
    eu(iEu).SpikeTimes = st;
    ISI{iEu} = isi;
end

for iEu = 1:length(eu)
    prcLowISI(iEu) = nnz(ISI{iEu} < p.ISIThreshold) ./ length(ISI{iEu});
end

histogram(prcLowISI, 0:0.01:1)
cat.isMultiUnit = prcLowISI > 0.05;
cat.isSingleUnit = prcLowISI <= 0.05;
eu = eu(cat.isSingleUnit);

eu = eu';

% 1.2. Load Video Tracking Data (vtd) and ArduinoConnection (ac), and group into experiments
exp = CompleteExperiment(eu);

% 1.3 Align video and ephys timestamps
exp.alignTimestamps();

% Find laterality
implantSide = repmat('L', 1, length(exp));
leverSide = repmat('R', 1, length(exp));
[lia, locb] = ismember(exp.animalName, {ai.name});
assert(all(lia))
ai = ai(locb);
implantSide([ai.ml] > 0) = 'R';
leverSide([ai.ml] > 0) = 'L';

clear iEu st isi ISI prcLowISI cat lia locb

%% 1.1. For CompleteExperiment data (per session), plot trial-aligned average (RMS) bodypart velocity traces
%% 1.1.1 Get trial aligned movement velocity data (slow)
p.velETAWindow = [-10, 3];
p.velETABinWidth = 0.025;
p.minTrialLength = 4;
p.maxTrialLength = Inf;

t = flip(p.velETAWindow(2):-p.velETABinWidth:p.velETAWindow(1));
f = cell(1, length(exp));
fnames = {'handL', 'footL', 'handR', 'footR', 'spine', 'nose', 'tongue'};
[kernels, ~, ~] = CompleteExperiment.makeConsineKernels(0, width=0.1); % Kernels for smoothing velocity traces
for iExp = 1:length(exp)
    clear trials
    f{iExp} = struct('press', [], 'lick', [], 'any', []);
    for trialType = {'press', 'lick', {'press', 'lick'}}
        trialType = trialType{1};
        trials = exp(iExp).eu(1).getTrials(trialType);

        if iscell(trialType)
            trialTypeName = 'any';
        else
            trialTypeName = trialType;
        end

        f{iExp}.(trialTypeName) = NaN(length(t), length(fnames), length(trials));
        for iTrial = 1:length(trials)
            if trials(iTrial).duration >= p.minTrialLength && trials(iTrial).duration <= p.maxTrialLength
                tGlobal = flip(trials(iTrial).Stop + p.velETAWindow(2):-p.velETABinWidth:trials(iTrial).Stop + p.velETAWindow(1));
                F = exp(iExp).getFeatures(timestamps=tGlobal, features=fnames, stats={'xPos'}, useGlobalNormalization=true);
                F = CompleteExperiment.convolveFeatures(F, kernels, kernelNames={'_smooth'}, ...
                    features=fnames, ...
                    stats={'xPos'}, ...
                    mode='replace', normalize='none');
%                 inTrial = F.inTrial;
                inTrial = F.t >= trials(iTrial).Start;
                if iTrial < length(trials)
                    inTrial = inTrial & F.t <= trials(iTrial + 1).Start;
                end
                F(:, {'t', 'inTrial'}) = [];
                thisData = table2array(F);
                thisData(~inTrial, :) = NaN;
                f{iExp}.(trialTypeName)(:, :, iTrial) = thisData;
            end
        end
    end
end

% Average by trial
fstats = cell(1, length(exp));
for iExp = 1:length(exp)
    statStruct = struct('mean', [], 'nTrials', [], 'sd', []);
    fstats{iExp} = struct('press', statStruct, 'lick', statStruct, 'any', statStruct);
    for trialTypeName = {'press', 'lick', 'any'}
        trialTypeName = trialTypeName{1};
        fstats{iExp}.(trialTypeName).mean = array2table(mean(f{iExp}.(trialTypeName), 3, 'omitnan'), VariableNames=fnames);
        fstats{iExp}.(trialTypeName).mean.t = t';
        fstats{iExp}.(trialTypeName).nTrials = nnz(any(~isnan(f{iExp}.(trialTypeName)), [1, 2]));
        fstats{iExp}.(trialTypeName).sd = array2table(std(f{iExp}.(trialTypeName), 0, 3, 'omitnan'), VariableNames=fnames);
    end
end

fall = struct('press', statStruct, 'lick', statStruct, 'any', statStruct);
iExp = iExp + 1;
for trialTypeName = {'press', 'lick', 'any'}
    trialTypeName = trialTypeName{1};
    ff = cellfun(@(f) f.(trialTypeName), f, UniformOutput=false);
    fall.(trialTypeName) = cat(3, ff{:});
    fstats{iExp}.(trialTypeName).mean = array2table(mean(fall.(trialTypeName), 3, 'omitnan'), VariableNames=fnames);
    fstats{iExp}.(trialTypeName).mean.t = t';
    fstats{iExp}.(trialTypeName).nTrials = nnz(any(~isnan(fall.(trialTypeName)), [1, 2]));
    fstats{iExp}.(trialTypeName).sd = array2table(std(fall.(trialTypeName), 0, 3, 'omitnan'), VariableNames=fnames);
end

iExp = iExp + 1;
fallL = struct('press', statStruct, 'lick', statStruct, 'any', statStruct);
for trialTypeName = {'press', 'lick', 'any'}
    trialTypeName = trialTypeName{1};
    ff = cellfun(@(f) f.(trialTypeName), f(implantSide=='L'), UniformOutput=false);
    fallL.(trialTypeName) = cat(3, ff{:});
    fstats{iExp}.(trialTypeName).mean = array2table(mean(fallL.(trialTypeName), 3, 'omitnan'), VariableNames=fnames);
    fstats{iExp}.(trialTypeName).mean.t = t';
    fstats{iExp}.(trialTypeName).nTrials = nnz(any(~isnan(fallL.(trialTypeName)), [1, 2]));
    fstats{iExp}.(trialTypeName).sd = array2table(std(fallL.(trialTypeName), 0, 3, 'omitnan'), VariableNames=fnames);
end

iExp = iExp + 1;
fallR = struct('press', statStruct, 'lick', statStruct, 'any', statStruct);
for trialTypeName = {'press', 'lick', 'any'}
    trialTypeName = trialTypeName{1};
    ff = cellfun(@(f) f.(trialTypeName), f(implantSide=='R'), UniformOutput=false);
    fallR.(trialTypeName) = cat(3, ff{:});
    fstats{iExp}.(trialTypeName).mean = array2table(mean(fallR.(trialTypeName), 3, 'omitnan'), VariableNames=fnames);
    fstats{iExp}.(trialTypeName).mean.t = t';
    fstats{iExp}.(trialTypeName).nTrials = nnz(any(~isnan(fallR.(trialTypeName)), [1, 2]));
    fstats{iExp}.(trialTypeName).sd = array2table(std(fallR.(trialTypeName), 0, 3, 'omitnan'), VariableNames=fnames);
end

clear iExp kernels trials trialTypeName trialType iTrial tGlobal inTrial thisData F mu sd n ff

% 1.1.2 Plot average traces
close all
for iExp = 1:length(fstats)
    fig = figure(Position=[0+(iExp-1)*300, 50, 300, 900]);
    axAll = gobjects(1, 3);
    iTrialType = 0;
    for trialTypeName = {'press', 'lick', 'diff'}
        trialTypeName = trialTypeName{1};
        iTrialType = iTrialType + 1;
        ax = subplot(3, 1, iTrialType);
        axAll(iTrialType) = ax;
        hold(ax, 'on')

        switch trialTypeName
            case {'press', 'lick'}
                mu = table2array(fstats{iExp}.(trialTypeName).mean(:, fnames));
                sd = table2array(fstats{iExp}.(trialTypeName).sd(:, fnames));
                n = fstats{iExp}.(trialTypeName).nTrials;
            case 'diff'
                mu = table2array(fstats{iExp}.press.mean(:, fnames)) - table2array(fstats{iExp}.lick.mean(:, fnames));
        end

        h = gobjects(1, length(fnames));
        for iVar = 1:length(fnames)
            col = hsl2rgb([0.8*(iVar-1)/(length(fnames)-1), 1, 0.5]);
            h(iVar) = plot(ax, t, mu(:, iVar), Color=col, LineWidth=1.5, DisplayName=fnames{iVar});
            if ~strcmp(trialTypeName, 'diff')
                sel = ~isnan(mu(:, iVar)+sd(:, iVar));
                patch(ax, [t(sel)'; flip(t(sel)')], [mu(sel, iVar)-sd(sel, iVar); flip(mu(sel, iVar)+sd(sel, iVar))], 'r', ...
                    LineStyle='none', FaceAlpha=0.1, FaceColor=col)
            end
        end
        switch trialTypeName
            case 'press'
                xlabel(ax, 'time to lever contact (s)')
                ylabel(ax, 'z-scored speed (a.u.)')
                trialTypeDispName = 'Press';
            case 'lick'
                xlabel(ax, 'time to spout contact (s)')
                ylabel(ax, 'z-scored speed (a.u.)')
                trialTypeDispName = 'Lick';
            case 'diff'
                xlabel(ax, 'Time to contact (s)')
                ylabel(ax, '\Deltaz-scored speed (a.u.)')
                trialTypeDispName = 'Diff';
                plot(ax, [-3, 3], [0, 0], 'k--')
        end
        plot(ax, [0, 0], [-100, 100], 'k--')
        legend(ax, h, Location='northwest')
        if iExp <= length(implantSide)
            leverSideText = implantSide(iExp);
        elseif iExp == length(implantSide) + 1
            leverSideText = 'both';
        elseif iExp == length(implantSide) + 2
            leverSideText = 'left lever';
        elseif iExp == length(implantSide) + 3
            leverSideText = 'right lever';
        end
        if strcmp(trialTypeName, 'diff')
            title(ax, sprintf('Press minus Lick (%s)', leverSideText));
        else
            title(ax, sprintf('%s trials (N=%d, %s)', trialTypeDispName, n, leverSideText));
        end
        hold(ax, 'off')
    end

%     yl = vertcat(axAll(1:2).YLim);
%     yl = [0, max(yl(:, 2))];
    set(axAll(1:2), YLim=[0, 15])
    set(axAll(3), YLim=[-2.5, 12.5])
    set(axAll, XLim=[-3, 3])
end
clear iTrialType iExp ax fig trialTypeName h iVar axAll n

%% 1.2. For CompleteExperiment data, do GLM (with ITI) and compare weights.
alpha = 0.05;
binEdges = 2:2:10;
nBins = length(binEdges) - 1;

clear bootBTA
for iEu = find(c.isPressResponsive)
    assert(length(iEu) == 1)
    fprintf(1, '%d/%d ', iEu, length(eu))
    if mod(iEu, 15) == 0
        fprintf(1, '\n')
    end
    [xx, tt] = eu(iEu).getTrialAlignedData('count', [-10, 0], 'press', allowedTrialDuration=[0, Inf], alignTo='stop', resolution=0.1, includeInvalid=false);
    dd = eu(iEu).getTrials('press').duration;
    assert(length(dd) == size(xx, 1))
    
    selTrials = dd >= binEdges(1) & dd <= binEdges(end);
    selTime = tt >= -4 & tt <= 0;
    xx = xx(selTrials, selTime);
    dd = dd(selTrials);
    tt = tt(selTime);
    
    [N, ~, bins] = histcounts(dd, binEdges);
    xxMean = NaN(nBins, nnz(selTime));
%     figure(), hold on
    for iBin = 1:length(binEdges) - 1
        xxMean(iBin, :) = mean(xx(bins == iBin, :), 1, 'omitnan');
%         plot(tt, mean(xx(bins == iBin, :), 1, 'omitnan'), Color=hsl2rgb([0.8*(iBin-1)/(nBins-1), 1, 0.5]))
    end
    
    pairs = nchoosek(1:nBins, 2);
    nPairs = size(pairs, 1);
    dist = NaN(nPairs, 1);
    for iPair = 1:nPairs
        x1 = xxMean(pairs(iPair, 1), :);
        x2 = xxMean(pairs(iPair, 2), :);
        sqrs = (x1 - x2).^2;
        valid = ~isnan(sqrs);
        n = nnz(valid);
        dist(iPair) = sqrt(sum(sqrs, 'omitnan')) ./ n;
    end
    
    nboot = 10000;
    bsample = zeros(sum(N), nboot);
    for iboot = 1:nboot
        bsample(:, iboot) = randperm(sum(N));
    end
    
    bootBinEdges = [0, cumsum(N)];
    xxMeanBoot = NaN(nBins, nnz(selTime), nboot);
    for iBin = 1:length(binEdges) - 1
        xxBoot = NaN(N(iBin), nnz(selTime), nboot);
        for iboot = 1:nboot
            selTrialsBoot = bsample(bootBinEdges(iBin)+1:bootBinEdges(iBin+1), iboot);
            xxMeanBoot(iBin, :, iboot) = mean(xx(selTrialsBoot, :), 1, 'omitnan');
            xxBoot(:, :, iboot) = xx(selTrialsBoot, :);
        end
    end
    
    distBoot = NaN(nPairs, nboot);
    distCI = NaN(nPairs, 2);
    distH = NaN(nPairs, 1);
    for iPair = 1:nPairs
        x1 = xxMeanBoot(pairs(iPair, 1), :, :);
        x2 = xxMeanBoot(pairs(iPair, 2), :, :);
        sqrs = (x1 - x2).^2;
        valid = ~isnan(sqrs);
        n = sum(valid, 2);
        distBoot(iPair, :) = sqrt(sum(sqrs, 2, 'omitnan')) ./ n;
        distCI(iPair, :) = prctile(distBoot(iPair, :), [alpha*50/nPairs, 100-alpha*50/nPairs]); % Bonferoni correction: alpha/nPairs
        distH(iPair) = dist(iPair) >= distCI(iPair, 2) || dist(iPair) <= distCI(iPair, 1);
    end
    bootBTA.h{iEu} = distH;
    bootBTA.distCI{iEu} = distCI;
    bootBTA.distObs{iEu} = dist;
end
bootBTA.h = cat(2, bootBTA.h{:});
bootBTA.distObs = cat(2, bootBTA.distObs{:});

clear iEu xx tt dd selTrials selTime N bins xxMean iBin pairs nPairs dist x1 x2 iPair sqrs valid n nboot bsample iboot bootBinEdges xxMeanBoot xxBoot iboot selTrialsBoot distBoot distCI distH

%% 1.2.1 Fit GLMs (very slow)
p.glmMinTrialLength = 0;
p.glmMaxTrialLength = 4;
p.glmITIWindowMargin = [1, 1];

MDL = cell(length(exp), 1);
SR = cell(length(exp), 1);
F = cell(length(exp), 1);

for iExp = 1:length(exp)
    fprintf(1, 'Fitting %g models...\n', length(exp(iExp).eu)); tTicAll = tic();
    thisF = exp(iExp).getFeatures(sampleRate=30, stats={'xVel', 'yVel'}, ...
        features={'handL', 'footL', 'handR', 'footR', 'spine', 'tongue'}, ...
        likelihoodThreshold=0.95);
    [kernels, ~, delays] = CompleteExperiment.makeConsineKernels(4, width=0.1, overlap=0.5, direction='both');
    delays = round(delays * 1000);
    thisF = CompleteExperiment.convolveFeatures(thisF, kernels, kernelNames=delays, ...
        features={'handL', 'footL', 'handR', 'footR', 'spine', 'tongue'}, ...
        stats={'xVel', 'yVel'}, ...
        mode='replace', normalize='maxabs');

    trials = exp(iExp).eu(1).getTrials({'press', 'lick'}, sorted=true);
    trials = trials(trials.duration > p.glmMinTrialLength & trials.duration < p.glmMaxTrialLength);
    start = [trials.Start];
    stop = [trials.Stop];
    iti = trials.iti();
    iti(end) = mean(iti(1:end-1));
    assert(all(isfinite(iti)));
    trials = Trial(stop + p.glmITIWindowMargin(1), stop + iti - p.glmITIWindowMargin(2))';

    tAll = thisF.t;
    inTrial = trials.inTrial(tAll);
    thisF = thisF(inTrial, :);
    thisF.constant = ones(height(thisF), 1);
    predictorVars = thisF.Properties.VariableNames(~ismember(thisF.Properties.VariableNames, {'t', 'inTrial'}));

    mdl = cell(length(exp(iExp).eu), 1);
    sr = cell(length(exp(iExp).eu), 1);
    warning('off','all')
    for iEu = 1:length(exp(iExp).eu)
        fprintf(1, '\tFitting model %g of %g...', iEu, length(exp(iExp).eu)); tTic = tic();
        eu = exp(iExp).eu(iEu);
        srAll = [0, eu.getSpikeRates('gaussian', 0.1, tAll)]';
        sr{iEu} = srAll(inTrial);

        FCopy = thisF;
        FCopy.SpikeRate = double(sr{iEu});

        mdl{iEu} = fitglm(FCopy, ResponseVar='SpikeRate', PredictorVars=predictorVars, ...
            Distribution='poisson');
        fprintf(1, '\tR^2 = %.2f\n', mdl{iEu}.Rsquared.Ordinary);
    end
    warning('on','all')
    warning('query','all')
    fprintf(1, 'Fitted %g units in %.2f seconds.\n', length(exp(iExp).eu), toc(tTicAll));

    MDL{iExp} = cellfun(@compact, mdl, UniformOutput=false);
    SR{iExp} = sr;
    F{iExp} = thisF;
end
mdl = cat(1, MDL{:});
sr = cat(1, SR{:});
expIndices = zeros(size(mdl, 1), 1);
i = 0;
for iExp = 1:length(exp)
    expIndices(i + 1:i + length(exp(iExp).eu)) = iExp;
    i = i + length(exp(iExp).eu);
end
clear MDL SR i iExp FCopy trials inTrial t iti stop thisF kernels delays tAll srAll

% 1.2.2 Plot predictor weights and fitted vs. observed spike rates for individual SNr units
close all

R2 = cellfun(@(mdl) mdl.Rsquared.ordinary, mdl);
figure; histogram(R2, 0:0.02:1);

try
    eu = [exp.eu];
catch
    eu = vertcat(exp.eu);
end

for iEu = 1:length(mdl)
    iExp = expIndices(iEu);
    thisF = F{iExp};

    fig = figure(Units='normalized', OuterPosition=[0, 0, 0.8, 0.5]);

    ax = subplot(2, 1, 1);
    hold(ax, 'on')
    plot((1:length(sr{iEu}))./30, sr{iEu}, 'k', LineWidth=1.2, DisplayName='Observed');
    yHat = predict(mdl{iEu}, F{expIndices(iEu)});
    plot((1:length(sr{iEu}))./30, yHat, 'r', LineWidth=1.2, DisplayName=sprintf('Fitted (R^2=%.2f)', mdl{iEu}.Rsquared.Ordinary));

    hold(ax, 'off')
    title(ax, eu(iEu).getName('_'), Interpreter='none');
    xlabel(ax, 'Time (s)')
    ylabel(ax, 'Spike rate (sp/s)')
    xlim(ax, [30, 300])
    legend(ax, Location='northwest');

    ax = subplot(2, 1, 2);
    ax.TickLabelInterpreter = 'none';
    hold(ax, 'on')
    x = 2:height(mdl{iEu}.Coefficients);
    errorbar(ax, x, mdl{iEu}.Coefficients.Estimate(2:end), mdl{iEu}.Coefficients.SE(2:end));
    plot(ax, [x(1), x(end)], [0, 0], 'k--', LineWidth=1.5)
    hold(ax, 'off')
    xticks(ax, x);
    xticklabels(ax, mdl{iEu}.CoefficientNames(x));
    xlim(ax, [x(1), x(end)])
    ylabel(ax, 'Coefficient +/- SE')
    title(ax, 'Coefficients')
    set(ax, FontSize=9)

    if ~isfolder('C:\\SERVER\\Figures\\GLM_ITI')
        mkdir('C:\\SERVER\\Figures\\GLM_ITI')
    end
    print(fig, sprintf('C:\\SERVER\\Figures\\GLM_ITI\\%s', eu(iEu).getName('_')), '-dpng');
    close(fig)
end
clear iEu iExp thisF fig ax yHat x

% 1.2.1 Show correlation matrix between predictors
% 1.2.2 Categorize SNr units by selectivity using GLM weights (orofacial vs.
% limb/trunk)
% 1.2.3 Confirm contralateral selectivity by comparing left vs. right limb
% weights

%% 2.1. For all animals (with lick tube servo), compare press trial ETA vs. lick trial ETA traces
%% 2.2. For all units with osci-lick spiking, compare press trial ETA vs. lick trial ETA traces. NULL: latter is a time shifted version of the former
close all
fig = figure(Units='inches', Position=[0 0 6 9], DefaultAxesFontSize=11);

SEL = { ...
        c.hasPress & c.hasLick; ...
        c.hasPress & c.hasLick & (c.isPressUp & c.isLickUp); ...
        c.hasPress & c.hasLick & (c.isPressDown & c.isLickDown); ...
        c.hasPress & c.hasLick & (c.isPressUp & c.isLickDown); ...
        c.hasPress & c.hasLick & (c.isPressDown & c.isLickUp); ...
    };

nRows = length(SEL);
ax = gobjects(nRows, 2);
for i = 1:nRows
    ax(i, 1) = subplot(nRows, 2, 2*(i-1)+1);
    plotETAComparison(ax(i, 1), eta.press, eta.lick, SEL{i}, 'press', 'lick', yUnit='a.u.', error='sd');
    ax(i, 2) = subplot(nRows, 2, 2*(i-1)+2);
    plotETAComparison(ax(i, 2), eta.press, eta.lick, SEL{i} & c.isLick, 'press', 'lick', yUnit='a.u.', error='sd');
end
xlabel(ax, '')
clear ax fig i



%% 3. Compare ETA (binned by trial length) (euclidean distance, bootstrap). To show there is no baseline/ramp slope/peak differences

[bta.pressUpRaw.X, bta.pressUpRaw.T, bta.pressUpRaw.N, bta.pressUpRaw.S, bta.pressUpRaw.B] = eu(c.isPressUp).getBinnedTrialAverage('rate', p.binnedTrialEdges, 'press', 'window', [-10, 1], 'normalize', false);
[bta.pressDownRaw.X, bta.pressDownRaw.T, bta.pressDownRaw.N, bta.pressDownRaw.S, bta.pressDownRaw.B] = eu(c.isPressDown).getBinnedTrialAverage('rate', p.binnedTrialEdges, 'press', 'window', [-10, 1], 'normalize', false);
% [bta.pressRaw.X, bta.pressRaw.T, bta.pressRaw.N, bta.pressRaw.S, bta.pressRaw.B] = eu(c.hasPress).getBinnedTrialAverage('rate', p.binnedTrialEdges, 'press', 'window', [-10, 1], 'normalize', false);

[bta.lickUpRaw.X, bta.lickUpRaw.T, bta.lickUpRaw.N, bta.lickUpRaw.S, bta.lickUpRaw.B] = eu(c.isLickUp).getBinnedTrialAverage('rate', p.binnedTrialEdges, 'lick', 'window', [-10, 1], 'normalize', false);
[bta.lickDownRaw.X, bta.lickDownRaw.T, bta.lickDownRaw.N, bta.lickDownRaw.S, bta.lickDownRaw.B] = eu(c.isLickDown).getBinnedTrialAverage('rate', p.binnedTrialEdges, 'lick', 'window', [-10, 1], 'normalize', false);
% [bta.lickRaw.X, bta.lickRaw.T, bta.lickRaw.N, bta.lickRaw.S, bta.lickRaw.B] = eu(c.hasLick).getBinnedTrialAverage('rate', p.binnedTrialEdges, 'lick', 'window', [-10, 1], 'normalize', false);



%% 4. Use bootstraping to find significant movement responses
clear boot
p.bootAlpha = 0.01;
boot.press = struct('h', NaN(length(eu), 1), 'muDiffCI', NaN(length(eu), 2), 'muDiffObs', NaN(length(eu), 1));
boot.lick = struct('h', NaN(length(eu), 1), 'muDiffCI', NaN(length(eu), 2), 'muDiffObs', NaN(length(eu), 1));
[boot.press.h(c.hasPress), boot.press.muDiffCI(c.hasPress, :), boot.press.muDiffObs(c.hasPress)] = bootstrapMoveResponse( ...
    eu(c.hasPress), 'press', alpha=p.bootAlpha, withReplacement=false, oneSided=false, ...
    responseWindow=[-0.5, -0.2]);
[boot.lick.h(c.hasLick), boot.lick.muDiffCI(c.hasLick, :), boot.lick.muDiffObs(c.hasLick)] = bootstrapMoveResponse( ...
    eu(c.hasLick), 'lick', alpha=p.bootAlpha, withReplacement=false, oneSided=false, ...
    responseWindow=[-0.3, 0]);
fprintf(1, '\nAll done\n')

%% Report bootstraped movement response direction
assert(nnz(isnan(boot.lick.h(c.hasLick))) == 0)
assert(nnz(isnan(boot.press.h(c.hasPress))) == 0)

figure, histogram(boot.press.h)
c.isPressUp = boot.press.h' == 1 & c.hasPress;
c.isPressDown = boot.press.h' == -1 & c.hasPress;
c.isPressResponsive = c.isPressUp | c.isPressDown;

figure, histogram(boot.lick.h)
c.isLickUp = boot.lick.h' == 1 & c.hasLick;
c.isLickDown = boot.lick.h' == -1 & c.hasLick;
c.isLickResponsive = c.isLickUp | c.isLickDown;

fprintf(1, ['%g total SNr units (baseline spike rate > %g):\n' ...
    '\t%g with %d+ press trials;\n' ...
    '\t%g with %d+ lick trials;\n' ...
    '\t%g with either (%g+ trials);\n' ...
    '\t%g with both (%g+ trials).\n'], ...
    length(eu), p.minSpikeRate, nnz(c.hasPress), p.minNumTrials, ...
    nnz(c.hasLick), p.minNumTrials, ...
    nnz(c.hasPress | c.hasLick), p.minNumTrials, ...
    nnz(c.hasPress & c.hasLick), p.minNumTrials)

fprintf(1, ['%g units with %g+ press trials (%gs or longer):\n' ...
    '\t%g (%.0f%%) are excited (p<%g);\n' ...
    '\t%g (%.0f%%) are inhibited (p<%g).\n'], ...
    nnz(c.hasPress), p.minNumTrials, p.minTrialDuration, ...
    nnz(c.isPressUp), 100*nnz(c.isPressUp)/nnz(c.isPressResponsive), p.bootAlpha, ...
    nnz(c.isPressDown), 100*nnz(c.isPressDown)/nnz(c.isPressResponsive), p.bootAlpha);

fprintf(1, ['%g units with %g+ lick trials (%gs or longer):\n' ...
    '\t%g (%.0f%%) are excited (p<%g);\n' ...
    '\t%g (%.0f%%) are inhibited (p<%g).\n'], ...
    nnz(c.hasLick), p.minNumTrials, p.minTrialDuration, ...
    nnz(c.isLickUp), 100*nnz(c.isLickUp)/nnz(c.isLickResponsive), p.bootAlpha, ...
    nnz(c.isLickDown), 100*nnz(c.isLickDown)/nnz(c.isLickResponsive), p.bootAlpha);

nTotal = nnz(c.isPressResponsive & c.isLickResponsive);
fprintf(1, ['%g units with %d+ press AND lick trials (%gs or longer):\n' ...
    '\t%g (%.0f%%) are press-excited AND lick-excited;\n' ...
    '\t%g (%.0f%%) are press-inhibited AND lick-inhibited;\n' ...
    '\t%g (%.0f%%) are press-excited AND lick-inhibited;\n' ...
    '\t%g (%.0f%%) are press-inhibited AND lick-excited;\n'], ...
    nnz(c.hasPress & c.hasLick), p.minNumTrials, p.minTrialDuration, ...
    nnz(c.isPressUp & c.isLickUp), 100*nnz(c.isPressUp & c.isLickUp)/nTotal, ...
    nnz(c.isPressDown & c.isLickDown), 100*nnz(c.isPressDown & c.isLickDown)/nTotal, ...
    nnz(c.isPressUp & c.isLickDown), 100*nnz(c.isPressUp & c.isLickDown)/nTotal, ...
    nnz(c.isPressDown & c.isLickUp), 100*nnz(c.isPressDown & c.isLickUp)/nTotal)   
% clear nTotal
%% 5. Redo movement onset
p.etaSortWindow = [-3, 0];
p.etaSignWindow = [-0.3, 0];
p.etaLatencyThresholdPos = 0.5;
p.etaLatencyThresholdNeg = 0.25;

close all
llatency = NaN(size(c.hasPress));
[ax, ~, ~, llatency(c.hasPress)] = EphysUnit.plotETA(eta.press, c.hasPress, xlim=[-4,0], clim=[-2, 2], ...
    sortWindow=p.etaSortWindow, signWindow=p.etaSignWindow, ...
    sortThreshold=p.etaLatencyThresholdPos, negativeSortThreshold=p.etaLatencyThresholdNeg); 
title('Lever-touch PETH')
xlabel('Time to lever-contact (s)')
ax.Parent.Position(3) = 0.25;

figure(DefaultAxesFontSize=11, Position=[200,200,600,300])
histogram(llatency(c.isPressResponsive), -2:0.1:0, Normalization='probability')
title('Neural activity onset latency')
legend(sprintf('%d units', nnz(~isnan(llatency(c.isPressResponsive & c.hasPress)))), Location='northwest')
xlabel('Time to lever-contact (s)')
ylabel('Probability')

% [~, ~, ~, llatency2] = EphysUnit.plotETA(eta.press, c.isPressResponsive & c.hasPress, xlim=[-4,0], clim=[-2, 2], ...
%     sortWindow=p.etaSortWindow, signWindow=p.etaSignWindow, ...
%     sortThreshold=p.etaLatencyThresholdPos, negativeSortThreshold=p.etaLatencyThresholdNeg); 
% 
% figure(DefaultAxesFontSize=11, Position=[200,200,600,300])
% histogram(llatency2, -2:0.1:0, Normalization='probability')
% title('Neural activity onset latency')
% legend(sprintf('%d units', nnz(~isnan(llatency2))), Location='northwest')
% xlabel('Time to lever-contact (s)')
% ylabel('Probability')

% Test differect between paw/neural onset latency

clear latency
latency.neural = llatency;
latency.contraPaw = tst;
latency.pRankSum = ranksum(tst, llatency, tail='right');

fprintf(1, 'Median neural onset latency = %.1f ms \n', median(latency.neural*1000, 'omitnan'))
fprintf(1, 'Median contralateral paw movement onset latency = %.3f ms \n', median(latency.contraPaw*1000, 'omitnan'))
fprintf(1, 'One-tailed ranksum test p = %g\n', latency.pRankSum)

clear llatency ax


%% 6.1 Redo GLM training on [-10, -0.75] window, predict [-10, 0] or [-10, -0.75]
%% 6.2 Compare predicted vs. fitted ramp onset/peak, show beside \deltaR^2

%% Functions
function [h, muDiffCI, muDiffObs] = bootstrapMoveResponse(eu, trialType, varargin)
    p = inputParser();
    p.addRequired('eu', @(x) length(x) >= 1 && isa(x, 'EphysUnit'));
    p.addRequired('trialType', @(x) ismember(x, {'press', 'lick'}));
    p.addParameter('nboot', 10000, @isnumeric)
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

function [ax, fig] = plotETAComparison(varargin) % NonHeatmap
    p = inputParser();
    if isgraphics(varargin{1}, 'Axes')
        p.addRequired('ax')
    end
    p.addRequired('eta1', @isstruct)
    p.addRequired('eta2', @isstruct)
    p.addOptional('sel', [], @(x) isnumeric(x) || islogical(x))
    p.addOptional('label1', '', @ischar)
    p.addOptional('label2', '', @ischar)
    p.addParameter('xlim', [], @(x) isnumeric(x) || (iscell(x) && all(cellfun(@isnumeric, x))))
    p.addParameter('averageAcrossUnits', true, @islogical) % False to generate individual plots for each unit
    p.addParameter('yUnit', 'a.u.', @(x) ismember(x, {'a.u.', 'sp/s'}))
    p.addParameter('error', 'none', @(x) ismember(x, {'sd', 'sem', 'none'}))
    p.parse(varargin{:});
    r = p.Results;
    eta1 = r.eta1;
    eta2 = r.eta2;
    sel = r.sel;
    label1 = r.label1;
    label2 = r.label2;

    if ~isfield(r, 'ax')
        fig = figure();
        ax = axes(fig);
    else
        ax = r.ax;
    end

    t1 = eta1.t;
    t2 = eta2.t;
    if ~isempty(sel)
        X1 = eta1.X(sel, :);
        X2 = eta2.X(sel, :);
        N = nnz(sel);
    else
        X1 = eta1.X;
        X2 = eta2.X;
        N = size(X1, 1);
    end

    if r.averageAcrossUnits
        x1 = mean(X1, 1, 'omitnan');
        x2 = mean(X2, 1, 'omitnan');

        switch r.error
            case 'none'
                s1 = [];
                s2 = [];
            case 'sd'
                s1 = std(X1, 0, 1, 'omitnan');
                s2 = std(X2, 0, 1, 'omitnan');
            case 'sem'
                s1 = std(X1, 0, 1, 'omitnan') ./ sqrt(size(X1, 1));
                s2 = std(X2, 0, 1, 'omitnan') ./ sqrt(size(X2, 1));
        end

        hold(ax, 'on')
        switch r.yUnit
            case 'a.u.'
                ylabel('Normalized spike rate (a.u.)')
            case 'sp/s'
                x1 = x1 ./ 0.1;
                x2 = x2 ./ 0.1;
                s1 = s1 ./ 0.1;
                s2 = s2 ./ 0.1;
                ylabel('Spike rate (sp/s)')
        end
        h(1) = plot(ax, t1, x1, 'r', DisplayName=label1);
        h(2) = plot(ax, t2, x2, 'b', DisplayName=label2);
        if ~isempty(s1)
            patch(ax, [t1, flip(t1)], [x1 - s1, flip(x1 + s1)], 'r', FaceAlpha=0.2, LineStyle='none', EdgeColor='r');
        end
        if ~isempty(s2)
            patch(ax, [t2, flip(t2)], [x2 - s2, flip(x2 + s2)], 'b', FaceAlpha=0.2, LineStyle='none', EdgeColor='b');
        end
        hold(ax, 'off')
        if ~isempty(r.xlim)
            xlim(ax, r.xlim)
        end
        legend(ax, h, Location='northwest')
        xlabel('Time to lever/spout contact (s)')
        title(sprintf('Average across %d units', N))
    end
end
