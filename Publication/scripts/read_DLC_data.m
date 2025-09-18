% Load acute EU objects (duplicates already removed)
euAcuteNames = dir('C:\SERVER\Units\acute_2cam\*.mat');
euAcuteNames = arrayfun(@(x) strrep(x, '.mat', ''), string({euAcuteNames.name}));

euAcute = eu(ismember(eu.getName(), euAcuteNames));

% This is the list of good single units, without duplicates, no drifting
% assert(length(eu) == 814);
% cat.isGoodUnit = ismember(euAcute.getName(), eu.getName());
% euAcute = euAcute(cat.isGoodUnit);
% 
% clear cat
clear euAcuteNames

%% 1.2. Load Video Tracking Data (vtd) and ArduinoConnection (ac), and group into experiments
expAcute = CompleteExperiment(euAcute);

% 1.3 Align video and ephys timestamps
expAcute.alignTimestamps();

%% Find laterality
implantSide = repmat('L', 1, length(expAcute));
leverSide = repmat('R', 1, length(expAcute));
[lia, locb] = ismember(expAcute.animalName, {ai.name});
assert(all(lia))
aiAcute = ai(locb);
implantSide([aiAcute.ml] > 0) = 'R';
leverSide([aiAcute.ml] > 0) = 'L';
%% 1.1.1 Get trial aligned movement velocity data (slow)
p.velETAWindow = [-10, 3];
p.velETABinWidth = 0.025;
p.minTrialLength = 2;

close all
% statName = {'spd'};
statName = {'xPos', 'yPos', 'xVel', 'yVel', 'spd'}; % Side camera, xVel along AP, yVel along DV.
statNameDisp = {'AP pos', 'DV pos', 'AP vel', 'DV vel', 'spd'};
t = flip(p.velETAWindow(2):-p.velETABinWidth:p.velETAWindow(1));
fCorrect = cell(1, length(expAcute));
fIncorrect = cell(1, length(expAcute));
fnamesL = {'handL', 'footL', 'handR', 'footR', 'spine', 'nose', 'tongue'};
fnamesR = {'handR', 'footR', 'handL', 'footL', 'spine', 'nose', 'tongue'};
fnames = {'handContra', 'footContra', 'handIpsi', 'footIpsi', 'spine', 'nose', 'tongue'};
assert(strcmpi(fnames{end}, 'tongue'))
fnames = cellfun(@(f) cellfun(@(stat) [f, '_', stat], statName, UniformOutput=false)', fnames(1:end-1), UniformOutput=false);
fnames = [cat(1, fnames{:})', {'tongue'}];
fnamesDisp = {'contra hand', 'contra foot', 'ipsi hand', 'ipsi foot', 'spine', 'nose', 'tongue'};
assert(strcmpi(fnamesDisp{end}, 'tongue'))
fnamesDisp = cellfun(@(f) cellfun(@(stat) [f, ' ', stat], statNameDisp, UniformOutput=false)', fnamesDisp(1:end-1), UniformOutput=false);
fnamesDisp = [cat(1, fnamesDisp{:})', {'tongue'}];
% [kernels, ~, ~] = CompleteExperiment.makeConsineKernels(0, width=0.1, sampleRate=1./p.velETABinWidth); % Kernels for smoothing velocity traces
for iExp = 1:length(expAcute)
    clear trials
    switch leverSide(iExp)
        case 'L'
            theseNames = fnamesL;
%             theseNamesSmooth = fnamesSmoothL;
        case 'R'
            theseNames = fnamesR;
%             theseNamesSmooth = fnamesSmoothR;
    end

    fCorrect{iExp} = struct('press', [], 'lick', []);
    fIncorrect{iExp} = struct('press', [], 'lick', []);
    for trialType = {'press', 'lick'}
        trialType = trialType{1};
        switch trialType
            case {'press', 'lick'}
                trials = expAcute(iExp).eu(1).getTrials(trialType);
                nCorrectTrials = nnz(trials.duration >= 4);
                nIncorrectTrials = nnz(trials.duration < 4 & trials.duration >= p.minTrialLength);
                fIncorrect{iExp}.(trialType) = NaN(length(t), length(fnames), nIncorrectTrials);
                fCorrect{iExp}.(trialType) = NaN(length(t), length(fnames), nCorrectTrials);

                iTrialCorrect = 0;
                iTrialIncorrect = 0;
                for iTrial = 1:length(trials)
                    if trials(iTrial).duration < p.minTrialLength
                        continue;
                    end
        
                    tGlobal = flip(trials(iTrial).Stop + p.velETAWindow(2):-p.velETABinWidth:trials(iTrial).Stop + p.velETAWindow(1));
                    F = expAcute(iExp).getFeatures(timestamps=tGlobal, features=theseNames, stats=statName, useGlobalNormalization=true);
                    % F = CompleteExperiment.convolveFeatures(F, kernels, kernelNames={'_smooth'}, ...
                    %     features=theseNames, ...
                    %     stats=statName, ...
                    %     mode='replace', normalize='none');
                    inTrial = F.t >= trials(iTrial).Start;
                    if iTrial < length(trials)
                        inTrial = inTrial & F.t <= trials(iTrial + 1).Start;
                    end
                    F(:, {'t', 'inTrial'}) = [];
                    thisData = table2array(F);
                    thisData(~inTrial, :) = NaN;
        
                    % Incorrect
                    if trials(iTrial).duration < 4
                        iTrialIncorrect = iTrialIncorrect + 1;
                        fIncorrect{iExp}.(trialType)(:, :, iTrialIncorrect) = thisData;
                    % Correct
                    else
                        iTrialCorrect = iTrialCorrect + 1;
                        fCorrect{iExp}.(trialType)(:, :, iTrialCorrect) = thisData;
                    end
                end
            case {'press_release', 'press_retract'}
                correctTrials = expAcute(iExp).eu(1).getTrials(sprintf('%s_correct', trialType));
                incorrectTrials = expAcute(iExp).eu(1).getTrials(sprintf('%s_incorrect', trialType));
                nCorrectTrials = length(correctTrials);
                nIncorrectTrials = length(incorrectTrials);
                fIncorrect{iExp}.(trialType) = NaN(length(t), length(fnames), nIncorrectTrials);
                fCorrect{iExp}.(trialType) = NaN(length(t), length(fnames), nCorrectTrials);


                for iTrial = 1:length(correctTrials)
                    tGlobal = flip(correctTrials(iTrial).Stop + p.velETAWindow(2):-p.velETABinWidth:correctTrials(iTrial).Stop + p.velETAWindow(1));
                    F = expAcute(iExp).getFeatures(timestamps=tGlobal, features=theseNames, stats=statName, useGlobalNormalization=true);
                    inTrial = F.t >= correctTrials(iTrial).Start;
                    if iTrial < length(correctTrials)
                        inTrial = inTrial & F.t <= correctTrials(iTrial + 1).Start;
                    end
                    F(:, {'t', 'inTrial'}) = [];
                    thisData = table2array(F);
                    thisData(~inTrial, :) = NaN;
                    fCorrect{iExp}.(trialType)(:, :, iTrial) = thisData;
                end
                for iTrial = 1:length(incorrectTrials)
                    tGlobal = flip(incorrectTrials(iTrial).Stop + p.velETAWindow(2):-p.velETABinWidth:incorrectTrials(iTrial).Stop + p.velETAWindow(1));
                    F = expAcute(iExp).getFeatures(timestamps=tGlobal, features=theseNames, stats=statName, useGlobalNormalization=true);
                    inTrial = F.t >= incorrectTrials(iTrial).Start;
                    if iTrial < length(incorrectTrials)
                        inTrial = inTrial & F.t <= incorrectTrials(iTrial + 1).Start;
                    end
                    F(:, {'t', 'inTrial'}) = [];
                    thisData = table2array(F);
                    thisData(~inTrial, :) = NaN;
                    fIncorrect{iExp}.(trialType)(:, :, iTrial) = thisData;
                end
        end
    end
end

% Average by trial
% fnames = F.Properties.VariableNames;
clear fstats
statStruct = struct('mean', [], 'nTrials', [], 'sd', []);
fallIncorrect = struct('press', statStruct, 'lick', statStruct);
for trialTypeName = {'press', 'lick'}
    trialTypeName = trialTypeName{1};
    ff = cellfun(@(f) f.(trialTypeName), fIncorrect, UniformOutput=false);
    fallIncorrect.(trialTypeName) = cat(3, ff{:});
    fstats{1}.(trialTypeName).mean = array2table(mean(fallIncorrect.(trialTypeName), 3, 'omitnan'), VariableNames=fnames);
    fstats{1}.(trialTypeName).mean.t = t';
    fstats{1}.(trialTypeName).nTrials = nnz(any(~isnan(fallIncorrect.(trialTypeName)), [1, 2]));
    fstats{1}.(trialTypeName).sd = array2table(std(fallIncorrect.(trialTypeName), 0, 3, 'omitnan'), VariableNames=fnames);
end

fallCorrect = struct('press', statStruct, 'lick', statStruct);
for trialTypeName = {'press', 'lick'}
    trialTypeName = trialTypeName{1};
    ff = cellfun(@(f) f.(trialTypeName), fCorrect, UniformOutput=false);
    fallCorrect.(trialTypeName) = cat(3, ff{:});
    fstats{2}.(trialTypeName).mean = array2table(mean(fallCorrect.(trialTypeName), 3, 'omitnan'), VariableNames=fnames);
    fstats{2}.(trialTypeName).mean.t = t';
    fstats{2}.(trialTypeName).nTrials = nnz(any(~isnan(fallCorrect.(trialTypeName)), [1, 2]));
    fstats{2}.(trialTypeName).sd = array2table(std(fallCorrect.(trialTypeName), 0, 3, 'omitnan'), VariableNames=fnames);
end

clear iExp kernels trials trialTypeName trialType iTrial tGlobal inTrial thisData F mu sd n ff

%% Determine movement onset time (contra lateral paw, 0.25SD [50ms baseline, 50ms deviation] at 1ms res)
onsetThreshold = 0.25;
% onsetPattern = [zeros(1, 5), ones(1, 5)];
onsetPattern = [0, 1, 1];
onsetOffset = find(onsetPattern, 1, 'first');
t = flip(p.velETAWindow(2):-p.velETABinWidth:p.velETAWindow(1));

assert(length(t) == size(fallIncorrect.press, 1))
clear fAll
fAll.pressArray = cat(3, fallIncorrect.press, fallCorrect.press);
fAll.lickArray = cat(3, fallIncorrect.lick, fallCorrect.lick);
fAll.press.t = t;
fAll.lick.t = t;
fAll.press.correct = vertcat(false(size(fallIncorrect.press, 3), 1), true(size(fallCorrect.press, 3), 1));
fAll.lick.correct = vertcat(false(size(fallIncorrect.lick, 3), 1), true(size(fallCorrect.lick, 3), 1));
for iFeat = 1:length(fnames)
    fAll.press.(fnames{iFeat}) = transpose(squeeze(fAll.pressArray(:, iFeat, :)));
    fAll.lick.(fnames{iFeat}) = transpose(squeeze(fAll.lickArray(:, iFeat, :)));
end

% Press onset
nTrials = length(fAll.press.correct);
fAll.press.onset = NaN(nTrials, 1);
t = fAll.press.t;
dist = zeros(nTrials, length(t));
for iTrial = 1:nTrials
    x = fAll.press.handContra_xPos(iTrial, :);
    y = fAll.press.handContra_yPos(iTrial, :);
    x(isnan(x)) = interp1(t(~isnan(x)), x(~isnan(x)), t(isnan(x)), 'linear');
    y(isnan(y)) = interp1(t(~isnan(y)), y(~isnan(y)), t(isnan(y)), 'linear');
    pos = [ ...
            x; ...
            y; ...
        ];
    dist(iTrial, :) = sqrt(sum(pos.^2, 1));
end
distNorm = dist;
% distNorm = (dist - median(dist(:), 'omitnan'))./(mad(dist(:), 1)./0.6745); % Normalize using session median/mad
selT = t <= 0 & t >= -3;
distNorm = distNorm(:, selT);
t = t(selT);
for iTrial = 1:nTrials
    % if any(isnan(distNorm(iTrial, :)))
    %     continue
    % end

    isAbove = distNorm(iTrial, :) >= onsetThreshold;
    iLastAbove = find(isAbove, 1, 'last'); % We don't do abs since dist is already positive, although z-scoring will generate negatives, true movement should be positive.
    if isempty(iLastAbove)
        continue
    end
    iOnset = strfind(isAbove(1:iLastAbove), onsetPattern) + onsetOffset - 1;
    if isempty(iOnset)
        continue
    end
    iOnset = iOnset(end);
    fAll.press.onset(iTrial) = t(iOnset);
end

figure(1)
histogram(fAll.press.onset, -4:0.1:0)
title(sprintf('%i//%i', nnz(~isnan(fAll.press.onset)), nTrials))