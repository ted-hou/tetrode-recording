% Load acute EU objects (duplicates already removed)
euAcute = EphysUnit.load('C:\SERVER\Units\acute_2cam'); 
%%
clear ISI isi st prcLowISI cat
% Remove multiunit detected by ISI test.
p.ISIThreshold = 0.0015;
for iEu = 1:length(euAcute)
    st = euAcute(iEu).SpikeTimes;
    isi = [NaN, diff(st)];
    st(isi == 0) = [];
    isi = [NaN, diff(st)];
    euAcute(iEu).SpikeTimes = st;
    ISI{iEu} = isi;
end

for iEu = 1:length(euAcute)
    prcLowISI(iEu) = nnz(ISI{iEu} < p.ISIThreshold) ./ length(ISI{iEu});
end

histogram(prcLowISI, 0:0.01:1)
cat.isMultiUnit = prcLowISI > 0.05;
cat.isSingleUnit = prcLowISI <= 0.05;
euAcute = euAcute(cat.isSingleUnit);

euAcute = euAcute';
clear cat

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
p.velETABinWidth = 0.05;
p.minTrialLength = 2;

close all
% statName = {'spd'};
statName = {'xVel', 'yVel', 'spd'}; % Side camera, xVel along AP, yVel along DV.
statNameDisp = {'AP', 'DV', 'spd'};
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
[kernels, ~, ~] = CompleteExperiment.makeConsineKernels(0, width=0.1); % Kernels for smoothing velocity traces
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
        trials = expAcute(iExp).eu(1).getTrials(trialType);
        if strcmp(trialType, 'press')
            disp(length(trials))
        end

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
%             F = CompleteExperiment.convolveFeatures(F, kernels, kernelNames={'_smooth'}, ...
%                 features=theseNames, ...
%                 stats={statName}, ...
%                 mode='replace', normalize='none');
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

