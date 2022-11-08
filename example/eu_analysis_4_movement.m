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
[lia, locb] = ismember(exp.animalName, {ai.name});
assert(all(lia))
ai = ai(locb);
implantSide([ai.ml] > 0) = 'R';

clear iEu st isi ISI prcLowISI cat lia locb

%% 1.1. For CompleteExperiment data (per session), plot trial-aligned average (RMS) bodypart velocity traces
%% 1.1.1 Get trial aligned movement velocity data (slow)
p.velETAWindow = [-10, 0];
p.velETABinWidth = 0.025;
p.minTrialLength = 2;

t = flip(p.velETAWindow(2):-p.velETABinWidth:p.velETAWindow(1));
f = cell(1, length(exp));
fnames = {'handL', 'footL', 'handR', 'footR', 'nose', 'spine'};
[velocityKernels, ~, ~] = CompleteExperiment.makeConsineKernels(0, width=0.2); % Kernels for smoothing velocity traces
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

        f{iExp}.(trialTypeName) = NaN(length(t), 6, length(trials));
        for iTrial = 1:length(trials)
            if trials(iTrial).duration >= p.minTrialLength
                tGlobal = flip(trials(iTrial).Stop + p.velETAWindow(2):-p.velETABinWidth:trials(iTrial).Stop + p.velETAWindow(1));
                F = exp(iExp).getFeatures(timestamps=tGlobal, features=fnames, stats={'spd'});
                F = CompleteExperiment.convolveFeatures(F, velocityKernels, kernelNames={'_smooth'}, ...
                    features=fnames, ...
                    stats={'spd'}, ...
                    mode='replace', normalize='none');
                inTrial = F.inTrial;
                F(:, {'t', 'inTrial'}) = [];
                thisData = table2array(F);
                thisData(~inTrial, :) = NaN;
                f{iExp}.(trialTypeName)(:, :, iTrial) = thisData;
            end
        end
    end
end

%% Average by trial
fstats = cell(1, length(exp));
for iExp = 1:length(exp)
    statStruct = struct('mean', [], 'nTrials', [], 'sd', []);
    fstats{iExp} = struct('press', statStruct, 'lick', statStruct, 'any', statStruct);
    for trialTypeName = {'press', 'lick', 'any'}
        trialTypeName = trialTypeName{1};
        fstats{iExp}.(trialTypeName).mean = array2table(mean(f{iExp}.(trialTypeName), 3, 'omitnan'), VariableNames=fnames);
        fstats{iExp}.(trialTypeName).mean.t = t';
        fstats{iExp}.(trialTypeName).nTrials = size(f{iExp}.(trialTypeName), 3);
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
    fstats{iExp}.(trialTypeName).nTrials = size(fall.(trialTypeName), 3);
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
    fstats{iExp}.(trialTypeName).nTrials = size(fallL.(trialTypeName), 3);
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
    fstats{iExp}.(trialTypeName).nTrials = size(fallR.(trialTypeName), 3);
    fstats{iExp}.(trialTypeName).sd = array2table(std(fallR.(trialTypeName), 0, 3, 'omitnan'), VariableNames=fnames);
end

clear iExp velocityKernels trials trialTypeName trialType iTrial tGlobal inTrial thisData F mu sd n ff

%% 1.1.2 Plot average traces
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
        xlim(ax, [-4, 0])
        switch trialTypeName
            case 'press'
                xlabel(ax, 'time to lever contact (s)')
                ylabel(ax, 'speed (a.u.)')
            case 'lick'
                xlabel(ax, 'time to spout contact (s)')
                ylabel(ax, 'speed (a.u.)')
            case 'diff'
                xlabel(ax, 'Time to contact (s)')
                ylabel(ax, '\Deltaspeed (a.u.)')
        end
        legend(ax, h, Location='northwest')
        if iExp <= length(implantSide)
            implantSideText = implantSide(iExp);
        elseif iExp == length(implantSide) + 1
            implantSideText = 'Both';
        elseif iExp == length(implantSide) + 2
            implantSideText = 'Left';
        elseif iExp == length(implantSide) + 3
            implantSideText = 'Right';
        end
        if ~strcmp(trialTypeName, 'diff')
            title(ax, sprintf('%s (%d trials, %s implant)', trialTypeName, n, implantSideText));
        else
            title(ax, sprintf('%s (%s implant)', trialTypeName, implantSideText));
        end
    end

    yl = vertcat(axAll.YLim);
    yl = [min(yl(:, 1)), max(yl(:, 2))];
    set(axAll, 'YLim', yl)
end
clear iTrialType iExp ax fig trialTypeName h iVar axAll n

%% 1.2. For CompleteExperiment data, do GLM (with ITI) and compare weights.
% 1.2.1 Show correlation matrix between predictors
% 1.2.2 Categorize SNr units by selectivity using GLM weights (orofacial vs.
% limb/trunk)
% 1.2.3 Confirm contralateral selectivity by comparing left vs. right limb
% weights

%% 2.1. For all animals (with lick tube servo), compare press trial ETA vs. lick trial ETA traces

%% 2.2. For all units with osci-lick spiking, compare press trial ETA vs. lick trial ETA traces. NULL: latter is a time shifted version of the former