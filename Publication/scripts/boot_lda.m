%%
load_ephysunits
read_reachDir_2tgt

%% LDA results
load('C:\SERVER\Units\lda_pressVsLick_20241212.mat')
load('C:\SERVER\Units\lda_reach2tgt_20241213.mat')


%% Extract peri-move responses for lick vs. reach
close all

clear pLDA;
pLDA.window = [-4, 2];
pLDA.baselineWindow = [-4, -2];
pLDA.responseWindow = [-0.3, 0];
pLDA.responseWindow2tgt = [-0.1, 0.2];
pLDA.minTrialDuration = 2;
pLDA.minNumUnits = 10;
pLDA.res = 0.1;
pLDA.nBoot = 10000;
pLDA.kFold = 5;

% Select sessions with enough press/lick trials, enough units
allUnitIndices = find(c.hasPress & c.hasLick);
[goodExpNames, ~, ic] = unique(string({eu(allUnitIndices).ExpName}));
nUnits = histcounts(ic, 1:max(ic)+1);
goodExpNames = goodExpNames(nUnits >= pLDA.minNumUnits);
nUnits = nUnits(nUnits >= pLDA.minNumUnits);

fprintf('Selected %i sessions with >%i units.\n', length(goodExpNames), pLDA.minNumUnits)
fprintf('\tUnits per session: %s\n', num2str(nUnits));

% Extract data
t = pLDA.window(1):pLDA.res:pLDA.window(2);
t = (t(1:end-1) + t(2:end))./2;
clear sr resp
sr(length(goodExpNames)) = struct(press=[], lick=[], t=[]);
resp(length(goodExpNames)) = struct(press=[], lick=[], baseline=[]);
for iExp = 1:length(goodExpNames)
    unitIndicesInExp = allUnitIndices(strcmpi(string({eu(allUnitIndices).ExpName}), goodExpNames{iExp}));
    unitIndicesInExp = unitIndicesInExp(:)';
    pressTrials = eu(unitIndicesInExp(1)).getTrials('press');
    pressTrials = pressTrials(pressTrials.duration() >= pLDA.minTrialDuration);
    lickTrials = eu(unitIndicesInExp(1)).getTrials('lick');
    lickTrials = lickTrials(lickTrials.duration() >= pLDA.minTrialDuration);
    sr(iExp).press = NaN(length(pressTrials), length(t), length(unitIndicesInExp));
    sr(iExp).lick = NaN(length(lickTrials), length(t), length(unitIndicesInExp));
    for i = 1:length(unitIndicesInExp)
        iEu = unitIndicesInExp(i);
        [sr(iExp).press(:, :, i), ~] = eu(iEu).getTrialAlignedData('count', pLDA.window, 'press', trials=pressTrials, alignTo='stop', resolution=pLDA.res, includeInvalid=false);
        [sr(iExp).lick(:, :, i), ~] = eu(iEu).getTrialAlignedData('count', pLDA.window, 'lick', trials=lickTrials, alignTo='stop', resolution=pLDA.res, includeInvalid=false);
    end
    sr(iExp).press = sr(iExp).press ./ pLDA.res;    
    sr(iExp).lick = sr(iExp).lick ./ pLDA.res;
    sr(iExp).t = t;

    selT = t>=pLDA.baselineWindow(1) & t<=pLDA.baselineWindow(2);
    
    mu = mean(sr(iExp).press(:, selT, :), [1, 2], 'omitnan');
    sd = std(sr(iExp).press(:, selT, :), 0, [1, 2], 'omitnan');
    sr(iExp).press = (sr(iExp).press - mu) ./ sd;

    mu = mean(sr(iExp).lick(:, selT, :), [1, 2], 'omitnan');
    sd = std(sr(iExp).lick(:, selT, :), 0, [1, 2], 'omitnan');
    sr(iExp).lick = (sr(iExp).lick - mu) ./ sd;


    resp(iExp).press = squeeze(mean(sr(iExp).press(:, t>=pLDA.responseWindow(1) & t<=pLDA.responseWindow(2), :), 2, 'omitnan'));
    resp(iExp).lick = squeeze(mean(sr(iExp).lick(:, t>=pLDA.responseWindow(1) & t<=pLDA.responseWindow(2), :), 2, 'omitnan'));
    resp(iExp).baseline = cat(1, ...
        squeeze(mean(sr(iExp).press(:, t>=pLDA.baselineWindow(1) & t<=pLDA.baselineWindow(2), :), 2, 'omitnan')), ...
        squeeze(mean(sr(iExp).lick(:, t>=pLDA.baselineWindow(1) & t<=pLDA.baselineWindow(2), :), 2, 'omitnan')) ...
        );

    % Remove NaNs
    selPress = all(~isnan(resp(iExp).press), 2);
    selLick = all(~isnan(resp(iExp).lick), 2);
    selBaseline = all(~isnan(resp(iExp).baseline), 2);

    resp(iExp).press = resp(iExp).press(selPress, :);
    resp(iExp).lick = resp(iExp).lick(selLick, :);
    resp(iExp).baseline = resp(iExp).baseline(selBaseline, :);

    sr(iExp).press = sr(iExp).press(selPress, :, :);
    sr(iExp).lick = sr(iExp).lick(selLick, :, :);
end

clear t iExp unitIndicesInExp pressTrials lickTrials i iEu

% Fit LDA
clear likelihood
likelihood(length(sr)) = struct(press=[], lick=[], baseline=[]);
t = sr(1).t;

for iExp = 1:length(sr)
    nPress = size(resp(iExp).press, 1);
    nLick = size(resp(iExp).lick, 1);
    nBaseline = size(resp(iExp).baseline, 1);
    nTrials = nPress + nLick;

    % Fit model using response window
    X = vertcat(resp(iExp).press, resp(iExp).lick, resp(iExp).baseline);
    Y = vertcat(repmat("press", [nPress, 1]), repmat("lick", [nLick, 1]), repmat("baseline", [nBaseline, 1]));
    mdl = fitcdiscr(X, Y, Prior='empirical', CrossVal='on', KFold=pLDA.kFold);

    % Predict full timecourse using fitted model
    likelihood(iExp).press = NaN(nTrials, length(t));
    likelihood(iExp).lick = NaN(nTrials, length(t));
    likelihood(iExp).baseline = NaN(nTrials, length(t));
    likelihood(iExp).trueLabel = Y;
    likelihood(iExp).t = t;

    assert(length(mdl.ClassNames) == 3)
    [~, iClass] = ismember(["press", "lick", "baseline"], mdl.ClassNames);
    for i = 1:length(t)
        XPress = squeeze(sr(iExp).press(:, i, :));
        XLick = squeeze(sr(iExp).lick(:, i, :));
        XAll = vertcat(XPress, XLick);

        score = NaN(nPress + nLick, 3);

        for iFold = 1:pLDA.kFold
            testIndices = mdl.Partition.test(iFold);
            testIndices = testIndices(1:nTrials)';
            [~, score(testIndices, 1:3)] = mdl.Trained{iFold}.predict(XAll(testIndices, :));
        end

        likelihood(iExp).press(:, i) = score(:, iClass(1));
        likelihood(iExp).lick(:, i) = score(:, iClass(2));
        likelihood(iExp).baseline(:, i) = score(:, iClass(3));
    end
    likelihood(iExp).df = likelihood(iExp).press - likelihood(iExp).lick; % df = press - lick
end
clear iExp nPress nLick nTrials X Y mdl i t XPress XLick score iFold testIndices trainIndices

% Plot results (individual sessions)
fig = figure(Units='inches', Position=[1 1 14 6]);
tl = tiledlayout(fig, 3, 5, TileSpacing='tight', Padding='tight', TileIndexing='rowmajor');

t = likelihood(1).t;

for iExp = 1:length(likelihood)
    isPress = likelihood(iExp).trueLabel == "press";
    isLick = likelihood(iExp).trueLabel == "lick";

    ax = nexttile(tl);

    hold(ax, 'on')
    h = gobjects(2, 1);
    h(1) = plot(ax, t, mean(likelihood(iExp).press(isPress, :) - likelihood(iExp).lick(isPress, :), 1, 'omitnan'), 'red', LineWidth=1.5, DisplayName=sprintf('press (%i trials)', nnz(isPress)));
    h(2) = plot(ax, t, mean(likelihood(iExp).press(isLick, :) - likelihood(iExp).lick(isLick, :), 1, 'omitnan'), 'blue', LineWidth=1.5, DisplayName=sprintf('lick (%i trials)', nnz(isLick)));
    hold(ax, 'off')
    title(ax, sprintf('%s (%i units)', goodExpNames(iExp), nUnits(iExp)), Interpreter='none')

    ylim(ax, [-1, 1])
    xticks(ax, -4:2:2)
    xline(ax, 0, 'k:')
    yline(ax, 0, 'k:')

    legend(ax, h, Location='northoutside', Orientation='horizontal')

    fontsize(ax, p.fontSize, 'points')
end
xlabel(tl, 'Time to contact (s)', FontSize=p.fontSize)
ylabel(tl, 'p(press) - p(lick)', FontSize=p.fontSize)

%% Extract peri-move responses for reachDir 2tgt
assert(all(strcmpi(targetNames, {'contra-out', 'contra-in'})))

[expNames2tgt, ~, euExpIndex] = unique(string({euReachDir2tgt.ExpName}));

% Exclude trials with ipsi paw movement
nUnits2tgt = arrayfun(@(traj) size(traj.eta(1).X, 1), traj2tgt);
goodExpIndices2tgt = find(nUnits2tgt >= pLDA.minNumUnits);
nUnits2tgt = nUnits2tgt(goodExpIndices2tgt);

% Select sessions with enough press/lick trials, enough units
fprintf('Selected %i sessions with >%i units.\n', length(goodExpIndices2tgt), pLDA.minNumUnits)
fprintf('\tUnits per session: %s\n', num2str(nUnits2tgt));

% Extract data
t = pLDA.window(1):pLDA.res:pLDA.window(2);
t = (t(1:end-1) + t(2:end))./2;

clear sr2tgt resp2tgt
sr2tgt(length(goodExpIndices2tgt)) = struct(contraOut=[], contraIn=[], t=[]);
resp2tgt(length(goodExpIndices2tgt)) = struct(contraOut=[], contraIn=[]);
for iExp = 1:length(goodExpIndices2tgt)
    iExpRaw = goodExpIndices2tgt(iExp);
    unitIndicesInExp = find(euExpIndex == iExpRaw);
    unitIndicesInExp = unitIndicesInExp(:)';

    contraOutTrials = traj2tgt(iExpRaw).trials{1}(~isnan(traj2tgt(iExpRaw).correction{1}));
    contraInTrials = traj2tgt(iExpRaw).trials{2}(~isnan(traj2tgt(iExpRaw).correction{2}));    

    sr2tgt(iExp).contraOut = NaN(length(contraOutTrials), length(t), length(unitIndicesInExp));
    sr2tgt(iExp).contraIn = NaN(length(contraInTrials), length(t), length(unitIndicesInExp));
    for i = 1:length(unitIndicesInExp)
        iEu = unitIndicesInExp(i);
        [sr2tgt(iExp).contraOut(:, :, i), ~] = eu(iEu).getTrialAlignedData('count', pLDA.window, 'press_spontaneous', trials=traj2tgt(iExpRaw).trials{1}, correction=traj2tgt(iExpRaw).correction{1}, alignTo='stop', resolution=pLDA.res, includeInvalid=true, allowedTrialDuration=[-Inf, Inf], correctionAdvancedValidation=false);
        [sr2tgt(iExp).contraIn(:, :, i), ~] = eu(iEu).getTrialAlignedData('count', pLDA.window, 'press_spontaneous', trials=traj2tgt(iExpRaw).trials{2}, correction=traj2tgt(iExpRaw).correction{2}, alignTo='stop', resolution=pLDA.res, includeInvalid=true, allowedTrialDuration=[-Inf, Inf], correctionAdvancedValidation=false);
    end
    sr2tgt(iExp).contraOut = sr2tgt(iExp).contraOut ./ pLDA.res;    
    sr2tgt(iExp).contraIn = sr2tgt(iExp).contraIn ./ pLDA.res;
    sr2tgt(iExp).t = t;

    selT = t>=pLDA.baselineWindow(1) & t<=pLDA.baselineWindow(2);

    mu = mean(sr2tgt(iExp).contraOut(:, selT, :), [1, 2], 'omitnan');
    sd = std(sr2tgt(iExp).contraOut(:, selT, :), 0, [1, 2], 'omitnan');
    sr2tgt(iExp).contraOut = (sr2tgt(iExp).contraOut - mu) ./ sd;

    mu = mean(sr2tgt(iExp).contraIn(:, selT, :), [1, 2], 'omitnan');
    sd = std(sr2tgt(iExp).contraIn(:, selT, :), 0, [1, 2], 'omitnan');
    sr2tgt(iExp).contraIn = (sr2tgt(iExp).contraIn - mu) ./ sd;


    resp2tgt(iExp).contraOut = squeeze(mean(sr2tgt(iExp).contraOut(:, t>=pLDA.responseWindow2tgt(1) & t<=pLDA.responseWindow2tgt(2), :), 2, 'omitnan'));
    resp2tgt(iExp).contraIn = squeeze(mean(sr2tgt(iExp).contraIn(:, t>=pLDA.responseWindow2tgt(1) & t<=pLDA.responseWindow2tgt(2), :), 2, 'omitnan'));
    resp2tgt(iExp).baseline = cat(1, ...
        squeeze(mean(sr2tgt(iExp).contraOut(:, t>=pLDA.baselineWindow(1) & t<=pLDA.baselineWindow(2), :), 2, 'omitnan')), ...
        squeeze(mean(sr2tgt(iExp).contraIn(:, t>=pLDA.baselineWindow(1) & t<=pLDA.baselineWindow(2), :), 2, 'omitnan')) ...
        );

    % Remove NaNs
    selContraOut = all(~isnan(resp2tgt(iExp).contraOut), 2);
    selContraIn = all(~isnan(resp2tgt(iExp).contraIn), 2);
    selBaseline = all(~isnan(resp2tgt(iExp).baseline), 2);

    resp2tgt(iExp).contraOut = resp2tgt(iExp).contraOut(selContraOut, :);
    resp2tgt(iExp).contraIn = resp2tgt(iExp).contraIn(selContraIn, :);
    resp2tgt(iExp).baseline = resp2tgt(iExp).baseline(selBaseline, :);

    sr2tgt(iExp).contraOut = sr2tgt(iExp).contraOut(selContraOut, :, :);
    sr2tgt(iExp).contraIn = sr2tgt(iExp).contraIn(selContraIn, :, :);
end

clear t iExp iExp unitIndicesInExp contraOutTrials contraInTrials i iEu selT mu sd

% Fit LDA
clear likelihood2tgt
likelihood2tgt(length(sr2tgt)) = struct(contraOut=[], contraIn=[], baseline=[]);
t = sr2tgt(1).t;

for iExp = 1:length(sr2tgt)
    nContraOut = size(resp2tgt(iExp).contraOut, 1);
    nContraIn = size(resp2tgt(iExp).contraIn, 1);
    nBaseline = size(resp2tgt(iExp).baseline, 1);
    nTrials = nContraOut + nContraIn;

    % Fit model using response window
    X = vertcat(resp2tgt(iExp).contraOut, resp2tgt(iExp).contraIn, resp2tgt(iExp).baseline);
    Y = vertcat(repmat("lateral", [nContraOut, 1]), repmat("medial", [nContraIn, 1]), repmat("baseline", [nBaseline, 1]));
    mdl = fitcdiscr(X, Y, Prior='empirical', CrossVal='on', KFold=pLDA.kFold);

    % Predict full timecourse using fitted model
    likelihood2tgt(iExp).contraOut = NaN(nTrials, length(t));
    likelihood2tgt(iExp).contraIn = NaN(nTrials, length(t));
    likelihood2tgt(iExp).baseline = NaN(nTrials, length(t));
    likelihood2tgt(iExp).trueLabel = Y;
    likelihood2tgt(iExp).t = t;

    assert(length(mdl.ClassNames) == 3)
    [~, iClass] = ismember(["lateral", "medial", "baseline"], mdl.ClassNames);
    for i = 1:length(t)
        XContraOut = squeeze(sr2tgt(iExp).contraOut(:, i, :));
        XContraIn = squeeze(sr2tgt(iExp).contraIn(:, i, :));
        XAll = vertcat(XContraOut, XContraIn);

        score = NaN(nContraOut + nContraIn, 3);

        for iFold = 1:pLDA.kFold
            testIndices = mdl.Partition.test(iFold);
            testIndices = testIndices(1:nTrials)';
            [~, score(testIndices, 1:3)] = mdl.Trained{iFold}.predict(XAll(testIndices, :));
        end

        likelihood2tgt(iExp).contraOut(:, i) = score(:, iClass(1));
        likelihood2tgt(iExp).contraIn(:, i) = score(:, iClass(2));
        likelihood2tgt(iExp).baseline(:, i) = score(:, iClass(3));
    end
    likelihood2tgt(iExp).df = likelihood2tgt(iExp).contraOut - likelihood2tgt(iExp).contraIn; % df = contraOut - contraIn
end
clear iExp nContraOut nContraIn nTrials X Y mdl i t XContraOut XContraIn score

% Plot results (individual sessions)
fig = figure(Units='inches', Position=[1 1 14 6]);
tl = tiledlayout(fig, 3, ceil(length(sr2tgt)./3), TileSpacing='tight', Padding='tight', TileIndexing='rowmajor');

t = likelihood2tgt(1).t;

for iExp = 1:length(likelihood2tgt)
    isContraOut = likelihood2tgt(iExp).trueLabel == "lateral";
    isContraIn = likelihood2tgt(iExp).trueLabel == "medial";

    ax = nexttile(tl);

    hold(ax, 'on')
    h = gobjects(2, 1);
    h(1) = plot(ax, t, mean(likelihood2tgt(iExp).contraOut(isContraOut, :) - likelihood2tgt(iExp).contraIn(isContraOut, :), 1, 'omitnan'), 'red', LineWidth=1.5, DisplayName=sprintf('lateral (%i trials)', nnz(isContraOut)));
    h(2) = plot(ax, t, mean(likelihood2tgt(iExp).contraOut(isContraIn, :) - likelihood2tgt(iExp).contraIn(isContraIn, :), 1, 'omitnan'), 'blue', LineWidth=1.5, DisplayName=sprintf('medial (%i trials)', nnz(isContraIn)));
    hold(ax, 'off')
    title(ax, sprintf('%s (%i units)', goodExpNames(iExp), nUnits(iExp)), Interpreter='none')

    ylim(ax, [-1, 1])
    xticks(ax, -4:2:2)
    xline(ax, 0, 'k:')
    yline(ax, 0, 'k:')

    legend(ax, h, Location='northoutside', Orientation='horizontal')

    fontsize(ax, p.fontSize, 'points')
end
xlabel(tl, 'Time to reach onset (s)', FontSize=p.fontSize)
ylabel(tl, 'p(lateral) - p(medial)', FontSize=p.fontSize)

clear fig tl iExp isContraOut isContraIn ax h

%% Bootstrap LDA (perm test) for lick vs reach

% Fit LDA
t = sr(1).t;

dfBoot = arrayfun(@(llh) llh.df, likelihood, UniformOutput=false);
dfBoot = cat(1, dfBoot{:});
dfBoot = NaN([size(dfBoot), pLDA.nBoot]);

rng(42); % Woah double-rainbow!

pool = parpool();
parfor iBoot = 1:pLDA.nBoot
    % fprintf('%i\n', iBoot);
    df = cell(length(sr), 1);
    for iExp = 1:length(sr)
        nPress = size(resp(iExp).press, 1);
        nLick = size(resp(iExp).lick, 1);
        nBaseline = size(resp(iExp).baseline, 1);
        nTrials = nPress + nLick;
    
        % Fit model using response window
        X = vertcat(resp(iExp).press, resp(iExp).lick);
        X = X(randperm(size(X, 1)), :);
        X = vertcat(X, resp(iExp).baseline);
        Y = vertcat(repmat("press", [nPress, 1]), repmat("lick", [nLick, 1]), repmat("baseline", [nBaseline, 1]));
        mdl = fitcdiscr(X, Y, Prior='empirical', CrossVal='on', KFold=pLDA.kFold);
    
        % Predict full timecourse using fitted model
        press = NaN(nTrials, length(t));
        lick = NaN(nTrials, length(t));
        assert(length(mdl.ClassNames) == 3)
        [~, iClass] = ismember(["press", "lick", "baseline"], mdl.ClassNames);
        for i = 1:length(t)
            XPress = squeeze(sr(iExp).press(:, i, :));
            XLick = squeeze(sr(iExp).lick(:, i, :));
            XAll = vertcat(XPress, XLick);

            score = NaN(nPress + nLick, 3);
    
            for iFold = 1:pLDA.kFold
                testIndices = mdl.Partition.test(iFold);
                testIndices = testIndices(1:nTrials)';
                [~, score(testIndices, 1:3)] = mdl.Trained{iFold}.predict(XAll(testIndices, :));
            end

            press(:, i) = score(:, iClass(1));
            lick(:, i) = score(:, iClass(2));
        end
        df{iExp} = press - lick; % df = press - lick
    end
    dfBoot(:, :, iBoot) = cat(1, df{:});
end

delete(pool)
clear iExp nPress nLick nTrials X Y mdl i press lick XPress XLick score likelihoodBoot df pool iBoot

fprintf('\nDone.\n')

% Save to disk this is a few GBs, but needed to calculate the CI
save('C:\SERVER\Units\lda_pressVsLick_fullBootData_20241216.mat', 'dfBoot', '-v7.3')

% Quick summary (99% CI, mean) of bootstrap for lick vs reach
clear dfBootStats;

Y = cell(length(sr), 1);
for iExp = 1:length(sr)
    nPress = size(resp(iExp).press, 1);
    nLick = size(resp(iExp).lick, 1);
    Y{iExp} = vertcat(repmat("press", [nPress, 1]), repmat("lick", [nLick, 1]));
end
Y = cat(1, Y{:});

isPress = Y == "press";
isLick = Y == "lick";

dfBootStats.press.X = transpose(squeeze(mean(dfBoot(isPress, :, :), 1, 'omitnan')));
dfBootStats.press.mu = mean(dfBootStats.press.X, 1, 'omitnan');
dfBootStats.press.ci = quantile(dfBootStats.press.X, [0.01, 0.99], 1);

dfBootStats.lick.X = transpose(squeeze(mean(dfBoot(isLick, :, :), 1, 'omitnan')));
dfBootStats.lick.mu = mean(dfBootStats.lick.X, 1, 'omitnan');
dfBootStats.lick.ci = quantile(dfBootStats.lick.X, [0.01, 0.99], 1);

dfBootStats.all.X = transpose(squeeze(mean(dfBoot, 1, 'omitnan')));
dfBootStats.all.mu = mean(dfBootStats.all.X, 1, 'omitnan');
dfBootStats.all.ci = quantile(dfBootStats.all.X, [0.01, 0.99], 1);

save('C:\SERVER\Units\lda_pressVsLick_20241216.mat', 'pLDA', 'likelihood', 'sr', 'resp', 't', 'goodExpNames', 'nUnits', 'dfBootStats')


% Bootstrap LDA (perm test) for reach 2tgt

% Fit LDA
t = sr2tgt(1).t;

dfBoot2tgt = arrayfun(@(llh) llh.df, likelihood2tgt, UniformOutput=false);
dfBoot2tgt = cat(1, dfBoot2tgt{:});
dfBoot2tgt = NaN([size(dfBoot2tgt), pLDA.nBoot]);

rng(42); % Woah double-rainbow!

pool = parpool();
parfor iBoot = 1:pLDA.nBoot
    df = cell(length(sr2tgt), 1);
    for iExp = 1:length(sr2tgt)
        nContraOut = size(resp2tgt(iExp).contraOut, 1);
        nContraIn = size(resp2tgt(iExp).contraIn, 1);
        nBaseline = size(resp2tgt(iExp).baseline, 1);
        nTrials = nContraOut + nContraIn;
    
        % Fit model using response window
        X = vertcat(resp2tgt(iExp).contraOut, resp2tgt(iExp).contraIn);
        X = X(randperm(size(X, 1)), :);
        X = vertcat(X, resp2tgt(iExp).baseline);
        Y = vertcat(repmat("lateral", [nContraOut, 1]), repmat("medial", [nContraIn, 1]), repmat("baseline", [nBaseline, 1]));
        mdl = fitcdiscr(X, Y, Prior='empirical', CrossVal='on', KFold=pLDA.kFold);
    
        % Predict full timecourse using fitted model
        contraOut = NaN(nTrials, length(t));
        contraIn = NaN(nTrials, length(t));
        assert(length(mdl.ClassNames) == 3)
        [~, iClass] = ismember(["lateral", "medial", "baseline"], mdl.ClassNames);
        for i = 1:length(t)
            XContraOut = squeeze(sr2tgt(iExp).contraOut(:, i, :));
            XContraIn = squeeze(sr2tgt(iExp).contraIn(:, i, :));
            XAll = vertcat(XContraOut, XContraIn);

            score = NaN(nContraOut + nContraIn, 3);

            for iFold = 1:pLDA.kFold
                testIndices = mdl.Partition.test(iFold);
                testIndices = testIndices(1:nTrials)';
                [~, score(testIndices, 1:3)] = mdl.Trained{iFold}.predict(XAll(testIndices, :));
            end

            contraOut(:, i) = score(:, iClass(1));
            contraIn(:, i) = score(:, iClass(2));
        end
        df{iExp} = contraOut - contraIn; % df = lateral - medial
    end
    dfBoot2tgt(:, :, iBoot) = cat(1, df{:});
end

delete(pool)
clear iExp nContraOut nContraIn nTrials X Y mdl i contraOut contraIn XContraOut XContraIn score df pool iBoot

fprintf('\nDone.\n')

% Save to disk this is a few GBs, but needed to calculate the CI
save('C:\SERVER\Units\lda_reach2tgt_fullBootData_20241216.mat', 'dfBoot2tgt', '-v7.3')

% Quick summary (99% CI, mean) of bootstrap for reach 2tgt
clear dfBootStats2tgt;

Y = cell(length(sr2tgt), 1);
for iExp = 1:length(sr2tgt)
    nContraOut = size(resp2tgt(iExp).contraOut, 1);
    nContraIn = size(resp2tgt(iExp).contraIn, 1);
    Y{iExp} = vertcat(repmat("lateral", [nContraOut, 1]), repmat("medial", [nContraIn, 1]));
end
Y = cat(1, Y{:});

isContraOut = Y == "lateral";
isContraIn = Y == "medial";

dfBootStats2tgt.contraOut.X = transpose(squeeze(mean(dfBoot2tgt(isContraOut, :, :), 1, 'omitnan')));
dfBootStats2tgt.contraOut.mu = mean(dfBootStats2tgt.contraOut.X, 1, 'omitnan');
dfBootStats2tgt.contraOut.ci = quantile(dfBootStats2tgt.contraOut.X, [0.01, 0.99], 1);

dfBootStats2tgt.contraIn.X = transpose(squeeze(mean(dfBoot2tgt(isContraIn, :, :), 1, 'omitnan')));
dfBootStats2tgt.contraIn.mu = mean(dfBootStats2tgt.contraIn.X, 1, 'omitnan');
dfBootStats2tgt.contraIn.ci = quantile(dfBootStats2tgt.contraIn.X, [0.01, 0.99], 1);

dfBootStats2tgt.all.X = transpose(squeeze(mean(dfBoot2tgt, 1, 'omitnan')));
dfBootStats2tgt.all.mu = mean(dfBootStats2tgt.all.X, 1, 'omitnan');
dfBootStats2tgt.all.ci = quantile(dfBootStats2tgt.all.X, [0.01, 0.99], 1);

save('C:\SERVER\Units\lda_reach2tgt_20241216.mat', 'pLDA', 'likelihood2tgt', 'sr2tgt', 'resp2tgt', 't', 'expNames2tgt', 'euExpIndex', 'goodExpIndices2tgt', 'nUnits2tgt', 'dfBootStats2tgt')

clear Y iExp nContraOut nContraIn Y isContraOut isContraIn


%% Plots
close all

fig = figure(Units='inches', Position=[1 1 2.25 4]);
tl = tiledlayout(fig, 2, 1);


% Plot results reach2tgt (average across sessions, same plot for both reach and lick, bootCI is for all trials)

t = likelihood2tgt(1).t;
latTrialLat = arrayfun(@(llh) llh.contraOut(llh.trueLabel=="lateral", :), likelihood2tgt, UniformOutput=false);
latTrialMed = arrayfun(@(llh) llh.contraIn(llh.trueLabel=="lateral", :), likelihood2tgt, UniformOutput=false);
medTrialLat = arrayfun(@(llh) llh.contraOut(llh.trueLabel=="medial", :), likelihood2tgt, UniformOutput=false);
medTrialMed = arrayfun(@(llh) llh.contraIn(llh.trueLabel=="medial", :), likelihood2tgt, UniformOutput=false);
latTrialLat = cat(1, latTrialLat{:});
latTrialMed = cat(1, latTrialMed{:});
medTrialLat = cat(1, medTrialLat{:});
medTrialMed = cat(1, medTrialMed{:});

DATA = { ...
    latTrialLat, latTrialMed; ...
    medTrialLat, medTrialMed ...
    };
COLOR = {getColor(1, 4, 0.8), getColor(3, 4, 0.8)};
LABEL = ["lateral", "medial"];

ax = nexttile(tl);
hold(ax, 'on')
h = gobjects(3, 1);
for iMove = 1:2
    mu = mean(DATA{iMove, 1} - DATA{iMove, 2}, 1, 'omitnan');
    h(iMove) = plot(ax, t, mu, Color=COLOR{iMove}, LineWidth=1.5, DisplayName=sprintf('%s', LABEL(iMove)));
end
h(3) = plot(ax, t, dfBootStats2tgt.all.mu, 'black', LineStyle='--', LineWidth=1.5, DisplayName='shuffle');
patch(ax, [t, flip(t)], [dfBootStats2tgt.all.ci(1, :), flip(dfBootStats2tgt.all.ci(2, :))], 'black', FaceAlpha=0.1, EdgeAlpha=0.5, DisplayName='99% CI');
patch(ax, [pLDA.responseWindow2tgt, flip(pLDA.responseWindow2tgt)], [-1, -1, 1, 1], 'yellow', FaceAlpha=0.1, EdgeAlpha=0, DisplayName='training');
ylim(ax, [-1, 1])
xline(ax, 0, ':')
yline(ax, 0, ':')
hold(ax, 'off')

xlabel(ax, 'Time to reach onset (s)')
ylabel(ax, 'p(lateral) - p(medial)')

fontsize(ax, p.fontSize, 'points')

lgd = legend(ax, h(1:2), Orientation='horizontal', Location='northoutside', FontSize=p.fontSize-1, AutoUpdate=false);

selUnits2tgt = ismember(string({euReachDir2tgt.ExpName}), expNames2tgt(goodExpIndices2tgt));
nAnimals2tgt = length(unique(euReachDir2tgt(selUnits2tgt).getAnimalName()));

fprintf('\nReach (2tgt) lateral vs. medial decoder:\n');
fprintf('\t1) Selected %i sessions (%i animals) with >=%i units (mean=%g, total=%i).\n', length(goodExpIndices2tgt), nAnimals2tgt, pLDA.minNumUnits, mean(nUnits2tgt), sum(nUnits2tgt))
fprintf('\t\t Units per session (sorted): %s\n', num2str(sort(nUnits2tgt)));
fprintf('\t2) %i lateral trials, %i medial trials on aggregate.\n', size(latTrialLat, 1), size(medTrialLat, 1))
fprintf('\t3) Shuffled trial labels %i times to retrain LDA, shuffled mean and 99%%CI is calculated by averaging across all %i trials for each shuffle, then averaging across all shuffles.\n', pLDA.nBoot, size(latTrialLat, 1)+size(medTrialLat, 1))
fprintf('\t4) Training window was chosen at [%g, %g] s.\n', pLDA.responseWindow2tgt(1), pLDA.responseWindow2tgt(2))

% Plot results lick vs reach (average across sessions, same plot for both reach and lick, bootCI is for all trials)
t = likelihood(1).t;
pressTrialPress = arrayfun(@(llh) llh.press(llh.trueLabel=="press", :), likelihood, UniformOutput=false); % Press likelihood for true-press trials
pressTrialLick = arrayfun(@(llh) llh.lick(llh.trueLabel=="press", :), likelihood, UniformOutput=false); % Lick likelihood for true-press trials
lickTrialPress = arrayfun(@(llh) llh.press(llh.trueLabel=="lick", :), likelihood, UniformOutput=false);
lickTrialLick = arrayfun(@(llh) llh.lick(llh.trueLabel=="lick", :), likelihood, UniformOutput=false);
pressTrialPress = cat(1, pressTrialPress{:});
pressTrialLick = cat(1, pressTrialLick{:});
lickTrialPress = cat(1, lickTrialPress{:});
lickTrialLick = cat(1, lickTrialLick{:});   

DATA = { ...
    pressTrialPress, pressTrialLick; ...
    lickTrialPress, lickTrialLick ...
    };
COLOR = ["red", "blue"];
LABEL = ["reach", "lick"];

ax = nexttile(tl);
hold(ax, 'on')
h = gobjects(3, 1);
for iMove = 1:2
    mu = mean(DATA{iMove, 1} - DATA{iMove, 2}, 1, 'omitnan');
    h(iMove) = plot(ax, t, mu, COLOR(iMove), LineWidth=1.5, DisplayName=sprintf('%s', LABEL(iMove)));
end
h(3) = plot(ax, t, dfBootStats.all.mu, 'black', LineStyle='--', LineWidth=1.5, DisplayName='shuffle');
patch(ax, [t, flip(t)], [dfBootStats.all.ci(1, :), flip(dfBootStats.all.ci(2, :))], 'black', FaceAlpha=0.1, EdgeAlpha=0.5, DisplayName='99% CI');
patch(ax, [pLDA.responseWindow, flip(pLDA.responseWindow)], [-1, -1, 1, 1], 'yellow', FaceAlpha=0.1, EdgeAlpha=0, DisplayName='training');
ylim(ax, [-1, 1])
xline(ax, 0, ':')
yline(ax, 0, ':')
hold(ax, 'off')

xlabel(ax, 'Time to bar/spout contact (s)')
ylabel(ax, 'p(reach) - p(lick)')

fontsize(ax, p.fontSize, 'points')
lgd = legend(ax, h(1:2), Orientation='horizontal', Location='northoutside', FontSize=p.fontSize-1, AutoUpdate=false);

nAnimals = length(unique(eu(ismember(string({eu.ExpName}), goodExpNames)).getAnimalName()));

fprintf('\nReach vs. lick decoder:\n');
fprintf('\t1) Selected %i sessions (%i animals) with >=%i units (mean=%g, total=%i).\n', length(goodExpNames), nAnimals, pLDA.minNumUnits, mean(nUnits), sum(nUnits))
fprintf('\t\t Units per session (sorted): %s\n', num2str(sort(nUnits)));
fprintf('\t2) %i reach trials, %i lick trials on aggregate.\n', size(pressTrialPress, 1), size(lickTrialPress, 1))
fprintf('\t3) Shuffled trial labels %i times to retrain LDA, shuffled mean and 99%%CI is calculated by averaging across all %i trials for each shuffle, then averaging across all shuffles.\n', pLDA.nBoot, size(pressTrialPress, 1)+size(lickTrialPress, 1))
fprintf('\t4) Training window was chosen at [%g, %g] s.\n', pLDA.responseWindow(1), pLDA.responseWindow(2))


copygraphics(fig, ContentType='vector', BackgroundColor='none')