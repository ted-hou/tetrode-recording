%%
load_ephysunits

%%
close all

clear pLDA;
pLDA.window = [-4, 2];
pLDA.baselineWindow = [-4, -2];
pLDA.responseWindow = [-0.3, 0];
pLDA.minTrialDuration = 2;
pLDA.minNumUnits = 15;
pLDA.res = 0.1;
pLDA.nBoot = 10000;

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
resp(length(goodExpNames)) = struct(press=[], lick=[]);
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
end

clear t iExp unitIndicesInExp pressTrials lickTrials i iEu

% Fit LDA
clear likelihood
likelihood(length(sr)) = struct(press=[], lick=[]);
t = sr(1).t;

for iExp = 1:length(sr)
    nPress = size(resp(iExp).press, 1);
    nLick = size(resp(iExp).lick, 1);
    nTrials = nPress + nLick;

    % Fit model using response window
    X = vertcat(resp(iExp).press, resp(iExp).lick);
    Y = vertcat(repmat("press", [nPress, 1]), repmat("lick", [nLick, 1]));
    mdl = fitcdiscr(X, Y, Prior='uniform');

    % Predict full timecourse using fitted model
    likelihood(iExp).press = NaN(nTrials, length(t));
    likelihood(iExp).lick = NaN(nTrials, length(t));
    likelihood(iExp).trueLabel = Y;
    likelihood(iExp).t = t;
    for i = 1:length(t)
        XPress = squeeze(sr(iExp).press(:, i, :));
        XLick = squeeze(sr(iExp).lick(:, i, :));
        [~, score] = mdl.predict(vertcat(XPress, XLick));
        likelihood(iExp).press(:, i) = score(:, strcmpi('press', mdl.ClassNames));
        likelihood(iExp).lick(:, i) = score(:, strcmpi('lick', mdl.ClassNames));
    end
    likelihood(iExp).df = likelihood(iExp).press - likelihood(iExp).lick; % df = press - lick
end
clear iExp nPress nLick nTrials X Y mdl i t XPress XLick score

%% Bootstrap LDA (perm test)

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
        nTrials = nPress + nLick;
    
        % Fit model using response window
        X = vertcat(resp(iExp).press, resp(iExp).lick);
        X = X(randperm(size(X, 1)), :);
        Y = vertcat(repmat("press", [nPress, 1]), repmat("lick", [nLick, 1]));
        mdl = fitcdiscr(X, Y, Prior='uniform');
    
        % Predict full timecourse using fitted model
        press = NaN(nTrials, length(t));
        lick = NaN(nTrials, length(t));
        for i = 1:length(t)
            XPress = squeeze(sr(iExp).press(:, i, :));
            XLick = squeeze(sr(iExp).lick(:, i, :));
            [~, score] = mdl.predict(vertcat(XPress, XLick));
            press(:, i) = score(:, strcmpi('press', mdl.ClassNames));
            lick(:, i) = score(:, strcmpi('lick', mdl.ClassNames));
        end
        df{iExp} = press - lick; % df = press - lick
    end
    dfBoot(:, :, iBoot) = cat(1, df{:});
end

clear iExp nPress nLick nTrials X Y mdl i press lick XPress XLick score likelihoodBoot df

fprintf('\nDone.\n')

delete(pool)

save('E:\DATA\Units\lda_pressVsLick_20241212.mat', 'pLDA', 'dfBoot', 'likelihood', 'sr', 'resp', 't', 'goodExpNames', 'nUnits', '-v7.3')

clear iBoot t likelihoodBoot

save('E:\DATA\Units\lda_pressVsLick_fullBootData_20241212.mat', 'dfBoot')

%%
clear dfBootStats;

Y = cell(iExp, 1);
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

save('E:\DATA\Units\lda_pressVsLick_20241212.mat', 'pLDA', 'likelihood', 'sr', 'resp', 't', 'goodExpNames', 'nUnits', 'dfBootStats', '-v7.3')

%% Plot results (individual sessions)
fig = figure(Units='normalized', Position=[0.1 0.1 0.8 0.8]);
tlp = tiledlayout(fig, 3, 5, TileSpacing='tight', Padding='tight', TileIndexing='rowmajor');

t = likelihood(1).t;

for iExp = 1:length(likelihood)
    isPress = likelihood(iExp).trueLabel == "press";
    isLick = likelihood(iExp).trueLabel == "lick";

    tl = tiledlayout(tlp, 2, 1, TileSpacing='tight', Padding='tight');
    tl.Layout.Tile = iExp;
    ax = gobjects(2, 1);

    ax(1) = nexttile(tl);
    hold(ax(1), 'on')
    plot(ax(1), t, mean(likelihood(iExp).press(isPress, :), 1, 'omitnan'), 'red', LineWidth=1.5, DisplayName='press')
    plot(ax(1), t, mean(likelihood(iExp).lick(isPress, :), 1, 'omitnan'), 'blue', LineWidth=1.5, DisplayName='lick')
    hold(ax(1), 'off')
    title(ax(1), sprintf('press trials (n=%i)', nnz(isPress)), Color='red')

    ax(2) = nexttile(tl);
    hold(ax(2), 'on')
    plot(ax(2), t, mean(likelihood(iExp).press(isLick, :), 1, 'omitnan'), 'red', LineWidth=1.5, DisplayName='press')
    plot(ax(2), t, mean(likelihood(iExp).lick(isLick, :), 1, 'omitnan'), 'blue', LineWidth=1.5, DisplayName='lick')
    hold(ax(2), 'off')
    title(ax(2), sprintf('lick trials (n=%i)', nnz(isLick)), Color='blue')

    ylim(ax, [0, 1])
    xticks(ax, -4:2:2)

    title(tl, sprintf('%s (%i units)', goodExpNames(iExp), nUnits(iExp)), Interpreter='none')
    legend(ax(1))
    legend(ax(2))
end
xlabel(tlp, 'Time to contact (s)')
ylabel(tlp, 'Probability')


%% Plot results (average across sessions)
fig = figure(Units='inches', Position=[1 1 3 3]);
tl = tiledlayout(fig, 2, 1);

t = likelihood(1).t;
pressTrialPress = arrayfun(@(llh) llh.press(llh.trueLabel=="press", :), likelihood, UniformOutput=false);
pressTrialLick = arrayfun(@(llh) llh.lick(llh.trueLabel=="press", :), likelihood, UniformOutput=false);
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
BOOT = {dfBootStats.press, dfBootStats.lick};
COLOR = ["red", "blue"];
LABEL = ["reach", "lick"];

ax = gobjects(2, 1);

h = gobjects(2, 3);
for iAx = 1:2
    ax(iAx) = nexttile(tl);
    hold(ax(iAx), 'on')
    h = gobjects(2, 1);
    mu = mean(DATA{iAx, 1} - DATA{iAx, 2}, 1, 'omitnan');
    err = std(DATA{iAx, 1} - DATA{iAx, 2}, 0, 1, 'omitnan');
    h(iAx, 1) = plot(ax(iAx), t, mu, 'black', LineWidth=1.5, DisplayName='observed');
    h(iAx, 2) = plot(ax(iAx), t, BOOT{iAx}.mu, 'black', LineStyle='--', LineWidth=1.5, DisplayName='shuffle');
    h(iAx, 3) = patch(ax(iAx), [t, flip(t)], [BOOT{iAx}.ci(1, :), flip(BOOT{iAx}.ci(2, :))], 'black', FaceAlpha=0.1, DisplayName='99% CI');
    hold(ax(iAx), 'off')
    title(ax(iAx), sprintf('%s (%i trials)', LABEL(iAx), size(DATA{iAx, 1}, 1)), Color=COLOR(iAx))
    fontsize(ax(iAx), p.fontSize, 'points')
end

lgd = legend(ax(1), Orientation='horizontal', FontSize=p.fontSize);
lgd.Layout.Tile = 'north';

ylim(ax, [-1, 1])

xlabel(tl, 'Time to contact (s)', fontSize=p.fontSize)
ylabel(tl, 'p(reach) - p(lick)', fontSize=p.fontSize)
title(tl, sprintf('%i sessions, %i units', length(nUnits), sum(nUnits)), FontWeight='bold', fontSize=p.fontSize)




%% Plot results (average across sessions, same plot for both reach and lick, bootCI is for all trials)
close all

fig = figure(Units='inches', Position=[1 1 3.5 2]);
tl = tiledlayout(fig, 1, 1);

t = likelihood(1).t;
pressTrialPress = arrayfun(@(llh) llh.press(llh.trueLabel=="press", :), likelihood, UniformOutput=false);
pressTrialLick = arrayfun(@(llh) llh.lick(llh.trueLabel=="press", :), likelihood, UniformOutput=false);
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

h = gobjects(2, 3);
ax = nexttile(tl);
hold(ax, 'on')
h = gobjects(4, 1);
for iMove = 1:2
    mu = mean(DATA{iMove, 1} - DATA{iMove, 2}, 1, 'omitnan');
    h(iMove) = plot(ax, t, mu, COLOR(iMove), LineWidth=1.5, DisplayName=sprintf('%s', LABEL(iMove)));
end
h(3) = plot(ax, t, dfBootStats.all.mu, 'black', LineStyle='--', LineWidth=1.5, DisplayName='shuffle');
h(4) = patch(ax, [t, flip(t)], [dfBootStats.all.ci(1, :), flip(dfBootStats.all.ci(2, :))], 'black', FaceAlpha=0.1, EdgeAlpha=0.5, DisplayName='99% CI');
h(5) = patch(ax, [pLDA.responseWindow, flip(pLDA.responseWindow)], [-1, -1, 1, 1], 'yellow', FaceAlpha=0.1, EdgeAlpha=0.5, DisplayName='training');
fontsize(ax, p.fontSize, 'points')

lgd = legend(ax, Orientation='vertical', Location='eastoutside', FontSize=p.fontSize, AutoUpdate=false);

ylim(ax, [-1, 1])
xline(ax, 0, ':')
yline(ax, 0, ':')

hold(ax, 'off')

xlabel(tl, 'Time to bar/spout contact (s)', fontSize=p.fontSize)
ylabel(tl, 'p(reach) - p(lick)', fontSize=p.fontSize)
% title(tl, sprintf('%i sessions, %i units', length(nUnits), sum(nUnits)), FontWeight='bold', fontSize=p.fontSize)

nAnimals = length(unique(eu(ismember(string({eu.ExpName}), goodExpNames)).getAnimalName()));

fprintf('\nReach vs. lick decoder:\n');
fprintf('\t1) Selected %i sessions (%i animals) with >=%i units (mean=%g, total=%i).\n', length(goodExpNames), nAnimals, pLDA.minNumUnits, mean(nUnits), sum(nUnits))
fprintf('\t\t Units per session (sorted): %s\n', num2str(sort(nUnits)));
fprintf('\t2) %i reach trials, %i lick trials on aggregate.\n', size(pressTrialPress, 1), size(lickTrialPress, 1))
fprintf('\t3) Shuffled trial labels %i times to retrain LDA, shuffled mean and 99%%CI is calculated by averaging across all %i trials for each shuffle, then averaging across all shuffles.\n', pLDA.nBoot, size(pressTrialPress, 1)+size(lickTrialPress, 1))
fprintf('\t4) Training window was chosen at [%g, %g] s.\n', pLDA.responseWindow(1), pLDA.responseWindow(2))


copygraphics(fig, ContentType='vector', BackgroundColor='none')