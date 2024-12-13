%%
load_ephysunits
read_reachDir_2tgt

%% LDA results
% load('E:\Data\Units\lda_pressVsLick_20241212.mat')
% load('E:\Data\Units\lda_pressVsLick_fullBootData_20241212.mat')


%% Extract peri-move responses for lick vs. reach
close all

clear pLDA;
pLDA.window = [-4, 2];
pLDA.baselineWindow = [-4, -2];
pLDA.responseWindow = [-0.3, 0];
pLDA.minTrialDuration = 2;
pLDA.minNumUnits = 10;
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
end

clear t iExp unitIndicesInExp pressTrials lickTrials i iEu

% Fit LDA
clear likelihood
likelihood(length(sr)) = struct(press=[], lick=[], baseline=[]);
t = sr(1).t;

for iExp = 1:length(sr)
    nPress = size(resp(iExp).press, 1);
    nLick = size(resp(iExp).lick, 1);
    nTrials = nPress + nLick;

    % Fit model using response window
    X = vertcat(resp(iExp).press, resp(iExp).lick, resp(iExp).baseline);
    Y = vertcat(repmat("press", [nPress, 1]), repmat("lick", [nLick, 1]), repmat("baseline", [nTrials, 1]));
    mdl = fitcdiscr(X, Y, Prior='empirical'); % try fitclinear (Eden says need to balance #n trials)

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
        likelihood(iExp).baseline(:, i) = score(:, strcmpi('baseline', mdl.ClassNames));
    end
    likelihood(iExp).df = likelihood(iExp).press - likelihood(iExp).lick; % df = press - lick
end
clear iExp nPress nLick nTrials X Y mdl i t XPress XLick score

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
        X = vertcat(X, resp(iExp).baseline);
        Y = vertcat(repmat("press", [nPress, 1]), repmat("lick", [nLick, 1]), repmat("baseline", [nTrials, 1]));
        mdl = fitcdiscr(X, Y, Prior='empirical');
    
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

delete(pool)
clear iExp nPress nLick nTrials X Y mdl i press lick XPress XLick score likelihoodBoot df pool iBoot

fprintf('\nDone.\n')

save('E:\DATA\Units\lda_pressVsLick_fullBootData_20241212.mat', 'dfBoot', '-v7.3')

%%
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

save('E:\DATA\Units\lda_pressVsLick_20241212.mat', 'pLDA', 'likelihood', 'sr', 'resp', 't', 'goodExpNames', 'nUnits', 'dfBootStats')



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
h(5) = patch(ax, [pLDA.responseWindow, flip(pLDA.responseWindow)], [-1, -1, 1, 1], 'yellow', FaceAlpha=0.1, EdgeAlpha=0, DisplayName='training');
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