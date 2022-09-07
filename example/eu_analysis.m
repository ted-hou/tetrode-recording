load('C:\SERVER\PETH_All_aggregate_20220201.mat')

%% Make eu objects (SLOW, takes ~60min)
eu = EphysUnit(PETHPress, 'cullITI', true, 'extendedWindow', [-1, 2], 'readWaveforms', true);

%% Alternatively, load eu objects from disk (SLOW)
eu = EphysUnit.load('C:\SERVER\Units');

%% Calculate eta/peth
[eta, t, ~] = eu.getPETH('count', 'press', [-4, 0], minTrialDuration=2, normalize=[-4, -2]);
eta = transpose(mean(eta(:, t >= -0.5), 2, 'omitnan'));

%% Calculate median spike rate vs. press response magnitude
stats = [eu.SpikeRateStats];
msr = [stats.medianITI]; % Median ITI spike rate
clear stats

% DA
isDA = msr < 15;

% SNr
isSNr = msr >= 15;

% SNr with press trials
isValid = isSNr & ~isnan(eta);

% Plot distribution of median spike rates during ITI
f = figure('Units', 'normalized', 'OuterPosition', [0, 0.33, 1, 0.5], 'DefaultAxesFontSize', 14);
ax = subplot(1, 3, 1);
histogram(ax, msr(isValid))
xlabel('Median spike rate (sp/s)')

% Plot distribution of response magnitude for SNr cells
ax = subplot(1, 3, 2);
histogram(ax, eta(isValid))
xlabel('Mean response magnitude (modified z-score)')

% Plot median spike rate vs response magnitude for SNr cells
ax = subplot(1, 3, 3);
hold(ax, 'on')
h1 = scatter(ax, msr(isValid), eta(isValid), '.', 'DisplayName', sprintf('%i units', nnz(isValid)));
h2 = plot(ax, repmat(median(msr(isValid)), [1, 2]), [min(eta(isValid)), max(eta(isValid))], 'k--', 'DisplayName', sprintf('Median spike rate = %.1f sp/s', median(msr(isValid))));
h3 = plot(ax, [min(msr(isValid)), max(msr(isValid))], [0, 0], 'k');
legend(ax, [h1, h2])
xlabel('Median spike rate (sp/s)')
ylabel('Mean response magnitude (modified z-score)')

clear ax f h1 h2 h3

%% H1: PressUp cells have a higher baseline firing rate than PressDown cells.
sigmaThreshold = 0.3;
isPressUp = eta >= sigmaThreshold & isValid;
isPressDown = eta <= -sigmaThreshold & isValid;

f = figure('Units', 'normalized', 'OuterPosition', [0, 0.33, 1, 0.5], 'DefaultAxesFontSize', 14);
ax = gobjects(1, 3);
ax(1) = subplot(1, 3, 1);
histogram(ax(1), msr(isPressUp), 'Normalization', 'probability');
title(ax(1), sprintf('Median spike rates, %g press-activated units', nnz(isPressUp)))

ax(2) = subplot(1, 3, 2);
histogram(ax(2), msr(isPressDown), 'Normalization', 'probability');
title(ax(2), sprintf('Median spike rates, %g press-inhibited units', nnz(isPressDown)))

ax(3) = subplot(1, 3, 3);
nBins = 10;
[NUp, edgesUp] = histcounts(msr(isPressUp), nBins, 'Normalization', 'probability');
[NDown, edgesDown] = histcounts(msr(isPressDown), nBins, 'Normalization', 'probability');
centersUp = 0.5*(edgesUp(2:end) + edgesUp(1:end-1));
centersDown = 0.5*(edgesDown(2:end) + edgesDown(1:end-1));
hold(ax(3), 'on')
plot(ax(3), centersUp, NUp, 'LineWidth', 2, 'Color', 'red', 'DisplayName', sprintf('Response>=%g (N=%g)', sigmaThreshold, nnz(isPressUp)));
plot(ax(3), centersDown, NDown, 'LineWidth', 2, 'Color', 'blue', 'DisplayName', sprintf('Response<=-%g (N=%g)', sigmaThreshold, nnz(isPressDown)))
hold(ax(3), 'off')
legend(ax(3), 'Location', 'northeast');

xlabel(ax, 'Median ITI spike rate (sp/s)')
ylabel(ax, 'Frequency')

nboot = 100;
ciUp = bootci(nboot, @median, msr(isPressUp))
ciDown = bootci(nboot, @median, msr(isPressDown))
ciAll = bootci(nboot, @median, msr)

clear sigmaThreshold isPressUp isPressDown f ax nBins NUp NDown edgesUp edgesDown centersUp centersDown


%% Separate press-units by direction
sigmaThreshold = 0.3;
etaValid = eta(isValid);
euPress = eu(isValid);
isPressUp = etaValid >= sigmaThreshold;
isPressDown = etaValid <= -sigmaThreshold;
isResponsive = isPressUp | isPressDown;
fprintf(1, 'With a zscore threshold of %.2f, found %i (%.f%%) press-up and %i(%.f%%) press-down units, %i non-responsive (%i total).\n', sigmaThreshold, sum(isPressUp), sum(isPressUp)/sum(isResponsive)*100, sum(isPressDown), sum(isPressDown)/sum(isResponsive)*100, sum(~isResponsive), sum(isValid))

%% Calculate binned average
edges = 0:1:10;
[PressUp.X, PressUp.T, PressUp.N, PressUp.S, PressUp.B] = euPress(isPressUp).getBinnedTrialAverage('rate', 1:1.5:10, 'press', 'window', [-10, 1], 'normalize', true);
[PressDown.X, PressDown.T, PressDown.N, PressDown.S, PressDown.B] = euPress(isPressDown).getBinnedTrialAverage('rate', 1:1.5:10, 'press', 'window', [-10, 1], 'normalize', true);
[PressUpRaw.X, PressUpRaw.T, PressUpRaw.N, PressUpRaw.S, PressUpRaw.B] = euPress(isPressUp).getBinnedTrialAverage('rate', 1:1.5:10, 'press', 'window', [-10, 1], 'normalize', false);
[PressDownRaw.X, PressDownRaw.T, PressDownRaw.N, PressDownRaw.S, PressDownRaw.B] = euPress(isPressDown).getBinnedTrialAverage('rate', 1:1.5:10, 'press', 'window', [-10, 1], 'normalize', false);

%% Plot binned average
figure();
axUp = subplot(2, 2, 1);
axUpRaw = subplot(2, 2, 3);
axDown = subplot(2, 2, 2);
axDownRaw = subplot(2, 2, 4);
ax = [axUp, axUpRaw, axDown, axDownRaw];
EphysUnit.plotBinnedTrialAverage(ax(1), PressUp, [1,12]);
EphysUnit.plotBinnedTrialAverage(ax(2), PressUpRaw, [1,12]);
EphysUnit.plotBinnedTrialAverage(ax(3), PressDown, [1,12]);
EphysUnit.plotBinnedTrialAverage(ax(4), PressDownRaw, [1,12]);
title([axUp, axUpRaw], sprintf('Press Up (%i units)', nnz(isPressUp)))
title([axDown, axDownRaw], sprintf('Press Down (%i units)', nnz(isPressDown)))
xlabel(ax, 'Time from cue (s)')
ylabel([axUp, axDown], 'Spike rate (modified z-score)')
ylabel([axUpRaw, axDownRaw], 'Spike rate (sp/s)')
clear axUp axDown axUpRaw axDownRaw ax

%% Calculate binned average using randomly partitioned partial data set
IUp = randperm(nnz(isPressUp), 25);
IDown = randperm(nnz(isPressDown), 25);
euPressUp = euPress(isPressUp);
euPressDown = euPress(isPressDown);

[PressUp.X, PressUp.T, PressUp.N, PressUp.S, PressUp.B] = euPressUp(IUp).getBinnedTrialAverage('rate', 1:1.5:10, 'press', 'window', [-10, 1], 'normalize', true);
[PressDown.X, PressDown.T, PressDown.N, PressDown.S, PressDown.B] = euPressDown(IDown).getBinnedTrialAverage('rate', 1:1.5:10, 'press', 'window', [-10, 1], 'normalize', true);
[PressUpRaw.X, PressUpRaw.T, PressUpRaw.N, PressUpRaw.S, PressUpRaw.B] = euPressUp(IUp).getBinnedTrialAverage('rate', 1:1.5:10, 'press', 'window', [-10, 1], 'normalize', false);
[PressDownRaw.X, PressDownRaw.T, PressDownRaw.N, PressDownRaw.S, PressDownRaw.B] = euPressDown(IDown).getBinnedTrialAverage('rate', 1:1.5:10, 'press', 'window', [-10, 1], 'normalize', false);

figure();
axUp = subplot(2, 2, 1);
axUpRaw = subplot(2, 2, 3);
axDown = subplot(2, 2, 2);
axDownRaw = subplot(2, 2, 4);
ax = [axUp, axUpRaw, axDown, axDownRaw];
EphysUnit.plotBinnedTrialAverage(ax(1), PressUp, [1,12]);
EphysUnit.plotBinnedTrialAverage(ax(2), PressUpRaw, [1,12]);
EphysUnit.plotBinnedTrialAverage(ax(3), PressDown, [1,12]);
EphysUnit.plotBinnedTrialAverage(ax(4), PressDownRaw, [1,12]);
title([axUp, axUpRaw], sprintf('Press Up (%i units)', nnz(isPressUp)))
title([axDown, axDownRaw], sprintf('Press Down (%i units)', nnz(isPressDown)))
xlabel(ax, 'Time from cue (s)')
ylabel([axUp, axDown], 'Spike rate (modified z-score)')
ylabel([axUpRaw, axDownRaw], 'Spike rate (sp/s)')
clear axUp axDown axUpRaw axDownRaw ax


%% Plot single units and save to disk
for i = 1:length(euPressUp)
    try
        [Sr.X, Sr.T, Sr.N, Sr.S, Sr.B] = euPressUp(i).getBinnedTrialAverage('rate', 1:1.5:10, 'press', 'window', [-10, 1], 'normalize', false);
        [Sn.X, Sn.T, Sn.N, Sn.S, Sn.B] = euPressUp(i).getBinnedTrialAverage('rate', 1:1.5:10, 'press', 'window', [-10, 1], 'normalize', true);
        
        fig = figure('Units', 'normalized', 'Position', [0, 0, 0.6, 0.9]);
        ax(1) = subplot(2, 1, 1);
        ax(2) = subplot(2, 1, 2);
        EphysUnit.plotBinnedTrialAverage(ax(1), Sr, [1, 12]);
        EphysUnit.plotBinnedTrialAverage(ax(2), Sn, [1, 12]);
        suptitle(euPressUp(1).getName('_'));
        
        print(fig, sprintf('C:\\SERVER\\Figures\\Single Units\\PressUp\\%s', euPressUp(i).getName('_')), '-dpng');
        
        close(fig)
    end
    clear Sr Sn
    close all
end

for i = 1:length(euPressDown)
    try
        [Sr.X, Sr.T, Sr.N, Sr.S, Sr.B] = euPressDown(i).getBinnedTrialAverage('rate', 1:1.5:10, 'press', 'window', [-10, 1], 'normalize', false);
        [Sn.X, Sn.T, Sn.N, Sn.S, Sn.B] = euPressDown(i).getBinnedTrialAverage('rate', 1:1.5:10, 'press', 'window', [-10, 1], 'normalize', true);
        
        fig = figure('Units', 'normalized', 'Position', [0, 0, 0.6, 0.9]);
        ax(1) = subplot(2, 1, 1);
        ax(2) = subplot(2, 1, 2);
        EphysUnit.plotBinnedTrialAverage(ax(1), Sr, [1, 12]);
        EphysUnit.plotBinnedTrialAverage(ax(2), Sn, [1, 12]);
        suptitle(euPressDown(1).getName('_'));
        
        print(fig, sprintf('C:\\SERVER\\Figures\\Single Units\\PressDown\\%s', euPressDown(i).getName('_')), '-dpng');
        
        close(fig)
    end
    clear Sr Sn
    close all
end

clear fig ax Sr Sn i


%% Acute recording to EphysUnit
ar = AcuteRecording.load();
%%
for i = 1:length(ar)
    clear eu
    try  
        eu = EphysUnit(ar(i), 'cullITI', true, 'extendedWindow', [-1, 2], 'readWaveforms', true);
    catch ME
        warning('Error while processing file %g (%s)', i, ar(i).expName);
    end
end

%% 
stats = [eu.SpikeRateStats];
msr = [stats.medianITI]; % Median ITI spike rate
clear stats

% Select SNr cells
minSpikeRate = 15;
minTrialLength = 2;
minNumTrials = 15;

isSNr = msr >= minSpikeRate;
hasPress = arrayfun(@(e) ~isempty(e.Trials.Press) && nnz(e.Trials.Press.duration() >= minTrialLength) >= minNumTrials, eu);
hasLick = arrayfun(@(e) ~isempty(e.Trials.Lick) && nnz(e.Trials.Lick.duration() >= minTrialLength) >= minNumTrials, eu);

[PETHPress.X, PETHPress.t, PETHPress.N] = eu(isSNr & hasPress).getPETH('count', 'press', [-6, 1], minTrialDuration=minTrialLength, normalize=[-4, -2]);
[PETHLick.X, PETHLick.t, PETHLick.N] = eu(isSNr & hasLick).getPETH('count', 'lick', [-6, 1], minTrialDuration=minTrialLength, normalize=[-4, -2]);
[PETHPressCompare.X, PETHPressCompare.t, PETHPressCompare.N] = eu(isSNr & hasPress & hasLick).getPETH('count', 'press', [-6, 1], minTrialDuration=minTrialLength, normalize=[-4, -2]);
[PETHLickCompare.X, PETHLickCompare.t, PETHLickCompare.N] = eu(isSNr & hasPress & hasLick).getPETH('count', 'lick', [-6, 1], minTrialDuration=minTrialLength, normalize=[-4, -2]);


%%
fprintf(1, 'Average press trial duration = %g sec\n', mean(arrayfun(@(e) mean(e.Trials.Press.duration(), 'omitnan'), eu), 'omitnan'));
fprintf(1, 'Average lick trial duration = %g sec\n', mean(arrayfun(@(e) mean(e.Trials.Lick.duration(), 'omitnan'), eu), 'omitnan'));
EphysUnit.plotPETH(PETHPress, xlim=[-3.5,0], clim=[-2, 2], sortWindow=[-3.5, 0], signWindow=[-0.2, 0], sortThreshold=0.6, negativeSortThreshold=0.3); title('Lever-press PETH')
EphysUnit.plotPETH(PETHLick, xlim=[-3.5,0], clim=[-2, 2], sortWindow=[-3.5, 0], signWindow=[-0.2, 0], sortThreshold=0.6, negativeSortThreshold=0.3); title('Lick PETH')



%%
theta = 0:0.1:2;
[eta, nUp, nDown, pUp, pDown] = plotResponseForThresholds(PETHPress, 0:0.1:2, [-0.25, 0]);


function [eta, nUp, nDown, pUp, pDown] = plotResponseForThresholds(PETH, theta, window)
    eta = transpose(mean(PETH.X(:, PETH.t >= window(1) & PETH.t <= window(2)), 2, 'omitnan'));
    nUp = zeros(length(theta), 1);
    nDown = zeros(length(theta), 1);
    for i = 1:length(theta)
        nUp(i) = nnz(eta >= theta(i));
        nDown(i) = nnz(eta <= -theta(i));
    end
    pUp = nUp ./ (nUp + nDown);
    pDown = nDown ./ (nUp + nDown);
    f = figure();
    ax(1) = subplot(1, 2, 1);
    xlabel(ax(1), 'Threshold (a.u.)');
    ylabel(ax(1), 'Number of units')
    hold(ax(1), 'on')
    plot(ax(1), theta, nUp, 'r', DisplayName='Up', LineWidth=2);
    plot(ax(1), theta, nDown, 'b', DisplayName='Down', LineWidth=2);
    hold(ax(1), 'off')
    ax(2) = subplot(1, 2, 2);
    xlabel(ax(2), 'Threshold (a.u.)');
    ylabel(ax(2), 'Percentage of responsive units')
    hold(ax(2), 'on')
    plot(ax(2), theta, pUp, 'r', DisplayName='Up', LineWidth=2);
    plot(ax(2), theta, pDown, 'b', DisplayName='Down', LineWidth=2);
    hold(ax(2), 'off')
    legend(ax(1))

    for i = 1:length(theta)
        fprintf(1, 'theta=%g, %g up (%.1f%%), %g down (%.1f%%)\n', theta(i), nUp(i), pUp(i), nDown(i), pDown(i));
    end

    f = figure();
    ax = axes(f);
    histogram(ax, eta, Normalization='probability');
    xlabel(ax, sprintf('Normalized move response [%g, %g]s', window(1), window(2)))
    ylabel(ax, 'Probability')
    title(ax, 'Distribution of lever-press response')
end

%% Lick
