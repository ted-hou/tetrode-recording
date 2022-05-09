load('C:\SERVER\PETH_All_aggregate_20220201.mat')

%% Make eu objects (SLOW, takes ~60min)
eu = EphysUnit(PETH, 'cullITI', true, 'extendedWindow', [-1, 2], 'readWaveforms', true);

%% Calculate eta/peth
[eta, ~, ~] = eu.getPETH('count', 'press', [-0.45, -0.05], 'minTrialDuration', 1, 'normalize', 'iti');
eta = transpose(nanmean(eta, 2));

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
f = figure('Units', 'normalized', 'OuterPosition', [0, 0.33, 1, 0.5]);
ax = subplot(1, 3, 1);
histogram(ax, msr(isValid))
xlabel('Median spike rate (sp/s)')

% Plot distribution of response magnitude for SNr cells
ax = subplot(1, 3, 2);
histogram(ax, eta(isValid))
xlabel('Mean response magnitude (modified z-score')

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

%% Separate press-units by direction
sigmaThreshold = 1.2;
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
    clear Sr Sn
end

for i = 1:length(euPressDown)
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
    clear Sr Sn
end

clear fig ax Sr Sn i