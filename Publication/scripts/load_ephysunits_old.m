%% 1.2Alt Or just load lite version, without non-SNr cells, without waveforms, spikecounts or spikerates.
eu = EphysUnit.load('C:\SERVER\Units\Lite_NonDuplicate');
ar = AcuteRecording.load('C:\SERVER\Acute\AcuteRecording');

% 1.3 AnimalInfo
animalInfo = { ...
%     'daisy1', 'wt', 'F', -3.2, -1.6, 'tetrode'; ...
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

% 2.1 Parameters
clear p
p.minSpikeRate = 15;
p.minTrialDuration = 2;
p.minNumTrials = 30;
p.etaNorm = [-4, -2];
p.etaWindow = [-4, 2];
p.cueEtaWindow = [-2, 4];
p.metaWindowPress = [-0.3, 0];
p.metaWindowLick = [-0.3, 0];
% p.posRespThreshold = 1;
% p.negRespThreshold = -0.5;
p.binnedTrialEdges = 2:2:10;

p.minStimDuration = 1e-2;
p.maxStimDuration = 5e-2;
p.errStimDuration = 1e-3;
p.allowAltStimDuration = true;
p.etaWindowStim = [-0.2, 0.5];
p.metaWindowStim = [0, 0.1];


% %% 2.2.1 Cull non-SNr units to save memory (only once)
% msr = arrayfun(@(stats) stats.medianITI, [eu.SpikeRateStats]);
% isSNr = msr >= p.minSpikeRate;
% eu = eu(isSNr);
% fprintf(1, 'Kept %g out of %g SNr units with spike rate >= %g.\n', nnz(isSNr), length(msr), p.minSpikeRate)
% clearvars -except eu p ai

% Multiunit detection by ISI.
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
c.isMultiUnit = prcLowISI > 0.05;
c.isSingleUnit = prcLowISI <= 0.05;
euAll = eu;
eu = euAll(c.isSingleUnit);


%% Remove drifting units
% This is cool because it actually requires spike rate to be above 15sp/s
% for 95% of the session
c.isDrifting = detectDriftingUnits(eu, smoothWindow=300, tolerance=0.05, spikeRateThreshold=15);

eu = eu(~c.isDrifting);


% Find location of units from AR
euPos = NaN(length(eu), 3); % ml dv ap
c.hasPos = false(1, length(eu));
for iEu = 1:length(eu)
    iAr = find(strcmpi(eu(iEu).ExpName, {ar.expName}));
    if ~isempty(iAr)
        euPos(iEu, :) = ar(iAr).getProbeCoords(eu(iEu).Channel);
        c.hasPos(iEu) = true;
    end
end

%%

% 2.3.1  Basic summaries
% Baseline (median) spike rates
msr = arrayfun(@(stats) stats.median, [eu.SpikeRateStats]);

% Lick/Press responses
eta.press = eu.getETA('count', 'press', p.etaWindow, minTrialDuration=p.minTrialDuration, normalize=p.etaNorm);
eta.lick = eu.getETA('count', 'lick', p.etaWindow, minTrialDuration=p.minTrialDuration, normalize=p.etaNorm);
eta.pressRaw = eu.getETA('count', 'press', p.etaWindow, minTrialDuration=p.minTrialDuration, normalize='none');
eta.lickRaw = eu.getETA('count', 'lick', p.etaWindow, minTrialDuration=p.minTrialDuration, normalize='none');
eta.pressCue = eu.getETA('count', 'press', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize=eta.press.stats, includeInvalid=true);
eta.lickCue = eu.getETA('count', 'lick', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize=eta.lick.stats, includeInvalid=true);
eta.pressCueRaw = eu.getETA('count', 'press', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize='none', includeInvalid=true);
eta.lickCueRaw = eu.getETA('count', 'lick', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize='none', includeInvalid=true);

meta.press = transpose(mean(eta.press.X(:, eta.press.t >= p.metaWindowPress(1) & eta.press.t <= p.metaWindowPress(2)), 2, 'omitnan'));
meta.lick = transpose(mean(eta.lick.X(:, eta.lick.t >= p.metaWindowLick(1) & eta.lick.t <= p.metaWindowLick(2)), 2, 'omitnan'));
meta.pressRaw = transpose(mean(eta.pressRaw.X(:, eta.pressRaw.t >= p.metaWindowPress(1) & eta.pressRaw.t <= p.metaWindowPress(2)), 2, 'omitnan'));
meta.lickRaw = transpose(mean(eta.lickRaw.X(:, eta.lickRaw.t >= p.metaWindowLick(1) & eta.lickRaw.t <= p.metaWindowLick(2)), 2, 'omitnan'));
meta.pressRawBaseline = transpose(mean(eta.pressRaw.X(:, eta.pressRaw.t >= p.etaNorm(1) & eta.pressRaw.t <= p.etaNorm(2)), 2, 'omitnan'));
meta.lickRawBaseline = transpose(mean(eta.lickRaw.X(:, eta.lickRaw.t >= p.etaNorm(1) & eta.lickRaw.t <= p.etaNorm(2)), 2, 'omitnan'));

p.metaWindowCue = [-0.7, 0.2];
meta.pressCue = transpose(mean(eta.pressCue.X(:, eta.pressCue.t >= p.metaWindowCue(1) & eta.pressCue.t <= p.metaWindowCue(2)), 2, 'omitnan'));
meta.lickCue = transpose(mean(eta.lickCue.X(:, eta.lickCue.t >= p.metaWindowCue(1) & eta.lickCue.t <= p.metaWindowCue(2)), 2, 'omitnan'));
meta.pressCueRaw = transpose(mean(eta.pressCueRaw.X(:, eta.pressCueRaw.t >= p.metaWindowCue(1) & eta.pressCueRaw.t <= p.metaWindowCue(2)), 2, 'omitnan'));
meta.lickCueRaw = transpose(mean(eta.lickCueRaw.X(:, eta.lickCueRaw.t >= p.metaWindowCue(1) & eta.lickCueRaw.t <= p.metaWindowCue(2)), 2, 'omitnan'));

% etaSmooth.pressRaw = eu.getETA('rate', 'press', p.etaWindow, minTrialDuration=p.minTrialDuration, normalize='none');
% etaSmooth.lickRaw = eu.getETA('rate', 'lick', p.etaWindow, minTrialDuration=p.minTrialDuration, normalize='none');
% etaSmooth.pressCueRaw = eu.getETA('rate', 'press', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize='none', includeInvalid=true);
% etaSmooth.lickCueRaw = eu.getETA('rate', 'lick', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize='none', includeInvalid=true);
etaSmooth.pressCue = eu.getETA('rate', 'press', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize=etaSmooth.press.stats, includeInvalid=true);
etaSmooth.lickCue = eu.getETA('rate', 'lick', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize=etaSmooth.lick.stats, includeInvalid=true);


% 2.3.2 Basic summaries (fast)
% hasPress/hasLick
c.hasPress = arrayfun(@(e) nnz(e.getTrials('press').duration() >= p.minTrialDuration) >= p.minNumTrials, eu);
c.hasLick = arrayfun(@(e) nnz(e.getTrials('lick').duration() >= p.minTrialDuration) >= p.minNumTrials, eu);

% press/lick x Up/Down
% c.isPressUp =         c.hasPress & meta.press >= p.posRespThreshold;
% c.isPressDown =       c.hasPress & meta.press <= p.negRespThreshold;
% c.isPressResponsive = c.isPressUp | c.isPressDown;
% c.isLickUp =          c.hasLick & meta.lick >= p.posRespThreshold;
% c.isLickDown =        c.hasLick & meta.lick <= p.negRespThreshold;
% c.isLickResponsive =  c.isLickUp | c.isLickDown;

% animal info
c.isWT = arrayfun(@(eu) strcmpi(getAnimalInfo(eu, ai, 'strain'), 'wt'), eu);
c.isD1 = arrayfun(@(eu) strcmpi(getAnimalInfo(eu, ai, 'strain'), 'd1-cre'), eu);
c.isA2A = arrayfun(@(eu) strcmpi(getAnimalInfo(eu, ai, 'strain'), 'a2a-cre'), eu);
c.isAi80 = arrayfun(@(eu) strcmpi(getAnimalInfo(eu, ai, 'strain'), 'd1-cre;dlx-flp;ai80'), eu);
c.isDAT = arrayfun(@(eu) strcmpi(getAnimalInfo(eu, ai, 'strain'), 'dat-cre'), eu);
c.isAcute = ismember(eu.getAnimalName, {'daisy14', 'daisy15', 'daisy16', 'desmond23', 'desmond24', 'desmond25', 'desmond26', 'desmond27'});

%% Calculate neural onset times
% Use smooth ETA for better temporal resolution
etaSmooth.press = eu.getETA('rate', 'press', p.etaWindow, minTrialDuration=p.minTrialDuration, normalize=p.etaNorm);
etaSmooth.lick = eu.getETA('rate', 'lick', p.etaWindow, minTrialDuration=p.minTrialDuration, normalize=p.etaNorm);

%% Calculate onset
p.etaOnsetThreshold = 0.25;
p.etaSortWindow = [-3, 0];
p.etaSignWindow = [-0.3, 0];
p.etaOnsetPattern = [zeros(1, 50), ones(1, 50)];

fig = figure(Units='normalized', Position=[0.1, 0.1, 0.8, 0.8]);
ax(1) = subplot(1, 2, 1);
ax(2) = subplot(1, 2, 2);
n = length(eu);
onset = struct(press=NaN(n, 1), lick=NaN(n, 1), pressOrder=NaN(n, 1), lickOrder=NaN(n, 1));
[~, onset.pressOrder(c.hasPress), ~, onset.press(c.hasPress)] = EphysUnit.plotETA(ax(1), etaSmooth.press, c.hasPress, sortWindow=p.etaSortWindow, signWindow=p.etaSignWindow, ...
    sortThreshold=p.etaOnsetThreshold, negativeSortThreshold=p.etaOnsetThreshold, clim=[-2, 2], onsetPattern=p.etaOnsetPattern, xlim=[-4, 0.5], ...
    onsetDirection='reverse');
hold(ax(1), 'on')
x = onset.press(c.hasPress);
I = onset.pressOrder(c.hasPress);
plot(ax(1), x(I), 1:length(I))

[~, onset.lickOrder(c.hasLick), ~, onset.lick(c.hasLick)] = EphysUnit.plotETA(ax(2), etaSmooth.lick, c.hasLick, sortWindow=p.etaSortWindow, signWindow=p.etaSignWindow, ...
    sortThreshold=p.etaOnsetThreshold, negativeSortThreshold=p.etaOnsetThreshold, clim=[-2, 2], onsetPattern=p.etaOnsetPattern, xlim=[-4, 0.5], ...
    onsetDirection='reverse');
hold(ax(2), 'on')
x = onset.lick(c.hasLick);
I = onset.lickOrder(c.hasLick);
plot(ax(2), x(I), 1:length(I))

clim(ax(1), [-1.5, 1.5])
clim(ax(2), [-1.5, 1.5])
clear n x I

%% Save metadata and units
eu.save('C:\SERVER\Units\Lite_NonDuplicate_NonDrift')
save('C:\SERVER\Units\meta_Lite_NonDuplicate_NonDrift.mat', 'p', 'c', 'eta', 'etaSmooth', 'euPos', 'meta', 'msr', 'onset')

%%
function info = getAnimalInfo(eu, ai, field)
    i = find(strcmpi({ai.name}, eu.getAnimalName()));
    assert(length(i) == 1, eu.getAnimalName())
    info = ai(i).(field);
end
