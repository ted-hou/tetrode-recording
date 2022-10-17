%% 1.1 Make eu objects (SLOW, takes ~60min)
load('C:\SERVER\PETH_All_aggregate_20220201.mat')
ar = AcuteRecording.load();

eu = EphysUnit(PETH, 'cullITI', true, 'extendedWindow', [-1, 2], 'readWaveforms', true);
for i = 1:length(ar)
    clear eu
    try  
        eu = EphysUnit(ar(i), 'cullITI', true, 'extendedWindow', [-1, 2], 'readWaveforms', true);
    catch ME
        warning('Error while processing file %g (%s)', i, ar(i).expName);
    end
end

%% 1.2 Alternatively, load eu objects from disk (SLOW, ~20min)
eu = EphysUnit.load('C:\SERVER\Units', waveforms=false, spikecounts=false, spikerates=false);
%% 1.2Alt Or just load lite version, without non-SNr cells, without waveforms, spikecounts or spikerates.
eu = EphysUnit.load('C:\SERVER\Units\Lite_NonDuplicate');
ar = AcuteRecording.load('C:\SERVER\Acute\AcuteRecording');

% 1.3 AnimalInfo
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

% 2.1 Parameters
clear p
p.minSpikeRate = 15;
p.minTrialDuration = 2;
p.minNumTrials = 30;
p.etaNorm = [-4, -2];
p.etaWindow = [-4, 2];
p.metaWindow = [-0.2, 0];
p.metaWindowLickVsPress = [-0.1, 0.1];
p.posRespThreshold = 1;
p.negRespThreshold = -0.5;
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




% 2.3.1  Basic summaries
% Baseline (median) spike rates
msr = arrayfun(@(stats) stats.medianITI, [eu.SpikeRateStats]);

% Lick/Press responses
eta.press = eu.getETA('count', 'press', p.etaWindow, minTrialDuration=p.minTrialDuration, normalize=p.etaNorm);
meta.press = transpose(mean(eta.press.X(:, eta.press.t >= p.metaWindow(1) & eta.press.t <= p.metaWindow(2)), 2, 'omitnan'));
eta.lick = eu.getETA('count', 'lick', p.etaWindow, minTrialDuration=p.minTrialDuration, normalize=p.etaNorm);
meta.lick = transpose(mean(eta.lick.X(:, eta.lick.t >= p.metaWindowLickVsPress(1) & eta.lick.t <= p.metaWindowLickVsPress(2)), 2, 'omitnan'));

eta.pressRaw = eu.getETA('count', 'press', p.etaWindow, minTrialDuration=p.minTrialDuration, normalize='none');
meta.pressRaw = transpose(mean(eta.pressRaw.X(:, eta.pressRaw.t >= p.metaWindow(1) & eta.pressRaw.t <= p.metaWindow(2)), 2, 'omitnan'));

eta.lickRaw = eu.getETA('count', 'lick', p.etaWindow, minTrialDuration=p.minTrialDuration, normalize='none');
meta.lickRaw = transpose(mean(eta.lickRaw.X(:, eta.lickRaw.t >= p.metaWindowLickVsPress(1) & eta.lickRaw.t <= p.metaWindowLickVsPress(2)), 2, 'omitnan'));


% 2.3.2 Basic summaries (fast)
% hasPress/hasLick
c.hasPress = arrayfun(@(e) nnz(e.getTrials('press').duration() >= p.minTrialDuration) >= p.minNumTrials, eu);
c.hasLick = arrayfun(@(e) nnz(e.getTrials('lick').duration() >= p.minTrialDuration) >= p.minNumTrials, eu);

% press/lick x Up/Down
c.isPressUp =         c.hasPress & meta.press >= p.posRespThreshold;
c.isPressDown =       c.hasPress & meta.press <= p.negRespThreshold;
c.isPressResponsive = c.isPressUp | c.isPressDown;
c.isLickUp =          c.hasLick & meta.lick >= p.posRespThreshold;
c.isLickDown =        c.hasLick & meta.lick <= p.negRespThreshold;
c.isLickResponsive =  c.isLickUp | c.isLickDown;

% animal info
c.isWT = arrayfun(@(eu) strcmpi(getAnimalInfo(eu, ai, 'strain'), 'wt'), eu);
c.isD1 = arrayfun(@(eu) strcmpi(getAnimalInfo(eu, ai, 'strain'), 'd1-cre'), eu);
c.isA2A = arrayfun(@(eu) strcmpi(getAnimalInfo(eu, ai, 'strain'), 'a2a-cre'), eu);
c.isAi80 = arrayfun(@(eu) strcmpi(getAnimalInfo(eu, ai, 'strain'), 'd1-cre;dlx-flp;ai80'), eu);
c.isDAT = arrayfun(@(eu) strcmpi(getAnimalInfo(eu, ai, 'strain'), 'dat-cre'), eu);
c.isAcute = ismember(eu.getAnimalName, {'daisy14', 'daisy15', 'daisy16', 'desmond23', 'desmond24', 'desmond25', 'desmond26', 'desmond27'});

fprintf(1, ['%g total SNr units (baseline spike rate > %g):\n' ...
    '\t%g with press trials;\n' ...
    '\t%g with lick trials;\n' ...
    '\t%g with both (at least %g trials).\n'], ...
    length(eu), p.minSpikeRate, nnz(c.hasPress), nnz(c.hasLick), nnz(c.hasPress & c.hasLick), p.minNumTrials)

fprintf(1, ['%g units with at least %g press trials:\n' ...
    '\t%g (%.0f%%) are excited (meta>=%g);\n' ...
    '\t%g (%.0f%%) are inhibited (meta<=%g).\n'], ...
    nnz(c.hasPress), p.minNumTrials, ...
    nnz(c.isPressUp), 100*nnz(c.isPressUp)/nnz(c.isPressResponsive), p.posRespThreshold, ...
    nnz(c.isPressDown), 100*nnz(c.isPressDown)/nnz(c.isPressResponsive), p.negRespThreshold);

fprintf(1, ['%g units with at least %g lick trials:\n' ...
    '\t%g (%.0f%%) are excited (meta>=%g);\n' ...
    '\t%g (%.0f%%) are inhibited (meta<=%g).\n'], ...
    nnz(c.hasLick), p.minNumTrials, ...
    nnz(c.isLickUp), 100*nnz(c.isLickUp)/nnz(c.isLickResponsive), p.posRespThreshold, ...
    nnz(c.isLickDown), 100*nnz(c.isLickDown)/nnz(c.isLickResponsive), p.negRespThreshold);

nTotal = nnz(c.isPressResponsive & c.isLickResponsive);
fprintf(1, ['%g units with press AND lick trials:\n' ...
    '\t%g (%.0f%%) are press-excited AND lick-excited;\n' ...
    '\t%g (%.0f%%) are press-inhibited AND lick-inhibited;\n' ...
    '\t%g (%.0f%%) are press-excited AND lick-inhibited;\n' ...
    '\t%g (%.0f%%) are press-inhibited AND lick-excited;\n'], ...
    nnz(c.hasPress & c.hasLick), ...
    nnz(c.isPressUp & c.isLickUp), 100*nnz(c.isPressUp & c.isLickUp)/nTotal, ...
    nnz(c.isPressDown & c.isLickDown), 100*nnz(c.isPressDown & c.isLickDown)/nTotal, ...
    nnz(c.isPressUp & c.isLickDown), 100*nnz(c.isPressUp & c.isLickDown)/nTotal, ...
    nnz(c.isPressDown & c.isLickUp), 100*nnz(c.isPressDown & c.isLickUp)/nTotal)   
clear nTotal

%% 4. Stim response
%

%% 4.1.1.2 Example stim responses (2dSPN 2iSPN)
close all
maxTrials = 200;
clear name
name{1} = 'daisy9_20211103_Electrode5_Unit1'; % ai80.inhibited
name{2} = 'daisy9_20211028_Electrode10_Unit1'; % ai80.excited

clear rrdd
for iUnit = 1:2
    iEu = find(strcmpi(name{iUnit}, eu.getName('_')));
%     rrdd(i) = rdStimFiltered(iEu);
    rrdd(iUnit) = eu(iEu).getRasterData('stim', window=[-0.1, 0.1], minTrialDuration=0.01, maxTrialDuration=0.01, sort=false, alignTo='start');
end


% 4.1.1.2 Select example units and plot stim raster (GALVO ACUTE UNITS)
clear selUnit
maxTrials = Inf;

selUnit(3).expName = 'daisy15_20220511'; %A2A
selUnit(3).channel = 37;
selUnit(3).unit = 1;

selUnit(4).expName = 'desmond24_20220510'; %A2A-suppressed, kind of weak
selUnit(4).channel = 52;
selUnit(4).unit = 1;

for iUnit = 3:4
    iEu = find(strcmpi(selUnit(iUnit).expName, {eu.ExpName}) & [eu.Channel] == selUnit(iUnit).channel & [eu.Unit] == selUnit(iUnit).unit);
    iAr = find(strcmpi(selUnit(iUnit).expName, {ar.expName}));
    iBsr = find([ar(iAr).bsr.channel] == selUnit(iUnit).channel & [ar(iAr).bsr.unit] == selUnit(iUnit).unit);
    [~, IPulse] = ar(iAr).selectStimResponse(Light=2, Duration=0.01, MLRank=2, DVRank=3);
    trials = Trial(ar(iAr).stim.tOn(IPulse), ar(iAr).stim.tOff(IPulse));
    assert(~isempty(iEu))
%     rrdd(i) = rdStimFiltered(iEu);
    rrdd(iUnit) = eu(iEu).getRasterData('stim', trials=trials, window=[-0.1, 0.1], minTrialDuration=0, maxTrialDuration=Inf, sort=false, alignTo='start');
end



ax = plotDoubleRaster(rrdd(1), rrdd(2), xlim=[-100, 100], timeUnit='ms', maxTrials=maxTrials);
fig = ax(1).Parent;
set(fig, Units='pixels', Position=[0, 0, 300, 600])
set(ax, FontSize=10)
xlabel(ax(1), '')
title(ax, '')
title(ax(1), 'dSPN stim')
ax(1).Legend.Orientation = 'horizontal';
ax(1).Legend.Location = 'southoutside';
ax(2).Legend.Visible = false;

% 4.1.1.2 Select example units and plot stim raster (GALVO ACUTE UNITS)
clear selUnit
maxTrials = 50;

selUnit(1).expName = 'daisy15_20220511'; %A2A
selUnit(1).channel = 37;
selUnit(1).unit = 1;

selUnit(2).expName = 'desmond24_20220510'; %A2A-suppressed, kind of weak
selUnit(2).channel = 52;
selUnit(2).unit = 1;


clear rrdd
for i = 1:2
    iEu = find(strcmpi(selUnit(i).expName, {eu.ExpName}) & [eu.Channel] == selUnit(i).channel & [eu.Unit] == selUnit(i).unit);
    iAr = find(strcmpi(selUnit(i).expName, {ar.expName}));
    iBsr = find([ar(iAr).bsr.channel] == selUnit(i).channel & [ar(iAr).bsr.unit] == selUnit(i).unit);
    [~, IPulse] = ar(iAr).selectStimResponse(Light=2, Duration=0.01, MLRank=2, DVRank=3);
    trials = Trial(ar(iAr).stim.tOn(IPulse), ar(iAr).stim.tOff(IPulse));
    assert(~isempty(iEu))
%     rrdd(i) = rdStimFiltered(iEu);
    rrdd(i) = eu(iEu).getRasterData('stim', trials=trials, window=[-0.1, 0.1], minTrialDuration=0, maxTrialDuration=Inf, sort=false, alignTo='start');
end

ax = plotDoubleRaster(rrdd(1), rrdd(2), xlim=[-100, 100], timeUnit='ms', maxTrials=maxTrials);
fig = ax(1).Parent;
set(fig, Units='pixels', Position=[0, 0, 300, 600])
set(ax, FontSize=10)
xlabel(ax(1), '')
title(ax, '')
title(ax(1), 'iSPN stim')
ax(1).Legend.Orientation = 'horizontal';
ax(1).Legend.Location = 'southoutside';
ax(2).Legend.Visible = false;


%% 4.2 ISI analysis to find latencies
clear isi rdStim

if ~exist('ar')
    ar = AcuteRecording.load('C:\SERVER\Acute\AcuteRecording');
end

c.hasAnyStimTrials = arrayfun(@(eu) ~isempty(eu.getTrials('stim')), eu);
c.excludeD1 = strcmpi('daisy4_20190429', {eu.ExpName}) | strcmpi('daisy4_20190404', {eu.ExpName});% Bad timing alignment?
c.exclude = ...
    (strcmpi('desmond26_20220531', {eu.ExpName})  & [eu.Channel] == 13 & [eu.Unit] == 1) | ...
    (strcmpi('daisy16_20220502', {eu.ExpName})  & [eu.Channel] == 76 & [eu.Unit] == 1);% No spikes in DLS stim

xl = [-.2, .1]; % Extend left by an extra 100ms to get accurate ISI curves. 

for iEu = find(c.hasAnyStimTrials & (c.isAi80 | c.isA2A) & ~c.exclude)
    % Attempt to select specific galvo trials to simplify conditions
    iAr = find(strcmpi(eu(iEu).ExpName, {ar.expName}));
    if isempty(iAr)
        rdStim(iEu) = eu(iEu).getRasterData('stim', window=xl, ...
            minTrialDuration=0.01, maxTrialDuration=0.01, sort=false, alignTo='start');
    else
        iBsr = find([ar(iAr).bsr.channel] == eu(iEu).Channel & [ar(iAr).bsr.unit] == eu(iEu).Unit);
        [~, IPulse] = ar(iAr).selectStimResponse(Light=2, Duration=0.01, MLRank=2, DVRank=3);
        trials = Trial(ar(iAr).stim.tOn(IPulse), ar(iAr).stim.tOff(IPulse));
        rdStim(iEu) = eu(iEu).getRasterData('stim', window=xl, ...
            minTrialDuration=0.01, maxTrialDuration=0.01, sort=false, alignTo='start', trials=trials);
    end
end

c.hasStim = arrayfun(@(rd) ~isempty(rd.I), rdStim);

% D1 version (100ms pulses, 100ms ITI)
xlD1 = [-.2, .1]; % Extend left by an extra 100ms to get accurate ISI curves. 

for iEu = find(c.hasAnyStimTrials & c.isD1 & ~c.excludeD1)
    rdStim(iEu) = eu(iEu).getRasterData('stimfirstpulse', window=xlD1, ...
        minTrialDuration=0.01, maxTrialDuration=0.101, sort=false, alignTo='start');
end

c.hasStim = arrayfun(@(rd) ~isempty(rd.I), rdStim);

fprintf(1, '%d units with requested stim conditions.\n', nnz(c.hasStim))

 
warning('off')
close all
clear isi sel
sel = c.hasStim;
isi(sel) = getISILatencies(rdStim(sel), xlim=[-100, 100], peakThreshold=0.75, posMultiplier=1, minProminence=2, onsetThreshold=0.25, ...
    showPlots=false, savePlots=false, maxTrials=Inf);
c.isStimUp = false(1, length(eu));
c.isStimDown = false(1, length(eu));
for iEu = 1:length(isi)
    if ~isempty(isi(iEu).peak) && ~isnan(isi(iEu).peak)
        c.isStimUp(iEu) = isi(iEu).peak < isi(iEu).baseline;
        c.isStimDown(iEu) = isi(iEu).peak > isi(iEu).baseline;
    end
end
c.hasStimResponse = c.isStimUp | c.isStimDown;
clear iEu sel

fprintf(1, '%d A2A units, %d up, %d down.\n', nnz(c.hasStim & c.isA2A), nnz(c.hasStim & c.isA2A & c.isStimUp), nnz(c.hasStim & c.isA2A & c.isStimDown))
fprintf(1, '%d Ai80 units, %d up, %d down.\n', nnz(c.hasStim & c.isAi80), nnz(c.hasStim & c.isAi80 & c.isStimUp), nnz(c.hasStim & c.isAi80 & c.isStimDown))
fprintf(1, '%d D1 units, %d up, %d down.\n', nnz(c.hasStim & c.isD1), nnz(c.hasStim & c.isD1 & c.isStimUp), nnz(c.hasStim & c.isD1 & c.isStimDown))
warning('on')

% Summarize responses
isiStim = isi(c.hasStim & c.hasStimResponse);
latencies = [isiStim.onsetLatency];
responses = [isiStim.peak] - [isiStim.isi0];

subplot(1, 2, 1)
histogram(latencies, 0:2:100)
title(sprintf('Distribution of latencies\nmin=%d, max=%d, mean=%.1f, median=%.1f', ...
    min(latencies), max(latencies), mean(latencies), median(latencies)))
subplot(1, 2, 2)
histogram(responses, 20)
title(sprintf('Distribution of responses\nMean=%.1f, %.1f, Median=%.1f, %.1f', mean(responses(responses<0)), mean(responses(responses>0)), median(responses(responses<0)), median(responses(responses>0))))

[~, pp] = ttest2(-responses(responses < 0), responses(responses > 0));
xlabel(sprintf('p=%.5f', pp))

%% Plot and save ISI analysis for manual checking
close all
sel = c.hasStimResponse;
getISILatencies(rdStim(sel), xlim=[-100, 100], peakThreshold=0.75, posMultiplier=1, minProminence=2, onsetThreshold=0.25, ...
    showPlots=true, savePlots=true, maxTrials=Inf);

%% Example plots, ISI analysis for 2 example units
close all
clear selUnits
clear selEu
selUnits(1).expName = 'daisy9_20211028';
selUnits(2).expName = 'daisy9_20211028';

selUnits(1).electrode = 4;
selUnits(2).electrode = 10;
selUnits(1).unit = 1;
selUnits(2).unit = 1;


for i = 1:length(selUnits)
    if ~isempty(selUnits(i).electrode)
        selEu(i) = find(strcmpi({eu.ExpName}, selUnits(i).expName) & [eu.Electrode] == selUnits(i).electrode & [eu.Unit] == selUnits(i).unit);
    else
        selEu(i) = find(strcmpi({eu.ExpName}, selUnits(i).expName) & [eu.Channel] == selUnits(i).channel & [eu.Unit] == selUnits(i).unit);
    end
end
I = find(c.hasStim);
selEu = [selEu, I([isi(c.hasStim).baseline] > 1000/15)];

getISILatencies(rdStim(selEu), xlim=[-100, 100], peakThreshold=0.75, posMultiplier=1, minProminence=2, onsetThreshold=0.25, ...
    showPlots=true, savePlots=false, maxTrials=Inf);
% clear selUnits i

%% Calculate ETAs
p.binWidthStim = 0.005;

eta.stim = eu.getETA('count', 'stim', [-0.5, 0.5], ...
    resolution=p.binWidthStim, ...
    minTrialDuration=0.01, maxTrialDuration=0.01, ...
    findSingleTrialDuration='min', normalize=[-0.2, 0], includeInvalid=false);

eta.stimD1 = eu.getETA('count', 'stimfirstpulse', [-0.5, 0.5], ...
    resolution=p.binWidthStim, ...
    minTrialDuration=0.01, maxTrialDuration=0.101, ...
    findSingleTrialDuration='max', normalize=[-0.2, 0], includeInvalid=false);
% For acute recording, pick DLS at 2.8DV
for iEu = 1:length(eu)
    iAr = find(strcmpi(eu(iEu).ExpName, {ar.expName}));
    if ~isempty(iAr)
        iBsr = find([ar(iAr).bsr.channel] == eu(iEu).Channel & [ar(iAr).bsr.unit] == eu(iEu).Unit);
        [~, IPulse] = ar(iAr).selectStimResponse(Light=2, Duration=0.01, MLRank=2, DVRank=3);
        trials = Trial(ar(iAr).stim.tOn(IPulse), ar(iAr).stim.tOff(IPulse));
        tempEta = eu(iEu).getETA('count', 'stim', [-0.5, 0.5], ...
            trials=trials, resolution=p.binWidthStim, ...
            minTrialDuration=0.01, maxTrialDuration=0.01, ...
            findSingleTrialDuration='min', normalize=[-0.2, 0], includeInvalid=false);
        eta.stim.X(iEu, :) = tempEta.X;
        eta.stim.N(iEu) = tempEta.N;
        eta.stim.D(iEu) = tempEta.D;
        eta.stim.stats(iEu) = tempEta.stats;
        clear tempEta
    end
end

%% Summarize and plot stim ETAs, latencies
close all
N = nnz(c.hasStimResponse & c.isA2A);
fprintf(1, 'Out of %g tested units, %g are responsive.\n\t%g (%.2f%%) were excited by A2A, %g(%.2f%%) were inhibited by A2A.\n', ...
    nnz(c.hasStim & c.isA2A), N, ...
    nnz(c.isStimUp & c.isA2A), nnz(c.isStimUp & c.isA2A) / N * 100, ...
    nnz(c.isStimDown & c.isA2A), nnz(c.isStimDown & c.isA2A) / N * 100);

N = nnz(c.hasStimResponse & c.isAi80);
fprintf(1, 'Out of %g tested units, %g are responsive.\n\t%g (%.2f%%) were excited by Ai80, %g(%.2f%%) were inhibited by Ai80.\n', ...
    nnz(c.hasStim & c.isAi80), N, ...
    nnz(c.isStimUp & c.isAi80), nnz(c.isStimUp & c.isAi80) / N * 100, ...
    nnz(c.isStimDown & c.isAi80), nnz(c.isStimDown & c.isAi80) / N * 100);

N = nnz(c.hasStimResponse & c.isD1);
fprintf(1, 'Out of %g tested units, %g are responsive.\n\t%g (%.2f%%) were excited by D1 (100ms), %g(%.2f%%) were inhibited by D1 (100ms).\n', ...
    nnz(c.hasStim & c.isD1), N, ...
    nnz(c.isStimUp & c.isD1), nnz(c.isStimUp & c.isD1) / N * 100, ...
    nnz(c.isStimDown & c.isD1), nnz(c.isStimDown & c.isD1) / N * 100);
clear N

% 4.2.2 Plot stim heatmap
SEL = { ...
    c.hasStimResponse & c.isA2A; ...
    c.hasStimResponse & c.isAi80; ...
    c.hasStimResponse & c.isD1};
ETA = {eta.stim, eta.stim, eta.stimD1};
NAME = {'iSPN', 'dSPN', 'Drd1-Cre'};

for i = 1:3
    sel = SEL{i};

    dmin = min(ETA{i}.D(sel));
    dmax = max(ETA{i}.D(sel));
    if dmin==dmax
        text = sprintf('%s stim (%g ms)', NAME{i}, 1000*dmin);
    else
        text = sprintf('%s stim (%g-%g ms)', NAME{i}, 1000*dmin, 1000*dmax);
    end
    
    latencies = [isi(sel).onsetLatency];
    responses = [isi(sel).peak] - [isi(sel).isi0];
    responseSigns = sign(responses);
    peakLatencies = [isi(sel).peakLatency];
    [~, ISort] = sort(responseSigns.*(latencies.*1e4 + peakLatencies*1e3 + abs(responses)), 'ascend');
    
    ax = EphysUnit.plotETA(ETA{i}, sel, xlim=[-100, 100], clim=[-1, 1], ...
        timeUnit='ms', order=ISort, ...
        event='opto onset'); 
    title(ax, text)
    ax.FontSize = 10;
    hold(ax, 'on')
    if dmin==dmax
        yl = ax.YLim;
        patch(ax, [0, dmin, dmin, 0].*1000, [yl(1), yl(1), yl(2), yl(2)], 'b', FaceAlpha=0.1, EdgeAlpha=0)
    end
    hold(ax, 'off')
    title(ax, text)
    ax.Parent.OuterPosition(3) = 0.2;
    ax.Parent.OuterPosition(4) = 0.5;
end

SEL = { ...
    c.isStimUp & c.isA2A; ...
    c.isStimUp & c.isAi80; ...
    c.isStimUp & c.isD1; ...
    c.isStimDown & c.isA2A; ...
    c.isStimDown & c.isAi80; ...
    c.isStimDown & c.isD1; ...
    };
NAME = { ...
    'iSPN-excited';
    'dSPN-excited';
    'Drd1-excited';
    'iSPN-inhibited';
    'dSPN-inhibited';
    'Drd1-inhibited';
    };
COLOR = {'red', 'red', 'red', 'default', 'default', 'default'};

figure(DefaultAxesFontSize=12, Units='normalized', OuterPosition=[0, 0, 0.55, 0.35])

latencies = cell(6, 1);
for i = 1:6
    subplot(2, 3, i)
    latencies{i} = [isi(SEL{i}).onsetLatency];
    histogram(latencies{i}, 1:4:100, DisplayName=sprintf('N=%g', length(latencies{i})), FaceColor=COLOR{i})
    title(sprintf('%s, median=%.0fms', NAME{i}, median(latencies{i})))
    legend()
    ylabel('Count')

    if i >= 4
        xlabel('SNr response latency (ms)')
    end
    if i == 1 || i == 4
        ylabel('Count')
    end
end

p1 = ranksum(latencies{1}, latencies{4}, tail='left')
p2 = ranksum(latencies{2}, latencies{5}, tail='right')
p3 = ranksum(latencies{3}, latencies{6}, tail='right')


%% TODO: Embed in EphysUnit as static method
function info = getAnimalInfo(eu, ai, field)
    i = find(strcmpi({ai.name}, eu.getAnimalName()));
    assert(length(i) == 1, eu.getAnimalName())
    info = ai(i).(field);
end

function ax = plotDoubleRaster(rd1, rd2, varargin)
    p = inputParser();
    p.addRequired('rd1', @isstruct)
    p.addRequired('rd2', @isstruct)
    p.addOptional('label1', '', @ischar)
    p.addOptional('label2', '', @ischar)
    p.addParameter('xlim', [-6, 1], @(x) isnumeric(x) && length(x) == 2 && x(2) > x(1));
    p.addParameter('iti', false, @islogical);
    p.addParameter('timeUnit', 's', @(x) ismember(x, {'s', 'ms'}))
    p.addParameter('maxTrials', Inf, @isnumeric)
    p.parse(rd1, rd2, varargin{:});
    r = p.Results;
    rd(1) = r.rd1;
    rd(2) = r.rd2;
    label{1} = r.label1;
    label{2} = r.label2;
    maxTrials = p.Results.maxTrials;


    f = figure(Units='normalized', OuterPosition=[0, 0, 0.5, 1], DefaultAxesFontSize=14);
    nTrials(1) = min(length(rd(1).duration), maxTrials);
    nTrials(2) = min(length(rd(2).duration), maxTrials);
    xmargin = 0.16;
    ymargin = 0.09;
    ax(1) = axes(f, Position=[xmargin, 2*ymargin+nTrials(2)/sum(nTrials)*(1-0.09*3), 0.7, nTrials(1)/sum(nTrials)*(1-ymargin*3)]);
    ax(2) = axes(f, Position=[xmargin, ymargin, 0.7, nTrials(2)/sum(nTrials)*(1-ymargin*3)]);

    for i = 1:2
        EphysUnit.plotRaster(ax(i), rd(i), xlim=r.xlim, iti=r.iti, ...
            timeUnit=p.Results.timeUnit, maxTrials=maxTrials);
        if ~isempty(label{i})
            title(ax(i), label{i})
        else
            switch lower(rd(i).trialType)
                case 'press'
                    name = 'Lever-press';
                case 'lick'
                    name = 'Lick';
                case {'stim', 'stimtrain', 'stimfirstpulse'}
                    name = 'Opto';
            end
            title(ax(i), name)
        end
%         suptitle(rd(1).name);
    end
end

function [isiStruct, ax] = getISILatencies(rd, varargin)
    p = inputParser();
    p.addRequired('rd', @isstruct);
    p.addParameter('xlim', [-100, 100], @isnumeric);
    p.addParameter('showPlots', false, @islogical);
    p.addParameter('savePlots', false, @islogical);
    p.addParameter('maxTrials', Inf, @isnumeric);
    p.addParameter('peakThreshold', 1, @isnumeric); % sd multiplier
    p.addParameter('posMultiplier', 0.5, @isnumeric); % sd multiplier
    p.addParameter('onsetThreshold', 0.1, @isnumeric); % sd multiplier
    p.addParameter('minProminence', 4, @isnumeric); % absolute unit, ms
    p.parse(rd, varargin{:});
    rd = p.Results.rd;
    xl = p.Results.xlim;
    onsetThreshold = p.Results.onsetThreshold;
    posMultiplier = p.Results.posMultiplier;
            

    isiStruct(length(rd), 1) = struct('t', [], 'isi', [], 'isi0', [], 'onsetLatency', [], 'peakLatency', [], 'peak', [], 'baseline', [], 'baselineSD', [], 'peakWidth', []);

    for iUnit = 1:length(rd)
        spikeTimes = rd(iUnit).t .* 1000;
        trialIndices = rd(iUnit).I;
    
        % Sort them
        [~, ISort] = sort(trialIndices*1000 + spikeTimes);
        trialIndices = trialIndices(ISort);
        spikeTimes = spikeTimes(ISort);

        if ~all(diff(trialIndices) >= 0)

            error('Well fuck you still.')
            
            isiStruct(iUnit).isi = [];
            isiStruct(iUnit).t = [];
            isiStruct(iUnit).onsetLatency = NaN;
            isiStruct(iUnit).peakLatency = NaN;
            isiStruct(iUnit).peak = NaN;
            isiStruct(iUnit).baseline = [];
            isiStruct(iUnit).baselineSD = [];
            isiStruct(iUnit).width = NaN;

            continue
        end

        [~, IFirstInTrial] = unique(trialIndices);
%         [~, ILastInTrial] = unique(flip(trialIndices));
%         ILastInTrial = length(trialIndices) + 1 - ILastInTrial;
    
        % Time to next spike
        isi = diff(spikeTimes);
        isi = [NaN, isi];
        isi(IFirstInTrial) = NaN;
%         isi(IFirstInTrial) = isiBaseline;
%         isiStd = mad(isi(spikeTimes < 0), 1, 'all') / 0.6745;

        % Always use 1ms bins
        t = xl(1):1:xl(2);
        isiContinuous = NaN(length(unique(rd(iUnit).I)), length(t));   
        iTrial = 0;
        isiStd = NaN(length(unique(trialIndices)), 1);
        for trialIndex = unique(trialIndices)
            iTrial = iTrial + 1;
            sel = trialIndices == trialIndex;
            if nnz(sel) >= 2 && length(unique(spikeTimes(sel))) == nnz(sel)
                isiContinuous(iTrial, :) = interp1(spikeTimes(sel), isi(sel), t, 'linear', NaN);
                isiStd(iTrial) = std(isiContinuous(iTrial, t<0 & t>xl(1)), 1, 'all', 'omitnan');
            end
        end

        x = mean(isiContinuous, 1, 'omitnan');
%         tPost = t(t>=0);
%         xPost = x(t>=0);
        isiBaseline = mean(x(t < 0 & t > xl(1)), 'omitnan');
        isiStd = mean(isiStd, 'omitnan');
        xStart = x(t == 0);
        x0 = x - xStart; % isi minus on baseline
        [posPeaks, posPeakTimes, posPeakWidths, posPeakProminences] = findpeaks(x0, t, MinPeakHeight=p.Results.peakThreshold*isiStd, MinPeakProminence=p.Results.minProminence);
        [negPeaks, negPeakTimes, negPeakWidths, negPeakProminences] = findpeaks(-x0, t, MinPeakHeight=posMultiplier * p.Results.peakThreshold*isiStd, MinPeakProminence=p.Results.minProminence);

        posPeaks = xStart + posPeaks;
        negPeaks = xStart - negPeaks;

        % Discard peaks before 0
        sel = posPeakTimes > 0;
        posPeaks = posPeaks(sel);
        posPeakTimes = posPeakTimes(sel);
        posPeakWidths = posPeakWidths(sel);
        posPeakProminences = posPeakProminences(sel);

        sel = negPeakTimes > 0;
        negPeaks = negPeaks(sel);
        negPeakTimes = negPeakTimes(sel);
        negPeakWidths = negPeakWidths(sel);
        negPeakProminences = negPeakProminences(sel);

        % Merge positive and negative peaks
        peaks = [posPeaks, negPeaks];
        peakSigns = [ones(size(posPeaks)), -1*ones(size(negPeaks))];

        % Skip the rest if no peaks detected
        peakTimes = [posPeakTimes, negPeakTimes];
        peakWidths = [posPeakWidths, negPeakWidths];
        peakProminences = [posPeakProminences, negPeakProminences];
    
        [peakTimes, ISort] = sort(peakTimes);
        ISort = ISort(peakTimes > 0);
        peakTimes = peakTimes(peakTimes > 0);
        peaks = peaks(ISort);
        peakSigns = peakSigns(ISort);
        peakWidths = peakWidths(ISort);
        peakProminences = peakProminences(ISort);
    
        if (isempty(peakTimes))
            isiStruct(iUnit).isi = x;
            isiStruct(iUnit).t = t;
            isiStruct(iUnit).isi0 = xStart;
            isiStruct(iUnit).onsetLatency = NaN;
            isiStruct(iUnit).peakLatency = NaN;
            isiStruct(iUnit).peak = NaN;
            isiStruct(iUnit).baseline = isiBaseline;
            isiStruct(iUnit).baselineSD = isiStd;
            isiStruct(iUnit).width = NaN;

            if p.Results.showPlots
                fig = figure(Position=[300 0 300, 600], DefaultAxesFontSize=12);
                ax = subplot(2, 1, 1);
                try
                    EphysUnit.plotRaster(ax, rd(iUnit), xlim=xl, ...
                        timeUnit='ms', maxTrials=p.Results.maxTrials);
                catch
                    disp()
                end
            
                ax = subplot(2, 1, 2);
                hold on
                h = gobjects(3, 1);
                h(1) = plot(ax, t, x, 'k', LineWidth=2, DisplayName='ISI');
                xlim(xl)
                h(2) = patch(ax, [xl, flip(xl)], [isiBaseline + isiStd, isiBaseline + isiStd, isiBaseline - isiStd, isiBaseline - isiStd], 'k', FaceAlpha=0.2);
                h(3) = plot(ax, xl, [isiBaseline, isiBaseline], 'k:', LineWidth=2, DisplayName='baseline');
                plot(ax, [0, 0], ax.YLim, 'k:')
                title('ISI (Time to next spike)')
                ylabel('ISI (ms)')
                yl = ax.YLim;
                plot(ax, [0, 0], ax.YLim, 'k:')
                ax.YLim = yl;
                xlabel('Time from opto on (ms)')
                legend(ax, h, Location='best')
        
                if p.Results.savePlots
                    if ~isfolder(sprintf('C:\\SERVER\\Figures\\Single Units\\Raster_%s_ISILatency', rd(iUnit).trialType))
                        mkdir(sprintf('C:\\SERVER\\Figures\\Single Units\\Raster_%s_ISILatency', rd(iUnit).trialType));
                    end
                    print(fig, sprintf('C:\\SERVER\\Figures\\Single Units\\Raster_%s_ISILatency\\%s (%s)', rd(iUnit).trialType, rd(iUnit).name, rd(iUnit).trialType), '-dpng')
                    close(fig)
                end     
            end
        else
            % Traceback from first isi peak to isi ramp onset 
            if peakSigns(1) > 0
                onsetSampleIndex = find(flip(x0 >= onsetThreshold * posMultiplier * isiStd & t >= 0), 1, 'last');
            else
                onsetSampleIndex = find(flip(x0 <= -onsetThreshold * isiStd & t >= 0), 1, 'last');
            end
            onsetSampleIndex = length(x) + 1 - onsetSampleIndex;
            tOnset = t(onsetSampleIndex);
            xOnset = x(onsetSampleIndex);
    
            % Store data
            isiStruct(iUnit).isi = x;
            isiStruct(iUnit).t = t;
            isiStruct(iUnit).isi0 = xStart;
            isiStruct(iUnit).onsetLatency = tOnset;
            isiStruct(iUnit).peakLatency = peakTimes(1);
            isiStruct(iUnit).peak = peaks(1);
            isiStruct(iUnit).baseline = isiBaseline;
            isiStruct(iUnit).baselineSD = isiStd;
            isiStruct(iUnit).width = peakWidths(1);
    
            if p.Results.showPlots
                fig = figure(Position=[300 0 300, 600], DefaultAxesFontSize=11);
                ax = subplot(2, 1, 1);
                EphysUnit.plotRaster(ax, rd(iUnit), xlim=xl, ...
                        timeUnit='ms', maxTrials=p.Results.maxTrials);
            
                ax = subplot(2, 1, 2);
                hold on
                clear h
                h(1) = plot(ax, t, x, 'k', LineWidth=2, DisplayName='ISI');
                xlim(xl)
                h(2) = patch(ax, [xl, flip(xl)], [isiBaseline + isiStd, isiBaseline + isiStd, isiBaseline - isiStd, isiBaseline - isiStd], 'k', FaceAlpha=0.2, DisplayName='std');
                h(3) = plot(ax, xl, [isiBaseline, isiBaseline], 'k:', LineWidth=2, DisplayName='baseline');
                if nnz(posPeakTimes>0) >= 1
                    h(end+1) = scatter(ax, posPeakTimes(posPeakTimes>0), posPeaks(posPeakTimes>0), 100, 'b', 'filled', DisplayName='peak');
                end
                if nnz(negPeakTimes>0) >= 1
                    h(end+1) = scatter(ax, negPeakTimes(negPeakTimes>0), negPeaks(negPeakTimes>0), 100, 'r', 'filled', DisplayName='peak');
                end
                plot(ax, [0, 0], ax.YLim, 'k:')
                title('ISI (Time to next spike)')
                ylabel('ISI (ms)')
                yl = ax.YLim;
                plot(ax, [0, 0], ax.YLim, 'k:')
                ax.YLim = yl;
                xlabel('Time from opto onset (ms)')
                text(peakTimes(1), peaks(1) - peakSigns(1) * 0.1 * diff(yl), arrayfun(@(x) sprintf('%.0fms', x), peakTimes(1), UniformOutput=false))
                h(end+1) = scatter(tOnset, xOnset, 100, 'yellow', 'filled', DisplayName='onset');
                text(tOnset, xOnset - peakSigns(1) * 0.1 * diff(yl), sprintf('%.0fms', tOnset));
                legend(ax, h, Location='best', FontSize=9)
        
                if p.Results.savePlots
                    if ~isfolder(sprintf('C:\\SERVER\\Figures\\Single Units\\Raster_%s_ISILatency', rd(iUnit).trialType))
                        mkdir(sprintf('C:\\SERVER\\Figures\\Single Units\\Raster_%s_ISILatency', rd(iUnit).trialType));
                    end
                    print(fig, sprintf('C:\\SERVER\\Figures\\Single Units\\Raster_%s_ISILatency\\%s (%s)', rd(iUnit).trialType, rd(iUnit).name, rd(iUnit).trialType), '-dpng')
                    close(fig)
                else
                    ax = subplot(2, 1, 1);
                    title(ax, 'Spike raster'); xlabel(ax, '');
                    ax = subplot(2, 1, 2);
                    title(ax, 'ISI (trial-average)')
                end
            end
        end
    end
end
