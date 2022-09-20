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
eu = EphysUnit.load('C:\SERVER\Units', waveforms=false, spikecounts=true, spikerates=true);

%% 2.1 Parameters
clear p
p.minSpikeRate = 15;
p.minTrialDuration = 2;
p.minNumTrials = 30;
p.etaNorm = [-4, -2];
p.etaWindow = [-4, 0];
p.metaWindow = [-0.2, 0];
p.posRespThreshold = 1;
p.negRespThreshold = -0.5;
p.binnedTrialEdges = 2:2:10;

%% 2.2 Cull non-SNr units to save memory (only once)
msr = arrayfun(@(stats) stats.medianITI, [eu.SpikeRateStats]);
isSNr = msr >= p.minSpikeRate;
eu = eu(isSNr);
fprintf(1, 'Kept %g out of %g SNr units with spike rate >= %g.\n', nnz(isSNr), length(msr), p.minSpikeRate)
clearvars -except eu p

%% 2.3.1  Basic summaries
% Baseline (median) spike rates
msr = arrayfun(@(stats) stats.medianITI, [eu.SpikeRateStats]);

% Lick/Press responses
eta.press = eu.getETA('count', 'press', p.etaWindow, minTrialDuration=p.minTrialDuration, normalize=p.etaNorm);
meta.press = transpose(mean(eta.press.X(:, eta.press.t >= p.metaWindow(1) & eta.press.t <= p.metaWindow(2)), 2, 'omitnan'));
eta.lick = eu.getETA('count', 'lick', p.etaWindow, minTrialDuration=2, normalize=p.etaNorm);
meta.lick = transpose(mean(eta.lick.X(:, eta.lick.t >= p.metaWindow(1) & eta.lick.t <= p.metaWindow(2)), 2, 'omitnan'));

%% 2.3.2 Basic summaries (fast)
clearvars -except eu p msr eta meta

% hasPress/hasLick
cat.hasPress = arrayfun(@(e) nnz(e.getTrials('press').duration() >= p.minTrialDuration) >= p.minNumTrials, eu);
cat.hasLick = arrayfun(@(e) nnz(e.getTrials('lick').duration() >= p.minTrialDuration) >= p.minNumTrials, eu);

% press/lick x Up/Down
cat.isPressUp =         cat.hasPress & meta.press >= p.posRespThreshold;
cat.isPressDown =       cat.hasPress & meta.press <= p.negRespThreshold;
cat.isPressResponsive = cat.isPressUp | cat.isPressDown;
cat.isLickUp =          cat.hasLick & meta.lick >= p.posRespThreshold;
cat.isLickDown =        cat.hasLick & meta.lick <= p.negRespThreshold;
cat.isLickResponsive =  cat.isLickUp | cat.isLickDown;

fprintf(1, ['%g total SNr units (baseline spike rate > %g):\n' ...
    '\t%g with press trials;\n' ...
    '\t%g with lick trials;\n' ...
    '\t%g with both (at least %g trials).\n'], ...
    length(eu), p.minSpikeRate, nnz(cat.hasPress), nnz(cat.hasLick), nnz(cat.hasPress & cat.hasLick), p.minNumTrials)

fprintf(1, ['%g units with at least %g press trials:\n' ...
    '\t%g (%.0f%%) are excited (meta>=%g);\n' ...
    '\t%g (%.0f%%) are inhibited (meta<=%g).\n'], ...
    nnz(cat.hasPress), p.minNumTrials, ...
    nnz(cat.isPressUp), 100*nnz(cat.isPressUp)/nnz(cat.isPressResponsive), p.posRespThreshold, ...
    nnz(cat.isPressDown), 100*nnz(cat.isPressDown)/nnz(cat.isPressResponsive), p.negRespThreshold);

fprintf(1, ['%g units with at least %g lick trials:\n' ...
    '\t%g (%.0f%%) are excited (meta>=%g);\n' ...
    '\t%g (%.0f%%) are inhibited (meta<=%g).\n'], ...
    nnz(cat.hasLick), p.minNumTrials, ...
    nnz(cat.isLickUp), 100*nnz(cat.isLickUp)/nnz(cat.isLickResponsive), p.posRespThreshold, ...
    nnz(cat.isLickDown), 100*nnz(cat.isLickDown)/nnz(cat.isLickResponsive), p.negRespThreshold);

nTotal = nnz(cat.isPressResponsive & cat.isLickResponsive);
fprintf(1, ['%g units with press AND lick trials:\n' ...
    '\t%g (%.0f%%) are press-excited AND lick-excited;\n' ...
    '\t%g (%.0f%%) are press-inhibited AND lick-inhibited;\n' ...
    '\t%g (%.0f%%) are press-excited AND lick-inhibited;\n' ...
    '\t%g (%.0f%%) are press-inhibited AND lick-excited;\n'], ...
    nnz(cat.hasPress & cat.hasLick), ...
    nnz(cat.isPressUp & cat.isLickUp), 100*nnz(cat.isPressUp & cat.isLickUp)/nTotal, ...
    nnz(cat.isPressDown & cat.isLickDown), 100*nnz(cat.isPressDown & cat.isLickDown)/nTotal, ...
    nnz(cat.isPressUp & cat.isLickDown), 100*nnz(cat.isPressUp & cat.isLickDown)/nTotal, ...
    nnz(cat.isPressDown & cat.isLickUp), 100*nnz(cat.isPressDown & cat.isLickUp)/nTotal)   
clear nTotal

%% 3.1 Plot basics
% Distribution of baseline, press response, lick response
% Baseline vs press, lick, press vs lick
close all
sz = 10;

fig = figure(Units='pixels', OuterPosition=[0, 0, 1920, 1080], DefaultAxesFontSize=14);
ax = subplot(2, 3, 1);
histogram(ax, msr, Normalization='probability');
xlabel(ax, 'Spike rate (sp/s)'), ylabel(ax, 'Probability'), title(ax, 'Baseline spike rate (median, ITI)')

ax = subplot(2, 3, 2);
histogram(ax, meta.press(cat.hasPress))
xlabel(ax, '\Deltaz (press, a.u.)'), ylabel(ax, 'Probability'), title(ax, sprintf('Normalized press response (%gs-%gs, %g units)', p.metaWindow(1), p.metaWindow(2), nnz(cat.hasPress)))

ax = subplot(2, 3, 3);
histogram(ax, meta.lick(cat.hasLick))
xlabel(ax, '\Deltaz (lick, a.u.)'), ylabel(ax, 'Probability'), title(ax, sprintf('Normalized lick response (%gs-%gs, %g units)', p.metaWindow(1), p.metaWindow(2), nnz(cat.hasLick)))

ax = subplot(2, 3, 4);
hold(ax, 'on')
isResp = cat.isPressResponsive;
isNonResp = cat.hasPress & ~cat.isPressResponsive;
h(1) = scatter(ax, msr(isResp), meta.press(isResp), sz, 'black', 'filled', DisplayName=sprintf('%g press-responsive units', nnz(isResp)));
h(2) = scatter(ax, msr(isNonResp), meta.press(isNonResp), sz, [0.5, 0.5, 0.5], 'filled', DisplayName=sprintf('%g non-responsive units', nnz(isNonResp)));
hold(ax, 'off')
xlabel(ax, 'Baseline spike rate (sp/s)'), ylabel(ax, '\Deltaz (press, a.u.)'), title(ax, 'Press response vs baseline')
legend(ax, h)

ax = subplot(2, 3, 5);
hold(ax, 'on')
isResp = cat.isLickResponsive;
isNonResp = cat.hasLick & ~cat.isLickResponsive;
h(1) = scatter(ax, msr(isResp), meta.lick(isResp), sz, 'black', 'filled', DisplayName=sprintf('%g lick-responsive units', nnz(isResp)));
h(2) = scatter(ax, msr(isNonResp), meta.lick(isNonResp), sz, [0.5, 0.5, 0.5], 'filled', DisplayName=sprintf('%g non-responsive units', nnz(isNonResp)));
hold(ax, 'off')
xlabel(ax, 'Baseline spike rate (sp/s)'), ylabel(ax, '\Deltaz (lick, a.u.)'), title(ax, 'Lick response vs baseline')
legend(ax, h)

ax = subplot(2, 3, 6);
hold(ax, 'on')
isResp = cat.isPressResponsive & cat.isLickResponsive;
isNonResp = cat.hasLick & cat.hasPress & (~cat.isLickResponsive | ~cat.isPressResponsive);
h(1) = scatter(ax, meta.lick(isResp), meta.press(isResp), sz, 'black', 'filled', DisplayName=sprintf('%g lick-responsive units', nnz(isResp)));
h(2) = scatter(ax, meta.lick(isNonResp), meta.press(isNonResp), sz, [0.5, 0.5, 0.5], 'filled', DisplayName=sprintf('%g non-responsive units', nnz(isNonResp)));
xl = ax.XLim; yl = ax.YLim;
ax.XLimMode = 'manual'; ax.YLimMode = 'manual'; 
plot(ax, xl, [0, 0], 'k:')
plot(ax, [0, 0], yl, 'k:')
hold(ax, 'off')
xlabel(ax, '\Deltaz (lick, a.u.)'), ylabel(ax, '\Deltaz (press, a.u.)'), title(ax, 'Press vs lick response')
legend(ax, h)

clear fig ax h isResp isNonResp xl yl sz

%% 3.2 Plot binned averaged (BTA)
% Calculate BTA
[bta.pressUp.X, bta.pressUp.T, bta.pressUp.N, bta.pressUp.S, bta.pressUp.B] = eu(cat.isPressUp).getBinnedTrialAverage('rate', p.binnedTrialEdges, 'press', 'window', [-10, 1], 'normalize', true);
[bta.pressDown.X, bta.pressDown.T, bta.pressDown.N, bta.pressDown.S, bta.pressDown.B] = eu(cat.isPressDown).getBinnedTrialAverage('rate', p.binnedTrialEdges, 'press', 'window', [-10, 1], 'normalize', true);
[bta.pressUpRaw.X, bta.pressUpRaw.T, bta.pressUpRaw.N, bta.pressUpRaw.S, bta.pressUpRaw.B] = eu(cat.isPressUp).getBinnedTrialAverage('rate', p.binnedTrialEdges, 'press', 'window', [-10, 1], 'normalize', false);
[bta.pressDownRaw.X, bta.pressDownRaw.T, bta.pressDownRaw.N, bta.pressDownRaw.S, bta.pressDownRaw.B] = eu(cat.isPressDown).getBinnedTrialAverage('rate', p.binnedTrialEdges, 'press', 'window', [-10, 1], 'normalize', false);

[bta.lickUp.X, bta.lickUp.T, bta.lickUp.N, bta.lickUp.S, bta.lickUp.B] = eu(cat.isLickUp).getBinnedTrialAverage('rate', p.binnedTrialEdges, 'lick', 'window', [-10, 1], 'normalize', true);
[bta.lickDown.X, bta.lickDown.T, bta.lickDown.N, bta.lickDown.S, bta.lickDown.B] = eu(cat.isLickDown).getBinnedTrialAverage('rate', p.binnedTrialEdges, 'lick', 'window', [-10, 1], 'normalize', true);
[bta.lickUpRaw.X, bta.lickUpRaw.T, bta.lickUpRaw.N, bta.lickUpRaw.S, bta.lickUpRaw.B] = eu(cat.isLickUp).getBinnedTrialAverage('rate', p.binnedTrialEdges, 'lick', 'window', [-10, 1], 'normalize', false);
[bta.lickDownRaw.X, bta.lickDownRaw.T, bta.lickDownRaw.N, bta.lickDownRaw.S, bta.lickDownRaw.B] = eu(cat.isLickDown).getBinnedTrialAverage('rate', p.binnedTrialEdges, 'lick', 'window', [-10, 1], 'normalize', false);

% Plot
close all
figure();
axUp = subplot(2, 2, 1);
axUpRaw = subplot(2, 2, 3);
axDown = subplot(2, 2, 2);
axDownRaw = subplot(2, 2, 4);
ax = [axUp, axUpRaw, axDown, axDownRaw];
EphysUnit.plotBinnedTrialAverage(ax(1), bta.pressUp, [1,12]);
EphysUnit.plotBinnedTrialAverage(ax(2), bta.pressUpRaw, [1,12]);
EphysUnit.plotBinnedTrialAverage(ax(3), bta.pressDown, [1,12]);
EphysUnit.plotBinnedTrialAverage(ax(4), bta.pressDownRaw, [1,12]);
title([axUp, axUpRaw], sprintf('Press-excited (%i units)', nnz(cat.isPressUp)))
title([axDown, axDownRaw], sprintf('Press-inhibited (%i units)', nnz(cat.isPressDown)))
xlabel(ax, 'Time from cue (s)')
ylabel([axUp, axDown], 'Normalized spike rate (modified z-score)')
ylabel([axUpRaw, axDownRaw], 'Spike rate (sp/s)')
clear axUp axDown axUpRaw axDownRaw ax

figure();
axUp = subplot(2, 2, 1);
axUpRaw = subplot(2, 2, 3);
axDown = subplot(2, 2, 2);
axDownRaw = subplot(2, 2, 4);
ax = [axUp, axUpRaw, axDown, axDownRaw];
EphysUnit.plotBinnedTrialAverage(ax(1), bta.lickUp, [1,12]);
EphysUnit.plotBinnedTrialAverage(ax(2), bta.lickUpRaw, [1,12]);
EphysUnit.plotBinnedTrialAverage(ax(3), bta.lickDown, [1,12]);
EphysUnit.plotBinnedTrialAverage(ax(4), bta.lickDownRaw, [1,12]);
title([axUp, axUpRaw], sprintf('Lick-excited (%i units)', nnz(cat.isLickUp)))
title([axDown, axDownRaw], sprintf('Lick-inhibited (%i units)', nnz(cat.isLickDown)))
xlabel(ax, 'Time from cue (s)')
ylabel([axUp, axDown], 'Normalized spike rate (modified z-score, a.u.)')
ylabel([axUpRaw, axDownRaw], 'Spike rate (sp/s)')
clear axUp axDown axUpRaw axDownRaw ax

%% 3.3 Plot ETA Heatmap
close all
EphysUnit.plotETA(eta.press, cat.hasPress, xlim=[-4,0], clim=[-2, 2], sortWindow=[-3, 0], signWindow=[-0.2, 0], sortThreshold=0.6, negativeSortThreshold=0.3); title('Press ETA')
EphysUnit.plotETA(eta.lick, cat.hasLick, xlim=[-4,0], clim=[-2, 2], sortWindow=[-3, 0], signWindow=[-0.2, 0], sortThreshold=0.6, negativeSortThreshold=0.3); title('Lick ETA')
[~, ~, ~, latency] = EphysUnit.plotDoubleETA(eta.press, eta.lick, cat.hasPress & cat.hasLick, 'Lever-press', 'Lick', xlim=[-4,0], clim=[-2, 2], sortWindow=[-3, 0], signWindow=[-0.2, 0], sortThreshold=0.6, negativeSortThreshold=0.3);
EphysUnit.plotDoubleETA(eta.lick, eta.press, cat.hasPress & cat.hasLick, 'Lick', 'Lever-press', xlim=[-4,0], clim=[-2, 2], sortWindow=[-3, 0], signWindow=[-0.2, 0], sortThreshold=0.6, negativeSortThreshold=0.3);

%% 3.4 Single units raster and BTAs (plot and save to disk)
plotRasterForSingleUnits(eu(cat.isLickUp), 'lick', 'Raster_LickUp')
plotRasterForSingleUnits(eu(cat.isLickDown), 'lick', 'Raster_LickDown')
plotRasterForSingleUnits(eu(cat.isPressUp), 'press', 'Raster_PressUp')
plotRasterForSingleUnits(eu(cat.isPressDown), 'press', 'Raster_PressDown')

%% 3.5 Single unit double rasters
plotDoubleRasterForSingleUnits(eu(cat.isPressDown & cat.isLickUp), 'DoubleRaster_PressDownLickUp')
plotDoubleRasterForSingleUnits(eu(cat.isPressUp & cat.isLickDown), 'DoubleRaster_PressUpLickDown')
plotDoubleRasterForSingleUnits(eu(cat.isPressDown & cat.isLickDown), 'DoubleRaster_PressDownLickDown')
plotDoubleRasterForSingleUnits(eu(cat.isPressUp & cat.isLickUp), 'DoubleRaster_PressUpLickUp')

%% 4. Stim response
%
%% 4.1 Single unit raster for stim responses (first pulse in each train)
clear rdStim
cat.hasStim = arrayfun(@(eu) ~isempty(eu.getTrials('stim')), eu);
rdStim(cat.hasStim) = eu(cat.hasStim).getRasterData('stimfirstpulse', window=[-0.5, 0.5], sort=false);

%% 4.1.1 Filter out specific stim durations (try 10ms, then the smallest above 10ms)
p.minStimDuration = 1e-2;
p.errStimDuration = 1e-3;
p.allowAltStimDuration = true;

rdStimFiltered = rdStim;

for i = 1:length(rdStim)
    if isempty(rdStim(i).t)
        continue
    end

    roundedDuration = round(rdStim(i).duration ./ p.errStimDuration) * p.errStimDuration;
    isMinDur = roundedDuration == p.minStimDuration;
    if any(isMinDur)
        sel = ismember(rdStim(i).I, find(isMinDur));
        rdStimFiltered(i).t = rdStim(i).t(sel);
        newI = rdStim(i).I(sel); newI = changem(newI, 1:length(unique(newI)), unique(newI));
        rdStimFiltered(i).I = newI;
        rdStimFiltered(i).duration = roundedDuration(isMinDur);
        rdStimFiltered(i).iti = rdStim(i).iti(isMinDur);
        fprintf('%g: found %g trials with length %g.\n', i, nnz(isMinDur), p.minStimDuration);
    else
        isAboveMinDur = roundedDuration > p.minStimDuration;
        if p.allowAltStimDuration && any(isAboveMinDur)
            altMinDur = min(roundedDuration(isAboveMinDur));
            isAltMinDur = roundedDuration == altMinDur;
            sel = ismember(rdStim(i).I, find(isAltMinDur));
            rdStimFiltered(i).t = rdStim(i).t(sel);
            newI = rdStim(i).I(sel); newI = changem(newI, 1:length(unique(newI)), unique(newI));
            rdStimFiltered(i).I = newI;
            rdStimFiltered(i).duration = roundedDuration(isAltMinDur);
            rdStimFiltered(i).iti = rdStim(i).iti(isAltMinDur);
            fprintf('%g: could not find requested duration, found %g trials with length %g instead.\n', i, nnz(isAltMinDur), altMinDur);
        else
            fprintf('%g: could not find requested duration.', i)
            rdStimFiltered(i).t = [];
            rdStimFiltered(i).I = [];
            rdStimFiltered(i).duration = [];
            rdStimFiltered(i).iti = [];
        end
    end
end

cat.hasStim = arrayfun(@(rd) ~isempty(rd.t), rdStimFiltered);
fprintf(1, '%g units with requested stim duration.\n', nnz(cat.hasStim))

clear i roundedDuration isMinDur isAboveMinDur altMinDur isAltMinDur sel

%% 4.1.2 Save single unit opto rasters to disk
if ~isfolder('C:\SERVER\Figures\Single Units\Raster_Stim')
    mkdir('C:\SERVER\Figures\Single Units\Raster_Stim')
end
for rrdd = rdStimFiltered(hasStim)
    try
        nTrials = length(rrdd.duration);
        fig = figure(Units='normalized', Position=[0, 0, 0.4, max(0.2, min(nTrials*0.004, 0.8))]);
        ax = axes(fig);
        EphysUnit.plotRaster(ax, rrdd, xlim=[-0.2, unique(rrdd.duration) + min(rrdd.iti)]);
        print(fig, sprintf('C:\\SERVER\\Figures\\Single Units\\Raster_Stim\\%s', rrdd.name), '-dpng');
        close(fig)
    catch
            fprintf(1, 'Error while processing %s.\n', rrdd.name);
    end
end
clear rrdd ax fig nTrials

%% 4.2 Stim ETA (PSTH) as heatmap


%%
function plotRasterForSingleUnits(eu, moveType, category)
    if ~isfolder(sprintf('C:\\SERVER\\Figures\\Single Units\\%s', category))
        mkdir(sprintf('C:\\SERVER\\Figures\\Single Units\\%s', category))
    end
    for e = eu
        try
            rd = e.getRasterData(moveType, window=[0, 0], sort=true);
            ax = EphysUnit.plotRaster(rd, xlim=[-4, 0]);
            fig = ax.Parent;
            print(fig, sprintf('C:\\SERVER\\Figures\\Single Units\\%s\\%s', category, e.getName('_')), '-dpng');            
            close(fig)
        catch ME
            fprintf(1, 'Error while processing %s.\n', e.getName('_'));
        end
        close all
    end
end

function plotDoubleRasterForSingleUnits(eu, category)
    if ~isfolder(sprintf('C:\\SERVER\\Figures\\Single Units\\%s', category))
        mkdir(sprintf('C:\\SERVER\\Figures\\Single Units\\%s', category))
    end
    for e = eu
        try
            rdpress = e.getRasterData('press', window=[0, 2], sort=true);
            rdlick = e.getRasterData('lick', window=[0, 2], sort=true);
            ax = plotDoubleRaster(rdpress, rdlick, xlim=[-4, 2], iti=false);
            fig = ax(1).Parent;
            print(fig, sprintf('C:\\SERVER\\Figures\\Single Units\\%s\\%s', category, e.getName('_')), '-dpng');            
            close(fig)
        catch ME
            fprintf(1, 'Error while processing %s.\n', e.getName('_'));
        end
        close all
    end
end

%% TODO: Embed in EphysUnit as static method
function ax = plotDoubleRaster(rd1, rd2, varargin)
    p = inputParser();
    p.addRequired('rd1', @isstruct)
    p.addRequired('rd2', @isstruct)
    p.addOptional('label1', '', @ischar)
    p.addOptional('label2', '', @ischar)
    p.addParameter('xlim', [-6, 1], @(x) isnumeric(x) && length(x) == 2 && x(2) > x(1));
    p.addParameter('iti', false, @islogical);
    p.parse(rd1, rd2, varargin{:});
    r = p.Results;
    rd(1) = r.rd1;
    rd(2) = r.rd2;
    label{1} = r.label1;
    label{2} = r.label2;


    f = figure(Units='normalized', OuterPosition=[0, 0, 0.5, 1], DefaultAxesFontSize=14);
    nTrials(1) = length(rd(1).duration);
    nTrials(2) = length(rd(2).duration);
    ax(1) = axes(f, OuterPosition=[0, nTrials(2)/sum(nTrials), 1, nTrials(1)/sum(nTrials)]);
    ax(2) = axes(f, OuterPosition=[0, 0, 1, nTrials(2)/sum(nTrials)]);

    for i = 1:2
        EphysUnit.plotRaster(ax(i), rd(i), xlim=r.xlim, iti=r.iti);
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
        suptitle(rd(1).name);
    end
end