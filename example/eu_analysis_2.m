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
eu = EphysUnit.load('C:\SERVER\Units\Lite');
%% 1.3 AnimalInfo
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

p.minStimDuration = 1e-2;
p.maxStimDuration = 5e-2;
p.errStimDuration = 1e-3;
p.allowAltStimDuration = true;
p.etaWindowStim = [-0.2, 0.5];
p.metaWindowStim = [0, 0.1];



%% 2.2.1 Cull non-SNr units to save memory (only once)
msr = arrayfun(@(stats) stats.medianITI, [eu.SpikeRateStats]);
isSNr = msr >= p.minSpikeRate;
eu = eu(isSNr);
fprintf(1, 'Kept %g out of %g SNr units with spike rate >= %g.\n', nnz(isSNr), length(msr), p.minSpikeRate)
clearvars -except eu p ai

%% 2.2.2 Cull duplicate units by thresholding cross-correlation of binned spike counts
clearvars -except eu ai p
[uniqueExpNames, ~, ic] = unique({eu.ExpName});
nUniqueExps = length(uniqueExpNames);
r = cell(nUniqueExps, 1);
lags = cell(nUniqueExps, 1);
pairIndices = cell(nUniqueExps, 1);
for iExp = 1:nUniqueExps
    iEus = find(ic == iExp);
    tTic = tic();
    if length(iEus) <= 1
        fprintf(1, 'Session %g (of %g) has only %g unit, no calculation needed.\n', iExp, nUniqueExps, length(iEus));
        continue;
    end
    fprintf(1, 'Session %g (of %g), calculating %g pair-wise correlations...', iExp, nUniqueExps, nchoosek(length(iEus), 2));
    [r{iExp}, lags{iExp}, pairIndices{iExp}] = eu(iEus).xcorr('count', resolution=0.005, maxlag=10, normalize=true);

    pairIndices{iExp} = iEus(pairIndices{iExp});

    fprintf(1, 'Done (%.2f sec).\n', toc(tTic))
end

assert(size(unique(cat(1, lags{:}), 'rows'), 1) == 1)
lags = unique(cat(1, lags{:}), 'rows');
pairIndices = cat(1, pairIndices{:});
r = cat(1, r{:});
clearvars -except eu ai p lags pairIndices r 


%% 2.2.2.1 Plot a distribution of all pairwise correlations, for all experiments
close all
ax = axes(figure());
hold(ax, 'on')
r0 = r(:, lags==0);
h(1) = histogram(ax, r(:, lags==0), 0:0.01:1, DisplayName='R(0)', Normalization='probability');
yrange = ax.YLim;
h(2) = plot(ax, repmat(mean(r0), [2, 1]), yrange, DisplayName=sprintf('mean=%.2f', mean(r0)), LineWidth=2);
h(3) = plot(ax, repmat(median(r0), [2, 1]), yrange, DisplayName=sprintf('median=%.2f', median(r0)), LineWidth=2);
h(4) = plot(ax, repmat(mean(r0)+2*std(r0, 0), [2, 1]), yrange, DisplayName=sprintf('mean+2*std=%.2f', mean(r0)++2*std(r0, 0)), LineWidth=2);
h(5) = plot(ax, repmat(median(r0)+2*mad(r0, 1)/0.6745, [2, 1]), yrange, DisplayName=sprintf('median+2*mad//0.6745=%.2f', median(r0)++2*mad(r0, 1)/0.6745), LineWidth=2);
h(6) = plot(ax, repmat(prctile(r0, 95), [2, 1]), yrange, DisplayName=sprintf('95 percentile=%.2f', prctile(r0, 95)), LineWidth=2);
ax.YLim = yrange;
hold(ax, 'off')
legend(ax, h)

clear ax h yrange r0

%%
rd.press = eu.getRasterData('press');
rd.lick = eu.getRasterData('lick');


%% Choose a threshold and plot double rasters to compare
clear ax figname fig iPair dirname r0 nPairs i j rrddi rrddj ifig ME

rTheta = 0.40;
dirname = sprintf('C:\\SERVER\\Figures\\duplicate_detect_xcorr\\rTheta=%.2f', rTheta);

nPairs = size(pairIndices, 1);
r0 = r(:, lags==0);
for iPair = find(r0 > rTheta)'
    assert(numel(iPair) == 1)
    try
        i = pairIndices(iPair, 1);
        j = pairIndices(iPair, 2);
    
        if ~isempty(rd.press(i).t) && ~isempty(rd.press(j).t)
            rrddi = rd.press(i);
            rrddj = rd.press(j);
        elseif ~isempty(rd.lick(i).t) && ~isempty(rd.lick(j).t)
            rrddi = rd.lick(i);
            rrddj = rd.lick(j);
        else
            warning('No lick or press raster data found for units %s and %s', eu(i).getName(), eu(j).getName())
            continue
        end
    
        figname{1} = sprintf('r=%.2f (%g=%g+%g)', r0(iPair), iPair, i, j);
        figname{2} = sprintf('r=%.2f (%g=%g+%g) (i=%g)', r0(iPair), iPair, i, j, i);
        figname{3} = sprintf('r=%.2f (%g=%g+%g) (j=%g)', r0(iPair), iPair, i, j, j);

        % Plot double raster
        ax = plotDoubleRaster(rrddi, rrddj, ...
            sprintf('%g_%s', i, eu(i).getName('_')), ...
            sprintf('%g_%s', j, eu(j).getName('_')));
        fig(1) = ax(1).Parent;
        suptitle(figname{1})

        % Plot single figures
        ax = EphysUnit.plotRaster(rrddi);
        fig(2) = ax.Parent;
        ax = EphysUnit.plotRaster(rrddj);
        fig(3) = ax.Parent;

        if ~isfolder(dirname)
            mkdir(dirname);
        end
        
        for ifig = 1:3
            print(fig(ifig), sprintf('%s\\%s.jpg', dirname, figname{ifig}), '-djpeg');
        end
    catch ME
        warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
    end
    close all
end

clear ax figname fig iPair dirname r0 nPairs i j rrddi rrddj ifig ME

% Observations:
% With 5ms binning:
% rTheta > 0.8: definitely duplicate units
% rTheta in [0.6, 0.8] could be multiunits, or one channel (singleunit) is 
% a subset of another (multiunit), consider keeping the one with lower 
% firing rate? Might need to check waveforms and ISI first to verify this.

%% Filter out duplicate units with high correlation. 
% Keep the one with the lower unit index. This ensures when there are more
% than 2 duplicates, all duplicates are removed. 
rTheta = 0.7;
duplicatePairs = pairIndices(r(:, lags==0) > rTheta, :);

isDuplicate = false(length(eu), 1);
for iPair = 1:size(duplicatePairs, 1)
    isDuplicate(duplicatePairs(iPair, 2)) = true;
end
disp(find(isDuplicate))

eu = eu(~isDuplicate);

fprintf(1, 'Removed %g putative duplicate units with R(0) > %.2f.\n', nnz(isDuplicate), rTheta);

clearvars -except eu p ai



%% 2.3.1  Basic summaries
% Baseline (median) spike rates
msr = arrayfun(@(stats) stats.medianITI, [eu.SpikeRateStats]);

% Lick/Press responses
eta.press = eu.getETA('count', 'press', p.etaWindow, minTrialDuration=p.minTrialDuration, normalize=p.etaNorm);
meta.press = transpose(mean(eta.press.X(:, eta.press.t >= p.metaWindow(1) & eta.press.t <= p.metaWindow(2)), 2, 'omitnan'));
eta.lick = eu.getETA('count', 'lick', p.etaWindow, minTrialDuration=2, normalize=p.etaNorm);
meta.lick = transpose(mean(eta.lick.X(:, eta.lick.t >= p.metaWindow(1) & eta.lick.t <= p.metaWindow(2)), 2, 'omitnan'));

%% 2.3.2 Basic summaries (fast)
clearvars -except eu p msr eta meta ai

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

% animal info
cat.isWT = arrayfun(@(eu) strcmpi(getAnimalInfo(eu, ai, 'strain'), 'wt'), eu);
cat.isD1 = arrayfun(@(eu) strcmpi(getAnimalInfo(eu, ai, 'strain'), 'd1-cre'), eu);
cat.isA2A = arrayfun(@(eu) strcmpi(getAnimalInfo(eu, ai, 'strain'), 'a2a-cre'), eu);
cat.isAi80 = arrayfun(@(eu) strcmpi(getAnimalInfo(eu, ai, 'strain'), 'd1-cre;dlx-flp;ai80'), eu);
cat.isDAT = arrayfun(@(eu) strcmpi(getAnimalInfo(eu, ai, 'strain'), 'dat-cre'), eu);
cat.isAcute = ismember(eu.getAnimalName, {'daisy14', 'daisy15', 'daisy16', 'desmond23', 'desmond24', 'desmond25', 'desmond26', 'desmond27'});

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

%% 4.2.1 Stim ETA (PSTH) as heatmap
eta.stim = eu.getETA('count', 'stimfirstpulse', [-0.5, 0.5], ...
    resolution=1e-2, ...
    minTrialDuration=0.01, maxTrialDuration=0.1, ...
    findSingleTrialDuration='min', normalize=[-0.5, 0], includeInvalid=false);
meta.stim = transpose(mean(eta.stim.X(:, eta.stim.t >= 0 & eta.stim.t <= 0.1), 2, 'omitnan'));

%% 4.2.2 Plot heatmap
close all
cat.hasStim = ~all(isnan(eta.stim.X), 2)';

sel = cat.hasStim & cat.isA2A;
dmin = min(eta.stim.D(sel));
dmax = max(eta.stim.D(sel));
if dmin==dmax
    text = sprintf('iSPN stim (%g ms)', 1000*dmin);
else
    text = sprintf('iSPN stim (%g-%g ms)', 1000*dmin, 1000*dmax);
end
ax = EphysUnit.plotETA(eta.stim, sel, xlim=[-0.2, 0.1], clim=[-1.5, 1.5], sortWindow=[0, 0.1], signWindow=[0, 0.1], sortThreshold=0.3, negativeSortThreshold=0.2, event='opto onset'); title(ax, text)

sel = cat.hasStim & cat.isAi80;
dmin = min(eta.stim.D(sel));
dmax = max(eta.stim.D(sel));
if dmin==dmax
    text = sprintf('dSPN stim (%g ms)', 1000*dmin);
else
    text = sprintf('dSPN stim (%g-%g ms)', 1000*dmin, 1000*dmax);
end
ax = EphysUnit.plotETA(eta.stim, sel, xlim=[-0.2, 0.1], clim=[-1.5, 1.5], sortWindow=[0, 0.1], signWindow=[0, 0.1], sortThreshold=0.3, negativeSortThreshold=0.2, event='opto onset'); title(ax, text)

sel = cat.hasStim & cat.isD1;
dmin = min(eta.stim.D(sel));
dmax = max(eta.stim.D(sel));
if dmin==dmax
    text = sprintf('D1-Cre stim (%g ms)', 1000*dmin);
else
    text = sprintf('D1-Cre stim (%g-%g ms)', 1000*dmin, 1000*dmax);
end
ax = EphysUnit.plotETA(eta.stim, sel, xlim=[-0.2, 0.1], clim=[-1.5, 1.5], sortWindow=[0, 0.1], signWindow=[0, 0.1], sortThreshold=0.3, negativeSortThreshold=0.2, event='opto onset'); title(ax, text)

clear sel dmin dmax ax text

%% Find acute units for video analysis
unitNames = eu(cat.isAcute).getName()';
expNames = cellfun(@(x) strsplit(x, ' '), unitNames, UniformOutput=false);
expNames = unique(cellfun(@(x) x{1}, expNames, UniformOutput=false));

fprintf(1, 'Found %g sessions for video analysis:\n', length(expNames));
cellfun(@(x) fprintf(1, '\t%s\n', x), expNames);


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

function info = getAnimalInfo(eu, ai, field)
    i = find(strcmpi({ai.name}, eu.getAnimalName()));
    assert(length(i) == 1, eu.getAnimalName())
    info = ai(i).(field);
end

