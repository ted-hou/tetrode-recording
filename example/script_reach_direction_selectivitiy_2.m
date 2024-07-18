%% (!!!ONLY DO ONCE PER DATASET!!!) Cull multi-units, duplicate-units, DA untis, and write the remainder to server
% clear p
% p.maxISI = 1.5e-3;
% p.maxFraction = 0.05;
% p.minSpikeRate = 15;
% 
% %
% eu = EphysUnit.load('C:\SERVER\Units\acute_3cam_reach_direction_2tgts');
% euAll = eu;
% % Multiunit detection by ISI.
% clear cAll
% [eu, cAll.isMultiUnit] = eu.removeMultiUnits(maxISI=p.maxISI, maxFraction=p.maxFraction, cullZeros=true);
% 
% % Recalculate SpikeRateStats after multiunit culling
% for iEu = 1:length(eu)
%     resolution_sc = 0.1;
%     resolution_sr = 1e-3;
%     
%     [sc, tsc] = eu(iEu).getSpikeCounts(resolution_sc);
%     [sr, tsr, kernel] = eu(iEu).getSpikeRates('gaussian', 0.1, resolution_sr, 'kernelWidth', 1);
%     
%     % Calculate 
%     [~, ~, scInTrial] = eu(iEu).cullITIData(tsc, sc, 'all'); % This is sus but we don't use ITI stats for this exp since there's no ITI anyway.
%     [~, ~, srInTrial] = eu(iEu).cullITIData(tsr, sr, 'all'); % This is sus but we don't use ITI stats for this exp since there's no ITI anyway.
%     
%     eu(iEu).SpikeCountStats = struct('median', median(double(sc)), 'mad', mad(double(sc), 1), 'medianITI', median(double(sc(~scInTrial))), 'madITI', mad(double(sc(~scInTrial)), 1), 'resolution', resolution_sc);
%     eu(iEu).SpikeRateStats = struct('median', median(sr), 'mad', mad(sr, 1), 'medianITI', median(sr(~srInTrial)), 'madITI', mad(sr(~srInTrial), 1), 'resolution', resolution_sr);
% end
% clear iEu resolution_sc resolution_sr sc tsc sr tsr kernel scInTrial
% 
% % 2.2.1 Cull non-SNr units to save memory (only once)
% msr = arrayfun(@(stats) stats.medianITI, [eu.SpikeRateStats]);
% figure(), histogram(msr, 0:5:max(msr)+5);
% isSNr = msr >= p.minSpikeRate;
% eu = eu(isSNr);
% fprintf(1, 'Kept %g out of %g SNr units with spike rate >= %g.\n', nnz(isSNr), length(msr), p.minSpikeRate)
% 
% % 2.2.2 Cull duplicate units by thresholding cross-correlation of binned spike counts
% [uniqueExpNames, ~, ic] = unique({eu.ExpName});
% nUniqueExps = length(uniqueExpNames);
% r = cell(nUniqueExps, 1);
% lags = cell(nUniqueExps, 1);
% pairIndices = cell(nUniqueExps, 1);
% for iExp = 1:nUniqueExps
%     iEus = find(ic == iExp);
%     tTic = tic();
%     if length(iEus) <= 1
%         fprintf(1, 'Session %g (of %g) has only %g unit, no calculation needed.\n', iExp, nUniqueExps, length(iEus));
%         continue;
%     end
%     fprintf(1, 'Session %g (of %g), calculating %g pair-wise correlations...', iExp, nUniqueExps, nchoosek(length(iEus), 2));
%     [r{iExp}, lags{iExp}, pairIndices{iExp}] = eu(iEus).xcorr('count', resolution=0.005, maxlag=10, normalize=true);
% 
%     pairIndices{iExp} = iEus(pairIndices{iExp});
% 
%     fprintf(1, 'Done (%.2f sec).\n', toc(tTic))
% end
% 
% assert(size(unique(cat(1, lags{:}), 'rows'), 1) == 1)
% lags = unique(cat(1, lags{:}), 'rows');
% pairIndices = cat(1, pairIndices{:});
% r = cat(1, r{:});
% 
% % 2.2.2.1 Plot a distribution of all pairwise correlations, for all experiments
% close all
% ax = axes(figure());
% hold(ax, 'on')
% r0 = r(:, lags==0);
% h(1) = histogram(ax, r(:, lags==0), 0:0.01:1, DisplayName='R(0)', Normalization='probability');
% yrange = ax.YLim;
% h(2) = plot(ax, repmat(mean(r0), [2, 1]), yrange, DisplayName=sprintf('mean=%.2f', mean(r0)), LineWidth=2);
% h(3) = plot(ax, repmat(median(r0), [2, 1]), yrange, DisplayName=sprintf('median=%.2f', median(r0)), LineWidth=2);
% h(4) = plot(ax, repmat(mean(r0)+2*std(r0, 0), [2, 1]), yrange, DisplayName=sprintf('mean+2*std=%.2f', mean(r0)++2*std(r0, 0)), LineWidth=2);
% h(5) = plot(ax, repmat(median(r0)+2*mad(r0, 1)/0.6745, [2, 1]), yrange, DisplayName=sprintf('median+2*mad//0.6745=%.2f', median(r0)++2*mad(r0, 1)/0.6745), LineWidth=2);
% h(6) = plot(ax, repmat(prctile(r0, 95), [2, 1]), yrange, DisplayName=sprintf('95 percentile=%.2f', prctile(r0, 95)), LineWidth=2);
% ax.YLim = yrange;
% hold(ax, 'off')
% legend(ax, h)
% 
% clear ax h yrange r0
% 
% 
% % Filter out duplicate units with high correlation. 
% % Keep the one with the lower unit index. This ensures when there are more
% % than 2 duplicates, all duplicates are removed. 
% 
% rTheta = 0.6;
% duplicatePairs = pairIndices(r(:, lags==0) > rTheta, :);
% 
% isDuplicate = false(length(euWithDup), 1);
% for iPair = 1:size(duplicatePairs, 1)
%     isDuplicate(duplicatePairs(iPair, 2)) = true;
% end
% disp(find(isDuplicate))
% 
% eu = euWithDup(~isDuplicate);
% 
% fprintf(1, 'Removed %g putative duplicate units with R(0) > %.2f.\n', nnz(isDuplicate), rTheta);
% 
% 
% % Save good units post-cull
% eu.save('C:\SERVER\Units\acute_3cam_reach_direction_2tgts\SingleUnits_NonDuplicate')

%% Plot distribution of spike rates (pre- vs. post-cull)
figure, histogram(msrAll, 0:10:300, Normalization='pdf', DisplayName='all units'), hold on, histogram(msr, 0:10:300, Normalization='pdf', DisplayName='single units, >15sp/s, non-duplicate')
xlabel('spike rate (sp/s)')
ylabel('prob')
legend

%%
for iEu = 1:length(eu)
    trials = eu(iEu).makeTrials('press_spontaneous2');
    goodTrials = trials(trials.duration() > 4);
    eu(iEu).Trials.PressSpontaneous = goodTrials;
end

%% Calculate ETA/PETH
eta = eu.getETA('count', 'press_spontaneous', window=[-4, 0], includeInvalid=true, normalize=[-4, -2]);
plot(eta.t, eta.X);

%% Calculate separate ETAs for medial vs. lateral reach
for iEu = 1:length(eu)
    pressTrials = eu(iEu).Trials.PressSpontaneous;
    pressTimes = [pressTrials.Stop];

    motorTrials = Trial(eu(1).EventTimes.Mot1LoOn, eu(1).EventTimes.Mot1LoOff, advancedValidation=false);
    isLateral = motorTrials.inTrial(pressTimes);

    fprintf(1, '%i med, %i lat\n', nnz(~isLateral), nnz(isLateral))

    eu(iEu).Trials.PressSpontaneousMedial = pressTrials(~isLateral);
    eu(iEu).Trials.PressSpontaneousLateral = pressTrials(isLateral);
end

etaMedial = eu.getETA('count', 'press_spontaneous_medial', window=[-4, 0], includeInvalid=true, normalize=[-4, -2]);
etaLateral = eu.getETA('count', 'press_spontaneous_lateral', window=[-4, 0], includeInvalid=true, normalize=[-4, -2]);

%% Plot ETA
% [~, I] = EphysUnit.plotETA(etaLateral, clim=[-1, 1], signWindow=[-0.5, -0.2], sortWindow=[-3, 0], sortThreshold=0.5, negativeSortThreshold=0.25);

fig = figure(Units='inches', Position=[0 0 8 8]);
tl = tiledlayout(fig, 1, 3);
ax = gobjects(1, 3);
ax(1) = nexttile(tl);
ax(2) = nexttile(tl);
ax(3) = nexttile(tl);
[~, I] = EphysUnit.plotETA(ax(1), eta, clim=[-1, 1], signWindow=[-0.5, -0.2], sortWindow=[-3, 0], sortThreshold=0.5, negativeSortThreshold=0.25, hideColorbar=true);

EphysUnit.plotETA(ax(2), etaMedial, clim=[-1, 1], order=I, hideColorbar=true);

EphysUnit.plotETA(ax(3), etaLateral, clim=[-1, 1], order=I);
ylabel(ax(2:3), '')
xlabel(ax(1:3), '')
title(ax(1), 'all reach trials')
title(ax(2), 'medial reach trials')
title(ax(3), 'lateral reach trials')
colorbar(ax(3), Position=[0.925683549861173,0.108072916666667,0.014453125,0.815104166666667])

xlabel(tl, 'Time from bar-contact (s)')

%%
t = etaLateral.t;
metaLat = mean(etaLateral.X(:, t < -0.2 & t > -0.5), 2);
metaMed = mean(etaMedial.X(:, t < -0.2 & t > -0.5), 2);
ax = axes(figure); 
scatter(ax, metaLat, metaMed, 25, 'k', 'o')
axis(ax, 'equal')
xl = ax.XLim;
yl = ax.YLim;
hold(ax, 'on')
plot(ax, xl, [0, 0], 'k')
plot(ax, [0, 0], yl, 'k')
plot(ax, xl, xl, 'k--')
xlabel(ax, 'Lateral [-0.2, -0.5]')
ylabel(ax, 'Medial [-0.5, -0.2]')
xlim(ax, xl)
ylim(ax, yl)