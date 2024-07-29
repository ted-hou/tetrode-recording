%% Load EphysUnit
eu = EphysUnit.load('\\research.files.med.harvard.edu\Neurobio\NEUROBIOLOGY SHARED\Assad Lab\Lingfeng\Data\Units\acute_3cam_reach_direction_2tgts\SingleUnits_NonDuplicate');

%% Pick spontaneous trials
for iEu = 1:length(eu)
    trials = eu(iEu).makeTrials('press_spontaneous2');
    goodTrials = trials(trials.duration() > 4);
    eu(iEu).Trials.PressSpontaneous = goodTrials;
end

% split trials by lateral/medial
for iEu = 1:length(eu)
    pressTrials = eu(iEu).Trials.PressSpontaneous;
    pressTimes = [pressTrials.Stop];

    motorTrials = Trial(eu(1).EventTimes.Mot1LoOn, eu(1).EventTimes.Mot1LoOff, advancedValidation=false);
    isLateral = motorTrials.inTrial(pressTimes);

    fprintf(1, '%i med, %i lat\n', nnz(~isLateral), nnz(isLateral))

    eu(iEu).Trials.PressSpontaneousMedial = pressTrials(~isLateral);
    eu(iEu).Trials.PressSpontaneousLateral = pressTrials(isLateral);
end

%% Make raster/PETH for two example units
% rd = eu.getRasterData('press_spontaneous', sort=true, alignTo='stop');
% %
% fig = figure();
% ax = axes(fig);
% 
% for iEu = 1:length(eu)
%     cla(ax)
%     EphysUnit.plotRaster(ax, rd(iEu), xlim=[-4, 0])
%     print(fig, sprintf('C:\\SERVER\\Figures\\Press_Spontaneous_2tgts\\%s.png', rd(iEu).name), '-dpng')
% end

%% Pick example unit, make raster
exampleUnitName = 'daisy25_20240628_Channel114_Unit1';
iEu = find(strcmp(eu.getName(), exampleUnitName));
rdMed = eu(iEu).getRasterData('press_spontaneous', sort=true, alignTo='stop', trials=eu(iEu).Trials.PressSpontaneousMedial);
rdLat = eu(iEu).getRasterData('press_spontaneous', sort=true, alignTo='stop', trials=eu(iEu).Trials.PressSpontaneousLateral);
% fig = figure();
% ax = gobjects(2, 1);
% ax(1) = subplot(2, 1, 1);
% ax(2) = subplot(2, 1, 2);
ax = EphysUnit.plotMultiRaster([rdMed, rdLat], label=["Medial reach", "Lateral reach"], xlim=[-4, 0], sz=3);
xlabel(ax(1), '')
xlabel(ax(2), 'Time to touch (s)')
delete(ax(1).Legend)
delete(ax(2).Legend)

xticks(ax(1), [])
ax(1).Parent.Units = 'inches';
ax(1).Parent.Position = [0.708333333333333,1.395833333333333,4.71875,5.78125];

%% Make PETH
ax = axes(figure(DefaultAxesFontSize=14));
etaMedExample = eu(iEu).getETA('count', 'press_spontaneous', window=[-4, 0], includeInvalid=true, normalize='none', trials=eu(iEu).Trials.PressSpontaneousMedial);
etaLatExample = eu(iEu).getETA('count', 'press_spontaneous', window=[-4, 0], includeInvalid=true, normalize='none', trials=eu(iEu).Trials.PressSpontaneousLateral);
hold(ax, 'on')
plot(ax, etaMedExample.t, etaMedExample.X./0.1, LineWidth=1.5, DisplayName='Medial');
plot(ax, etaLatExample.t, etaLatExample.X./0.1, LineWidth=1.5, DisplayName='Lateral');
legend(ax, Location='northwest')
xlabel(ax, 'Time to touch (s)')
ylabel(ax, 'Spike rate (sp/s)')

%% Make heatmap
eta = eu.getETA('count', 'press_spontaneous', window=[-4, 0], includeInvalid=true, normalize=[-4, -2]);
etaMedial = eu.getETA('count', 'press_spontaneous_medial', window=[-4, 0], includeInvalid=true, normalize=[-4, -2]);
etaLateral = eu.getETA('count', 'press_spontaneous_lateral', window=[-4, 0], includeInvalid=true, normalize=[-4, -2]);

%%
fig = figure(Units='inches', Position=[0 0 9 7.5], DefaultAxesFontSize=14);
tl = tiledlayout(fig, 1, 3);
ax = gobjects(1, 3);
ax(1) = nexttile(tl);
ax(2) = nexttile(tl);
ax(3) = nexttile(tl);
[~, I] = EphysUnit.plotETA(ax(1), eta, clim=[-1, 1], signWindow=[-0.5, -0.2], sortWindow=[-3, 0], sortThreshold=0.5, negativeSortThreshold=0.25, hideColorbar=true);

EphysUnit.plotETA(ax(2), etaMedial, clim=[-1, 1], order=I, hideColorbar=true);

EphysUnit.plotETA(ax(3), etaLateral, clim=[-1, 1], order=I);
ylabel(ax(2:3), '')
ylabel(ax(1), 'Unit', FontSize=14)
xlabel(ax(1:3), '')
title(ax(1), 'all reach trials', FontSize=14)
title(ax(2), 'medial reach trials', FontSize=14)
title(ax(3), 'lateral reach trials', FontSize=14)
colorbar(ax(3), Position=[0.910,0.108,0.014,0.815])

ax(3).Colorbar.Label.String = 'Normalized spike rate (a.u.)';

xlabel(tl, 'Time from bar-contact (s)', FontSize=14)
