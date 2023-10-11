files = uipickfiles();

% 'C:\SERVER' is 
% '\\research.files.med.harvard.edu\neurobio\NEUROBIOLOGY SHARED\Assad Lab\Lingfeng\Data'


%% Create EphysUnit objects from TetrodeRecording.
for iTr = 1:length(files)
    clear S tr ar eu
    try
        S = load(files{iTr});
        tr = S.tr;
        assert(ismember(tr.GetAnimalName, {'daisy17', 'daisy18', 'desmond31'}))
        ar = AcuteRecording(tr, 'N/A');
        ar.binMoveResponse(tr, 'press_spontaneous', Window=[-1, 0], Store=true);
        eu = EphysUnit(ar, tr=tr, savepath='\\research.files.med.harvard.edu\neurobio\NEUROBIOLOGY SHARED\Assad Lab\Lingfeng\Data\Units\acute_DA\spontaneous', ...
            cullITI=false, readWaveforms=true);
    catch ME
        warning('Error while processing file %g (%s)', iTr, files{iTr});
        warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
    end
end

clear S tr ar eu iTr

%% Load EhpysUnit
eu = EphysUnit.load('\\research.files.med.harvard.edu\neurobio\NEUROBIOLOGY SHARED\Assad Lab\Lingfeng\Data\Units\acute_DA\spontaneous');

%% Calculate mean spike rate (across whole session) for each unit
meanSpikeRates = arrayfun(@(eu) mean(eu.SpikeRates), eu);

baselineSpikeRates = zeros(1, length(eu));
for iEu = 1:length(eu)
    eta = eu(iEu).getETA('count', 'press', resolution=0.1, window=[-6, -2], minTrialDuration=2, alignTo='stop', includeInvalid=false);
    baselineSpikeRates(iEu) = mean(eta.X)./0.1;
end
clear iEu eta

%% Cull units with spike rate>12 sp/s
eu(baselineSpikeRates >= 12 | baselineSpikeRates <= 1) = [];
%%
baselineSpikeRates = zeros(1, length(eu));
for iEu = 1:length(eu)
    eta = eu(iEu).getETA('count', 'press', resolution=0.1, window=[-6, -2], minTrialDuration=2, alignTo='stop', includeInvalid=false);
    baselineSpikeRates(iEu) = mean(eta.X)./0.1;
end
%% Calculate mean spike rate (across whole session) for each unit
eta = eu.getETA('count', 'press', resolution=0.1, window=[-6, 2], minTrialDuration=0, maxTrialDuration=Inf, alignTo='stop', includeInvalid=true);
eta.X = eta.X ./ 0.1;

etaNormalized = eu.getETA('count', 'press', resolution=0.1, window=[-6, 2], minTrialDuration=0, maxTrialDuration=Inf, alignTo='stop', includeInvalid=true, normalize=[-4, -2]);


%% Plot ETAs for each unit
close all

ax = axes(figure);
hold(ax, 'on')
plot(ax, eta.t, eta.X);
hold(ax, 'off')

ax = axes(figure);
EphysUnit.plotETA(ax, etaNormalized, clim=[-1.5, 1.5], sortWindow=[-2, 0], signWindow=[-0.5, 0], sortThreshold=0.5)
title(ax, 'Spontaneous')

%% Plot Rasters for each unit
close all
rd = eu.getRasterData('press', [-1, 2], alignTo='stop');
for iEu = 1:length(rd)
    f = figure(iEu);
    ax = axes(f);
    EphysUnit.plotRaster(ax, rd(iEu), xlim=[-6, 2], sz=2.5);
end