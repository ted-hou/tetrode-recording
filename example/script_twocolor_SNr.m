%% Load TetrodeRecording objects and convert to EphysUnits
sessionInfo(1) = struct('name', 'daisy21_20240506',     'ml', +1300, 'ap', -3280, 'dv', -4600, 'duraOffset', -200, 'facing', 'front', 'model', '128DN');
sessionInfo(2) = struct('name', 'daisy21_20240510',     'ml', +1300, 'ap', -3280, 'dv', -4600, 'duraOffset', -200, 'facing', 'front', 'model', '128DN');
sessionInfo(3) = struct('name', 'daisy21_20240516',     'ml', +1300, 'ap', -3280, 'dv', -4600, 'duraOffset', -200, 'facing', 'front', 'model', '128DN');
sessionInfo(4) = struct('name', 'daisy21_20240523',     'ml', +1300, 'ap', -3280, 'dv', -4600, 'duraOffset', -200, 'facing', 'front', 'model', '128DN');
sessionInfo(5) = struct('name', 'daisy22_20240530',     'ml', -1300, 'ap', -3280, 'dv', -4600, 'duraOffset', -200, 'facing', 'front', 'model', '128DN');
sessionInfo(6) = struct('name', 'daisy22_20240531',     'ml', -1300, 'ap', -3280, 'dv', -4600, 'duraOffset', -200, 'facing', 'front', 'model', '128DN');

for iExp = 1:length(sessionInfo)
    try
        expName = sessionInfo(iExp).name;
        animalName = strsplit(expName, '_');
        animalName = animalName{1};
        files = dir(sprintf('C:\\SERVER\\%s\\SpikeSort\\tr_sorted_%s*.mat', animalName, expName));
        assert(~isempty(files), 'Cannot find files under C:\\SERVER\\%s\\SpikeSort\\tr_sorted_%s*.mat', animalName, expName)
        files = arrayfun(@(x) [x.folder, '\', x.name], files, UniformOutput=false);
        tr = TetrodeRecording.BatchLoadSimple(files);
        ar = AcuteRecording(tr, 'N/A');
        ar.binMoveResponse(tr, 'none', Window=[-1, 0], Store=true);
        eu = EphysUnit(ar, readWaveforms=true, cullITI=false, savepath='C:\SERVER\Units\TwoColor_SNr\BatchOne', tr=tr, sessionInfo=sessionInfo(iExp));
        clear tr ar eu
    catch
        warning('Failed to process session %i', iExp)
    end
end
clear iExp
%% Load EU
eu = EphysUnit.load('C:\SERVER\Units\TwoColor_SNr\BatchOne\SingleUnits_NonDuplicate');

%% Make rasters and save to disk
rd = eu.getRasterData('stimtwocolor', window=[-0.1, 0.4], durErr=1e-2, shutterDelay=0.01);
%%
fig = figure(Units='inches', Position=[0, 0, 6, 8]);
ax = axes(fig);
if ~exist('C:\SERVER\Figures\TwoColor_SNr\BatchOne', 'dir')
    mkdir('C:\SERVER\Figures\TwoColor_SNr\BatchOne')
end
for iEu = 1:length(eu)
    cla(ax)
    EphysUnit.plotRaster(ax, rd(iEu), xlim=[-0.1, 0.3], sz=2);
    print(fig, sprintf('C:\\SERVER\\Figures\\TwoColor_SNr\\BatchOne\\%s.png', eu(iEu).getName()), '-dpng');
end

%% Make PE-ISI and save to disk
clear p
p.isiWindow = [-0.5, 0.5];
p.isiRes = 1e-3;
p.xlim = [-0.1, 0.1];
p.path = 'C:\SERVER\Figures\TwoColor_SNr\BatchOne\PEISI\wavelength_power_duration';
if ~exist(p.path, 'dir')
    mkdir(p.path)
end
for iEu = 1:length(eu)
% for iEu = 1:length(eu)
    ax = axes(figure(Units='inches', Position=[0, 0, 10, 8]));
    groups = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'duration'});
    isi = NaN(length(groups), length(p.isiWindow(1):p.isiRes:p.isiWindow(2)));
    normSR = isi;
    for iGrp = 1:length(groups)
        [isi(iGrp, :), t] = eu(iEu).getMeanPEISI('stimtwocolor', groups(iGrp).trials, window=p.isiWindow, resolution=p.isiRes);
        normSR(iGrp, :) = 1./isi(iGrp, :) - mean(1./isi(iGrp, t<0), 'omitnan');
    end
    imagesc(ax, 1e3*t, [], normSR)
    xlim(ax, p.xlim*1e3)
    clim(ax, [-25, 25])
    ax.YAxisLocation = 'right';
    colormap(ax, 'jet')
    h = colorbar(ax, 'westoutside');
    h.Label.String = '\Deltasp/s';
    yticks(ax, 1:length(groups));
    yticklabels(ax, {groups.label})
    xlabel(ax, 'Time from opto onset (ms)')
    title(ax, eu(iEu).getName, Interpreter="none")
    print(ax.Parent, sprintf('%s\\%s.png', p.path, eu(iEu).getName), '-dpng')
    close(ax.Parent)
end
%% Create plot across units
p.path2 = 'C:\SERVER\Figures\TwoColor_SNr\BatchOne\PEISI\wavelength_power_location';
if ~exist(p.path2, 'dir')
    mkdir(p.path2)
end
for iEu = 1:length(eu)
    ax = axes(figure(Units='inches', Position=[0, 0, 6, 8]));
    groups = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'location'});
    isi = NaN(length(groups), length(p.isiWindow(1):p.isiRes:p.isiWindow(2)));
    normSR = isi;
    for iGrp = 1:length(groups)
        [isi(iGrp, :), t] = eu(iEu).getMeanPEISI('stimtwocolor', groups(iGrp).trials, window=p.isiWindow, resolution=p.isiRes);
        normSR(iGrp, :) = 1./isi(iGrp, :) - mean(1./isi(iGrp, t<0), 'omitnan');
    end
    imagesc(ax, 1e3*t, [], normSR)
    xlim(ax, p.xlim*1e3)
    clim(ax, [-25, 25])
    ax.YAxisLocation = 'right';
    colormap(ax, 'jet')
    h = colorbar(ax, 'westoutside');
    h.Label.String = '\Deltasp/s';
    yticks(ax, 1:length(groups));
    yticklabels(ax, {groups.label})
    xlabel(ax, 'Time from opto onset (ms)')
    title(ax, eu(iEu).getName, Interpreter="none")
    print(ax.Parent, sprintf('%s\\%s.png', p.path2, eu(iEu).getName), '-dpng')
    close(ax.Parent)
end

%% Create plot across units
close all
p.path3 = 'C:\SERVER\Figures\TwoColor_SNr\BatchOne\PEISI\wavelength_power';
p.xlim = [-0.05, 0.1];
if ~exist(p.path3, 'dir')
    mkdir(p.path3)
end
for iEu = 1:length(eu)
    ax = axes(figure(Units='inches', Position=[0, 0, 6, 4]));
    groups = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power'});
    isi = NaN(length(groups), length(p.isiWindow(1):p.isiRes:p.isiWindow(2)));
    normSR = isi;
    for iGrp = 1:length(groups)
        [isi(iGrp, :), t] = eu(iEu).getMeanPEISI('stimtwocolor', groups(iGrp).trials, window=p.isiWindow, resolution=p.isiRes);
        normSR(iGrp, :) = (1./isi(iGrp, :) - mean(1./isi(iGrp, t<=0), 'omitnan'))./std(1./isi(iGrp, t<0), 0, 'omitnan');
    end
    imagesc(ax, 1e3*t, [], normSR)
    xlim(ax, p.xlim*1e3)
    clim(ax, [-3, 3])
    ax.YAxisLocation = 'right';
    colormap(ax, 'jet')
    h = colorbar(ax, 'westoutside');
    h.Label.String = 'Normalized spike rate (a.u.)';
    yticks(ax, 1:length(groups));
    yticklabels(ax, {groups.label})
    xlabel(ax, 'Time from opto onset (ms)')
    title(ax, eu(iEu).getName, Interpreter="none")
    print(ax.Parent, sprintf('%s\\%s.png', p.path3, eu(iEu).getName), '-dpng')
    close(ax.Parent)
end
