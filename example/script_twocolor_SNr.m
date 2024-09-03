%% Load TetrodeRecording objects and convert to EphysUnits
sessionInfo(1) = struct(name='daisy21_20240506', ml=+1300, ap=-3280, dv=-4600, duraOffset=-200, facing='front', model='128DN', firstShutterControlTrain=NaN);
sessionInfo(2) = struct(name='daisy21_20240510', ml=+1300, ap=-3280, dv=-4600, duraOffset=-200, facing='front', model='128DN', firstShutterControlTrain=NaN);
sessionInfo(3) = struct(name='daisy21_20240516', ml=+1300, ap=-3280, dv=-4600, duraOffset=-200, facing='front', model='128DN', firstShutterControlTrain=NaN);
sessionInfo(4) = struct(name='daisy21_20240523', ml=+1300, ap=-3280, dv=-4600, duraOffset=-200, facing='front', model='128DN', firstShutterControlTrain=NaN);
sessionInfo(5) = struct(name='daisy22_20240530', ml=-1300, ap=-3280, dv=-4600, duraOffset=-200, facing='front', model='128DN', firstShutterControlTrain=121);
sessionInfo(6) = struct(name='daisy22_20240531', ml=-1300, ap=-3280, dv=-4600, duraOffset=-200, facing='front', model='128DN', firstShutterControlTrain=122);
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
p.isiWindow = [-0.98, 0.98];
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
%% Make PE-ISI and save to disk
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

%% Make PE-ISI and save to disk
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

%% Make PE-ISI and save to disk
close all
p.path4 = 'C:\SERVER\Figures\TwoColor_SNr\BatchOne\PEISI\ReceptiveField';
p.xlim = [-0.05, 0.1];
if ~exist(p.path4, 'dir')
    mkdir(p.path4)
end
for iEu = 1:length(eu)
    ax = axes(figure(Units='inches', Position=[0, 0, 6, 4]));
    groups = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'duration', 'location'}, selectBy=struct(power=25e-6, duration=0.010, location=[], wavelength=[]));

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
    print(ax.Parent, sprintf('%s\\%s.png', p.path4, eu(iEu).getName), '-dpng')
    close(ax.Parent)
end

%% Make ETA and save to disk
close all
p.path = 'C:\SERVER\Figures\TwoColor_SNr\BatchOne\ETA\ReceptiveField_25uW';
if ~exist(p.path, 'dir')
    mkdir(p.path)
end
for iEu = 1:length(eu)
    ax = axes(figure(Units='inches', Position=[2, 2, 8, 1.5]));
    groups = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'duration', 'location'}, selectBy=struct(power=25e-6, duration=0.010, location=[], wavelength=[]));

    isi = NaN(length(groups), length(p.isiWindow(1):p.isiRes:p.isiWindow(2)));
    clear eta
    eta(length(groups)) = struct(X=[], t=[], N=[], D=[], stats=[]);
    normSR = isi;
    for iGrp = 1:length(groups)
        eta(iGrp) = eu(iEu).getETA('count', 'stimtwocolor', window=[-0.5, 1], resolution=0.010, normalize=[-0.5, 0], trials=groups(iGrp).trials);
    end
    hold(ax, 'on')
    t = eta(1).t;
    X = vertcat(eta.X);
    imagesc(ax, 1e3*t, [], X);
    plot(ax, [0, 0], 0.5+[0, length(groups)], 'k-', LineWidth=2)
    ylim(ax, 0.5+[0, length(groups)])
    hold(ax, 'off')
    clim(ax, [-3, 3])
    xlim(ax, [-100, 500])
    colormap(ax, 'jet')
    h = colorbar(ax, 'westoutside');
    h.Label.String = 'Normalized spike rate (a.u.)';
    yticks(ax, 1:length(groups));
    yticklabels(ax, {groups.label})
    xlabel(ax, 'Time from opto onset (ms)')
    ax.YAxisLocation = 'right';
    title(ax, eu(iEu).getName, Interpreter="none")
    print(ax.Parent, sprintf('%s\\%s.png', p.path, eu(iEu).getName), '-dpng')
    close(ax.Parent)
end

close all
p.path = 'C:\SERVER\Figures\TwoColor_SNr\BatchOne\ETA\ReceptiveField_50uW';
if ~exist(p.path, 'dir')
    mkdir(p.path)
end
for iEu = 1:length(eu)
    ax = axes(figure(Units='inches', Position=[2, 2, 8, 1.5]));
    groups = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'duration', 'location'}, selectBy=struct(power=50e-6, duration=0.010, location=[], wavelength=[]));

    isi = NaN(length(groups), length(p.isiWindow(1):p.isiRes:p.isiWindow(2)));
    clear eta
    eta(length(groups)) = struct(X=[], t=[], N=[], D=[], stats=[]);
    normSR = isi;
    for iGrp = 1:length(groups)
        eta(iGrp) = eu(iEu).getETA('count', 'stimtwocolor', window=[-0.5, 1], resolution=0.010, normalize=[-0.5, 0], trials=groups(iGrp).trials);
    end
    hold(ax, 'on')
    t = eta(1).t;
    X = vertcat(eta.X);
    imagesc(ax, 1e3*t, [], X);
    plot(ax, [0, 0], 0.5+[0, length(groups)], 'k-', LineWidth=2)
    ylim(ax, 0.5+[0, length(groups)])
    hold(ax, 'off')
    clim(ax, [-3, 3])
    xlim(ax, [-100, 500])
    colormap(ax, 'jet')
    h = colorbar(ax, 'westoutside');
    h.Label.String = 'Normalized spike rate (a.u.)';
    yticks(ax, 1:length(groups));
    yticklabels(ax, {groups.label})
    xlabel(ax, 'Time from opto onset (ms)')
    ax.YAxisLocation = 'right';
    title(ax, eu(iEu).getName, Interpreter="none")
    print(ax.Parent, sprintf('%s\\%s.png', p.path, eu(iEu).getName), '-dpng')
    close(ax.Parent)
end

close all
p.path = 'C:\SERVER\Figures\TwoColor_SNr\BatchOne\ETA\ReceptiveField_100uW';
if ~exist(p.path, 'dir')
    mkdir(p.path)
end
for iEu = 1:length(eu)
    ax = axes(figure(Units='inches', Position=[2, 2, 8, 1.5]));
    groups = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'duration', 'location'}, selectBy=struct(power=100e-6, duration=0.010, location=[], wavelength=[]));

    isi = NaN(length(groups), length(p.isiWindow(1):p.isiRes:p.isiWindow(2)));
    clear eta
    eta(length(groups)) = struct(X=[], t=[], N=[], D=[], stats=[]);
    normSR = isi;
    for iGrp = 1:length(groups)
        eta(iGrp) = eu(iEu).getETA('count', 'stimtwocolor', window=[-0.5, 1], resolution=0.010, normalize=[-0.5, 0], trials=groups(iGrp).trials);
    end
    hold(ax, 'on')
    t = eta(1).t;
    X = vertcat(eta.X);
    imagesc(ax, 1e3*t, [], X);
    plot(ax, [0, 0], 0.5+[0, length(groups)], 'k-', LineWidth=2)
    ylim(ax, 0.5+[0, length(groups)])
    hold(ax, 'off')
    clim(ax, [-3, 3])
    xlim(ax, [-100, 500])
    colormap(ax, 'jet')
    h = colorbar(ax, 'westoutside');
    h.Label.String = 'Normalized spike rate (a.u.)';
    yticks(ax, 1:length(groups));
    yticklabels(ax, {groups.label})
    xlabel(ax, 'Time from opto onset (ms)')
    ax.YAxisLocation = 'right';
    title(ax, eu(iEu).getName, Interpreter="none")
    print(ax.Parent, sprintf('%s\\%s.png', p.path, eu(iEu).getName), '-dpng')
    close(ax.Parent)
end

close all
p.path = 'C:\SERVER\Figures\TwoColor_SNr\BatchOne\ETA\ReceptiveField_2mW';
if ~exist(p.path, 'dir')
    mkdir(p.path)
end
for iEu = 1:length(eu)
    ax = axes(figure(Units='inches', Position=[2, 2, 8, 1.5]));
    groups = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'duration', 'location'}, selectBy=struct(power=2e-3, duration=0.010, location=[], wavelength=[]));

    isi = NaN(length(groups), length(p.isiWindow(1):p.isiRes:p.isiWindow(2)));
    clear eta
    eta(length(groups)) = struct(X=[], t=[], N=[], D=[], stats=[]);
    normSR = isi;    
    for iGrp = 1:length(groups)
        eta(iGrp) = eu(iEu).getETA('count', 'stimtwocolor', window=[-0.5, 1], resolution=0.010, normalize=[-0.5, 0], trials=groups(iGrp).trials);
    end
    hold(ax, 'on')
    t = eta(1).t;
    X = vertcat(eta.X);
    imagesc(ax, 1e3*t, [], X);
    plot(ax, [0, 0], 0.5+[0, length(groups)], 'k-', LineWidth=2)
    ylim(ax, 0.5+[0, length(groups)])
    hold(ax, 'off')
    clim(ax, [-3, 3])
    xlim(ax, [-100, 500])
    colormap(ax, 'jet')
    h = colorbar(ax, 'westoutside');
    h.Label.String = 'Normalized spike rate (a.u.)';
    yticks(ax, 1:length(groups));
    yticklabels(ax, {groups.label})
    xlabel(ax, 'Time from opto onset (ms)')
    ax.YAxisLocation = 'right';
    title(ax, eu(iEu).getName, Interpreter="none")
    print(ax.Parent, sprintf('%s\\%s.png', p.path, eu(iEu).getName), '-dpng')
    close(ax.Parent)
end

%% All powers, 10ms duration
close all
p.path = 'C:\SERVER\Figures\TwoColor_SNr\BatchOne\ETA\ReceptiveField_AllPowers_10ms_withControl';
if ~exist(p.path, 'dir')
    mkdir(p.path)
end
for iEu = 1:length(eu)
    ax = axes(figure(Units='inches', Position=[2, 2, 8, 5]));
    iSession = find(strcmpi(eu(iEu).ExpName, {sessionInfo.name}));
    assert(~isempty(iSession))
    groups = eu(iEu).groupTwoColorStimTrials({'duration', 'power', 'wavelength', 'location'}, selectBy=struct(power=[], duration=0.010, location=[], wavelength=[]), firstShutterControlTrain=sessionInfo(iSession).firstShutterControlTrain);

    isi = NaN(length(groups), length(p.isiWindow(1):p.isiRes:p.isiWindow(2)));
    clear eta
    eta(length(groups)) = struct(X=[], t=[], N=[], D=[], stats=[]);
    normSR = isi;    
    for iGrp = 1:length(groups)
        eta(iGrp) = eu(iEu).getETA('count', 'stimtwocolor', window=[-0.5, 1], resolution=0.010, normalize=[-0.5, 0], trials=groups(iGrp).trials);
    end
    hold(ax, 'on')
    t = eta(1).t;
    X = vertcat(eta.X);
    imagesc(ax, 1e3*t, [], X);
    plot(ax, [0, 0], 0.5+[0, length(groups)], 'k-', LineWidth=2)
    ylim(ax, 0.5+[0, length(groups)])
    minorTicks = 4.5:8:length(groups)-3.5;
    majorTicks = 0.5:8:length(groups)+0.5;
    if ~isnan(sessionInfo(iSession).firstShutterControlTrain)
        minorTicks = minorTicks + 1;
        majorTicks = majorTicks + 1;
    end
    for i = majorTicks
        plot(ax, [-100, 500], [i, i], 'k-', LineWidth=1.5)
    end
    for i = minorTicks
        plot(ax, [-100, 500], [i, i], 'k--', LineWidth=1.5)
    end
    hold(ax, 'off')
    clim(ax, [-1.5, 1.5])
    xlim(ax, [-100, 500])
    colormap(ax, 'jet')
    h = colorbar(ax, 'westoutside');
    h.Label.String = 'Normalized spike rate (a.u.)';
    yticks(ax, 1:length(groups));
    yticklabels(ax, {groups.label})
    xlabel(ax, 'Time from opto onset (ms)')
    ax.YAxisLocation = 'right';
    title(ax, eu(iEu).getName, Interpreter="none")
    print(ax.Parent, sprintf('%s\\%s.png', p.path, eu(iEu).getName), '-dpng')
    close(ax.Parent)
end



%% All powers, first 10 ms
groups = cell(length(eu), 1);
meanX = zeros(32, length(eu));
for iEu = 1:length(eu)
    groups{iEu} = eu(iEu).groupTwoColorStimTrials({'duration', 'power', 'wavelength', 'location'}, selectBy=struct(power=[25e-6, 50e-6, 100e-6, 2e-3], duration=0.010, location=[], wavelength=[]));
%     isi = NaN(length(groups), length(p.isiWindow(1):p.isiRes:p.isiWindow(2)));
%     clear eta
%     eta(length(groups)) = struct(X=[], t=[], N=[], D=[], stats=[]);
%     normSR = isi;    
%     for iGrp = 1:length(groups)
%         eta(iGrp) = eu(iEu).getETA('count', 'stimtwocolor', window=[-0.5, 1], resolution=0.010, normalize=[-0.5, 0], trials=groups(iGrp).trials);
%     end
    if iEu > 1
        assert(isequal({groups{iEu}.label}, {groups{1}.label}))
    end
    clear eta
    for iGrp = 1:length(groups{iEu})
        eta(iGrp) = eu(iEu).getETA('count', 'stimtwocolor', window=[-0.5, 1], resolution=0.010, normalize=[-0.5, 0], trials=groups{iEu}(iGrp).trials);
    end
    t = eta(1).t;
    X = vertcat(eta.X);
    meanX(:, iEu) = mean(X(:, t<0.1 & t>0), 2, 'omitnan');
end
%%
ax = axes(figure);
imagesc(X')
clim(ax, [-1.5, 1.5])
colormap(ax, 'jet')
xticks(ax, 1:length(groups{1}));
xticklabels(ax, {groups{1}.label})