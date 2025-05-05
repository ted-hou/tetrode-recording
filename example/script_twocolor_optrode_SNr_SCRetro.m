
%%
eu = EphysUnit.load('\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\Units\TwoColor_SNr_SCRetro', waveforms=false, spikecounts=false, spikerates=false, animalNames={'daisy26'});
%%
euAll = eu;

eu = eu.removeMultiUnits(cullZeros=true);

%% Make Raster
rd = eu.getRasterData('stimtwocolor', window=[-0.1, 0.4], durErr=1e-2, shutterDelay=0);
%%
fig = figure(Units='inches', Position=[0, 0, 6, 8]);
ax = axes(fig);
if ~exist('C:\SERVER\Figures\TwoColor_SNr_SCRetro', 'dir')
    mkdir('C:\SERVER\Figures\TwoColor_SNr_SCRetro')
end
for iEu = 8
    cla(ax)
    EphysUnit.plotRaster(ax, rd(iEu), xlim=[-0.1, 0.3], sz=2);
    print(fig, sprintf('C:\\SERVER\\Figures\\TwoColor_SNr_SCRetro\\%s.png', eu(iEu).getName()), '-dpng');
end

%% Make PE-ISI
% eu = EphysUnit.load();
p.isiWindow = [-0.98, 0.98];
p.isiRes = 1e-3;
p.xlim = [-0.1, 0.1];

for iEu = 2
    ax = axes(figure(Units='inches', Position=[0, 0, 10, 8]));
    groups = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power'});
    isi = NaN(length(groups), length(p.isiWindow(1):p.isiRes:p.isiWindow(2)));
    normSR = isi;
    for iGrp = 1:length(groups)
        [isi(iGrp, :), t] = eu(iEu).getMeanPEISI('stimtwocolor', groups(iGrp).trials, window=p.isiWindow, resolution=p.isiRes);
        normSR(iGrp, :) = 1./isi(iGrp, :) - mean(1./isi(iGrp, t<0), 'omitnan');
    end
    imagesc(ax, 1e3*t, [], normSR)
    xlim(ax, p.xlim*1e3)
    clim(ax, [0, 100])
    ax.YAxisLocation = 'right';
    colormap(ax, 'turbo')
    h = colorbar(ax, 'westoutside');
    h.Label.String = '\Deltasp/s';
    yticks(ax, 1:length(groups));
    yticklabels(ax, {groups.label})
    xlabel(ax, 'Time from opto onset (ms)')
    title(ax, eu(iEu).getName, Interpreter="none")
%     print(ax.Parent, sprintf('%s\\%s.png', p.path, eu(iEu).getName), '-dpng')
%     close(ax.Parent)
end

%% All powers, 10ms duration
close all
% p.path = 'C:\SERVER\Figures\TwoColor_SNr\BatchOne\ETA\ReceptiveField_AllPowers_10ms_withControl';
% if ~exist(p.path, 'dir')
%     mkdir(p.path)
% end
for iEu = 1:length(eu)
    ax = axes(figure(Units='inches', Position=[2, 2, 8, 5]));
    groups = eu(iEu).groupTwoColorStimTrials({'duration', 'power', 'wavelength', 'location'}, selectBy=struct(power=[], duration=0.020, location=[], wavelength=[]));

    isi = NaN(length(groups), length(p.isiWindow(1):p.isiRes:p.isiWindow(2)));
    clear eta
    eta(length(groups)) = struct(X=[], t=[], N=[], D=[]);
    normSR = isi;    
    for iGrp = 1:length(groups)
        eta(iGrp) = eu(iEu).getETA('count', 'stimtwocolor', window=[-0.5, 0.5], resolution=0.010, normalize='none', trials=groups(iGrp).trials);
    end
    hold(ax, 'on')
    t = eta(1).t;
    X = vertcat(eta.X);
    imagesc(ax, 1e3*t, [], X);
    plot(ax, [0, 0], 0.5+[0, length(groups)], 'k-', LineWidth=2)
    ylim(ax, 0.5+[0, length(groups)])
    minorTicks = 4.5:8:length(groups)-3.5;
    majorTicks = 0.5:8:length(groups)+0.5;
    for i = majorTicks
        plot(ax, [-100, 500], [i, i], 'k-', LineWidth=1.5)
    end
    for i = minorTicks
        plot(ax, [-100, 500], [i, i], 'k--', LineWidth=1.5)
    end
    hold(ax, 'off')
    clim(ax, [-1.5, 1.5])
    xlim(ax, [-500, 500])
    colormap(ax, 'jet')
    h = colorbar(ax, 'westoutside');
    h.Label.String = 'Normalized spike rate (a.u.)';
    yticks(ax, 1:length(groups));
    yticklabels(ax, {groups.label})
    xlabel(ax, 'Time from opto onset (ms)')
    ax.YAxisLocation = 'right';
    title(ax, eu(iEu).getName, Interpreter="none")
%     print(ax.Parent, sprintf('%s\\%s.png', p.path, eu(iEu).getName), '-dpng')
%     close(ax.Parent)
end


%% Make ETA All powers, first 20 ms
groups = cell(length(eu), 1);
meanX = zeros(32, length(eu));
for iEu = 1:length(eu)
    groups{iEu} = eu(iEu).groupTwoColorStimTrials({'duration', 'power', 'wavelength', 'location'}, selectBy=struct(power=[25e-6, 50e-6, 100e-6, 500e-6, 2e-3], duration=0.020, location=[], wavelength=[]));
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