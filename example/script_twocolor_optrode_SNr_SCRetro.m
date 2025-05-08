

%% Remove duplicates (Slow)
eu = EphysUnit.load('\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\Units\TwoColor_SNr_SCRetro', waveforms=false, spikecounts=false, spikerates=false, animalNames={'desmond38'});
%
euAll = eu;

eu = eu.removeMultiUnits(cullZeros=true);
[eu, isDuplicate] = eu.removeDuplicates(0.7);
eu.save('\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\Units\TwoColor_SNr_SCRetro\SingleUnit_NonDuplicate')

%% Remove drift, low spike rate units
clear c
eu = EphysUnit.load('\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\Units\TwoColor_SNr_SCRetro\SingleUnit_NonDuplicate', waveforms=false, spikecounts=false, spikerates=false);

% Remove drift
c.isDrifting = detectDriftingUnits(eu, smoothWindow=300, tolerance=0.05, spikeRateThreshold=15, includeITI=true);


% Filter by spike rate
msr = arrayfun(@(eu) eu.SpikeRateStats.median, eu);
p.minSpikeRate = 15;
c.isSNr = msr >= p.minSpikeRate;

eu = eu(c.isSNr & ~c.isDrifting);

eu.save('\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\Units\TwoColor_SNr_SCRetro\SingleUnit_NonDuplicate_NonDrift_SNr')

%% Fix extra red stim event for bad session
for iEu = find(string({eu.ExpName}) == "desmond38_20250403")
    assert(isscalar(iEu))
    assert(eu(iEu).EventTimes.LaserModRedOn(1) - eu(iEu).EventTimes.LaserModRedOff(1) < 1e-6)
    iBadStimOn = find(eu(iEu).EventTimes.StimOn == eu(iEu).EventTimes.LaserModRedOn(1));
    iBadStimOff = find(eu(iEu).EventTimes.StimOff == eu(iEu).EventTimes.LaserModRedOff(1));
    eu(iEu).EventTimes.StimOn(iBadStimOn) = [];
    eu(iEu).EventTimes.StimOff(iBadStimOff) = [];
    eu(iEu).EventTimes.LaserModRedOn(1) = [];
    eu(iEu).EventTimes.LaserModRedOff(1) = [];
    eu(iEu).Trials.Stim = eu(iEu).makeTrials('stim');
end


eu.save('\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\Units\TwoColor_SNr_SCRetro\SingleUnit_NonDuplicate_NonDrift_SNr')

%% Make spontaneous reach/lick trials
p.minTrialLength = 4;
for iEu = 1:length(eu)
    eu(iEu).Trials.Press = eu(iEu).makeTrials('press_spontaneous_clean', minSpontaneousTrialDuration=p.minTrialLength);
    eu(iEu).Trials.Lick = eu(iEu).makeTrials('lick_spontaneous_clean', minSpontaneousTrialDuration=p.minTrialLength);
end

eu.save('\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\Units\TwoColor_SNr_SCRetro\SingleUnit_NonDuplicate_NonDrift_SNr')

%%
eu = EphysUnit.load('\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\Units\TwoColor_SNr_SCRetro\SingleUnit_NonDuplicate_SNr', waveforms=false, spikecounts=false, spikerates=false);

% Make Raster
clear rd
rd.stim = eu.getRasterData('stimtwocolor', window=[-0.1, 0.4], durErr=1e-3, shutterDelay=0, photoelectricBlankDuration=0.5e-3);
rd.press = eu.getRasterData('press', window=[-4, 0], alignTo='stop');
rd.lick = eu.getRasterData('lick', window=[-4, 0], alignTo='stop');

% Make ETA
clear eta
eta.pressRaw = eu.getETA('count', 'press', [-4, 0], resolution=0.1, alignTo='stop', includeInvalid=false, normalize='none');
eta.lickRaw = eu.getETA('count', 'lick', [-4, 0], resolution=0.1, alignTo='stop', includeInvalid=false, normalize='none');
eta.press = eu.getETA('count', 'press', [-4, 0], resolution=0.1, alignTo='stop', includeInvalid=false, normalize=[-4, -2]);
eta.lick = eu.getETA('count', 'lick', [-4, 0], resolution=0.1, alignTo='stop', includeInvalid=false, normalize=[-4, -2]);
eta.pressRaw.X = eta.pressRaw.X ./ 0.1;
eta.lickRaw.X = eta.lickRaw.X ./ 0.1;

% Calculate META
clear meta
p.metaWindow = [-0.3, 0];
p.posRespThreshold = 1;
p.negRespThreshold = -0.5;
t = eta.press.t;
meta.press = mean(eta.press.X(:, t>=p.metaWindow(1) & t<=p.metaWindow(2)), 2, 'omitnan');
meta.lick = mean(eta.lick.X(:, t>=p.metaWindow(1) & t<=p.metaWindow(2)), 2, 'omitnan');
t = eta.pressRaw.t;
meta.pressRaw = mean(eta.press.X(:, t>=p.metaWindow(1) & t<=p.metaWindow(2)), 2, 'omitnan');
meta.lickRaw = mean(eta.lick.X(:, t>=p.metaWindow(1) & t<=p.metaWindow(2)), 2, 'omitnan');
clear t

c.isPressUp =         meta.press >= p.posRespThreshold;
c.isPressDown =       meta.press <= p.negRespThreshold;
c.isPressResponsive = c.isPressUp | c.isPressDown;
c.isLickUp =          meta.lick >= p.posRespThreshold;
c.isLickDown =        meta.lick <= p.negRespThreshold;
c.isLickResponsive =  c.isLickUp | c.isLickDown;

c.isLickUnresponsiveButUp = ~c.isLickResponsive & meta.lick > 0;
c.isLickUnresponsiveButDown = ~c.isLickResponsive & meta.lick < 0;
c.isPressUnresponsiveButUp = ~c.isPressResponsive & meta.press > 0;
c.isPressUnresponsiveButDown = ~c.isPressResponsive & meta.press < 0;


%% Combined PEISI and Stim Raster
% Stim Rasters
p.isiWindow = [-0.4, 0.4];
p.isiRes = 1e-3;
p.xlim.stim = [-0.1, 0.3];
p.xlim.move = [-4, 0];
p.path = 'C:\SERVER\Figures\TwoColor_SNr_SCRetro';
p.rasterSzStim = 1;
p.rasterSzMove = 1;

fig = figure(Units='inches', Position=[0, 0, 12, 8]);
clear layout

layout.w = [3, 4]; % Stim, reach/lick
layout.h = [3, 3, 4]; % Raster, peisi/peth

layout.tl = tiledlayout(fig, sum(layout.h), sum(layout.w), TileSpacing='tight', Padding='tight', TileIndexing='columnmajor');
layout.ax = gobjects(length(layout.h), length(layout.w));
layout.ax(1, 1) = nexttile(layout.tl, [sum(layout.h(1:2)), layout.w(1)]);
layout.ax(2, 1) = nexttile(layout.tl, [layout.h(3), layout.w(1)]);
layout.ax(1, 2) = nexttile(layout.tl, [layout.h(1), layout.w(2)]);
layout.ax(2, 2) = nexttile(layout.tl, [layout.h(2), layout.w(2)]);
layout.ax(3, 2) = nexttile(layout.tl, [layout.h(3), layout.w(2)]);

if ~exist(p.path, 'dir')
    mkdir(p.path)
end

for iEu = 1:length(eu)
    % try
        cla(layout.ax(1:2, 1))
        cla(layout.ax(:, 2))
        % Make Raster
        ax = layout.ax(1, 1);
        EphysUnit.plotRaster(ax, rd.stim(iEu), xlim=p.xlim.stim*1e3, sz=p.rasterSzStim, timeUnit='ms');
        ax.Legend.Location = 'northeast';
        ax.Legend.FontSize = 6;
        title(ax, eu(iEu).getName(), Interpreter="none")

        % Make PE-ISI
        ax = layout.ax(2, 1);
        groups = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power'});
        isi = NaN(length(groups), length(p.isiWindow(1):p.isiRes:p.isiWindow(2)));
        deltaSR = isi;
        for iGrp = 1:length(groups)
            [isi(iGrp, :), t] = eu(iEu).getMeanPEISI('stimtwocolor', groups(iGrp).trials, window=p.isiWindow, resolution=p.isiRes, ...
                photoelectricBlankDuration=0.5e-3);
            deltaSR(iGrp, :) = 1./isi(iGrp, :) - mean(1./isi(iGrp, t<0), 'omitnan');
        end
        imagesc(ax, 1e3*t, [], deltaSR)
        xline(ax, 0)
        xlim(ax, p.xlim.stim*1e3)
        clim(ax, [-50, 50])
        ax.YAxisLocation = 'right';
        colormap(ax, 'turbo')
        h = colorbar(ax, 'westoutside');
        h.Label.String = '\Deltasp/s';
        yticks(ax, 1:length(groups));
        yticklabels(ax, {groups.label})
        xlabel(ax, 'Time from opto onset (ms)')
        title(ax, 'Stim PETH (from ISI)')

        % Make Raster Press
        trialTypes = ["press", "lick"];
        colors = ["red", "blue"];
        for iTrialType = 1:2
            trialType = trialTypes(iTrialType);
            ax = layout.ax(iTrialType, 2);
            EphysUnit.plotRaster(ax, rd.(trialType)(iEu), xlim=p.xlim.move, sz=p.rasterSzMove);
            title(ax, trialType, Interpreter="none")

            ax = layout.ax(3, 2);
            hold(ax, 'on')
            fieldname = sprintf('%sRaw', trialType);
            plot(ax, eta.(fieldname).t, eta.(fieldname).X(iEu, :), LineWidth=1.5, Color=colors(iTrialType), DisplayName=trialTypes(iTrialType))
            xlim(ax, p.xlim.move)
            title(ax, 'Move PETH')
            xlabel(ax, 'Time to bar/spout contact (s)')
            ylabel(ax, 'sp/s')
            legend(ax, Location='northwest')
        end

        print(fig, sprintf('%s\\%s.png', p.path, eu(iEu).getName()), '-dpng')
    % end
end

clear fig ax iEu

clear iEu ax groups isi deltaSR iGrp h


%% ETA Stim
p.isiBaselineWindow = [-0.2, 0];
close all
XBlue = cell(length(eu), 1);
XRed = cell(length(eu), 1);
iEuBad = [];
for iEu = 1:length(eu)
    groupsBlue = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'duration'}, selectBy=struct(power=[25e-6, 50e-6, 100e-6, 500e-6], duration=[0.010, 0.020], location=[], wavelength=[470, 473]));
    groupsRed = eu(iEu).groupTwoColorStimTrials({'wavelength', 'duration'}, selectBy=struct(power=[], duration=[0.010, 0.020], location=[], wavelength=635));

    [isi, t] = eu(iEu).getMeanPEISI('stimtwocolor', [groupsBlue.trials], window=p.isiWindow, resolution=p.isiRes, photoelectricBlankDuration=0.5e-3);
    selBaseline = t>p.isiBaselineWindow(1) & t<p.isiBaselineWindow(2);
    normSR = (1./isi - mean(1./isi(:, selBaseline), 'omitnan')) ./ std(1./isi(:, selBaseline), 0, 2, 'omitnan');
    XBlue{iEu} = normSR;
    if nnz(isnan(isi)) == 1001
        iEuBad = [iEuBad, iEu];
        fprintf('%s\n', eu(iEu).ExpName)
    end

    [isi, t] = eu(iEu).getMeanPEISI('stimtwocolor', [groupsRed.trials], window=p.isiWindow, resolution=p.isiRes, ...
        photoelectricBlankDuration=0.5e-3);
    selBaseline = t>p.isiBaselineWindow(1) & t<p.isiBaselineWindow(2);
    normSR = (1./isi - mean(1./isi(:, selBaseline), 'omitnan')) ./ std(1./isi(:, selBaseline), 0, 2, 'omitnan');
    XRed{iEu} = normSR;

    fprintf('nan=%i, nan=%i\n', nnz(isnan(XBlue{iEu})), nnz(isnan(XRed{iEu})))
end

eta.stimBlue = struct(X=cat(1, XBlue{:}), t=t, N=[], D=[], stats=[]);
eta.stimRed = struct(X=cat(1, XRed{:}), t=t, N=[], D=[], stats=[]);


% clear groups isi normSR iGrp selBaseline t XBlue XRed


%% Plot ETA Heatmaps Reach Lick Reach Lick StimBlue StimRed
groupVar = NaN(length(eu), 1);
groupVar(c.isPressUnresponsiveButDown & c.isLickUp) = 0;
groupVar(c.isPressUnresponsiveButDown & c.isLickUnresponsiveButUp) = 0;
groupVar(c.isPressDown & c.isLickUp) = 0;
groupVar(c.isPressDown & c.isLickUnresponsiveButUp) = 0;
groupVar(c.isPressUp & c.isLickUnresponsiveButDown) = 1;
groupVar(c.isPressUp & c.isLickDown) = 1;
groupVar(c.isPressUnresponsiveButUp & c.isLickUnresponsiveButDown) = 1;
groupVar(c.isPressUnresponsiveButUp & c.isLickDown) = 1;

groupVar(c.isPressUnresponsiveButDown & c.isLickUnresponsiveButDown) = 2;
groupVar(c.isPressUnresponsiveButDown & c.isLickDown) = 2;
groupVar(c.isPressDown & c.isLickUnresponsiveButDown) = 2;
groupVar(c.isPressDown & c.isLickDown) = 2;
groupVar(c.isPressUp & c.isLickUp) = 3;
groupVar(c.isPressUp & c.isLickUnresponsiveButUp) = 3;
groupVar(c.isPressUnresponsiveButUp & c.isLickUp) = 3;
groupVar(c.isPressUnresponsiveButUp & c.isLickUnresponsiveButUp) = 3;

fig = figure;
ax = arrayfun(@(i) subplot(1, 6, i), 1:6);

EphysUnit.plotETA(ax(1), eta.press, ...
    clim=[-1.5, 1.5], xlim=[-4, 0], sortWindow=[-3, 0], signWindow=[-0.3, 0], ...
    sortThreshold=0.25, negativeSortThreshold=0.25, hidecolorbar=true);
EphysUnit.plotETA(ax(2), eta.lick, ...
    clim=[-1.5, 1.5], xlim=[-4, 0], sortWindow=[-3, 0], signWindow=[-0.3, 0], ...
    sortThreshold=0.25, negativeSortThreshold=0.25, hidecolorbar=true);

[~, order] = EphysUnit.plotETA(ax(3), eta.press, sortGroup=groupVar, ...
    clim=[-1.5, 1.5], xlim=[-4, 0], sortWindow=[-3, 0], signWindow=[-0.3, 0], ...
    sortThreshold=0.25, negativeSortThreshold=0.25, hidecolorbar=true);
EphysUnit.plotETA(ax(4), eta.lick, order=order, ...
    clim=[-1.5, 1.5], xlim=[-4, 0], hidecolorbar=true);

EphysUnit.plotETA(ax(5), eta.stimBlue, order=order, ...
    clim=[-3, 3], xlim=[-0.1, 0.3], hidecolorbar=true);
EphysUnit.plotETA(ax(6), eta.stimRed, order=order, ...
    clim=[-3, 3], xlim=[-0.1, 0.3], hidecolorbar=true);

title(ax([1, 3]), 'Reach')
title(ax([2, 4]), 'Lick')
title(ax(5), 'CoChR2 (470nm)')
title(ax(6), 'ChrimsonR (635nm)')
for iAx = 1:4
    applyCustomColormap(ax(iAx), [-1.5, 1.5], hlim=[0.375, 0, 0, -0.375], llim=[0.125, 0.5, 0.5, 0.25], hpwr=.5, lpwr=1, h0=0.33);
end

N = histcounts(groupVar, -0.5:2:3.5);
yline(ax(3), cumsum(N(1:end-1)) + 1, 'k--');
yline(ax(4), cumsum(N(1:end-1)) + 1, 'k--');
