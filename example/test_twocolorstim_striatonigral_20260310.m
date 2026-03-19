
folders = { ...
    'C:\SERVER\daisy34\daisy34_20260309',
    };

chunkSize = 32; % NumChannelsPerChunk
for iSession = 1:length(folders)
    try
        tr = TetrodeRecording();
        tr.SelectFiles(NeuropixelPath=folders{iSession});
        tr.LoadNeuropixelIO();
        tr.ParseNeuropixelIO(DigitalChannels={'Sync', 0; 'Lick', 1; 'Press', 2; 'Reward', 3; 'Timeout', 4; 'Mot2Busy', 5; 'CueLeft', 6; 'CueRight', 7});

        % IterativeArtifactRemoval (OnionPeeling): Load spikes
        for iChunk = 1:(384/chunkSize)
            channels = (iChunk-1)*chunkSize + 1 : iChunk*chunkSize;
            tr.LoadSpikes(channels, Path='Spikes_AutoSortedIterative');

            if isempty(tr.Spikes) || isempty([tr.Spikes.Channel])
                tr.Spikes = [];
                continue
            end

            channels = [tr.Spikes.Channel];

            ar = AcuteRecording(tr, 'N/A');
            ar.binMoveResponse(tr, 'none', Window=[-1, 0], Store=true);
            eu = EphysUnit(ar, readWaveforms=false, cullITI=false, savepath='E:\Data\Units\Test', tr=tr);

            tr.Spikes = [];
            clear eu
        end
    catch ME
        warning('Could not process folder: %s', folders{iSession})
        warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
    end
end
%%
eu = EphysUnit.load('E:\Data\Units\Test', waveforms=false, spikecounts=false, spikerates=false);

%% Do psth
rd.stim = eu.getRasterData('stimtwocolor', window=[-0.1, 0.4], durErr=1e-3, shutterDelay=0, photoelectricBlankDuration=1.5e-3);

%%
% ETA Stim
p.isiBaselineWindow = [-0.2, 0];
p.stimBluePowers = [500, 2000]*1e-6; 
p.stimRedPowers = [500, 2000]*1e-6;
p.stimBlueDurations = [10]*1e-3;
p.stimRedDurations = [10]*1e-3;

p.isiWindow = [-0.4, 0.4];
p.isiRes = 1e-3;
p.xlim.stim = [-0.05, 0.1];
p.xlim.move = [-4, 2];
p.path = 'E:\Data\Figures\Test';
p.rasterSzStim = 1;
p.rasterSzMove = 1;

close all
XBlue = cell(length(eu), 1);
XRed = cell(length(eu), 1);
for iEu = 1:length(eu)
    groupsBlue = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'duration'}, selectBy=struct(power=p.stimBluePowers, duration=p.stimBlueDurations, location=[], wavelength=[470, 473]));
    groupsRed = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'duration'}, selectBy=struct(power=p.stimRedPowers, duration=p.stimRedDurations, location=[], wavelength=593));

    % rd = eu(iEu).getRasterData('stimtwocolor', p.isiWindow, trials=[groupsRed.trials], alignTo='start', shutterDelay=0, sort=false, photoelectricBlankDuration=0.5e-3);
    % EphysUnit.plotRaster(rd)

    [isi, t] = eu(iEu).getMeanPEISI('stimtwocolor', [groupsBlue.trials], window=p.isiWindow, resolution=p.isiRes, photoelectricBlankDuration=1.5e-3);
    selBaseline = t>p.isiBaselineWindow(1) & t<p.isiBaselineWindow(2);
    normSR = (1./isi - mean(1./isi(:, selBaseline), 'omitnan')) ./ std(1./isi(:, selBaseline), 0, 2, 'omitnan');
    XBlue{iEu} = normSR;

    [isi, t] = eu(iEu).getMeanPEISI('stimtwocolor', [groupsRed.trials], window=p.isiWindow, resolution=p.isiRes, ...
        photoelectricBlankDuration=1.5e-3);
    selBaseline = t>p.isiBaselineWindow(1) & t<p.isiBaselineWindow(2);
    normSR = (1./isi - mean(1./isi(:, selBaseline), 'omitnan')) ./ std(1./isi(:, selBaseline), 0, 2, 'omitnan');
    XRed{iEu} = normSR;

    % fprintf('nan=%i, nan=%i\n', nnz(isnan(XBlue{iEu})), nnz(isnan(XRed{iEu})))
end

eta.stimBlue = struct(X=cat(1, XBlue{:}), t=t, N=[], D=[], stats=[]);
eta.stimRed = struct(X=cat(1, XRed{:}), t=t, N=[], D=[], stats=[]);

clear XBlue XRed iEu groupsBlue groupsRed isi t selBaseline normSR

% Calculate META
clear meta

p.metaWindowStim = [0.005, 0.019];
p.posRespThresholdStim = 2;
p.negRespThresholdStim = -1;

t = eta.stimBlue.t;
meta.stimBlue = mean(eta.stimBlue.X(:, t>=p.metaWindowStim(1) & t<=p.metaWindowStim(2)), 2, 'omitnan');
t = eta.stimRed.t;
meta.stimRed = mean(eta.stimRed.X(:, t>=p.metaWindowStim(1) & t<=p.metaWindowStim(2)), 2, 'omitnan');
clear t

c.isStimBlueUp = meta.stimBlue >= p.posRespThresholdStim;
c.isStimBlueDown = meta.stimBlue <= p.negRespThresholdStim;
c.isStimRedUp = meta.stimRed >= p.posRespThresholdStim;
c.isStimRedDown = meta.stimRed <= p.negRespThresholdStim;

c.isStimBlueUpRedUpThereforeChrimsonMaybe = c.isStimBlueUp & c.isStimRedUp;
c.isStimBlueUpRedNotUpThereforeCoChrMaybe = c.isStimBlueUp & ~c.isStimRedUp;
c.isStimBlueNotUpRedUpThereforeChrimsonMaybe = ~c.isStimBlueUp & c.isStimRedUp;
c.isStimBlueNotUpRedNotUp = ~c.isStimBlueUp & ~c.isStimRedUp;

fprintf('ChrimsonR %i, CoChR %i\n', nnz(c.isStimRedUp), nnz(c.isStimBlueUpRedNotUpThereforeCoChrMaybe))

%% Combined PEISI and Stim Raster
% Stim Rasters

fig = figure(Units='inches', Position=[0, 0, 3.5, 5]);
clear layout

layout.w = [3]; % Stim, reach/lick
layout.h = [3, 3, 4]; % Raster, peisi/peth

layout.tl = tiledlayout(fig, sum(layout.h), sum(layout.w), TileSpacing='tight', Padding='tight', TileIndexing='columnmajor');
layout.ax = gobjects(length(layout.h), length(layout.w));
layout.ax(1, 1) = nexttile(layout.tl, [sum(layout.h(1:2)), layout.w(1)]);
layout.ax(2, 1) = nexttile(layout.tl, [layout.h(3), layout.w(1)]);
% layout.ax(1, 2) = nexttile(layout.tl, [layout.h(1), layout.w(2)]);
% layout.ax(2, 2) = nexttile(layout.tl, [layout.h(2), layout.w(2)]);
% layout.ax(3, 2) = nexttile(layout.tl, [layout.h(3), layout.w(2)]);

p.path = 'E:\Data\Figures\Test';
selUnits = find(c.isStimBlueUpRedNotUpThereforeCoChrMaybe);

% p.path = 'C:\SERVER\Figures\TwoColor_SNr_SCRetro\ReverseInjection\ChrimsonR';
% selUnits = find(c.isStimRedUp);

if ~exist(p.path, 'dir')
    mkdir(p.path)
end

for iEu = 1:length(eu)
    % try
        % Make Raster
        ax = layout.ax(1, 1);
        cla(ax)
        EphysUnit.plotRaster(ax, rd.stim(iEu), xlim=p.xlim.stim*1e3, sz=p.rasterSzStim, timeUnit='ms');
        delete(ax.Legend)
        % ax.Legend.Location = 'northeast';
        % ax.Legend.FontSize = 6;
        title(ax, eu(iEu).getName(), Interpreter="none")

        % Make ETA
        ax = layout.ax(2, 1);
        cla(ax)
        groups = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power'});
        % isi = NaN(length(groups), length(p.isiWindow(1):p.isiRes:p.isiWindow(2)));
        % deltaSR = isi;
        clear deltaSR
        for iGrp = 1:length(groups)
            % [isi(iGrp, :), t] = eu(iEu).getMeanPEISI('stimtwocolor', groups(iGrp).trials, window=p.isiWindow, resolution=p.isiRes, ...
            %     photoelectricBlankDuration=0.5e-3);
            % deltaSR(iGrp, :) = 1./isi(iGrp, :) - mean(1./isi(iGrp, t<0), 'omitnan');
            % etaTemp = eu(iEu).getETA('count', 'stimtwocolor', [-0.5, 0.5], normalize=[-0.5, -0.3], resolution=0.025, alignTo='start', ...
            %     trials = groups(iGrp).trials);
            etaTemp = eu(iEu).getETA('count', 'stimtwocolor', [-0.5, 0.5], normalize=struct(mean=mean(), sd=[]), resolution=0.025, alignTo='start', ...
                trials = groups(iGrp).trials);
            t = etaTemp.t;
            deltaSR(iGrp, :) = etaTemp.X;
        end
        imagesc(ax, 1e3*t, [], deltaSR)
        xline(ax, 0)
        xlim(ax, p.xlim.stim*1e3)
        % clim(ax, [-1.5, 1.5])
        ax.YAxisLocation = 'right';
        % colormap(ax, 'turbo')
        applyCustomColormap(ax, [-1.5, 1.5], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);
        h = colorbar(ax, 'westoutside');
        h.Label.String = '\Deltasp/s';
        yticks(ax, 1:length(groups));
        yticklabels(ax, {groups.label})
        xlabel(ax, 'Time from opto onset (ms)')
        title(ax, 'Stim PETH (from ISI)')

        % % Make Raster Press
        % trialTypes = ["press", "lick"];
        % colors = ["red", "blue"];
        % cla(layout.ax(1, 2))
        % cla(layout.ax(2, 2))
        % cla(layout.ax(3, 2))
        % for iTrialType = 1:2
        %     trialType = trialTypes(iTrialType);
        %     ax = layout.ax(iTrialType, 2);
        %     EphysUnit.plotRaster(ax, rd.(trialType)(iEu), xlim=p.xlim.move, sz=p.rasterSzMove);
        %     title(ax, trialType, Interpreter="none")
        %     xline(ax, 0, '--', DisplayName=trialType)
        % 
        %     ax = layout.ax(3, 2);
        %     hold(ax, 'on')
        %     fieldname = sprintf('%sRaw', trialType);
        %     plot(ax, eta.(fieldname).t, eta.(fieldname).X(iEu, :), LineWidth=1.5, Color=colors(iTrialType), DisplayName=trialTypes(iTrialType))
        %     xlim(ax, p.xlim.move)
        %     title(ax, 'Move PETH')
        %     xlabel(ax, 'Time to bar/spout contact (s)')
        %     ylabel(ax, 'sp/s')
        %     legend(ax, Location='northwest')
        %     if iTrialType == 1
        %         xline(ax, 0, '--', DisplayName='bar/spout contact')
        %     end
        % end

        print(fig, sprintf('%s\\%s.png', p.path, eu(iEu).getName()), '-dpng')
    % end
end

clear fig layout iEu ax groups isi deltaSR iGrp h trialTypes colors iTrialType trialType
