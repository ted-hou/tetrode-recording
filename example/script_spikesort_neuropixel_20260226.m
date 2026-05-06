%% Spike detection
pause(40*60);

clear, clc
folders = { ...
    ... 'C:\SERVER\daisy33\daisy33_20260226', ...
    ... 'C:\SERVER\daisy33\daisy33_20260227', ...
    ... 'C:\SERVER\daisy33\daisy33_20260302', ...
    ... 'C:\SERVER\daisy33\daisy33_20260303', ... 16 channels only
    ... 'C:\SERVER\daisy34\daisy34_20260303', ...
    ... 'C:\SERVER\daisy33\daisy33_20260304', ...
    ... 'C:\SERVER\daisy34\daisy34_20260304', ...
    ... 'C:\SERVER\daisy34\daisy34_20260305', ...
    ... 'C:\SERVER\daisy34\daisy34_20260306', ...
    ... 'C:\SERVER\daisy34\daisy34_20260309', ...
    ... 'C:\SERVER\desmond43\desmond43_20260309'...
    ... 'C:\SERVER\desmond43\desmond43_20260310'...
    ... 'C:\SERVER\daisy34\daisy34_20260325_incompleteStim' ...
    'C:\SERVER\daisy33\daisy33_20260320', ... HIGHERPOWER, SORTED
    'C:\SERVER\daisy34\daisy34_20260320', ... HIGHERPOWER, SORTED
    'C:\SERVER\daisy33\daisy33_20260323', ... HIGHERPOWER, SORTED
    ... 'C:\SERVER\daisy34\daisy34_20260323', ... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch1)
    'C:\SERVER\daisy33\daisy33_20260325', ... HIGHERPOWER, SORTED
    ... 'C:\SERVER\daisy33\daisy33_20260326', ... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch1)
    ... 'C:\SERVER\daisy34\daisy34_20260326', ... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch1)
    ... 'C:\SERVER\daisy33\daisy33_20260327', ... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch1)
    ... 'C:\SERVER\daisy34\daisy34_20260327', ... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch1)
    ... 'C:\SERVER\desmond43\desmond43_20260311'... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch1)
    ... 'C:\SERVER\desmond43\desmond43_20260312'... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch1)
    ... 'C:\SERVER\desmond43\desmond43_20260313'... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch1)
    };

for iSession = 1:length(folders)
    try
        tr = TetrodeRecording;
        tr.SelectFiles(NeuropixelPath=folders{iSession});
        tr.ReadFiles(Duration=240, NumSigmas=4, NumSigmasReturn=1.5, NumSigmasReject=40, WaveformWindow=[-0.5, 1])
        tr.SaveNeuropixelIO()

        % Read detected spikes and NIDQ digital/analog channels
        % tr.LoadNeuropixelIO();
        tr.ParseNeuropixelIO(DigitalChannels={'Sync', 0; 'Lick', 1; 'Press', 2; 'Reward', 3; 'Timeout', 4; 'Mot2Busy', 5; 'CueLeft', 6; 'CueRight', 7});

        % Spike sort
        tr.LoadSpikes(1:384);

        tr.IterativeArtifactRemoval(1:384, MinSpikeRate=0.5, KIterative=4, KFinal=2, MaxIters=5, ...
            DimensionIterative=3, DimensionFinal=10, FeatureMethod='PCA', ClusterMethod='kmeans', ...
            WaveformWindow=[-0.5, 0.5]);
        tr.SaveSpikes(Path='Spikes_AutoSortedIterative');
    catch ME
        warning('Could not process folder: %s', folders{iSession})
        warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
    end
end


%% Load sorted data on a different PC

clear, clc
tr = TetrodeRecording();
tr.SelectFiles(NeuropixelPath='C:\SERVER\daisy34\daisy34_20260320')
tr.LoadNeuropixelIO();
tr.ParseNeuropixelIO(DigitalChannels={'Sync', 0; 'Lick', 1; 'Press', 2; 'Reward', 3; 'Timeout', 4; 'Mot2Busy', 5; 'CueLeft', 6; 'CueRight', 7});

channels = 1:128;
tr.LoadSpikes(channels, Path='Spikes_AutoSortedIterative');
% tr.LoadSpikes(channels, Path='Spikes_AutoSortedIterative', ExpName='daisy29_20251120');
tr.PlotAllChannels(Channels=channels, plotMethod='mean')

%%
tr.SaveSpikes(Channels=channels, Path='Spikes_Sorted')
channels = 129:256;
tr.Spikes = [];
tr.LoadSpikes(channels, Path='Spikes_AutoSortedIterative');
tr.PlotAllChannels(Channels=channels, plotMethod='mean')
%%
tr.SaveSpikes(Channels=channels, Path='Spikes_Sorted')
channels = 257:384;
tr.Spikes = [];
tr.LoadSpikes(channels, Path='Spikes_AutoSortedIterative');
tr.PlotAllChannels(Channels=channels, plotMethod='mean')
%%
tr.SaveSpikes(Channels=channels, Path='Spikes_Sorted')

%% Convert to EphysUnits
folders = { ...
    ... 'C:\SERVER\daisy33\daisy33_20260226', ...
    ... 'C:\SERVER\daisy33\daisy33_20260227', ...
    ... 'C:\SERVER\daisy33\daisy33_20260302', ...
    ... 'C:\SERVER\daisy33\daisy33_20260303', ... 16 channels only
    ... 'C:\SERVER\daisy34\daisy34_20260303', ...
    ... 'C:\SERVER\daisy33\daisy33_20260304', ...
    ... 'C:\SERVER\daisy34\daisy34_20260304', ...
    ... 'C:\SERVER\daisy34\daisy34_20260305', ...
    ... 'C:\SERVER\daisy34\daisy34_20260306', ...
    ... 'C:\SERVER\daisy34\daisy34_20260309', ...
    ... 'C:\SERVER\desmond43\desmond43_20260309'...
    ... 'C:\SERVER\desmond43\desmond43_20260310'...
    ... 'C:\SERVER\daisy34\daisy34_20260325_incompleteStim' ...
    'C:\SERVER\daisy33\daisy33_20260320', ... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch2)
    'C:\SERVER\daisy34\daisy34_20260320', ... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch2)
    'C:\SERVER\daisy33\daisy33_20260323', ... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch2)
    ... 'C:\SERVER\daisy34\daisy34_20260323', ... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch1)
    'C:\SERVER\daisy33\daisy33_20260325', ... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch2)
    ... 'C:\SERVER\daisy33\daisy33_20260326', ... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch1)
    ... 'C:\SERVER\daisy34\daisy34_20260326', ... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch1)
    ... 'C:\SERVER\daisy33\daisy33_20260327', ... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch1)
    ... 'C:\SERVER\daisy34\daisy34_20260327', ... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch1)
    ... 'C:\SERVER\desmond43\desmond43_20260311'... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch1)
    ... 'C:\SERVER\desmond43\desmond43_20260312'... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch1)
    ... 'C:\SERVER\desmond43\desmond43_20260313'... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch1)
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
            tr.LoadSpikes(channels, Path='Spikes_Sorted');

            if isempty(tr.Spikes) || isempty([tr.Spikes.Channel])
                tr.Spikes = [];
                continue
            end

            channels = [tr.Spikes.Channel];

            ar = AcuteRecording(tr, 'N/A');
            ar.binMoveResponse(tr, 'none', Window=[-1, 0], Store=true);
            eu = EphysUnit(ar, readWaveforms=false, cullITI=false, savepath='C:\SERVER\Units\TwoColor_Striatonigral\Batch2', tr=tr);

            tr.Spikes = [];
            clear eu
        end
    catch ME
        warning('Could not process folder: %s', folders{iSession})
        warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
    end
end

%%
clear
eu = EphysUnit.load('C:\SERVER\Units\TwoColor_Striatonigral\Batch2', waveforms=false, spikecounts=false, spikerates=false);
% Remove multiunits, fast (ISS test)
eu = eu.removeMultiUnits(cullZeros=true);

% Remove drift, low spike rate units, fast
clear c
% Remove drift
c.isDrifting = detectDriftingUnits(eu, smoothWindow=300, tolerance=0.05, spikeRateThreshold=5, includeITI=true);

% Filter by spike rate
msr = arrayfun(@(eu) eu.SpikeRateStats.median, eu);
p.minSpikeRate = 15;
c.isSNr = msr >= p.minSpikeRate;

eu = eu(c.isSNr & ~c.isDrifting);

% Remove duplicates (slow, pairwise comparisons)
[eu, isDuplicate] = eu.removeDuplicates(0.7);

eu.save('C:\SERVER\Units\TwoColor_Striatonigral\SingleUnit_NonDuplicate_NonDrift_SNr')

%%

clear clc
eu = EphysUnit.load('C:\SERVER\Units\TwoColor_Striatonigral\SingleUnit_NonDuplicate_NonDrift_SNr', waveforms=false, spikecounts=false, spikerates=false, ...
    animalNames={'daisy33', 'daisy34', 'desmond43'});

% Do psth
rd.stim = eu.getRasterData('stimtwocolor', window=[-0.1, 0.4], durErr=1e-3, shutterDelay=0, photoelectricBlankDuration=1.5e-3, photoelectricOffsetBlankWindow=[10e-3, 10.5e-3]);

%%
% ETA Stim
p.isiBaselineWindow = [-0.2, 0];
p.stimBluePowers = [20000]*1e-6; 
p.stimRedPowers = [20000]*1e-6;
p.stimBlueDurations = [10]*1e-3;
p.stimRedDurations = [10]*1e-3;

p.isiWindow = [-0.4, 0.4];
p.isiRes = 1e-3;
p.xlim.stim = [-0.05, 0.1];
p.xlim.move = [-4, 2];
p.rasterSzStim = 3;
p.rasterSzMove = 1;

p.artifacts = struct(event=[], length=[], lengthUnit=[], direction=[]);
p.artifacts(1) = struct(event='StimOn', length=1, lengthUnit='ms', direction='right');
p.artifacts(2) = struct(event='StimOff', length=1, lengthUnit='ms', direction='right');

close all
XBlue = cell(length(eu), 1);
XRed = cell(length(eu), 1);
% for iEu = 1:length(eu)
%     groupsBlue = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'duration'}, selectBy=struct(power=p.stimBluePowers, duration=p.stimBlueDurations, location=[], wavelength=[470, 473]));
%     groupsRed = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'duration'}, selectBy=struct(power=p.stimRedPowers, duration=p.stimRedDurations, location=[], wavelength=593));
% 
%     % rd = eu(iEu).getRasterData('stimtwocolor', p.isiWindow, trials=[groupsRed.trials], alignTo='start', shutterDelay=0, sort=false, photoelectricBlankDuration=0.5e-3);
%     % EphysUnit.plotRaster(rd)
% 
%     [isi, t] = eu(iEu).getMeanPEISI('stimtwocolor', [groupsBlue.trials], window=p.isiWindow, resolution=p.isiRes, photoelectricBlankDuration=1.5e-3);
%     selBaseline = t>p.isiBaselineWindow(1) & t<p.isiBaselineWindow(2);
%     normSR = (1./isi - mean(1./isi(:, selBaseline), 'omitnan')) ./ std(1./isi(:, selBaseline), 0, 2, 'omitnan');
%     XBlue{iEu} = normSR;
% 
%     [isi, t] = eu(iEu).getMeanPEISI('stimtwocolor', [groupsRed.trials], window=p.isiWindow, resolution=p.isiRes, ...
%         photoelectricBlankDuration=1.5e-3);
%     selBaseline = t>p.isiBaselineWindow(1) & t<p.isiBaselineWindow(2);
%     normSR = (1./isi - mean(1./isi(:, selBaseline), 'omitnan')) ./ std(1./isi(:, selBaseline), 0, 2, 'omitnan');
%     XRed{iEu} = normSR;
% 
%     % fprintf('nan=%i, nan=%i\n', nnz(isnan(XBlue{iEu})), nnz(isnan(XRed{iEu})))
% end

for iEu = 1:length(eu)
    groupsBlue = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'duration'}, selectBy=struct(power=p.stimBluePowers, duration=p.stimBlueDurations, location=[], wavelength=[470, 473]));
    groupsRed = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'duration'}, selectBy=struct(power=p.stimRedPowers, duration=p.stimRedDurations, location=[], wavelength=593));

    % rd = eu(iEu).getRasterData('stimtwocolor', p.isiWindow, trials=[groupsRed.trials], alignTo='start', shutterDelay=0, sort=false, photoelectricBlankDuration=0.5e-3);
    % EphysUnit.plotRaster(rd)

    etaTemp = eu(iEu).getETA('count', 'stimtwocolor', [-0.5, 0.5], normalize=[-0.5, -0.1], resolution=0.010, alignTo='start', ...
        trials = [groupsBlue.trials], artifacts=p.artifacts);
    XBlue{iEu} = etaTemp.X;

    etaTemp = eu(iEu).getETA('count', 'stimtwocolor', [-0.5, 0.5], normalize=[-0.5, -0.1], resolution=0.010, alignTo='start', ...
        trials = [groupsRed.trials], artifacts=p.artifacts);
    XRed{iEu} = etaTemp.X;
end

eta.stimBlue = struct(X=cat(1, XBlue{:}), t=etaTemp.t, N=[], D=[], stats=[]);
eta.stimRed = struct(X=cat(1, XRed{:}), t=etaTemp.t, N=[], D=[], stats=[]);

clear XBlue XRed iEu groupsBlue groupsRed isi t selBaseline normSR etaTemp

% Calculate META
clear meta

p.metaWindowStim = [0, 0.050];
p.posRespThresholdStim = 1;
p.negRespThresholdStim = -0.5;

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

c.isStimResponsive = c.isStimBlueUp | c.isStimBlueDown | c.isStimRedUp | c.isStimRedDown;

fprintf('ChrimsonR %i, CoChR %i\n', nnz(c.isStimRedUp), nnz(c.isStimBlueUpRedNotUpThereforeCoChrMaybe))

% Combined PEISI and Stim Raster
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

p.path = 'C:\SERVER\Figures\TwoColor_SNr_Striatonigral\ErinsBatch';

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
            % etaTemp = eu(iEu).getETA('count', 'stimtwocolor', [-0.5, 0.5], normalize=[-0.5, -0.3], resolution=0.010, alignTo='start', ...
            %     trials = groups(iGrp).trials);
            etaTemp = eu(iEu).getETA('count', 'stimtwocolor', [-0.5, 0.5], normalize=[-0.5, -0.1], resolution=0.010, alignTo='start', ...
                trials = groups(iGrp).trials, artifacts=p.artifacts);
            t = etaTemp.t;
            deltaSR(iGrp, :) = etaTemp.X;
        end
        imagesc(ax, 1e3*t, [], deltaSR)
        xline(ax, 0)
        xlim(ax, p.xlim.stim*1e3)
        % clim(ax, [-1.5, 1.5])
        ax.YAxisLocation = 'right';
        % colormap(ax, 'turbo')
        applyCustomColormap(ax, [-1, 1], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);
        h = colorbar(ax, 'westoutside');
        h.Label.String = '\DeltaSR (a.u.)';
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

        print(fig, sprintf('%s\\%s.png', p.path, eu(iEu).getName()), '-dpng', '-r0')
    % end
end

clear fig layout iEu ax groups isi deltaSR iGrp h trialTypes colors iTrialType trialType