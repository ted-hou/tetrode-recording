%%

    % 'C:\SERVER\daisy37\daisy37_20260331', ... Striatum A2A Chrimson Optrode Control, AutoSortedIterative
    % 'C:\SERVER\daisy37\daisy37_20260401', ... Striatum A2A Chrimson Optrode Control, AutoSortedIterative
    % 'C:\SERVER\daisy37\daisy37_20260402', ... Striatum A2A Chrimson Optrode Control, AutoSortedIterative
    % 'C:\SERVER\daisy37\daisy37_20260403', ... Striatum A2A Chrimson Optrode Control, AutoSortedIterative
    % 'C:\SERVER\daisy37\daisy37_20260406', ... Striatum A2A Chrimson Optrode Control, AutoSortedIterative
    % 'C:\SERVER\daisy37\daisy37_20260407', ... Striatum A2A Chrimson Optrode Control, AutoSortedIterative
    % 'C:\SERVER\daisy37\daisy37_20260408', ... Striatum A2A Chrimson Optrode Control, Not processed yet


%% Load first 128 channels

clear, clc
fprintf("Loading data...\n");
tr = TetrodeRecording();
tr.SelectFiles(NeuropixelPath='C:\SERVER\daisy37\daisy37_20260331')
tr.LoadNeuropixelIO();
tr.ParseNeuropixelIO(DigitalChannels={'Sync', 0; 'Lick', 1; 'Press', 2; 'Reward', 3; 'Timeout', 4; 'Mot2Busy', 5; 'CueLeft', 6; 'CueRight', 7});

channels = 1:128;
fprintf("Loading channels %i-%i...\n", channels(1), channels(end));
tr.LoadSpikes(channels, Path='Spikes_AutoSortedIterative');
tr.PlotAllChannels(Channels=channels, plotMethod='mean')

%% Save 1:128
fprintf("Saving channels %i-%i...\n", channels(1), channels(end));
tr.SaveSpikes(Channels=channels, Path='Spikes_Sorted')

%% Load 129:256
channels = 129:256;
fprintf("Loading channels %i-%i...\n", channels(1), channels(end));
tr.Spikes = [];
tr.LoadSpikes(channels, Path='Spikes_AutoSortedIterative');
tr.PlotAllChannels(Channels=channels, plotMethod='mean')
%% Save 129:256
fprintf("Saving channels %i-%i...\n", channels(1), channels(end));
tr.SaveSpikes(Channels=channels, Path='Spikes_Sorted')

%% Load 257:384
channels = 257:384;
fprintf("Loading channels %i-%i...\n", channels(1), channels(end));
tr.Spikes = [];
tr.LoadSpikes(channels, Path='Spikes_AutoSortedIterative');
tr.PlotAllChannels(Channels=channels, plotMethod='mean')
%% Save 257:384
fprintf("Saving channels %i-%i...\n", channels(1), channels(end));
tr.SaveSpikes(Channels=channels, Path='Spikes_Sorted')



%% Convert to EphysUnits
folders = { ...
    'C:\SERVER\daisy37\daisy37_20260331', ... Striatum A2A Chrimson Optrode Control, AutoSortedIterative
    'C:\SERVER\daisy37\daisy37_20260401', ... Striatum A2A Chrimson Optrode Control, AutoSortedIterative
    'C:\SERVER\daisy37\daisy37_20260402', ... Striatum A2A Chrimson Optrode Control, AutoSortedIterative
    'C:\SERVER\daisy37\daisy37_20260403', ... Striatum A2A Chrimson Optrode Control, AutoSortedIterative
    'C:\SERVER\daisy37\daisy37_20260406', ... Striatum A2A Chrimson Optrode Control, AutoSortedIterative
    'C:\SERVER\daisy37\daisy37_20260407', ... Striatum A2A Chrimson Optrode Control, AutoSortedIterative
    % 'C:\SERVER\daisy37\daisy37_20260408', ... Striatum A2A Chrimson Optrode Control, Not processed yet
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
            eu = EphysUnit(ar, readWaveforms=false, cullITI=false, savepath='C:\SERVER\Units\Striatum_A2A_ChromsonR', tr=tr);

            tr.Spikes = [];
            clear eu
        end
    catch ME
        warning('Could not process folder: %s', folders{iSession})
        warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
    end
end

%
clear
eu = EphysUnit.load('C:\SERVER\Units\Striatum_A2A_ChromsonR', waveforms=false, spikecounts=false, spikerates=false);
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

eu.save('C:\SERVER\Units\Striatum_A2A_ChromsonR\SingleUnit_NonDuplicate_NonDrift_SNr')

% %% Load ephysUnits
% eu = EphysUnit.load('C:\SERVER\Units\SNr_CoChR_Optrode\SingleUnit_NonDuplicate_NonDrift_SNr');
% 
% rd.stim = eu.getRasterData('stim', window=[-2, 2], durErr=1e-3, shutterDelay=0, photoelectricBlankDuration=0.5e-3, photoelectricOffsetBlankWindow=[10e-3, 10.5e-3]);
% 
% 
% for iEu = 1%:length(eu)
%     EphysUnit.plotRaster(axes(figure), rd.stim(iEu), xlim=[-2, 2])
% end
% 
% %% Calculate ETA
% eta = eu.getETA('count', 'stim', [-1, 1], resolution=0.020, alignTo='start');