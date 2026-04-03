%% Spike detection
% pause(40*60);

clear, clc
folders = { ...
    'C:\SERVER\daisy37\daisy37_20260331', ...
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
tr.SelectFiles(NeuropixelPath='C:\SERVER\daisy33\daisy33_20260226')
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
% 
% %% Convert to EphysUnits
% folders = { ...
%     ... 'C:\SERVER\daisy29\daisy29_20251023', ... SORTED, EU
%     ... 'C:\SERVER\desmond41\desmond41_20251028', ... SORTED, EU
%     ... 'C:\SERVER\desmond42\desmond42_20251030', ... SORTED, EU
%     ... 'C:\SERVER\daisy29\daisy29_20251031', ... SORTED, EU
%     ... 'C:\SERVER\desmond41\desmond41_20251104', ... SORTED, EU
%     ... 'C:\SERVER\desmond41\desmond41_20251117', ... SORTED, EU
%     ... 'C:\SERVER\daisy29\daisy29_20251118', ... SORTED, EU
%     ... 'C:\SERVER\daisy30\daisy30_20251119', ... SORTED, EU
%     ... 'C:\SERVER\daisy29\daisy29_20251024', ... SORTED, EU
%     ... 'C:\SERVER\daisy29\daisy29_20251025', ... SORTED, EU
%     ... 'C:\SERVER\daisy29\daisy29_20251027', ... TCE ERROR (952 pulses?)
%     ... 'C:\SERVER\daisy29\daisy29_20251028', ... SORTED, EU
%     ... 'C:\SERVER\desmond41\desmond41_20251029', ... SORTED, EU
%     ... 'C:\SERVER\daisy30\daisy30_20251029', ... SORTED, EU
%     ... 'C:\SERVER\daisy30\daisy30_20251030', ... SORTED, EU
%     ... 'C:\SERVER\desmond42\desmond42_20251031', ... SORTED, EU
%     ... 'C:\SERVER\daisy29\daisy29_20251103', ... SORTED, EU
%     ... 'C:\SERVER\daisy29\daisy29_20251120', ... SORTED, EU
%     ... 'C:\SERVER\daisy30\daisy30_20251121', ... SORTED, EU
%     ... 'C:\SERVER\desmond42\desmond42_20251121', ... SORTED, EU
%     ... 'C:\SERVER\desmond41\desmond41_20251124', ... SORTED, EU
%     'C:\SERVER\daisy31\daisy31_20251208', ... This is control (ChrimsonR) striatal optrode recording
%     'C:\SERVER\daisy31\daisy31_20251209', ... This is control (ChrimsonR) striatal optrode recording
%     'C:\SERVER\daisy32\daisy32_20251210', ... This is control (ChrimsonR) striatal optrode recording
%     'C:\SERVER\daisy32\daisy32_20251211', ... This is control (ChrimsonR) striatal optrode recording
%     'C:\SERVER\daisy31\daisy31_20251211', ... This is control (ChrimsonR) striatal optrode recording    
%     };
% 
% chunkSize = 32; % NumChannelsPerChunk
% for iSession = 1:length(folders)
%     try
%         tr = TetrodeRecording();
%         tr.SelectFiles(NeuropixelPath=folders{iSession});
%         tr.LoadNeuropixelIO();
%         tr.ParseNeuropixelIO(DigitalChannels={'Sync', 0; 'Lick', 1; 'Press', 2; 'Reward', 3; 'Timeout', 4; 'Mot2Busy', 5; 'CueLeft', 6; 'CueRight', 7});
% 
%         % IterativeArtifactRemoval (OnionPeeling): Load spikes
%         for iChunk = 1:(384/chunkSize)
%             channels = (iChunk-1)*chunkSize + 1 : iChunk*chunkSize;
%             tr.LoadSpikes(channels, Path='Spikes_Sorted');
% 
%             if isempty(tr.Spikes) || isempty([tr.Spikes.Channel])
%                 tr.Spikes = [];
%                 continue
%             end
% 
%             channels = [tr.Spikes.Channel];
% 
%             ar = AcuteRecording(tr, 'N/A');
%             ar.binMoveResponse(tr, 'none', Window=[-1, 0], Store=true);
%             eu = EphysUnit(ar, readWaveforms=false, cullITI=false, savepath='C:\SERVER\Units\TwoColor_Striatonigral', tr=tr);
% 
%             tr.Spikes = [];
%             clear eu
%         end
%     catch ME
%         warning('Could not process folder: %s', folders{iSession})
%         warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
%     end
% end
% 
% %%
% clear
% eu = EphysUnit.load('C:\SERVER\Units\TwoColor_Striatonigral', waveforms=false, spikecounts=false, spikerates=false);
% % Remove multiunits, fast (ISS test)
% eu = eu.removeMultiUnits(cullZeros=true);
% 
% % Remove drift, low spike rate units, fast
% clear c
% % Remove drift
% c.isDrifting = detectDriftingUnits(eu, smoothWindow=300, tolerance=0.05, spikeRateThreshold=5, includeITI=true);
% 
% % Filter by spike rate
% msr = arrayfun(@(eu) eu.SpikeRateStats.median, eu);
% p.minSpikeRate = 15;
% c.isSNr = msr >= p.minSpikeRate;
% 
% eu = eu(c.isSNr & ~c.isDrifting);
% 
% % Remove duplicates (slow, pairwise comparisons)
% [eu, isDuplicate] = eu.removeDuplicates(0.7);
% 
% eu.save('C:\SERVER\Units\TwoColor_Striatonigral\SingleUnit_NonDuplicate_NonDrift_SNr')
