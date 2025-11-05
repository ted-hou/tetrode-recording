%% Spike detection
pause(60*60);

clear, clc
folders = { ...
    % 'C:\SERVER\daisy29\daisy29_20251023', ...
    % 'C:\SERVER\daisy29\daisy29_20251024', ...
    % 'C:\SERVER\daisy29\daisy29_20251025', ...
    % 'C:\SERVER\daisy29\daisy29_20251027', ...
    % 'C:\SERVER\desmond41\desmond41_20251028', ...
    % 'C:\SERVER\daisy29\daisy29_20251028', ...
    % 'C:\SERVER\desmond41\desmond41_20251029', ...
    % 'C:\SERVER\daisy30\daisy30_20251029', ...
    % 'C:\SERVER\daisy30\daisy30_20251030', ...
    % 'C:\SERVER\desmond42\desmond42_20251030', ...
    'C:\SERVER\daisy29\daisy29_20251031', ...
    'C:\SERVER\desmond42\desmond42_20251031', ...
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

        % % Remove spikes after train 80, i.e. only keep 10ms pulses
        % [tce, stimOn, stimOff, ~] = tr.LoadTwoColorExperiment();
        % selPulses = abs((stimOff - stimOn) - 0.010) < 0.001;
        % iLastPulse = find(selPulses, 1, 'last');
        % for iChn = 1:384
        %     toDiscard = tr.Spikes(iChn).Timestamps > stimOff(iLastPulse) + 10;
        %     tr.Spikes(iChn).SampleIndex(toDiscard) = [];
        %     tr.Spikes(iChn).Timestamps(toDiscard) = [];
        %     tr.Spikes(iChn).Waveforms(toDiscard, :) = [];
        % end
        % clear tce stimOn stimOff selPulses iLastPulse iChn toDiscard

        tr.IterativeArtifactRemoval(1:384, MinSpikeRate=0.5, KIterative=4, KFinal=2, MaxIters=5, ...
            DimensionIterative=3, DimensionFinal=10, FeatureMethod='PCA', ClusterMethod='kmeans', ...
            WaveformWindow=[-0.5, 0.5]);
        tr.SaveSpikes(Path='Spikes_AutoSortedIterative');
    catch ME
        warning('Could not process folder: %s', folders{iSession})
        warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
    end
end

%% Remove spikes after train 80 (i.e. only keep 10ms pulses) and redo spike sorting, this helps remove drift
folders = { ...
    'C:\SERVER\daisy29\daisy29_20251023', ...
    'C:\SERVER\daisy29\daisy29_20251024', ...
    'C:\SERVER\daisy29\daisy29_20251027', ...
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
            tr.LoadSpikes(channels, Path='Spikes');
            channels = [tr.Spikes.Channel];

            if isempty(channels)
                tr.Spikes = [];
                continue
            end

            % Remove spikes after train 80, i.e. only keep 10ms pulses
            [tce, stimOn, stimOff, ~] = tr.LoadTwoColorExperiment();
            selPulses = abs((stimOff - stimOn) - 0.010) < 0.001;
            iLastPulse = find(selPulses, 1, 'last');
            for iChn = channels(:)'
                toDiscard = tr.Spikes(iChn).Timestamps > stimOff(iLastPulse) + 10;
                tr.Spikes(iChn).SampleIndex(toDiscard) = [];
                tr.Spikes(iChn).Timestamps(toDiscard) = [];
                tr.Spikes(iChn).Waveforms(toDiscard, :) = [];
            end
            clear tce stimOn stimOff selPulses iLastPulse iChn toDiscard

            % IterativeArtifactRemoval (OnionPeeling): kmeans, pca
            tr.IterativeArtifactRemoval(channels, MinSpikeRate=0.5, KIterative=4, KFinal=2, MaxIters=5, ...
                DimensionIterative=3, DimensionFinal=10, FeatureMethod='PCA', ClusterMethod='kmeans', ...
                WaveformWindow=[-0.5, 0.5]);

            tr.SaveSpikes(Channels=channels, Path='Spikes_AutoSortedIterative')
            tr.Spikes = [];
        end
    catch ME
        warning('Could not process folder: %s', folders{iSession})
        warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
    end
end


%% Load sorted data on a different PC

clear, clc
tr = TetrodeRecording();
tr.SelectFiles(NeuropixelPath='C:\SERVER\daisy29\daisy29_20251023')
tr.LoadNeuropixelIO();
tr.ParseNeuropixelIO(DigitalChannels={'Sync', 0; 'Lick', 1; 'Press', 2; 'Reward', 3; 'Timeout', 4; 'Mot2Busy', 5; 'CueLeft', 6; 'CueRight', 7});

channels = 1:128;
tr.LoadSpikes(channels, Path='Spikes_AutoSortedIterative');
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
    %'C:\SERVER\daisy26\daisy26_20250521', ...
    %'C:\SERVER\desmond38\desmond38_20250520', ...
    %'C:\SERVER\desmond39\desmond39_20250522', ...
   % 'C:\SERVER\daisy27\daisy27_20250717', ...
   % 'C:\SERVER\daisy27\daisy27_20250721', ...
 % 'C:\SERVER\daisy28\daisy28_20250729'...
% 'C:\SERVER\daisy28\daisy28_20250728' ...
'C:\SERVER\daisy28\daisy28_20250718' ...
'C:\SERVER\daisy28\daisy28_20250716' ...
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
            eu = EphysUnit(ar, readWaveforms=false, cullITI=false, savepath='C:\SERVER\Units\TwoColor_SNr_SCRetro\Batch3', tr=tr);

            tr.Spikes = [];
            clear eu
        end
    catch ME
        warning('Could not process folder: %s', folders{iSession})
        warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
    end
end

%% Batch do stuff to SelectedChannels
tr.SelectedChannels = [8:11, 13, 15:35, 37, 39:128];

% for iChannel = tr.SelectedChannels
%     tr.ClusterRemove(iChannel, 1);
% end

tr.FeatureExtract(tr.SelectedChannels, 'Method', 'PCA', 'Dimension', 10, 'WaveformWindow', [-0.5, 0.5]);
tr.Cluster(tr.SelectedChannels, 'Clusters', [], 'Method', 'kmeans', 'NumClusters', 2);
tr.SpikeClusterAutoReorder(tr.SelectedChannels, verbose=false)

tr.SelectedChannels = [];
tr.PlotAllChannels(Channels=channels, plotMethod='mean')


%% Batch do stuff to SelectedChannels
tr.SelectedChannels = [257:384];

for iChannel = tr.SelectedChannels
    tr.ClusterRemove(iChannel, 1);
end

tr.FeatureExtract(tr.SelectedChannels, 'Method', 'PCA', 'Dimension', 10, 'WaveformWindow', [-0.5, 0.5]);
tr.Cluster(tr.SelectedChannels, 'Clusters', [], 'Method', 'kmeans', 'NumClusters', 2);
tr.SpikeClusterAutoReorder(tr.SelectedChannels, verbose=false)

tr.SelectedChannels = [];
tr.PlotAllChannels(Channels=channels, plotMethod='mean')
%%
for iChannel = tr.SelectedChannels
	for field = fieldnames(tr.Spikes)'
		tr.Spikes(iChannel).(field{1}) = [];
	end				
end