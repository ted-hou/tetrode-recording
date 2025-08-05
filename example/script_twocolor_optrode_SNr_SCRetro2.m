
%% Spike detection
clear, clc
folders = { ...
    % 'C:\SERVER\daisy27\daisy27_20250624', ...
    % 'C:\SERVER\daisy27\daisy27_20250626', ...
    % 'C:\SERVER\daisy28\daisy28_20250701', ...
    % 'C:\SERVER\daisy28\daisy28_20250702', ...
    % 'C:\SERVER\daisy27\daisy27_20250624', ...
    % 'C:\SERVER\daisy27\daisy27_20250626', ...
    % 'C:\SERVER\daisy27\daisy27_20250717', ...
    % 'C:\SERVER\daisy27\daisy27_20250721', ...
    %'C:\SERVER\daisy28\daisy28_20250716', ...
    %'C:\SERVER\daisy28\daisy28_20250718', ...
     %'C:\SERVER\daisy27\daisy27_20250717', ...
   % 'C:\SERVER\daisy27\daisy27_20250707', ...
   % 'C:\SERVER\daisy28\daisy28_20250729', ...
    'C:\SERVER\daisy27\daisy27_20250715' ...
    'C:\SERVER\daisy28\daisy28_20250714' ...
    };

for iSession = 1:length(folders)
    try
        tr = TetrodeRecording;
        tr.SelectFiles(NeuropixelPath=folders{iSession});
        tr.ReadFiles(Duration=240, NumSigmas=4, NumSigmasReturn=1.5, NumSigmasReject=40, WaveformWindow=[-0.5, 1])
        tr.SaveNeuropixelIO()
        
        % Read detected spikes and NIDQ digital/analog channels
        % tr.LoadNeuropixelIO();
        tr.ParseNeuropixelIO();
        
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

%% IterativeArtifactRemoval (OnionPeeling)
folders = { ...
   % 'C:\SERVER\daisy28\daisy28_20250716', ...
    %'C:\SERVER\daisy28\daisy28_20250718', ...
    };

chunkSize = 32; % NumChannelsPerChunk
for iSession = 1:length(folders)
    try
        tr = TetrodeRecording();
        tr.SelectFiles(NeuropixelPath=folders{iSession});
        tr.LoadNeuropixelIO();
        tr.ParseNeuropixelIO();
        
        % IterativeArtifactRemoval (OnionPeeling): Load spikes
        for iChunk = 1:(384/chunkSize)
            channels = (iChunk-1)*chunkSize + 1 : iChunk*chunkSize;
            tr.LoadSpikes(channels, Path='Spikes');
            channels = [tr.Spikes.Channel];

            if isempty(channels)
                tr.Spikes = [];
                continue
            end
            % tr.PlotAllChannels(Channels=channels, plotMethod='mean')

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
    % 'C:\SERVER\daisy27\daisy27_20250717', ...
    % 'C:\SERVER\daisy27\daisy27_20250721', ...

    %daisy28_20250728  18 16
    % 02

clear, clc
tr = TetrodeRecording();
tr.SelectFiles(NeuropixelPath='C:\SERVER\daisy28\daisy28_20250702')
tr.LoadNeuropixelIO();
tr.ParseNeuropixelIO();

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
        tr.ParseNeuropixelIO();
        
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