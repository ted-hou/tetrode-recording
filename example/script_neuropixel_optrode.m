%% Spike detection
tr = TetrodeRecording;
tr.SelectFiles();
tr.ReadFiles(Duration=240, NumSigmas=4, NumSigmasReturn=1.5, NumSigmasReject=20, WaveformWindow=[-0.5, 1])
tr.SaveNeuropixelIO()

% Read detected spikes and NIDQ digital/analog channels
% tr.LoadNeuropixelIO();
tr.ParseNeuropixelIO();

% Spike sort
tr.LoadSpikes(1:384);
tr.IterativeArtifactRemoval(1:384, MinSpikeRate=4, KIterative=4, KFinal=2, MaxIters=5, ...
    DimensionIterative=3, DimensionFinal=10, FeatureMethod='PCA', ClusterMethod='kmeans', ...
    WaveformWindow=[-0.5, 0.5]);
tr.SaveSpikes(Path='Spikes_AutoSortedIterative');

%% Spike detection
tr = TetrodeRecording;
tr.SelectFiles();
tr.LoadNeuropixelIO();
tr.ParseNeuropixelIO();

% Spike sort
tr.LoadSpikes(1:384, Path='Spikes');
tr.SpikeSort(1:384, Dimension=10, FeatureMethod='PCA', WaveformWindow=[-0.5, 0.5], ClusterMethod='kmeans', NumClusters=3);
tr.SpikeCullLowSpikeRateClusters(1:384, MinSpikeRate=1);
tr.SpikeSort(1:384, Dimension=10, FeatureMethod='PCA', WaveformWindow=[-0.5, 0.5], ClusterMethod='kmeans', NumClusters=3);
tr.SaveSpikes(Path='Spikes_AutoSorted');


%% Load sorted data on a different PC
clear, clc
tr = TetrodeRecording();
tr.SelectFiles();
tr.LoadNeuropixelIO();
tr.ParseNeuropixelIO();

%
% channels = 1:128;
% tr.LoadSpikes(channels, Path='Spikes_AutoSortedIterative');
% tr.PlotAllChannels(Channels=channels, plotMethod='mean')
% %%
% tr.SaveSpikes(Channels=channels, Path='Spikes_Sorted')
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


%% IterativeArtifactRemoval (OnionPeeling)
folders = { ...
    % 'C:\SERVER\desmond38\desmond38_20250401', ...
    % 'C:\SERVER\desmond38\desmond38_20250402', ...
    % 'C:\SERVER\desmond38\desmond38_20250403', ...
%     'C:\SERVER\desmond38\desmond38_20250407', ...
    % 'C:\SERVER\desmond38\desmond38_20250417', ...
    'C:\SERVER\desmond39\desmond39_20250404', ...
%     'C:\SERVER\desmond39\desmond39_20250407', ...
%     'C:\SERVER\desmond39\desmond39_20250408', ...
%     'C:\SERVER\desmond39\desmond39_20250423', ...
    };

chunkSize = 48; % NumChannelsPerChunk
for iSession = 1:length(folders)
    try
        tr = TetrodeRecording();
        tr.SelectFiles(NeuropixelPath=folders{iSession});
        tr.LoadNeuropixelIO();
        tr.ParseNeuropixelIO();
        
        % IterativeArtifactRemoval (OnionPeeling): Load spikes
        for iChunk = 1:(384/chunkSize)
            channels = (iChunk-1)*chunkSize + 1 : iChunk*chunkSize;
            tr.LoadSpikes(channels, Path='Spikes_AutoSorted');
            channels = [tr.Spikes.Channel];

            if isempty(channels)
                tr.Spikes = [];
                continue
            end
            % tr.PlotAllChannels(Channels=channels, plotMethod='mean')

            % IterativeArtifactRemoval (OnionPeeling): kmeans, pca
            tr.IterativeArtifactRemoval(channels, MinSpikeRate=4, KIterative=4, KFinal=2, MaxIters=5, ...
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

%% Convert to EphysUnits
folders = { ...
    'C:\SERVER\desmond38\desmond38_20250402', ...
    'C:\SERVER\desmond38\desmond38_20250407', ...
    % 'C:\SERVER\desmond38\desmond38_20250417', ...
    };

chunkSize = 48; % NumChannelsPerChunk
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
            eu = EphysUnit(ar, readWaveforms=false, cullITI=false, savepath='C:\SERVER\Units\TwoColor_SNr_SCRetro\Batch2', tr=tr);

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
tr.SelectedChannels = [257:265, 267:2:273, 274:315, 317, 319:384];

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