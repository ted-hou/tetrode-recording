%% Spike detection
tr = TetrodeRecording;
tr.SelectFiles();
tr.ReadFiles(Duration=120, NumSigmas=4, NumSigmasReturn=1.5, NumSigmasReject=20, WaveformWindow=[-0.5, 1])
tr.SaveNeuropixelIO()

% Read detected spikes and NIDQ digital/analog channels
% tr.LoadNeuropixelIO();
tr.ParseNeuropixelIO();

% Spike sort
tr.LoadSpikes(1:384);
tr.SpikeSort(1:384, Dimension=10, FeatureMethod='PCA', WaveformWindow=[-0.5, 0.5], ClusterMethod='kmeans', NumClusters=3);
tr.SpikeCullLowSpikeRateClusters(1:384, MinSpikeRate=1);
tr.SpikeSort(1:384, Dimension=10, FeatureMethod='PCA', WaveformWindow=[-0.5, 0.5], ClusterMethod='kmeans', NumClusters=3);
tr.SaveSpikes(Path='Spikes_AutoSorted');

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
tr = TetrodeRecording();
tr.SelectFiles();
tr.LoadNeuropixelIO();
tr.ParseNeuropixelIO();
tr.LoadSpikes(1:128, Path='Spikes_AutoSortedIterative');
tr.PlotAllChannels(Channels=1:128, plotMethod='mean')
% tr.SaveSpikes(Channels=1:128, Path='Spikes_Sorted')


%% IterativeArtifactRemoval (OnionPeeling)
folders = { ...
    'C:\SERVER\desmond38\desmond38_20250401', ...
    'C:\SERVER\desmond38\desmond38_20250402', ...
    'C:\SERVER\desmond38\desmond38_20250403', ...
%     'C:\SERVER\desmond38\desmond38_20250407', ...
%     'C:\SERVER\desmond38\desmond38_20250417', ...
%     'C:\SERVER\desmond39\desmond39_20250404', ...
%     'C:\SERVER\desmond39\desmond39_20250407', ...
%     'C:\SERVER\desmond39\desmond39_20250408', ...
%     'C:\SERVER\desmond39\desmond39_20250423', ...
    };

chunkSize = 16; % NumChannelsPerChunk
for iSession = 1:length(folders)
    try
        tr = TetrodeRecording();
        tr.SelectFiles(NeuropixelPath=folders{iSession});
%         tr.LoadNeuropixelIO();
%         tr.ParseNeuropixelIO();
        
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
            tr.IterativeArtifactRemoval(channels, MinSpikeRate=2, KIterative=4, KFinal=2, MaxIters=5, ...
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



%%
for iChannel = [obj.Spikes.Channel]
    clusters = unique(obj.Spikes(iChannel).Cluster.Classes);
    if length(clusters) == 1
        obj.SpikeSort(iChannel, Dimension=10, FeatureMethod='PCA', WaveformWindow=[-0.5, 0.5], ClusterMethod='kmeans', NumClusters=2)
    end
end