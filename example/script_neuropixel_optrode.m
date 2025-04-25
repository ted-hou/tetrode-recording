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
tr.SpikeSort(1:384, Dimension=10, FeatureMethod='PCA', WaveformWindow=[-0.5, 1], ClusterMethod='kmeans', NumClusters=3);
tr.SaveSpikes(Path='Spikes_AutoSorted');
%%
tr.LoadSpikes(1:384, Path='Spikes_AutoSorted');

%% Load sorted data on a different PC
tr = TetrodeRecording();
tr.SelectFiles();
tr.LoadNeuropixelIO();
tr.ParseNeuropixelIO();
tr.LoadSpikes(1:384, Path='Spikes');

%%
tr.PlotAllChannels(Channels=129:256, plotMethod='mean')