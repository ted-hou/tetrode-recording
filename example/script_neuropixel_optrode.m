%% Spike detection
tr = TetrodeRecording;
tr.SelectFiles();
tr.ReadFiles(Duration=120, NumSigmas=4, NumSigmasReturn=1.25, NumSigmasReject=20, WaveformWindow=[-1, 1])
tr.SaveNeuropixelIO()
%% Read detected spikes and NIDQ digital/analog channels
tr.LoadNeuropixelIO();
tr.ParseNeuropixelIO();
tr.LoadSpikes(1:384);

%% Spike sort
tr.SpikeSort(1:384, Dimension=10, FeatureMethod='PCA', WaveformWindow=[-1, 1], ClusterMethod='kmeans', NumClusters=3);