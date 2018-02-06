% Prototype
tr = TetrodeRecording();
chunkSize = 10;
tr.ReadFiles(chunkSize, 'DigitalDetect', true);
channels = tr.MapChannelID([tr.Amplifier.Channels.NativeOrder]);
tr.SpikeDetect(channels, 'NumSigmas', 4, 'WaveformWindow', [-1, 1]);
tr.SpikeSort(channels, 'ClusterMethod', 'kmeans', 'FeatureMethod', 'PCA', 'Dimension', 3);
clear chunkSize channels
TetrodeRecording.RandomWords();

for iChannel = 1:16
	tr.PlotChannelRaw(iChannel, [100, 101], [-300, 300]);
end

for iChannel = 17:32
	tr.PlotChannelRaw(iChannel, [100, 101], [-300, 300]);
end

% 1101
% No spikes