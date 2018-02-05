% Prototype
tr = TetrodeRecording();
chunkSize = 10;
tr.ReadFiles(chunkSize);
channels = [tr.Amplifier.Channels.NativeOrder] + 1;
tr.SpikeDetect(channels, 'NumSigmas', 4, 'WaveformWindow', [-1, 1]);
tr.SpikeSort(channels, 'ClusterMethod', 'kmeans', 'WaveformWindow', [], 'FeatureMethod', 'WaveletTransform', 'Dimension', 10);
clear chunkSize channels


for iChannel = 1:16
	tr.PlotChannelRaw(iChannel, [100, 101], [-300, 300]);
end

for iChannel = 17:32
	tr.PlotChannelRaw(iChannel, [100, 101], [-300, 300]);
end



% 1101
tr = TetrodeRecording();
chunkSize = 10;
tr.ReadFiles(chunkSize);
tr.SpikeDetect(tr.MapChannelID([4, 7, 11]), 'NumSigmas', 4, 'WaveformWindow', [-1, 1]);
tr.ReadFiles(chunkSize, 'Chunks', 'remaining', 'SpikeDetect', true, 'DigitalDetect', true);
tr.SpikeSort(10, 'ClusterMethod', 'kmeans', 'WaveformWindow', [], 'FeatureMethod', 'WaveletTransform', 'Dimension', 10);
[clu, tree] = tr.SPC(10);
TetrodeRecording.RandomWords();
clear chunkSize
