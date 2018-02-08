% Prototype
tr = TetrodeRecording(); tr.Preview();

% 1113
% No spikes
tr = TetrodeRecording();
chunkSize = 10;
channels = [5 6 10 12 31 32];
tr.ReadFiles(chunkSize, 'DigitalDetect', true);
tr.SpikeDetect(channels, 'NumSigmas', 4, 'WaveformWindow', [-1, 1], 'ThresholdMode', 'MinPeakProminence');
tr.ReadFiles(chunkSize, 'Chunks', 'remaining', 'SpikeDetect', true);
tr.SpikeSort(channels, 'ClusterMethod', 'SPC', 'FeatureMethod', 'WaveletTransform', 'Dimension', 10);
clear chunkSize channels
TetrodeRecording.RandomWords();