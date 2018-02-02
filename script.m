for iChannel = 1:16
	tr.PlotChannelRaw(iChannel, [100, 101], [-300, 300]);
end

for iChannel = 17:32
	tr.PlotChannelRaw(iChannel, [100, 101], [-300, 300]);
end

% 1106
% 4, 10, 12, 19
tr = TetrodeRecording();

chunkSize = 10;

tr.ReadFiles(chunkSize);

tr.SpikeDetect(4, -25, 'ExclusionThreshold', -500, 'ExitThreshold', 15, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);
tr.SpikeDetect(10, -25, 'ExclusionThreshold', -500, 'ExitThreshold', 15, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);
tr.SpikeDetect(12, -25, 'ExclusionThreshold', -500, 'ExitThreshold', 15, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);
tr.SpikeDetect(19, -25, 'ExclusionThreshold', -500, 'ExitThreshold', 15, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);

tr.ReadFiles(chunkSize, 'Chunks', 'remaining', 'SpikeDetect', true, 'DigitalDetect', true);

TetrodeRecording.TTS(['All Done.\n',], true)

clear chunkSize

% 1107
% 12, 18, 21, 23, 30
tr = TetrodeRecording();

chunkSize = 10;

tr.ReadFiles(chunkSize);

tr.SpikeDetect(12, 25, 'ExclusionThreshold', 300, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);
tr.SpikeDetect(18, -30, 'ExclusionThreshold', 300, 'ExitThreshold', 15, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);
tr.SpikeDetect(21, -35, 'ExclusionThreshold', 300, 'ExitThreshold', 15, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);
tr.SpikeDetect(23, -35, 'ExclusionThreshold', 300, 'ExitThreshold', 15, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);
tr.SpikeDetect(30, -25, 'ExclusionThreshold', 300, 'ExitThreshold', 25, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);

tr.ReadFiles(chunkSize, 'Chunks', 'remaining', 'SpikeDetect', true, 'DigitalDetect', true);

TetrodeRecording.TTS(['All Done.\n',], true)

clear chunkSize

% 1108
% 3, 6, 8, 21, 23, 29
tr = TetrodeRecording();

chunkSize = 10;

tr.ReadFiles(chunkSize);

tr.SpikeDetect(3, -30, 'ExclusionThreshold', 400, 'ExitThreshold', 15, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);
tr.SpikeDetect(6, -30, 'ExclusionThreshold', 300, 'ExitThreshold', 15, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);
tr.SpikeDetect(8, -30, 'ExclusionThreshold', 300, 'ExitThreshold', 15, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);
tr.SpikeDetect(21, -30, 'ExclusionThreshold', 300, 'ExitThreshold', 15, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);
tr.SpikeDetect(23, -25, 'ExclusionThreshold', 300, 'ExitThreshold', 15, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);
tr.SpikeDetect(29, -25, 'ExclusionThreshold', 500, 'ExitThreshold', 15, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);

tr.ReadFiles(chunkSize, 'Chunks', 'remaining', 'SpikeDetect', true, 'DigitalDetect', true);

TetrodeRecording.TTS(['All Done.\n',], true)

clear chunkSize

% 1113
% 2 5 6 16 32
tr = TetrodeRecording();

chunkSize = 10;

tr.ReadFiles(chunkSize);

tr.SpikeDetect(2, 'NumSigmas', 4, 'WaveformWindow', [-1, 1]);
tr.SpikeDetect(5, 'NumSigmas', 4, 'WaveformWindow', [-1, 1]);
tr.SpikeDetect(6, 'NumSigmas', 4, 'WaveformWindow', [-1, 1]);
tr.SpikeDetect(16, 'NumSigmas', 4, 'WaveformWindow', [-1, 1]);
tr.SpikeDetect(32, 'NumSigmas', 4, 'WaveformWindow', [-1, 1]);


% tr.SpikeDetect(2, -25, 'ExclusionThreshold', 400, 'ExitThreshold', 15, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);
% tr.SpikeDetect(5, -25, 'ExclusionThreshold', 300, 'ExitThreshold', 15, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);
% tr.SpikeDetect(6, -25, 'ExclusionThreshold', 300, 'ExitThreshold', 15, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);
% tr.SpikeDetect(16, -25, 'ExclusionThreshold', 300, 'ExitThreshold', 15, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);
% tr.SpikeDetect(32, -25, 'ExclusionThreshold', 300, 'ExitThreshold', 15, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);

tr.ReadFiles(chunkSize, 'Chunks', 'remaining', 'SpikeDetect', true, 'DigitalDetect', true);

TetrodeRecording.TTS(['All Done.\n',], true)

clear chunkSize

% 1115
% 10 20
tr = TetrodeRecording();

chunkSize = 10;

tr.ReadFiles(chunkSize);

tr.SpikeDetect(10, 25, 'ExclusionThreshold', 500, 'ExitThreshold', -15, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);

tr.ReadFiles(chunkSize, 'Chunks', 'remaining', 'SpikeDetect', true, 'DigitalDetect', true);

TetrodeRecording.TTS(['All Done.\n',], true)

clear chunkSize

% 1116
tr = TetrodeRecording();
chunkSize = 10;
tr.ReadFiles(chunkSize);
tr.SpikeDetect(tr.MapChannelID([4, 7, 11]), 'NumSigmas', 4, 'WaveformWindow', [-1, 1]);
tr.ReadFiles(chunkSize, 'Chunks', 'remaining', 'SpikeDetect', true, 'DigitalDetect', true);
tr.SpikeSort(10, 'ClusterMethod', 'kmeans', 'WaveformWindow', [], 'FeatureMethod', 'WaveletTransform', 'Dimension', 10, 'HideResults', true);
[clu, tree] = tr.SPC(10);
TetrodeRecording.TTS(['All Done.\n',], true)
clear chunkSize

% tr.SpikeDetect(tr.MapChannelID(4), 40, 'ExclusionThreshold', 600, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);
% tr.SpikeDetect(tr.MapChannelID(7), -50, 'ExclusionThreshold', -500, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);
% tr.SpikeDetect(tr.MapChannelID(11), -30, 'ExitThreshold', 10, 'ExclusionThreshold', 1000, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);

% 1122
for iChannel = [1 4 5 24 31]
	tr.PlotChannelRaw(tr.MapChannelID(iChannel), [100, 101], [-300, 300]);
end

tr = TetrodeRecording();

chunkSize = 10;

tr.ReadFiles(chunkSize);

for iChannel = [1 4 10 12 15];
	tr.SpikeDetect(iChannel, -30, 'ExitThreshold', 15, 'ExclusionThreshold', 500, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);
end

tr.ReadFiles(chunkSize, 'Chunks', 'remaining', 'SpikeDetect', true, 'DigitalDetect', true);

TetrodeRecording.TTS(['All Done.\n',], true)

clear chunkSize