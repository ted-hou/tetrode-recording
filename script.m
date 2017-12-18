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

tr.SpikeDetect(2, -25, 'ExclusionThreshold', 400, 'ExitThreshold', 15, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);
tr.SpikeDetect(5, -25, 'ExclusionThreshold', 300, 'ExitThreshold', 15, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);
tr.SpikeDetect(6, -25, 'ExclusionThreshold', 300, 'ExitThreshold', 15, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);
tr.SpikeDetect(16, -25, 'ExclusionThreshold', 300, 'ExitThreshold', 15, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);
tr.SpikeDetect(32, -25, 'ExclusionThreshold', 300, 'ExitThreshold', 15, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);

tr.ReadFiles(chunkSize, 'Chunks', 'remaining', 'SpikeDetect', true, 'DigitalDetect', true);

TetrodeRecording.TTS(['All Done.\n',], true)

clear chunkSize

% 1116
tr = TetrodeRecording();

chunkSize = 10;

tr.ReadFiles(chunkSize);

tr.SpikeDetect(tr.MapChannelID(4), 40, 'ExclusionThreshold', 600, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);
tr.SpikeDetect(tr.MapChannelID(7), -50, 'ExclusionThreshold', -500, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);
% tr.SpikeDetect(tr.MapChannelID(11), -50, 'ExclusionThreshold', -400, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);

tr.ReadFiles(chunkSize, 'Chunks', 'remaining', 'SpikeDetect', true, 'DigitalDetect', true);

TetrodeRecording.TTS(['All Done.\n',], true)

clear chunkSize

