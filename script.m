tr = TetrodeRecording();

chunkSize = 10;

tr.ReadFiles(chunkSize);

tr.SpikeDetect(tr.MapChannelID(4), 40, 'ExclusionThreshold', 500, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);
tr.SpikeDetect(tr.MapChannelID(7), -50, 'ExclusionThreshold', -400, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);
tr.SpikeDetect(tr.MapChannelID(11), -50, 'ExclusionThreshold', -400, 'WaveformWindow', [-0.35, 0.35], 'WaveformWindowExtended', [-0.75, 1]);

tr.ReadFiles(chunkSize, 'Chunks', 'remaining', 'SpikeDetect', true);

TetrodeRecording.TTS(['All Done.\n',], true)

clear chunkSize