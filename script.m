tr = TetrodeRecording();

tr.ReadFiles(10, 'Chunks', 'first', 'SpikeDetect', false);

tr.SpikeDetect(tr.MapChannelID(4), 40, 'ExclusionThreshold', 500, 'WaveformWindow', [-0.3, 0.3], 'WaveformWindowExtended', [-1, 2]);
tr.SpikeDetect(tr.MapChannelID(7), -50, 'ExclusionThreshold', -400, 'WaveformWindow', [-0.3, 0.3], 'WaveformWindowExtended', [-1, 2]);
tr.SpikeDetect(tr.MapChannelID(11), -40, 'ExclusionThreshold', -400, 'WaveformWindow', [-0.3, 0.3], 'WaveformWindowExtended', [-1, 2]);