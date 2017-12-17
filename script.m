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