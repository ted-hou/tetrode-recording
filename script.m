% Prototype
tr = TetrodeRecording(); tr.Preview();

% Batch preview
tr = TetrodeRecording();
numFiles = 2;
dirs = uipickfiles('Prompt', 'Select (multiple) folders...');
dirs = dirs(isfolder(dirs));
for iDir = 1:length(dirs)
	files = dir([dirs{iDir}, '\*.rhd']);
	files = {files(1:numFiles).name};
	tr(iDir) = TetrodeRecording();
	tr(iDir).Path = [dirs{iDir}, '\'];
	tr(iDir).Files = files;
	tr(iDir).Preview();
end
TetrodeRecording.RandomWords();
clear dirs iDir files numFiles

% Batch process
chunkSize = 10;
selectedChannels = {tr.SelectedChannels};
allPaths = {tr.Path};
for iDir = 1:length(selectedChannels)
	channels = selectedChannels{iDir};
	if ~isempty(channels)
		tr = TetrodeRecording();
		tr.Path = allPaths{iDir};
		files = dir([tr.Path, '*.rhd']);
		tr.Files = {files.name};
		tr.ReadFiles(chunkSize, 'DigitalDetect', true);
		tr.SpikeDetect(channels, 'NumSigmas', 4, 'WaveformWindow', [-1, 1]);
		tr.ReadFiles(chunkSize, 'Chunks', 'remaining', 'SpikeDetect', true, 'DigitalDetect', true);
		tr.SpikeSort(channels, 'ClusterMethod', 'kmeans', 'FeatureMethod', 'WaveletTransform', 'Dimension', 10);
		tr.PlotAllChannels();
		tr.ClearCache();
		save([tr.Path, 'tr_', datestr(datetime, 'yyyymmdd_HHMM'), '.mat'], tr)
	end
end
TetrodeRecording.RandomWords();
clear iDir chunkSize channels allPaths