% Prototype
tr = TetrodeRecording(); tr.Preview();

% Batch preview
tr = TetrodeRecording();
dirs = uipickfiles('Prompt', 'Select (multiple) folders...');
dirs = dirs(isfolder(dirs));
for iDir = 1:length(dirs)
	files = dir([dirs{iDir}, '\*.rhd']);
	files = {files(unique([1:round(length(files)/4):length(files), length(files)])).name};
	tr(iDir) = TetrodeRecording();
	tr(iDir).Path = [dirs{iDir}, '\'];
	tr(iDir).Files = files;
	tr(iDir).Preview();
end
TetrodeRecording.RandomWords();
clear dirs iDir files numFiles

% Batch process
close all
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
		expName = strsplit(tr.Path, '\');
		expName = expName{end - 1};
		save([tr.Path, '..\SpikeSort\tr_', expName, '.mat'], 'tr')
	end
end
TetrodeRecording.RandomWords();
clear iDir chunkSize channels allPaths expName

% Batch load
files = uipickfiles('Prompt', 'Select .mat files containing TetrodeRecording objects to load...', 'Type', {'*.mat', 'MAT-files'});
for iFile = 1:length(files)
	S(iFile) = load(files{iFile}, 'tr');
end
tr = [S.tr];
clear iFile files S
