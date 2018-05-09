% Batch preview
ptr = TetrodeRecording.BatchPreview();

% Save/load batch preview files
TetrodeRecording.BatchSave(ptr, 'Prefix', 'ptr_', 'DiscardData', true);
ptr = TetrodeRecording.BatchLoad();

% Batch process
TetrodeRecording.BatchProcess(ptr, 'NumSigmas', 3, 'WaveformWindow', [-0.5, 0.5], 'FeatureMethod', 'WaveletTransform', 'ClusterMethod', 'kmeans', 'Dimension', 10);

% Batch load
tr = TetrodeRecording.BatchLoad();

% Batch plot
for iTr = 1:length(tr)
	for iChannel = [tr(iTr).Spikes.Channel]
		for iCluster = 1:(max(tr(iTr).Spikes(iChannel).Cluster.Classes) - 1)
			hFigure = tr(iTr).PlotChannel(iChannel, 'PrintMode', true, 'Clusters', iCluster, 'Reference', 'CueOn', 'Event', 'LickOn', 'Exclude', 'PressOn');
			tr(iTr).GUISavePlot([], [], hFigure, 'Filename', ['C:\MATLAB\DATA\daisy_1\SpikeSort\', tr(iTr).GetExpName(), '_Chn', num2str(iChannel), 'Clu', num2str(iCluster), '_FirstLick'])
			close(hFigure)
			hFigure = tr(iTr).PlotChannel(iChannel, 'PrintMode', true, 'Clusters', iCluster, 'Reference', 'CueOn', 'Event', 'PressOn', 'Exclude', 'LickOn');
			tr(iTr).GUISavePlot([], [], hFigure, 'Filename', ['C:\MATLAB\DATA\daisy_1\SpikeSort\', tr(iTr).GetExpName(), '_Chn', num2str(iChannel), 'Clu', num2str(iCluster), '_FirstPress'])
			close(hFigure)
		end
	end
end
clear iTr iChannel iCluster


% Batchplot based on date/channel, and paste to clipboard
batchPlotList = [...
	20171117, 10, 1;...
	20171122, 12, 1;...
	20171122, 15, 1;...
	20171128, 19, 1;...
	20171128, 28, 1;...
	20171121, 24, 1;...
	20171121, 28, 1;...
	20171114, 32, 1;...
	20171117, 7, 1;...
	20171121, 1, 1;...
	20171130, 19, 1 ...
	];

unique(batchPlotList(:, 1))

for iPlot = 1:size(batchPlotList, 1)
	thisDate = batchPlotList(iPlot, 1);
	thisChannel = batchPlotList(iPlot, 2);
	thisCluster = batchPlotList(iPlot, 3);

	for iTr = 1:length(tr)
		if ~isempty(strfind(tr(iTr).Path, num2str(thisDate)))
			thisRefCluster = max(tr(iTr).Spikes(thisChannel).Cluster.Classes);
			hFigure = tr(iTr).PlotChannel(thisChannel, 'PrintMode', true, 'Clusters', thisCluster, 'ReferenceCluster', thisRefCluster, 'Reference', 'CueOn', 'Event', 'PressOn', 'Exclude', 'LickOn', 'WaveformYLim', 'auto', 'RasterXLim', [-6, 0]);
			tr(iTr).GUISavePlot([], [], hFigure)
			input('Type anything to continue...\n');
			close(hFigure)
			break
		end
	end
end

clear iPlot thisDate thisChannel thisCluster iTr hFigure