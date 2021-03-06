% Batch preview
ptr = TetrodeRecording.BatchPreview();

% Save/load batch preview files
TetrodeRecording.BatchSave(ptr, 'Prefix', 'ptr_', 'DiscardData', true);
ptr = TetrodeRecording.BatchLoad();

% Batch process
TetrodeRecording.BatchProcess(ptr, 'NumSigmas', 3, 'WaveformWindow', [-0.5, 0.5], 'FeatureMethod', 'WaveletTransform', 'ClusterMethod', 'kmeans', 'Dimension', 10);
TetrodeRecording.BatchProcess(ptr, 'NumSigmas', 3, 'NumSigmasReturn', 1.5, 'NumSigmasReject', 40, 'WaveformWindow', [-0.35, 0.35], 'FeatureMethod', 'WaveletTransform', 'ClusterMethod', 'kmeans', 'Dimension', 10);
% Batch load
tr = TetrodeRecording.BatchLoad();

% Plot Channel
tr(iTr).PlotChannel([], 'Reference', 'CueOn', 'Event', 'PressOn', 'Exclude', 'LickOn', 'Event2', 'LickOn', 'Exclude2', 'PressOn', 'BinMethod', 'percentile', 'Bins', 3, 'RasterXLim', [-5, 1], 'ExtendedWindow', [-1, 1], 'WaveformYLim', 'auto');
tr(iTr).PlotChannel([], 'Reference', 'CueOn', 'Event', 'PressOn', 'Exclude', 'LickOn', 'Event2', '', 'Exclude2', '', 'RasterXLim', [-5, 1], 'ExtendedWindow', [-1, 1], 'WaveformYLim', 'auto', 'Clusters', 1);
tr(iTr).PlotChannel([], 'RasterXLim', [-5, 1], 'ExtendedWindow', [-1, 1], 'WaveformYLim', [-200, 200]);

iTrLastSaved = 0;
TetrodeRecording.BatchSave(tr(iTrLastSaved+1:iTr), 'Prefix', 'tr_sorted_', 'DiscardData', false); iTrLastSaved = iTr;

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
% Ref/noise cluster is last cluster

batchPlotList = {...
	'Daisy1', 20171114, 32, 1;...
	'Daisy1', 20171117, 10, 1;...
	'Daisy1', 20171117, 7, 1;...
	'Daisy1', 20171121, 24, 1;...
	'Daisy1', 20171121, 28, 1;...
	'Daisy1', 20171121, 1, 1;...
	'Daisy1', 20171122, 12, 1;...
	'Daisy1', 20171122, 15, 1;...
	'Daisy1', 20171128, 19, 1;...
	'Daisy1', 20171128, 28, 1;...
	'Daisy1', 20171130, 19, 1;...
	'Daisy3', 20180429, 30, 1 ...
	};

TetrodeRecording.BatchPlot(tr, batchPlotList, 'Reformat', 'RasterAndPETHAndWaveform', 'WaveformYLim', 'auto', 'RasterXLim', [-7, 0], 'ExtendedWindow', [-1, 0], 'CopyLegend', true);
TetrodeRecording.BatchPlot(tr, batchPlotList, 'Reformat', 'PETH', 'RasterXLim', [-7, 0], 'ExtendedWindow', [-1, 0], 'CopyLegend', false, 'CopyLabel', false);

clear batchPlotList