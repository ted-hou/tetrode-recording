batchPlotList = {...
	'desmond10', 20180910, 6, 1;... % DA up
	'desmond10', 20180910, 9, 1;... % DA up
	'desmond10', 20180910, 14, 1;... % DA
	'desmond10', 20180913, 5, 1;... % DA
	'desmond10', 20180915, 1, 1;... % DA
	'desmond10', 20180915, 1, 2;... % DA
	'desmond10', 20180917, 2, 1;... % DA up
	'desmond10', 20180917, 13, 1;... % DA
	'desmond10', 20180918, 2, 1;... % DA up
	'desmond10', 20180919, 2, 1;... % DA up
	'desmond10', 20180919, 5, 1;... % DA up
	'desmond10', 20180919, 10, 1;... % DA
	'desmond10', 20180919, 13, 1;... % DA
	'desmond10', 20180920, 4, 1;... % DA reward
	'desmond10', 20180920, 13, 1;... % DA
	'desmond10', 20180920, 16, 1;... % DA
	'desmond10', 20180922, 4, 1;... % DA
	'desmond10', 20180924, 25, 1;... % DA
	'desmond10', 20180925, 10, 1;... % DA
	'desmond10', 20180925, 20, 2;... % DA
	'desmond10', 20181017, 4, 1;... % DA
	'desmond10', 20181018, 7, 1;... % DA
	'desmond10', 20181018, 8, 2;... % DA
	'desmond10', 20181019, 1, 2;... % DA
	'desmond10', 20181019, 7, 1;... % DA
	'desmond10', 20181022, 1, 2;... % DA
	'desmond11', 20180911, 27, 1;... % DA
	'desmond11', 20180912, 28, 1;... % DA ?
	'desmond11', 20180914, 27, 1;... % DA
	'desmond11', 20180916, 20, 1;... % DA
	'desmond11', 20180922, 28, 1;... % DA
	'desmond11', 20180924, 28, 1;... % DA
	'desmond22', 20210624, 5, 2;... % DA, slow on
	'desmond22', 20210624, 7, 2;... % DA flat
	'desmond22', 20210630, 1, 1;... % DA flat
	'desmond22', 20210707, 12, 1;... % DA slow on
	'daisy8', 20210623, 4, 2;... % DA
	'daisy8', 20210625, 1, 1;... % DA inhibited by ChR2
	'daisy8', 20210707, 9, 1;... % DA slow off
	'daisy8', 20210709, 12, 1;... % DA fast off
};

expNames = cell(size(batchPlotList, 1), 1);
for iExp = 1:size(batchPlotList, 1)
	expNames{iExp} = [batchPlotList{iExp, 1}, '_', num2str(batchPlotList{iExp, 2})];
end

expNamesUnique = unique(expNames);

for iTr = 1:length(expNamesUnique)
    AHprogressBar(iTr, length(expNamesUnique), false, 1)
	tr = TetrodeRecording.BatchLoad(expNamesUnique(iTr));
    if contains(expNamesUnique(iTr), 'desmond22')
        tr.ReadDigitalEvents();
    end
	try
		TetrodeRecording.BatchPlot(tr, batchPlotList, 'RasterXLim', [-5, 2], 'ExtendedWindow', [0, 2], 'PlotStim', false, 'PlotLick', false, 'Reformat', 'RasterAndPETHAndWaveform');
	catch ME
		warning(['Error when processing iTr = ', num2str(iTr), ' - this one will be skipped.'])
		warning(sprintf('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message))
	end
end

