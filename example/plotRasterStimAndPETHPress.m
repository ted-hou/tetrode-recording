

bpl = {...
	'daisy4', 20190409, 4, 1;...
	'daisy4', 20190424, 6, 1;...
	'daisy4', 20190417, 1, 1;...
	'daisy4', 20190508, 8, 1;...
	'daisy4', 20190424, 5, 1;...
	'daisy4', 20190514, 6, 1 ...
};
expNames = cell(size(bpl, 1), 1);
for iExp = 1:size(bpl, 1)
	expNames{iExp} = [bpl{iExp, 1}, '_', num2str(bpl{iExp, 2})];
end

expNamesUnique = unique(expNames);

for iTr = 1:length(expNamesUnique)
	tr = TetrodeRecording.BatchLoad(expNamesUnique(iTr));
	try
		TetrodeRecording.BatchPlot(tr, bpl, 'PlotStim', true, 'copyLegend', true, 'Reformat', 'PETHAndRasterStimAndWaveform', 'ExtendedWindow', [-0.5, 1]);
	catch ME
		warning(['Error when processing iTr = ', num2str(iTr), ' - this one will be skipped.'])
		warning(sprintf('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message))
	end
end


% TetrodeRecording.BatchPlot(tr, bpl, 'PlotStim', true, 'copyLegend', true, 'Reformat', 'PETHAndRasterStimAndWaveform', 'ExtendedWindow', [-0.5, 0]);