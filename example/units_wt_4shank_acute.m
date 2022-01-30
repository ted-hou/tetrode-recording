%% Select files
files = uipickfiles('Prompt', 'Select .mat files containing TetrodeRecording objects to load...', 'Type', {'*.mat', 'MAT-files'});

%% Generate batchPlotList and PETH from file (SLOW)
batchPlotList = {};

for iTr = 1:length(files)
% 	try
		disp(['Loading file: ', files{iTr}, '...'])
		S = load(files{iTr}, 'tr');
		tr = S.tr;

		expName = tr.GetExpName();
		animalName = strsplit(expName, '_');
		date = str2num(animalName{2});
		animalName = animalName{1};

		channels = [tr.Spikes.Channel];
		thisBatchPlotList = {};
		for iChn = channels
			for iUnit = 1:max(tr.Spikes(iChn).Cluster.Classes) - 1
				thisBatchPlotList = vertcat(thisBatchPlotList, {animalName, date, iChn, iChn, iUnit});
			end
		end
		batchPlotList = vertcat(batchPlotList, thisBatchPlotList);

		thisPETH = TetrodeRecording.BatchPETHistCounts(tr, thisBatchPlotList, 'TrialLength', 6, 'ExtendedWindow', 1, 'SpikeRateWindow', 100, 'Press', true, 'Lick', true, 'Stim', false);
		if iTr == 1
			PETH = thisPETH;
		else
			PETH = [PETH, thisPETH];
		end

		clear S expName animalName date iChn thisPETH thisBatchPlotList channels tr iTr iUnit
% 	catch ME
% 		warning(['Error when processing iTr = ', num2str(iTr), ' - this one will be skipped.'])
% 		warning(sprintf('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message))
% 	end
end


%%
load('C:\SERVER\PETH_WT_4shank_acute.mat')
scatter_press_vs_lick(PETH, 'Daisy12/13 Press vs Lick', 2, [-2, 0], true)
scatter_press_vs_lick(PETH, 'Daisy12/13 Press vs Lick', 2, [-0, 0.5], true)