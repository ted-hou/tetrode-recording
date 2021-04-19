PETHStim = PETHStimD1;
clear pethStim
% Sorted stim PETHStim
for i = 1:length(PETHStim)
	Stim = PETHStim(i).Stim;
	sel = Stim(1).Timestamps <= 0.2;
	t = Stim(1).Timestamps(sel);
	thisPeth = zeros(1, sum(sel));
	for j = 1:length(Stim)
		thisPeth = thisPeth + Stim(j).SpikeRate(sel) * Stim(j).NumTrains;
	end
	pethStim(i, 1:length(t)) = thisPeth ./ sum([Stim.NumTrains]);
end
clear i sel Stim
%%
[pethSorted, Istim, whenDidFiringRateChange] =  TetrodeRecording.SortPETH(pethStim, 'Method', 'latency', 'LatencyThreshold', 0.675);
pethSortedNorm = TetrodeRecording.NormalizePETH(pethSorted, 'Method', 'zscore', 'BaselineSamples', t < 0);


figure()
axes()
hold on
for i = 1:size(pethSortedNorm, 1)
	sel = t >= -0.2;
	isDec = any(pethSortedNorm(i, t>=0 & t < 0.03) < -2);
    isInc = any(pethSortedNorm(i, t>=0 & t < 0.03) > 2);
	if isDec
		plot(t(sel)*1000, pethSortedNorm(i, sel), 'color', [.2, .8, .2, 0.67])
    elseif isInc
		plot(t(sel)*1000, pethSortedNorm(i, sel), 'color', [.8, .2, .2, 0.67])
    else
		plot(t(sel)*1000, pethSortedNorm(i, sel), 'color', [.2, .2, .2, 0.67])
    end
end
plot([0, 100], [55, 55], 'b', 'LineWidth', 3)
hold off
title(char(sprintf("Drd1-ChR2 stim response of %d SNr units", size(pethStim, 1))))
xlabel("Time (ms)")
ylabel("Spike rate (z-score)")
xlim([-100, 100])

%% Contingency table
% Normalize rates
pethStimNorm = TetrodeRecording.NormalizePETH(pethStim, 'Method', 'zscore', 'BaselineSamples', t < 0);


tPress = PETHStim(1).Time;
pethPress = transpose(reshape([PETHStim.Press], [length(tPress), length(PETHStim)]));
pethPressNorm = TetrodeRecording.NormalizePETH(pethPress, 'Method', 'zscore', 'BaselineSamples', tPress < -2 & tPress > -4);

% Find isDec/isInc
for i = 1:length(PETHStim)
	PETHStim(i).IsDecStim = any(pethStimNorm(i, t>=0 & t <= 0.02) < -2);
	PETHStim(i).IsIncStim = any(pethStimNorm(i, t>=0 & t <= 0.02) > 2);
	PETHStim(i).IsDecPress = any(pethPressNorm(i, tPress>=-2 & tPress <= 0) < -4);
	PETHStim(i).IsIncPress = any(pethPressNorm(i, tPress>=-2 & tPress <= 0) > 4);
end


sum([PETHStim.IsDecPress])
sum([PETHStim.IsIncPress])

sum([PETHStim.IsDecStim])
sum([PETHStim.IsIncStim])

contTable = zeros(3);
contTable(1, 1) = sum([PETHStim.IsDecPress] & [PETHStim.IsDecStim]);
contTable(1, 2) = sum([PETHStim.IsDecPress] & [PETHStim.IsIncStim]);
contTable(1, 3) = sum([PETHStim.IsDecPress] & ~([PETHStim.IsDecStim] | [PETHStim.IsIncStim]));
contTable(2, 1) = sum([PETHStim.IsIncPress] & [PETHStim.IsDecStim]);
contTable(2, 2) = sum([PETHStim.IsIncPress] & [PETHStim.IsIncStim]);
contTable(2, 3) = sum([PETHStim.IsIncPress] & ~([PETHStim.IsDecStim] | [PETHStim.IsIncStim]));
contTable(3, 1) = sum(~([PETHStim.IsDecPress] | [PETHStim.IsIncPress]) & [PETHStim.IsDecStim]);
contTable(3, 2) = sum(~([PETHStim.IsDecPress] | [PETHStim.IsIncPress]) & [PETHStim.IsIncStim]);
contTable(3, 3) = sum(~([PETHStim.IsDecPress] | [PETHStim.IsIncPress]) & ~([PETHStim.IsDecStim] | [PETHStim.IsIncStim]));





% Desmond 21

batchPlotListStim = {...
	'desmond21', 20210406, 2, 1;... % Up, stim flat
	'desmond21', 20210406, 4, 1;... % Down, stim flat
	'desmond21', 20210406, 5, 1;... % Down maybe, stim flat
	'desmond21', 20210406, 6, 1;... % Up, stim flat
	% 'desmond21', 20210406, 7, 1;... % DA
	'desmond21', 20210406, 9, 1;... % Up, stim flat
	'desmond21', 20210407, 1, 1;... % Up, stim flat
	'desmond21', 20210407, 2, 1;... % Down, stim flat
	'desmond21', 20210407, 6, 1;... % Up, stim flat
	'desmond21', 20210407, 7, 1;... % Up, stim flat
	'desmond21', 20210407, 8, 1;... % Up, stim flat
	'desmond21', 20210407, 9, 1;... % DA

	'desmond21', 20210408, 2, 1;... % Down, stim flat
	'desmond21', 20210408, 7, 1;... % Up, stim flat
	'desmond21', 20210409, 1, 1;... % Up, stim flat
	'desmond21', 20210409, 2, 1;... % Down, stim off maybe?
	'desmond21', 20210409, 5, 1;... % Up, stim flat
	'desmond21', 20210411, 1, 1;... % Up, stim flat
	'desmond21', 20210411, 3, 1;... % Up, stim flat
	'desmond21', 20210411, 5, 1;... % Up, stim flat
	};


expNames = cell(length(batchPlotListStim), 1);
for iExp = 1:length(batchPlotListStim)
	expNames{iExp} = [batchPlotListStim{iExp, 1}, '_', num2str(batchPlotListStim{iExp, 2})];
end

expNamesUnique = unique(expNames);

for iTr = 1:length(expNamesUnique)
	tr = TetrodeRecording.BatchLoad(expNamesUnique(iTr));
	try
		if iTr == 1;
			PETH = TetrodeRecording.BatchPETHistCounts(tr, batchPlotListStim, 'TrialLength', 6, 'ExtendedWindow', 1, 'SpikeRateWindow', 100, 'ExtendedWindowStim', [-1, 1], 'SpikeRateWindowStim', 10, 'Press', true, 'Lick', false, 'Stim', true);
		else
			PETH = [PETH, TetrodeRecording.BatchPETHistCounts(tr, batchPlotListStim, 'TrialLength', 6, 'ExtendedWindow', 1, 'SpikeRateWindow', 100, 'ExtendedWindowStim', [-1, 1], 'SpikeRateWindowStim', 10, 'Press', true, 'Lick', false, 'Stim', true)];
		end
	catch ME
		warning(['Error when processing iTr = ', num2str(iTr), ' - this one will be skipped.'])
		warning(sprintf('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message))
	end
end

PETHStim = PETH;

clear pethStim
% Sorted stim PETHStim
for i = 1:length(PETHStim)
	Stim = PETHStim(i).Stim;
	sel = Stim(1).Timestamps <= 0.2;
	t = Stim(1).Timestamps(sel);
	thisPeth = zeros(1, sum(sel));
	for j = 1:length(Stim)
		thisPeth = thisPeth + Stim(j).SpikeRate(sel) * Stim(j).NumTrains;
	end
	pethStim(i, 1:length(t)) = thisPeth ./ sum([Stim.NumTrains]);
end
clear i sel Stim
%%
[pethSorted, Istim, whenDidFiringRateChange] =  TetrodeRecording.SortPETH(pethStim, 'Method', 'latency', 'LatencyThreshold', 0.675);
pethSortedNorm = TetrodeRecording.NormalizePETH(pethSorted, 'Method', 'zscore', 'BaselineSamples', t < 0);



figure()
axes()
hold on
for i = 1:size(pethSortedNorm, 1)
	sel = t >= -0.2;
	isDec = any(pethSortedNorm(i, t>0 & t <= 0.1) < -2.5);
    isInc = any(pethSortedNorm(i, t>0 & t <= 0.1) > 2.5);
	if isDec
		h = plot(t(sel)*1000, pethSortedNorm(i, sel), 'color', [.2, .8, .2, 0.67]);
    elseif isInc
		h = plot(t(sel)*1000, pethSortedNorm(i, sel), 'color', [.8, .2, .2, 0.67]);
    else
		h = plot(t(sel)*1000, pethSortedNorm(i, sel), 'color', [.2, .2, .2, 0.67]);
    end
    h.DisplayName = num2str(i);
end
plot([0, 20], [0, 0], 'b', 'LineWidth', 3)
hold off
title(char(sprintf("Drd1-ChR2 stim response of %d SNr units", size(pethStim, 1))))
xlabel("Time (ms)")
ylabel("Spike rate (z-score)")
xlim([-100, 100])

figure()
axes()
hold on
for i = 1:size(pethSortedNorm, 1)
	sel = t >= -0.2;
	isDec = any(pethSortedNorm(i, t>0 & t <= 0.1) < -2.5);
    isInc = any(pethSortedNorm(i, t>0 & t <= 0.1) > 2.5);
	if isDec
		h = plot(t(sel)*1000, pethSorted(i, sel), 'color', [.2, .8, .2, 0.67]);
    elseif isInc
		h = plot(t(sel)*1000, pethSorted(i, sel), 'color', [.8, .2, .2, 0.67]);
    else
		h = plot(t(sel)*1000, pethSorted(i, sel), 'color', [.2, .2, .2, 0.67]);
    end
    h.DisplayName = num2str(i);
end
plot([0, 20], [0, 0], 'b', 'LineWidth', 3)
hold off
title(char(sprintf("Drd1-ChR2 stim response of %d SNr units", size(pethStim, 1))))
xlabel("Time (ms)")
ylabel("Spike rate (z-score)")
xlim([-100, 100])
