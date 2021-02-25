PETHStim = PETHStimA2A;
clear pethStim pethPress
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
[pethStimSorted, Istim, whenDidFiringRateChange] =  TetrodeRecording.SortPETH(pethStim, 'Method', 'latency', 'LatencyThreshold', 0.675);
pethStimSortedNorm = TetrodeRecording.NormalizePETH(pethStimSorted, 'Method', 'zscore', 'BaselineSamples', t < 0);


figure()
axes()
hold on
for i = 1:size(pethStimSortedNorm, 1)
	sel = t >= -0.2;
	isDec = any(pethStimSortedNorm(i, t>=0 & t < 0.03) < -2);
    isInc = any(pethStimSortedNorm(i, t>=0 & t < 0.03) > 2);
	if isDec
		plot(t(sel)*1000, pethStimSortedNorm(i, sel), 'color', [.2, .8, .2, 0.67])
    elseif isInc
		plot(t(sel)*1000, pethStimSortedNorm(i, sel), 'color', [.8, .2, .2, 0.67])
    else
		plot(t(sel)*1000, pethStimSortedNorm(i, sel), 'color', [.2, .2, .2, 0.67])
    end
end
plot([0, 10], [55, 55], 'b', 'LineWidth', 3)
hold off
title(char(sprintf("Adora2a-ChR2 stim response of %d SNr units", size(pethStim, 1))))
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
	PETHStim(i).IsDecStim = any(pethStimNorm(i, t>=0 & t <= 0.02) < -4);
	PETHStim(i).IsIncStim = any(pethStimNorm(i, t>=0 & t <= 0.02) > 4);
	PETHStim(i).IsDecPress = any(pethPressNorm(i, tPress>=-2 & tPress <= 0) < -4);
	PETHStim(i).IsIncPress = any(pethPressNorm(i, tPress>=-2 & tPress <= 0) > 4);
end


% sum([PETHStim.IsDecPress])
% sum([PETHStim.IsIncPress])

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
contTable(4, 1:3) = sum(contTable, 1);
contTable(1:3, 4) = sum(contTable(1:3, 1:3), 2);