function scatter_stim_vs_press(PETHStim, titletext, sigma_threshold)% Sorted stim PETHStim
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
	[pethStimSorted, Istim, whenDidFiringRateChange] =  TetrodeRecording.SortPETH(pethStim, 'Method', 'latency', 'LatencyThreshold', 0.675);
	pethStimSortedNorm = TetrodeRecording.NormalizePETH(pethStimSorted, 'Method', 'zscore', 'BaselineSamples', t < 0);

	% Normalize rates
	pethStimNorm = TetrodeRecording.NormalizePETH(pethStim, 'Method', 'zscore', 'BaselineSamples', t < 0);


	tPress = PETHStim(1).Time;
	pethPress = transpose(reshape([PETHStim.Press], [length(tPress), length(PETHStim)]));
	pethPressNorm = TetrodeRecording.NormalizePETH(pethPress, 'Method', 'zscore', 'BaselineSamples', tPress < -2 & tPress > -4);

	% Find isDec/isInc
	for i = 1:length(PETHStim)
		PETHStim(i).IsDecStim = any(pethStimNorm(i, t>=0 & t <= 0.025) < -2);
		PETHStim(i).IsIncStim = any(pethStimNorm(i, t>=0 & t <= 0.025) > 2);
		PETHStim(i).IsDecPress = any(pethPressNorm(i, tPress>=-2 & tPress <= 0) < -3);
		PETHStim(i).IsIncPress = any(pethPressNorm(i, tPress>=-2 & tPress <= 0) > 3);

		normStim = pethStimNorm(i, t>=0 & t <= 0.02);
		[~, iMaxStimEffect] = max(abs(normStim));
		PETHStim(i).MaxStimEffect = normStim(iMaxStimEffect);

		normPress = pethPressNorm(i, tPress>=-2 & tPress <= 0);
		[~, iMaxPressEffect] = max(abs(normPress));
		PETHStim(i).MaxPressEffect = normPress(iMaxPressEffect);
	end

	maxPressEffect = [PETHStim.MaxPressEffect];
	maxStimEffect = [PETHStim.MaxStimEffect];
	sel = abs(maxPressEffect) > sigma_threshold & abs(maxStimEffect) > sigma_threshold;
	f = figure;
	ax = axes(f);
	hold on
	hLine1 = plot(ax, maxPressEffect(~sel), maxStimEffect(~sel), 'o', 'MarkerEdgeColor', '#808080', 'DisplayName', sprintf('%d unresponsive units', sum(~sel)));
	hLine2 = plot(ax, maxPressEffect(sel), maxStimEffect(sel), 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'black', 'DisplayName', sprintf('%d responsive units', sum(sel)));
	plot(ax, ax.XLim, [0, 0], 'k:')
	plot(ax, [0, 0], ax.YLim, 'k:')
	hold off
	xlabel('Peak MOVE effect (\sigma) ([-2s, 0])', 'FontSize', 14)
	ylabel('Peak STIM effect (\sigma) ([0, 25ms]', 'FontSize', 14)
	title(sprintf('%s', titletext), 'FontSize', 14)
	legend([hLine1, hLine2], 'FontSize', 14)
