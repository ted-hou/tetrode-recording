
numBins = 10;
numSigmas = 1.5;


load PETH
% Select units with stim effects
for iPETH = 1:length(PETHStimD1)
	Stim = PETHStimD1(iPETH).Stim;
	Stim = 	Stim([Stim.TrainType] == 1910);
	inWindowBaseline = Stim.Timestamps < 0 & Stim.Timestamps >= -1;
	inWindowPulse = arrayfun(@(x,y) Stim.Timestamps > x & Stim.Timestamps <= y, Stim.PulseOn, Stim.PulseOff, 'UniformOutput', false);
	inWindowPulse = logical(sum(cell2mat(inWindowPulse')));
	inWindowTrain =  Stim.Timestamps >= 0 & Stim.Timestamps < 2;

	spikeRateBaseline = Stim.SpikeRate(inWindowBaseline);
	spikeRatePulse = Stim.SpikeRate(inWindowPulse);
	spikeRateTrain = Stim.SpikeRate(inWindowTrain);

	meanBaseline = mean(spikeRateBaseline);
	meanPulse = mean(spikeRatePulse);
	meanTrain = mean(spikeRateTrain);

	stdBaseline = std(spikeRateBaseline);
	stdPulse = std(spikeRatePulse);
	stdTrain = std(spikeRateTrain);

	timestamps = Stim.Timestamps(inWindowTrain);
	turningPoint = find(abs(spikeRateTrain - meanBaseline) >= numSigmas*stdBaseline, 1);
	if isempty(turningPoint)
		latency(iPETH) = timestamps(end);
	else
		latency(iPETH) = timestamps(turningPoint);
	end

	% h(iPETH) = logical(ttest2(spikeRateBaseline, spikeRatePulse));
	h(iPETH) = abs(meanPulse - meanBaseline) > numSigmas * stdBaseline;
	diffMean(iPETH) = meanPulse - meanBaseline;
end

TetrodeRecording.HeatMapStim(PETHStimD1(h), 'SortBy', sign(diffMean(h)).*latency(h), 'Window', [-4, 0.5], 'NormalizationBaselineWindow', [-6, 0], 'WindowStim', [-0.5, -1], 'NormalizationBaselineWindowStim', [-1, 0], 'CLimStim', [-3,3], 'StimTrainTypes', [1910], 'Sorting', 'latency', 'LatencyThreshold', numSigmas);	
I = TetrodeRecording.HeatMapStim(PETHStimD1(h), 'SortBy', sign(diffMean(h)).*latency(h), 'Window', [-4, 0.5], 'NormalizationBaselineWindow', [-6, 0], 'WindowStim', [-0.5, 0.5], 'NormalizationBaselineWindowStim', [-1, 0], 'CLimStim', [-3,3], 'StimTrainTypes', [1910, 1604, 4002], 'Sorting', 'latency', 'LatencyThreshold', numSigmas);	
% TetrodeRecording.HeatMapStim(PETHStimD1(h), 'SortBy', [], 'LatencyThreshold', 6, 'Window', [-4, 0.5], 'NormalizationBaselineWindow', [-6, 0], 'WindowStim', [-0.5, 0.5], 'NormalizationBaselineWindowStim', [-1, 0], 'CLimStim', [-3,3], 'StimTrainTypes', [1910, 1604, 4002], 'Sorting', 'latency');	


PETHStimD1H = PETHStimD1(h);
names = cellfun(@(x) x{1}, cellfun(@(x) strsplit(x, '_'), {PETHStimD1H(I).ExpName}', 'UniformOutput', false), 'UniformOutput', false);
dates = cellfun(@(x) str2num(x{2}), cellfun(@(x) strsplit(x, '_'), {PETHStimD1H(I).ExpName}', 'UniformOutput', false), 'UniformOutput', false);
bpl = [names, dates, {PETHStimD1H(I).Channel}', {PETHStimD1H(I).Cluster}'];


expNames = cell(length(bpl), 1);
for iExp = 1:length(bpl)
	expNames{iExp} = [bpl{iExp, 1}, '_', num2str(bpl{iExp, 2})];
end
for iTr = 1:length(expNames)
	tr = TetrodeRecording.BatchLoad(expNames(iTr));
	try
		TetrodeRecording.BatchPlot(tr, bpl, 'PlotStim', true, 'Reformat', 'PETHAndRasterStimAndWaveform');
	catch ME
		warning(['Error when processing iTr = ', num2str(iTr), ' - this one will be skipped.'])
		warning(sprintf('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message))
	end
end


lat = latency(h);
edges = linspace(min(lat), max(lat), numBins);
% histLatUp = 100*histcounts(lat(isUp), edges)/nnz(isUp);
% histLatDown = 100*histcounts(lat(~isUp), edges)/nnz(~isUp);
% histLatUp = histcounts(lat(isUp), edges);
% histLatDown = histcounts(lat(~isUp), edges);
% centers = (edges(1:end-1) + edges(2:end))/2;


f = figure;
set(f, 'DefaultAxesFontSize', 18)
ax = axes(f);
histogram(ax, lat*1000);
xlabel(ax, 'Stim Latency (ms)')
ylabel(ax, 'Occurance')

% h2 = ttest2(lat(isUp), lat(~isUp));

% edges = linspace(min(lat), max(lat), numBins);
% % histLatUp = 100*histcounts(lat(isUp), edges)/nnz(isUp);
% % histLatDown = 100*histcounts(lat(~isUp), edges)/nnz(~isUp);
% histLatUp = histcounts(lat(isUp), edges);
% histLatDown = histcounts(lat(~isUp), edges);
% centers = (edges(1:end-1) + edges(2:end))/2;

% f = figure;
% set(f, 'DefaultAxesFontSize', 18)
% ax = axes(f);
% hold(ax, 'on')
% plot(ax, centers*1000, histLatUp, 'r', 'DisplayName', 'Excited', 'LineWidth', 2);
% plot(ax, centers*1000, histLatDown, 'g', 'DisplayName', 'Inhibited', 'LineWidth', 2);
% hold (ax, 'off')
% legend(ax)
% xlabel(ax, 'Stim Latency (ms)')
% % ylabel(ax, 'Frequency (%)')
% ylabel(ax, 'Occurance')


numSigmas = 0.1;

% Select units with stim effects
for iPETH = 1:length(PETHStimDAT)
	Stim = PETHStimDAT(iPETH).Stim;
	Stim = 	Stim([Stim.TrainType] == 1910);
	inWindowBaseline = Stim.Timestamps < 0 & Stim.Timestamps >= -1;
	inWindowPulse = arrayfun(@(x,y) Stim.Timestamps > x & Stim.Timestamps <= y, Stim.PulseOn, Stim.PulseOff, 'UniformOutput', false);
	inWindowPulse = logical(sum(cell2mat(inWindowPulse')));
	inWindowTrain =  Stim.Timestamps >= 0 & Stim.Timestamps < 2;

	spikeRateBaseline = Stim.SpikeRate(inWindowBaseline);
	spikeRatePulse = Stim.SpikeRate(inWindowPulse);
	spikeRateTrain = Stim.SpikeRate(inWindowTrain);

	meanBaseline = mean(spikeRateBaseline);
	meanPulse = mean(spikeRatePulse);
	meanTrain = mean(spikeRateTrain);

	stdBaseline = std(spikeRateBaseline);
	stdPulse = std(spikeRatePulse);
	stdTrain = std(spikeRateTrain);

	timestamps = Stim.Timestamps(inWindowTrain);
	turningPoint = find(abs(spikeRateTrain - meanBaseline) >= numSigmas*stdBaseline, 1);
	if isempty(turningPoint)
		latency(iPETH) = timestamps(end);
	else
		latency(iPETH) = timestamps(turningPoint);
	end

	% h(iPETH) = logical(ttest2(spikeRateBaseline, spikeRatePulse));
	h(iPETH) = abs(meanPulse - meanBaseline) > numSigmas * stdBaseline;
	diffMean(iPETH) = meanPulse - meanBaseline;
end

I = TetrodeRecording.HeatMapStim(PETHStimDAT(h), 'SortBy', sign(diffMean(h)).*latency(h), 'Window', [-4, 0.5], 'NormalizationBaselineWindow', [-6, 0], 'WindowStim', [-0.5, 0.5], 'NormalizationBaselineWindowStim', [-1, 0], 'CLimStim', [-3,3], 'StimTrainTypes', [1910], 'Sorting', 'latency');	


lat = latency(h);
isUp = diffMean(h) > 0;


% Average isUp and isDown traces separately
PETHStimDATSig = PETHStimDAT(h);

PETHStimDATUp = PETHStimDAT(isUp);

for iPETH = 1:length(PETHStimDATUp)
	Stim = PETHStimDATUp(iPETH).Stim;
	Stim = Stim([Stim.TrainType] == 1910);
	spikeRate(iPETH, 1:length(Stim.SpikeRate)) = Stim.SpikeRate;
	timestamps = Stim.Timestamps;
end

figure, plot(timestamps, mean(spikeRate, 1));
xlabel('Time')
ylabel('Firing rate, population average, increasing (sp/s)')


PETHStimDATDown = PETHStimDAT(~isUp);

for iPETH = 1:length(PETHStimDATDown)
	Stim = PETHStimDATDown(iPETH).Stim;
	Stim = Stim([Stim.TrainType] == 1910);
	spikeRate(iPETH, 1:length(Stim.SpikeRate)) = Stim.SpikeRate;
	timestamps = Stim.Timestamps;
end

figure, plot(timestamps, mean(spikeRate, 1));
xlabel('Time')
ylabel('Firing rate, population average, decreasing (sp/s)')





h2 = ttest2(lat(isUp), lat(~isUp));

edges = linspace(min(lat), max(lat), numBins);
% histLatUp = 100*histcounts(lat(isUp), edges)/nnz(isUp);
% histLatDown = 100*histcounts(lat(~isUp), edges)/nnz(~isUp);
histLatUp = histcounts(lat(isUp), edges);
histLatDown = histcounts(lat(~isUp), edges);
centers = (edges(1:end-1) + edges(2:end))/2;

f = figure;
set(f, 'DefaultAxesFontSize', 18)
ax = axes(f);
hold(ax, 'on')
plot(ax, centers*1000, histLatUp, 'r', 'DisplayName', 'Excited', 'LineWidth', 2);
plot(ax, centers*1000, histLatDown, 'g', 'DisplayName', 'Inhibited', 'LineWidth', 2);
hold (ax, 'off')
legend(ax)
xlabel(ax, 'Stim Latency (ms)')
% ylabel(ax, 'Frequency (%)')
ylabel(ax, 'Occurance')




figure; 
hold on;
plot(P(1).Stim.Timestamps, P(1).Stim.SpikeRate)
plot(P(2).Stim(2).Timestamps, P(2).Stim(2).SpikeRate)
plot(P(3).Stim(2).Timestamps, P(3).Stim(2).SpikeRate)
plot(P(4).Stim(2).Timestamps, P(4).Stim(2).SpikeRate)
hold off;
xlim([-1, 2.5])
xlabel('Time to stim (s)')
ylabel('Firing rate (sp/s)')