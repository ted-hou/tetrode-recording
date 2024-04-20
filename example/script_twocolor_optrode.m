 tr = TetrodeRecording.BatchLoadSimple();
%% Load two color experiment and validate stim pulse timing
file = dir(sprintf('%s\\..\\%s.mat', tr.Path, tr.GetExpName(includeSuffix=false)));
assert(~isempty(file))

tce = load(sprintf('%s\\%s', file.folder, file.name));
tce = tce.obj;
assert(isa(tce, 'TwoColorExperiment'))
clear file

tOn = tr.DigitalEvents.StimOn;
tOff = tr.DigitalEvents.StimOff;

nPulsesPerTrain = arrayfun(@(log) log.params.nPulses, tce.Log);
pulseWidth = arrayfun(@(log) repmat(log.params.pulseWidth, [1, log.params.nPulses]), tce.Log, UniformOutput=false);
pulseWidth = cat(2, pulseWidth{:});

% Verify congruence between Intan shutter events and TwoColorExperiment.Log
assert(sum(nPulsesPerTrain) == nnz(tOn), 'Intan recorded %i shutterOn events while TwoColorExperiment.Log has %i.', nnz(tOn), sum(nPulsesPerTrain));
assert(sum(nPulsesPerTrain) == nnz(tOff), 'Intan recorded %i shutterOn events while TwoColorExperiment.Log has %i.', nnz(tOn), sum(nPulsesPerTrain));
nPulses = nnz(tOn);
fprintf('Intan & TwoColorExperiment.Log both recorded %i shutter pulses.\n', nnz(tOn))

% Pulse durations should match between Intan and TCE as well.
pulseWidthErrorMargin = 1e-3;
assert(all(abs(tOff - tOn - pulseWidth) < pulseWidthErrorMargin), 'Only %i/%i pulseWidths agree (df<%gs).', nnz(abs(tOff - tOn - pulseWidth) < pulseWidthErrorMargin), nPulses, pulseWidthErrorMargin)
fprintf('All %i pulseWidths agree (df<%gs).\n', nPulses, pulseWidthErrorMargin)

% Generate pulse->train map
iPulse = 0;
pulseToTrain = NaN(1, nPulses);
for iTrain = 1:length(tce.Log)
    pulseToTrain(iPulse + 1:iPulse + nPulsesPerTrain(iTrain)) = iTrain;
    iPulse = iPulse + nPulsesPerTrain(iTrain);
end


clear iPulse iTrain pulseWidth


% Make one raster 
% Stim condition has
% 1000000*wavelength (1, 2) + 100000*location (1, 2, 3, 4) + 1000*duration (1, 2, 4, 6, 10, 20, 40) + 1*nominalpower (1, 2, 4, 20, 80, 320)

