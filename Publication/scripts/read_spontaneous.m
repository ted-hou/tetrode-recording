
%% 1.1. Load acute EU objects (duplicates already removed)
euSpontaneous = EphysUnit.load('C:\SERVER\Units\acute_spontaneous_reach'); 
%%
clear isi ISI prcLowISI

% Remove multiunit detected by ISI test.
p.ISIThreshold = 0.0015;
ISI = cell(length(iEu), 1);
for iEu = 1:length(euSpontaneous)
    st = euSpontaneous(iEu).SpikeTimes;
    isi = [NaN, diff(st)];
    st(isi == 0) = [];
    isi = [NaN, diff(st)];
    euSpontaneous(iEu).SpikeTimes = st;
    ISI{iEu} = isi;
end

prcLowISI = NaN(length(iEu), 1);
for iEu = 1:length(euSpontaneous)
    prcLowISI(iEu) = nnz(ISI{iEu} < p.ISIThreshold) ./ length(ISI{iEu});
end
histogram(prcLowISI, 0:0.01:1)
cSpontaneous.isMultiUnit = prcLowISI > 0.05;
cSpontaneous.isSingleUnit = prcLowISI <= 0.05;
euSpontaneous = euSpontaneous(cSpontaneous.isSingleUnit);

euSpontaneous = euSpontaneous';

% 1.2. Load Video Tracking Data (vtd) and ArduinoConnection (ac), and group into experiments
expSpontaneous = CompleteExperiment2(euSpontaneous);

% 1.3 Align video and ephys timestamps
expSpontaneous.alignTimestamps();

% Calculate ETA
etaSpontaneous = euSpontaneous.getETA('count', 'press', window=[-4, 0], alignTo='stop', includeInvalid=true, ...
        normalize=[-4, -2], MinTrialDuration=8);

%%
cSpontaneous.hasPress = arrayfun(@(e) nnz(e.getTrials('press').duration() >= p.minTrialDuration) >= p.minNumTrials, euSpontaneous)';
cSpontaneous.hasLick = arrayfun(@(e) nnz(e.getTrials('lick').duration() >= p.minTrialDuration) >= p.minNumTrials, euSpontaneous)';


%% Bootstrap movement response
p.bootAlpha = 0.01;
bootSpontaneous.press = struct('h', NaN(length(euSpontaneous), 1), 'muDiffCI', NaN(length(euSpontaneous), 2), 'muDiffObs', NaN(length(euSpontaneous), 1));
[bootSpontaneous.press.h(cSpontaneous.hasPress), bootSpontaneous.press.muDiffCI(cSpontaneous.hasPress, :), bootSpontaneous.press.muDiffObs(cSpontaneous.hasPress)] = bootstrapMoveResponse( ...
    euSpontaneous(cSpontaneous.hasPress), 'press', alpha=p.bootAlpha, withReplacement=false, oneSided=false, ...
    responseWindow=[-0.5, -0.2]);
fprintf(1, '\nAll done\n')

% Report bootstraped movement response direction
assert(nnz(isnan(bootSpontaneous.press.h(cSpontaneous.hasPress))) == 0)

figure, histogram(bootSpontaneous.press.h)
cSpontaneous.isPressUp = bootSpontaneous.press.h' == 1 & cSpontaneous.hasPress;
cSpontaneous.isPressDown = bootSpontaneous.press.h' == -1 & cSpontaneous.hasPress;
cSpontaneous.isPressResponsive = cSpontaneous.isPressUp | cSpontaneous.isPressDown;

fprintf(1, ['%g total SNr units (baseline spike rate > %g):\n' ...
    '\t%g with %d+ press trials;\n'], ...
    length(euSpontaneous), p.minSpikeRate, nnz(cSpontaneous.hasPress), p.minNumTrials)

fprintf(1, ['%g units with %g+ press trials (%gs or longer):\n' ...
    '\t%g (%.0f%%) are excited (p<%g);\n' ...
    '\t%g (%.0f%%) are inhibited (p<%g).\n'], ...
    nnz(cSpontaneous.hasPress), p.minNumTrials, p.minTrialDuration, ...
    nnz(cSpontaneous.isPressUp), 100*nnz(cSpontaneous.isPressUp)/nnz(cSpontaneous.isPressResponsive), p.bootAlpha, ...
    nnz(cSpontaneous.isPressDown), 100*nnz(cSpontaneous.isPressDown)/nnz(cSpontaneous.isPressResponsive), p.bootAlpha);

%% Calculate bta's, split by press up or down

[btaSpontaneous.pressUpRaw.X, btaSpontaneous.pressUpRaw.T, btaSpontaneous.pressUpRaw.N, btaSpontaneous.pressUpRaw.S, btaSpontaneous.pressUpRaw.B] = euSpontaneous(cSpontaneous.isPressUp).getBinnedTrialAverage('rate', [8, 16, 32, 64], 'press', 'window', [-4, 0], 'normalize', false);
[btaSpontaneous.pressDownRaw.X, btaSpontaneous.pressDownRaw.T, btaSpontaneous.pressDownRaw.N, btaSpontaneous.pressDownRaw.S, btaSpontaneous.pressDownRaw.B] = euSpontaneous(cSpontaneous.isPressDown).getBinnedTrialAverage('rate', [8, 16, 32, 64], 'press', 'window', [-4, 0], 'normalize', false);

