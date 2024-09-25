
%% 1.1. Load acute EU objects (duplicates already removed)
euSpontaneous = EphysUnit.load('C:\SERVER\Units\acute_spontaneous_reach');
% Remove drifting units
cSpontaneous.isDrifting = detectDriftingUnits(euSpontaneous, smoothWindow=300, tolerance=0.05, spikeRateThreshold=15);
euSpontaneous = euSpontaneous(~cSpontaneous.isDrifting);
euSpontaneous = euSpontaneous.removeDuplicates();
%%
clear isi ISI prcLowISI

% Remove multiunit detected by ISI test.
p.ISIThreshold = 0.0015;
ISI = cell(length(euSpontaneous), 1);
for iEu = 1:length(euSpontaneous)
    st = euSpontaneous(iEu).SpikeTimes;
    isi = [NaN, diff(st)];
    st(isi == 0) = [];
    isi = [NaN, diff(st)];
    euSpontaneous(iEu).SpikeTimes = st;
    ISI{iEu} = isi;
end

prcLowISI = NaN(length(euSpontaneous), 1);
for iEu = 1:length(euSpontaneous)
    prcLowISI(iEu) = nnz(ISI{iEu} < p.ISIThreshold) ./ length(ISI{iEu});
end
histogram(prcLowISI, 0:0.01:1)
cSpontaneous.isMultiUnit = prcLowISI > 0.05;
cSpontaneous.isSingleUnit = prcLowISI <= 0.05;
euSpontaneous = euSpontaneous(cSpontaneous.isSingleUnit);

euSpontaneous = euSpontaneous';
%%
% 1.2. Load Video Tracking Data (vtd) and ArduinoConnection (ac), and group into experiments
expSpontaneous = CompleteExperiment2(euSpontaneous);

% 1.3 Align video and ephys timestamps
expSpontaneous.alignTimestamps();

% Calculate ETA
etaSpontaneous = euSpontaneous.getETA('count', 'press', window=[-4, 2], alignTo='stop', includeInvalid=true, ...
        normalize=[-4, -2], MinTrialDuration=6);

etaSpontaneousRaw = euSpontaneous.getETA('count', 'press', window=[-4, 2], alignTo='stop', includeInvalid=true, ...
        normalize='none', MinTrialDuration=6);

%%
cSpontaneous.hasPress = arrayfun(@(e) nnz(e.getTrials('press').duration() >= 6) >= p.minNumTrials, euSpontaneous)';
cSpontaneous.hasLick = arrayfun(@(e) nnz(e.getTrials('lick').duration() >= 6) >= p.minNumTrials, euSpontaneous)';


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
cSpontaneous.isPressUp = bootSpontaneous.press.h == 1 & cSpontaneous.hasPress';
cSpontaneous.isPressDown = bootSpontaneous.press.h == -1 & cSpontaneous.hasPress';
cSpontaneous.isPressResponsive = cSpontaneous.isPressUp | cSpontaneous.isPressDown;

fprintf(1, ['%g total SNr units (baseline spike rate > %g):\n' ...
    '\t%g with %d+ press trials;\n'], ...
    length(euSpontaneous), p.minSpikeRate, nnz(cSpontaneous.hasPress), p.minNumTrials)

fprintf(1, ['%g units with %g+ press trials (6s or longer):\n' ...
    '\t%g (%.0f%%) are excited (p<%g);\n' ...
    '\t%g (%.0f%%) are inhibited (p<%g).\n'], ...
    nnz(cSpontaneous.hasPress), p.minNumTrials, ...
    nnz(cSpontaneous.isPressUp), 100*nnz(cSpontaneous.isPressUp)./nnz(cSpontaneous.isPressResponsive), p.bootAlpha, ...
    nnz(cSpontaneous.isPressDown), 100*nnz(cSpontaneous.isPressDown)./nnz(cSpontaneous.isPressResponsive), p.bootAlpha);

%% Get paw trajectories and movement onset time

% Correct feature names
for iExp = 1:length(expSpontaneous)
    switch expReachDir(iExp).animalName
        case {'desmond28', 'desmond30'}
            expReachDir(iExp).vtdL = renamevars(expReachDir(iExp).vtdL, ...
                {'handIpsiCam_X', 'handIpsiCam_Y', 'handIpsiCam_Likelihood', 'footIpsiCam_X', 'footIpsiCam_Y', 'footIpsiCam_Likelihood'}, ...
                {'handContra_X', 'handContra_Y', 'handContra_Likelihood', 'footContra_X', 'footContra_Y', 'footContra_Likelihood'});
            expReachDir(iExp).vtdR = renamevars(expReachDir(iExp).vtdR, ...
                {'handIpsiCam_X', 'handIpsiCam_Y', 'handIpsiCam_Likelihood', 'footIpsiCam_X', 'footIpsiCam_Y', 'footIpsiCam_Likelihood'}, ...
                {'handIpsi_X', 'handIpsi_Y', 'handIpsi_Likelihood', 'footIpsi_X', 'footIpsi_Y', 'footIpsi_Likelihood'});
            posOrder(iExp, :) = [4, 3, 2, 1];

        case 'desmond29'
            expReachDir(iExp).vtdL = renamevars(expReachDir(iExp).vtdL, ...
                {'handIpsiCam_X', 'handIpsiCam_Y', 'handIpsiCam_Likelihood', 'footIpsiCam_X', 'footIpsiCam_Y', 'footIpsiCam_Likelihood'}, ...
                {'handIpsi_X', 'handIpsi_Y', 'handIpsi_Likelihood', 'footIpsi_X', 'footIpsi_Y', 'footIpsi_Likelihood'});
            expReachDir(iExp).vtdR = renamevars(expReachDir(iExp).vtdR, ...
                {'handIpsiCam_X', 'handIpsiCam_Y', 'handIpsiCam_Likelihood', 'footIpsiCam_X', 'footIpsiCam_Y', 'footIpsiCam_Likelihood'}, ...
                {'handContra_X', 'handContra_Y', 'handContra_Likelihood', 'footContra_X', 'footContra_Y', 'footContra_Likelihood'});
            posOrder(iExp, :) = [1, 2, 3, 4];            
    end
end

%% Save results (only once)
% euSpontaneous.save('C:\SERVER\Units\acute_spontaneous_reach\SNr_SingleUnit_NonDuplicate_NonDrift')
% if ~exist('C:\SERVER\Units\acute_spontaneous_reach\meta', 'dir')
%     mkdir('C:\SERVER\Units\acute_spontaneous_reach\meta')
% end
% save('C:\SERVER\Units\acute_spontaneous_reach\meta\SNr_SingleUnit_NonDuplicate_NonDrift.mat', 'bootSpontaneous', 'cSpontaneous', 'etaSpontaneous', 'etaSpontaneousRaw', '-v7.3')

%% Calculate bta's, split by press up or down
% 
% [btaSpontaneous.pressUpRaw.X, btaSpontaneous.pressUpRaw.T, btaSpontaneous.pressUpRaw.N, btaSpontaneous.pressUpRaw.S, btaSpontaneous.pressUpRaw.B] = euSpontaneous(cSpontaneous.isPressUp).getBinnedTrialAverage('rate', [8, 16, 32, 64], 'press', 'window', [-4, 0], 'normalize', false);
% [btaSpontaneous.pressDownRaw.X, btaSpontaneous.pressDownRaw.T, btaSpontaneous.pressDownRaw.N, btaSpontaneous.pressDownRaw.S, btaSpontaneous.pressDownRaw.B] = euSpontaneous(cSpontaneous.isPressDown).getBinnedTrialAverage('rate', [8, 16, 32, 64], 'press', 'window', [-4, 0], 'normalize', false);
% 
