
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

%% Correct feature names
for iExp = 1:length(expSpontaneous)
    switch expSpontaneous(iExp).animalName
        case 'daisy18'
            expSpontaneous(iExp).vtdL = renamevars(expSpontaneous(iExp).vtdL, ...
                {'handIpsiCam_X', 'handIpsiCam_Y', 'handIpsiCam_Likelihood'}, ...
                {'handContra_X', 'handContra_Y', 'handContra_Likelihood'});
            expSpontaneous(iExp).vtdR = renamevars(expSpontaneous(iExp).vtdR, ...
                {'handIpsiCam_X', 'handIpsiCam_Y', 'handIpsiCam_Likelihood'}, ...
                {'handIpsi_X', 'handIpsi_Y', 'handIpsi_Likelihood'});

        case {'daisy17', 'desmond31'}
            expSpontaneous(iExp).vtdL = renamevars(expSpontaneous(iExp).vtdL, ...
                {'handIpsiCam_X', 'handIpsiCam_Y', 'handIpsiCam_Likelihood'}, ...
                {'handIpsi_X', 'handIpsi_Y', 'handIpsi_Likelihood'});
            expSpontaneous(iExp).vtdR = renamevars(expSpontaneous(iExp).vtdR, ...
                {'handIpsiCam_X', 'handIpsiCam_Y', 'handIpsiCam_Likelihood'}, ...
                {'handContra_X', 'handContra_Y', 'handContra_Likelihood'});
        otherwise
            error();
    end
end

%% Extract trajectories of forepaws
% Whole session paw position
likelihoodThreshold = 0;
clear pos
pos(length(expSpontaneous)) = struct(contra=[], ipsi=[]);
for iExp = 1:length(expSpontaneous)
    switch expSpontaneous(iExp).animalName
        case 'daisy18'
            vtdContra = expSpontaneous(iExp).vtdL;
            vtdIpsi = expSpontaneous(iExp).vtdR;
        case {'daisy17', 'desmond31'}
            vtdContra = expSpontaneous(iExp).vtdR;
            vtdIpsi = expSpontaneous(iExp).vtdL;
        otherwise
            error();
    end
    pos(iExp).contra = struct(t=vtdContra.Timestamp, x=normalize(vtdContra.handContra_X, 'zscore'), y=normalize(vtdContra.handContra_Y, 'zscore'));
    sel = vtdContra.handContra_Likelihood < likelihoodThreshold;
    pos(iExp).contra.x(sel) = NaN;
    pos(iExp).contra.y(sel) = NaN;
    pos(iExp).ipsi = struct(t=vtdIpsi.Timestamp, x=vtdIpsi.handIpsi_X, y=vtdIpsi.handIpsi_Y);
    sel = vtdIpsi.handIpsi_Likelihood < likelihoodThreshold;
    pos(iExp).ipsi.x(sel) = NaN;
    pos(iExp).ipsi.y(sel) = NaN;
end

% Break into bar-contact-aligned trials
trajWindow = [-2, 0];
trajResolution = 1/100;
t = trajWindow(1):trajResolution:trajWindow(2);
clear traj
traj(length(expSpontaneous)) = struct(t=[], T=[], contra=[], ipsi=[]);
for iExp = 1:length(expSpontaneous)
    trials = expSpontaneous(iExp).eu(1).Trials.Press;
    traj(iExp).t = t;
    traj(iExp).T = zeros(length(trials), length(t));
    traj(iExp).contra = struct(x=NaN(length(trials), length(t)), y=NaN(length(trials), length(t)));
    for iTrial = 1:length(trials)
        traj(iExp).T(iTrial, :) = t + trials(iTrial).Stop;
        traj(iExp).contra.x(iTrial, :) = interp1(pos(iExp).contra.t, pos(iExp).contra.x, traj(iExp).T(iTrial, :), 'linear');
        traj(iExp).contra.y(iTrial, :) = interp1(pos(iExp).contra.t, pos(iExp).contra.y, traj(iExp).T(iTrial, :), 'linear');
        traj(iExp).ipsi.x(iTrial, :) = interp1(pos(iExp).ipsi.t, pos(iExp).ipsi.x, traj(iExp).T(iTrial, :), 'linear');
        traj(iExp).ipsi.y(iTrial, :) = interp1(pos(iExp).ipsi.t, pos(iExp).ipsi.y, traj(iExp).T(iTrial, :), 'linear');
    end
end

clear pos iExp sel
%% Plot trajectories
% ax = axes(figure);
% hold(ax, 'on')
% for iExp = 1:length(expSpontaneous)
%     plot(ax, mean(traj(iExp).contra.x - traj(iExp).contra.x(:, end-90), 1, 'omitnan'), mean(traj(iExp).contra.y - traj(iExp).contra.y(:, end-90), 1, 'omitnan'))
% end
%% Detect onset times
onsetThreshold = 0.25;
onsetDir = 'reverse'; % reverse, forward
onsetPattern = [0 1 1];
onsetOffset = find(onsetPattern, 1, 'first');

clear onset
onset(length(expSpontaneous)) = struct(contra=[], ipsi=[], either=[], both=[]);

for iExp = 1:length(expSpontaneous)
    for pawName = ["contra", "ipsi"]
        nTrials = size(traj(iExp).T, 1);
        onset(iExp).(pawName) = NaN(nTrials, 1);
        t = traj(iExp).t;
        dist = zeros(nTrials, length(t));
        for iTrial = 1:nTrials
            pos = [ ...
                traj(iExp).(pawName).x(iTrial, :); ...
                traj(iExp).(pawName).y(iTrial, :); ...
                ];
            dist(iTrial, :) = sqrt(sum(pos.^2, 1));
        end
        distNorm = (dist - median(dist(:), 'omitnan'))./(mad(dist(:), 1)./0.6745); % Normalize using session median/mad
        for iTrial = 1:nTrials
            if any(isnan(distNorm(iTrial, :)))
                continue
            end
    
            isAbove = distNorm(iTrial, :) >= onsetThreshold;
            iLastAbove = find(isAbove, 1, 'last'); % We don't do abs since dist is already positive, although z-scoring will generate negatives, true movement should be positive.
            if isempty(iLastAbove)
                continue
            end
            iOnset = strfind(isAbove(1:iLastAbove), onsetPattern) + onsetOffset; % Woah, offset ~= ~onset!
            if isempty(iOnset)
                continue
            end
            switch onsetDir
                case 'reverse'
                    iOnset = iOnset(end);
                case 'forward'
                    iOnset = iOnset(1);
            end
            onset(iExp).(pawName)(iTrial) = t(iOnset) - t(end);
        end
    end
    onset(iExp).either = min(onset(iExp).contra, onset(iExp).ipsi);
    onset(iExp).both = max(onset(iExp).contra, onset(iExp).ipsi);
end
clear iExp pawName nTrials t dist distNorm iTrial pos isAbove iLastAbove iOnset

% Plot onset times
figure, histogram(cat(1, onset.contra), -2:0.05:0), title(sprintf('median=%gs', median(cat(1, onset.contra), 'omitnan')))
xlabel('contra paw onset time before bar contact (s)')


% Calculate ETA with movement onset time correction
clear etaByExp etaRawByExp
onsetThreshold = -Inf;
etaByExp(length(expSpontaneous)) = struct(X=[], t=[], N=[], D=[], stats=[]);
etaRawByExp(length(expSpontaneous)) = struct(X=[], t=[], N=[], D=[]);
for iExp = 1:length(expSpontaneous)
    trials = expSpontaneous(iExp).eu(1).getTrials('press');
    selTrials = onset(iExp).contra >= onsetThreshold;
    etaByExp(iExp) = expSpontaneous(iExp).eu.getETA('count', 'press', window=[-4, 2], ...
        alignTo='stop', includeInvalid=true, normalize=[-4, -2], MinTrialDuration=6, ...
        correction=onset(iExp).contra(selTrials), trials=trials(selTrials));
    etaRawByExp(iExp) = expSpontaneous(iExp).eu.getETA('count', 'press', window=[-4, 2], ...
        alignTo='stop', includeInvalid=true, normalize='none', MinTrialDuration=6, ...
        correction=onset(iExp).contra(selTrials), trials=trials(selTrials));
end


X = arrayfun(@(eta) eta.X, etaByExp, 'UniformOutput', false);
N = arrayfun(@(eta) eta.N, etaByExp, 'UniformOutput', false);
D = arrayfun(@(eta) eta.D, etaByExp, 'UniformOutput', false);
stats = arrayfun(@(eta) eta.stats, etaByExp, 'UniformOutput', false);
etaSpontaneous.X = cat(1, X{:});
etaSpontaneous.t = etaByExp(1).t;
etaSpontaneous.N = cat(1, N{:});
etaSpontaneous.D = cat(1, D{:});
etaSpontaneous.stats = cat(2, stats{:});

X = arrayfun(@(eta) eta.X, etaRawByExp, 'UniformOutput', false);
N = arrayfun(@(eta) eta.N, etaRawByExp, 'UniformOutput', false);
D = arrayfun(@(eta) eta.D, etaRawByExp, 'UniformOutput', false);
etaSpontaneousRaw.X = cat(1, X{:});
etaSpontaneousRaw.t = etaRawByExp(1).t;
etaSpontaneousRaw.N = cat(1, N{:});
etaSpontaneousRaw.D = cat(1, D{:});

% etaSpontaneous = euSpontaneous.getETA('count', 'press', window=[-4, 2], alignTo='stop', includeInvalid=true, ...
%         normalize=[-4, -2], MinTrialDuration=6);
% 
% etaSpontaneousRaw = euSpontaneous.getETA('count', 'press', window=[-4, 2], alignTo='stop', includeInvalid=true, ...
%         normalize='none', MinTrialDuration=6);

%%
cSpontaneous.hasPress = arrayfun(@(e) nnz(e.getTrials('press').duration() >= 6) >= p.minNumTrials, euSpontaneous)';
cSpontaneous.hasLick = arrayfun(@(e) nnz(e.getTrials('lick').duration() >= 6) >= p.minNumTrials, euSpontaneous)';


%% Bootstrap movement response
[~, b] = ismember({euSpontaneous.ExpName}, {expSpontaneous.name});
correction = {onset.contra};
correction = correction(b);

p.bootAlpha = 0.01;
bootSpontaneous.press = struct('h', NaN(length(euSpontaneous), 1), 'muDiffCI', NaN(length(euSpontaneous), 2), 'muDiffObs', NaN(length(euSpontaneous), 1));
[bootSpontaneous.press.h(cSpontaneous.hasPress), bootSpontaneous.press.muDiffCI(cSpontaneous.hasPress, :), bootSpontaneous.press.muDiffObs(cSpontaneous.hasPress)] = bootstrapMoveResponse( ...
    euSpontaneous(cSpontaneous.hasPress), 'press', alpha=p.bootAlpha, withReplacement=false, oneSided=false, ...
    responseWindow=[-0.3, 0], correction=correction);
fprintf(1, '\nAll done\n')

% Report bootstraped movement response direction
assert(nnz(isnan(bootSpontaneous.press.h(cSpontaneous.hasPress))) == 0)

figure, histogram(bootSpontaneous.press.h)
cSpontaneous.isPressUp = bootSpontaneous.press.h == 1 & cSpontaneous.hasPress;
cSpontaneous.isPressDown = bootSpontaneous.press.h == -1 & cSpontaneous.hasPress;
cSpontaneous.isPressResponsive = cSpontaneous.isPressUp | cSpontaneous.isPressDown;

fprintf(1, ['%g total SNr units:\n' ...
    '\t%g with 30+ press trials;\n'], ...
    length(euSpontaneous), nnz(cSpontaneous.hasPress))

fprintf(1, ['%g units with 30+ press trials (6s or longer):\n' ...
    '\t%g (%.0f%%) are excited (p<%g);\n' ...
    '\t%g (%.0f%%) are inhibited (p<%g).\n'], ...
    nnz(cSpontaneous.hasPress), ...
    nnz(cSpontaneous.isPressUp), 100*nnz(cSpontaneous.isPressUp)./nnz(cSpontaneous.isPressResponsive), p.bootAlpha, ...
    nnz(cSpontaneous.isPressDown), 100*nnz(cSpontaneous.isPressDown)./nnz(cSpontaneous.isPressResponsive), p.bootAlpha);

%% Save results (only once)
% euSpontaneous.save('C:\SERVER\Units\acute_spontaneous_reach\SNr_SingleUnit_NonDuplicate_NonDrift')
% if ~exist('C:\SERVER\Units\acute_spontaneous_reach\meta', 'dir')
%     mkdir('C:\SERVER\Units\acute_spontaneous_reach\meta')
% end
% save('C:\SERVER\Units\acute_spontaneous_reach\meta\SNr_SingleUnit_NonDuplicate_NonDrift.mat', ...
%     'bootSpontaneous', 'cSpontaneous', 'etaSpontaneous', 'etaSpontaneousRaw', 'onset', 'traj', 'correction', ...
%     '-v7.3')

%% Calculate bta's, split by press up or down
% 
% [btaSpontaneous.pressUpRaw.X, btaSpontaneous.pressUpRaw.T, btaSpontaneous.pressUpRaw.N, btaSpontaneous.pressUpRaw.S, btaSpontaneous.pressUpRaw.B] = euSpontaneous(cSpontaneous.isPressUp).getBinnedTrialAverage('rate', [8, 16, 32, 64], 'press', 'window', [-4, 0], 'normalize', false);
% [btaSpontaneous.pressDownRaw.X, btaSpontaneous.pressDownRaw.T, btaSpontaneous.pressDownRaw.N, btaSpontaneous.pressDownRaw.S, btaSpontaneous.pressDownRaw.B] = euSpontaneous(cSpontaneous.isPressDown).getBinnedTrialAverage('rate', [8, 16, 32, 64], 'press', 'window', [-4, 0], 'normalize', false);
% 
