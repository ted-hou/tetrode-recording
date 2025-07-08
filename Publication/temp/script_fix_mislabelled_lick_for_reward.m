eu = EphysUnit.load('C:\SERVER\Units\Lite_NonDuplicate_NonDrift', waveforms=false, spikecounts=false, spikerates=false);
euNames = lower(eu.getName());
files = dir('C:\SERVER\Units\NonLite_PressVsLick\*.mat');
sel = ismember(cellfun(@(n) lower(strrep(n, '.mat', '')), {files.name}, UniformOutput=false), euNames);
files = files(sel);
cd('C:\SERVER\Units\NonLite_PressVsLick\')
euComplete = EphysUnit.load({files.name}, waveforms=false, spikecounts=false, spikerates=false);
% euComplete.save('C:\SERVER\Units\NonLite_PressVsLick_NonDuplicate_NonDrift');
[lia, locb] = ismember(eu.getName(), euComplete.getName());
eu(lia) = euComplete(locb(lia));
clear euComplete
load('C:\SERVER\Units\meta_Lite_NonDuplicate_NonDrift.mat')
%% Align ArduinoConnection events to ephys time using a common event (try CueOn)
clear ac
[ac.ac, ac.euIndicesFirstInSession, ac.expIndices, ac.uniqueExpNames] = eu.loadArduinoConnection();
%% Find retract times
tEU = eu.alignTimestamps(["LEVER_RELEASED", "LEVER_RETRACT_START", "LEVER_RETRACTED", "LEVER_RETRACT_END", "TUBE_RETRACT_START", "TUBE_RETRACT_END"], ac={ac.ac, ac.euIndicesFirstInSession, ac.expIndices, ac.uniqueExpNames});
for iEu = 1:length(eu)
    try
        leverReleaseTimes = eu(iEu).EventTimes.LEVER_RELEASED;
        leverRetractTimes = eu(iEu).EventTimes.LEVER_RETRACT_START;
        if isempty(leverRetractTimes)
            leverRetractTimes = eu(iEu).EventTimes.LEVER_RETRACTED;
        end
        % eu(iEu).Trials.Press = Trial(eu(iEu).EventTimes.Cue, eu(iEu).EventTimes.Press, 'first');
        eu(iEu).Trials.PressIncorrect = eu(iEu).Trials.Press(eu(iEu).Trials.Press.duration() < 4 & eu(iEu).Trials.Press.duration() >= 2);
        eu(iEu).Trials.PressCorrect = eu(iEu).Trials.Press(eu(iEu).Trials.Press.duration() >= 4);
        [~, leverRetractTimesIncorrect] = eu(iEu).Trials.PressIncorrect.inTrial(leverRetractTimes, [0, 2], windowMode='stop');
        assert(nnz(leverRetractTimesIncorrect) > 0)
        [~, leverRetractTimesCorrect] = eu(iEu).Trials.PressCorrect.inTrial(leverRetractTimes, [-0.1, 8], windowMode='stop');
        eu(iEu).Trials.RetractReleaseIncorrect = Trial([leverRetractTimesIncorrect, Inf], leverReleaseTimes, 'first');
        eu(iEu).Trials.RetractReleaseCorrect = Trial([leverRetractTimesCorrect, Inf], leverReleaseTimes, 'first');
        eu(iEu).Trials.PressReleaseIncorrect = Trial([eu(iEu).Trials.PressIncorrect.Start, Inf], [eu(iEu).Trials.RetractReleaseIncorrect.Stop], 'first');
        eu(iEu).Trials.PressReleaseCorrect = Trial([eu(iEu).Trials.PressCorrect.Stop, Inf], [eu(iEu).Trials.RetractReleaseCorrect.Stop], 'first');
        eu(iEu).Trials.PressRetractIncorrect = Trial([eu(iEu).Trials.PressIncorrect.Stop], leverRetractTimesIncorrect, 'first');
        eu(iEu).Trials.PressRetractCorrect = Trial([eu(iEu).Trials.PressCorrect.Stop], leverRetractTimesCorrect, 'first');

    catch
        warning('Counld not process iEu = %i, "%s"', iEu, eu(iEu).getName())
    end
end
%% Find and correct bad lick labels
for iEu = 1:length(eu)
    eu(iEu).EventTimes.LICK = [];
    eu(iEu).EventTimes.LICK_OFF = [];
end
eu.alignTimestamps(["LICK", "LICK_OFF", "REWARD_ON"], ac={ac.ac, ac.euIndicesFirstInSession, ac.expIndices, ac.uniqueExpNames});
clear nLicks
nLicks.EUvAC = arrayfun(@(eu) [length(eu.EventTimes.Lick), length(eu.EventTimes.LICK)], eu, UniformOutput=false);
nLicks.EUvAC = cat(1, nLicks.EUvAC{:});
nLicks.EUvReward = arrayfun(@(eu) [length(eu.EventTimes.Lick), length(eu.EventTimes.RewardTimes)], eu, UniformOutput=false);
nLicks.EUvReward = cat(1, nLicks.EUvReward{:});
lickIsReward = diff(nLicks.EUvReward, 1, 2) == 0; lickIsReward = lickIsReward(:)';
for iEu = find(lickIsReward)
    eu(iEu).Trials.Lick = Trial(eu(iEu).EventTimes.Cue, eu(iEu).EventTimes.LICK, 'first', eu(iEu).EventTimes.Press);
end
%%
% 2.3.1  Basic summaries
% Baseline (median) spike rates
msr = arrayfun(@(stats) stats.median, [eu.SpikeRateStats]);
% Lick/Press responses
eta.press = eu.getETA('count', 'press', p.etaWindow, minTrialDuration=p.minTrialDuration, normalize=p.etaNorm);
eta.lick = eu.getETA('count', 'lick', p.etaWindow, minTrialDuration=p.minTrialDuration, normalize=p.etaNorm);
eta.pressIncorrect = eu.getETA('count', 'press', p.etaWindow, minTrialDuration=p.minTrialDuration, maxTrialDuration=4, normalize=p.etaNorm);
eta.pressRaw = eu.getETA('count', 'press', p.etaWindow, minTrialDuration=p.minTrialDuration, normalize='none');
eta.lickRaw = eu.getETA('count', 'lick', p.etaWindow, minTrialDuration=p.minTrialDuration, normalize='none');
eta.pressCue = eu.getETA('count', 'press', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize=eta.press.stats, includeInvalid=true);
eta.lickCue = eu.getETA('count', 'lick', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize=eta.lick.stats, includeInvalid=true);
eta.pressCueRaw = eu.getETA('count', 'press', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize='none', includeInvalid=true);
eta.lickCueRaw = eu.getETA('count', 'lick', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize='none', includeInvalid=true);
meta.press = transpose(mean(eta.press.X(:, eta.press.t >= p.metaWindowPress(1) & eta.press.t <= p.metaWindowPress(2)), 2, 'omitnan'));
meta.lick = transpose(mean(eta.lick.X(:, eta.lick.t >= p.metaWindowLick(1) & eta.lick.t <= p.metaWindowLick(2)), 2, 'omitnan'));
meta.pressRaw = transpose(mean(eta.pressRaw.X(:, eta.pressRaw.t >= p.metaWindowPress(1) & eta.pressRaw.t <= p.metaWindowPress(2)), 2, 'omitnan'));
meta.lickRaw = transpose(mean(eta.lickRaw.X(:, eta.lickRaw.t >= p.metaWindowLick(1) & eta.lickRaw.t <= p.metaWindowLick(2)), 2, 'omitnan'));
meta.pressRawBaseline = transpose(mean(eta.pressRaw.X(:, eta.pressRaw.t >= p.etaNorm(1) & eta.pressRaw.t <= p.etaNorm(2)), 2, 'omitnan'));
meta.lickRawBaseline = transpose(mean(eta.lickRaw.X(:, eta.lickRaw.t >= p.etaNorm(1) & eta.lickRaw.t <= p.etaNorm(2)), 2, 'omitnan'));
p.metaWindowCue = [-0.7, 0.2];
meta.pressCue = transpose(mean(eta.pressCue.X(:, eta.pressCue.t >= p.metaWindowCue(1) & eta.pressCue.t <= p.metaWindowCue(2)), 2, 'omitnan'));
meta.lickCue = transpose(mean(eta.lickCue.X(:, eta.lickCue.t >= p.metaWindowCue(1) & eta.lickCue.t <= p.metaWindowCue(2)), 2, 'omitnan'));
meta.pressCueRaw = transpose(mean(eta.pressCueRaw.X(:, eta.pressCueRaw.t >= p.metaWindowCue(1) & eta.pressCueRaw.t <= p.metaWindowCue(2)), 2, 'omitnan'));
meta.lickCueRaw = transpose(mean(eta.lickCueRaw.X(:, eta.lickCueRaw.t >= p.metaWindowCue(1) & eta.lickCueRaw.t <= p.metaWindowCue(2)), 2, 'omitnan'));
% etaSmooth.pressRaw = eu.getETA('rate', 'press', p.etaWindow, minTrialDuration=p.minTrialDuration, normalize='none');
% etaSmooth.lickRaw = eu.getETA('rate', 'lick', p.etaWindow, minTrialDuration=p.minTrialDuration, normalize='none');
% etaSmooth.pressCueRaw = eu.getETA('rate', 'press', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize='none', includeInvalid=true);
% etaSmooth.lickCueRaw = eu.getETA('rate', 'lick', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize='none', includeInvalid=true);
etaSmooth.pressCue = eu.getETA('rate', 'press', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize=etaSmooth.press.stats, includeInvalid=true);
etaSmooth.lickCue = eu.getETA('rate', 'lick', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize=etaSmooth.lick.stats, includeInvalid=true);
% 2.3.2 Basic summaries (fast)
% hasPress/hasLick
c.hasPress = arrayfun(@(e) nnz(e.getTrials('press').duration() >= p.minTrialDuration) >= p.minNumTrials, eu);
c.hasLick = arrayfun(@(e) nnz(e.getTrials('lick').duration() >= p.minTrialDuration) >= p.minNumTrials, eu);
% press/lick x Up/Down
% c.isPressUp =         c.hasPress & meta.press >= p.posRespThreshold;
% c.isPressDown =       c.hasPress & meta.press <= p.negRespThreshold;
% c.isPressResponsive = c.isPressUp | c.isPressDown;
% c.isLickUp =          c.hasLick & meta.lick >= p.posRespThreshold;
% c.isLickDown =        c.hasLick & meta.lick <= p.negRespThreshold;
% c.isLickResponsive =  c.isLickUp | c.isLickDown;
% animal info
c.isWT = arrayfun(@(eu) strcmpi(getAnimalInfo(eu, ai, 'strain'), 'wt'), eu);
c.isD1 = arrayfun(@(eu) strcmpi(getAnimalInfo(eu, ai, 'strain'), 'd1-cre'), eu);
c.isA2A = arrayfun(@(eu) strcmpi(getAnimalInfo(eu, ai, 'strain'), 'a2a-cre'), eu);
c.isAi80 = arrayfun(@(eu) strcmpi(getAnimalInfo(eu, ai, 'strain'), 'd1-cre;dlx-flp;ai80'), eu);
c.isDAT = arrayfun(@(eu) strcmpi(getAnimalInfo(eu, ai, 'strain'), 'dat-cre'), eu);
c.isAcute = ismember(eu.getAnimalName, {'daisy14', 'daisy15', 'daisy16', 'desmond23', 'desmond24', 'desmond25', 'desmond26', 'desmond27'});
%% Calculate neural onset times
% Use smooth ETA for better temporal resolution
etaSmooth.press = eu.getETA('rate', 'press', p.etaWindow, minTrialDuration=p.minTrialDuration, normalize=p.etaNorm);
etaSmooth.lick = eu.getETA('rate', 'lick', p.etaWindow, minTrialDuration=p.minTrialDuration, normalize=p.etaNorm);
%%
p.kernel = struct(type='gaussian', params= struct(sigma=0.05, window=[-0.05, 0.05], resolution=1e-3, width=0.1));
etaSmooth.press2 = eu.getETA('rate', 'press', p.etaWindow, minTrialDuration=p.minTrialDuration, normalize=p.etaNorm, kernel=p.kernel);
etaSmooth.lick2 = eu.getETA('rate', 'lick', p.etaWindow, minTrialDuration=p.minTrialDuration, normalize=p.etaNorm, kernel=p.kernel);
etaSmooth.kernel2 = p.kernel;
%%
etaFine.press = eu.getETA('count', 'press', p.etaWindow, minTrialDuration=p.minTrialDuration, normalize=p.etaNorm, resolution=0.025);
etaFine.lick = eu.getETA('count', 'lick', p.etaWindow, minTrialDuration=p.minTrialDuration, normalize=p.etaNorm, resolution=0.025);
etaFine.p.minTrialDuration = p.minTrialDuration;
etaFine.p.norm = p.etaNorm;
etaFine.p.resolution = 0.025;
%% Calculate onset
p.etaOnsetThreshold = 0.25;
p.etaSortWindow = [-3, 0];
p.etaSignWindow = [-0.3, 0];
p.etaOnsetPattern = [zeros(1, 25), ones(1, 50)];
fig = figure(Units='normalized', Position=[0.1, 0.1, 0.8, 0.8]);
ax(1) = subplot(1, 2, 1);
ax(2) = subplot(1, 2, 2);
n = length(eu);
onset = struct(press=NaN(n, 1), lick=NaN(n, 1), pressOrder=NaN(n, 1), lickOrder=NaN(n, 1));
[~, onset.pressOrder(c.hasPress), ~, onset.press(c.hasPress)] = EphysUnit.plotETA(ax(1), etaSmooth.press, c.hasPress, sortWindow=p.etaSortWindow, signWindow=p.etaSignWindow, ...
    sortThreshold=p.etaOnsetThreshold, negativeSortThreshold=p.etaOnsetThreshold, clim=[-2, 2], onsetPattern=p.etaOnsetPattern, xlim=[-4, 0.5], ...
    onsetDirection='reverse');
hold(ax(1), 'on')
x = onset.press(c.hasPress);
I = onset.pressOrder(c.hasPress);
plot(ax(1), x(I), 1:length(I))
[~, onset.lickOrder(c.hasLick), ~, onset.lick(c.hasLick)] = EphysUnit.plotETA(ax(2), etaSmooth.lick, c.hasLick, sortWindow=p.etaSortWindow, signWindow=p.etaSignWindow, ...
    sortThreshold=p.etaOnsetThreshold, negativeSortThreshold=p.etaOnsetThreshold, clim=[-2, 2], onsetPattern=p.etaOnsetPattern, xlim=[-4, 0.5], ...
    onsetDirection='reverse');
hold(ax(2), 'on')
x = onset.lick(c.hasLick);
I = onset.lickOrder(c.hasLick);
plot(ax(2), x(I), 1:length(I))
clim(ax(1), [-1.5, 1.5])
clim(ax(2), [-1.5, 1.5])
title(ax(1), sprintf('%i NaN', nnz(isnan(onset.press(c.isPressResponsive)))))
title(ax(2), sprintf('%i NaN (%i NaN)', nnz(isnan(onset.lick(c.isLickResponsive))), nnz(isnan(onset.lick(c.isLickResponsive & c.isPressResponsive)))))
clear n x I
p.bootAlpha = 0.01;
p.nboot = 100000;
p.responseWindowPress = [-0.3, 0];
p.responseWindowLick = [-0.3, 0];
assert(isequal(p.responseWindowPress, [-0.3, 0]))
assert(isequal(p.responseWindowLick, [-0.3, 0]))
boot.press = struct('h', NaN(length(eu), 1), 'muDiffCI', NaN(length(eu), 2), 'muDiffObs', NaN(length(eu), 1));
boot.lick = struct('h', NaN(length(eu), 1), 'muDiffCI', NaN(length(eu), 2), 'muDiffObs', NaN(length(eu), 1));
[boot.press.h(c.hasPress), boot.press.muDiffCI(c.hasPress, :), boot.press.muDiffObs(c.hasPress)] = bootstrapMoveResponse( ...
    eu(c.hasPress), 'press', nboot=p.nboot, alpha=p.bootAlpha, withReplacement=false, oneSided=false, ...
    responseWindow=p.responseWindowPress);
[boot.lick.h(c.hasLick), boot.lick.muDiffCI(c.hasLick, :), boot.lick.muDiffObs(c.hasLick)] = bootstrapMoveResponse( ...
    eu(c.hasLick), 'lick', nboot=p.nboot, alpha=p.bootAlpha, withReplacement=false, oneSided=false, ...
    responseWindow=p.responseWindowLick);
fprintf(1, '\nAll done\n')
%% Report bootstraped movement response direction
assert(nnz(isnan(boot.lick.h(c.hasLick))) == 0)
assert(nnz(isnan(boot.press.h(c.hasPress))) == 0)
figure, histogram(boot.press.h)
c.isPressUp = boot.press.h' == 1 & c.hasPress;
c.isPressDown = boot.press.h' == -1 & c.hasPress;
c.isPressResponsive = c.isPressUp | c.isPressDown;
figure, histogram(boot.lick.h)
c.isLickUp = boot.lick.h' == 1 & c.hasLick;
c.isLickDown = boot.lick.h' == -1 & c.hasLick;
c.isLickResponsive = c.isLickUp | c.isLickDown;
fprintf(1, ['%g total SNr units (baseline spike rate > %g):\n' ...
    '\t%g with %d+ press trials;\n' ...
    '\t%g with %d+ lick trials;\n' ...
    '\t%g with either (%g+ trials);\n' ...
    '\t%g with both (%g+ trials).\n'], ...
    length(eu), p.minSpikeRate, nnz(c.hasPress), p.minNumTrials, ...
    nnz(c.hasLick), p.minNumTrials, ...
    nnz(c.hasPress | c.hasLick), p.minNumTrials, ...
    nnz(c.hasPress & c.hasLick), p.minNumTrials)
fprintf(1, ['%g units with %g+ press trials (%gs or longer):\n' ...
    '\t%g (%.0f%%) are excited (p<%g);\n' ...
    '\t%g (%.0f%%) are inhibited (p<%g).\n'], ...
    nnz(c.hasPress), p.minNumTrials, p.minTrialDuration, ...
    nnz(c.isPressUp), 100*nnz(c.isPressUp)/nnz(c.isPressResponsive), p.bootAlpha, ...
    nnz(c.isPressDown), 100*nnz(c.isPressDown)/nnz(c.isPressResponsive), p.bootAlpha);
fprintf(1, ['%g units with %g+ lick trials (%gs or longer):\n' ...
    '\t%g (%.0f%%) are excited (p<%g);\n' ...
    '\t%g (%.0f%%) are inhibited (p<%g).\n'], ...
    nnz(c.hasLick), p.minNumTrials, p.minTrialDuration, ...
    nnz(c.isLickUp), 100*nnz(c.isLickUp)/nnz(c.isLickResponsive), p.bootAlpha, ...
    nnz(c.isLickDown), 100*nnz(c.isLickDown)/nnz(c.isLickResponsive), p.bootAlpha);
nTotal = nnz(c.isPressResponsive & c.isLickResponsive);
fprintf(1, ['%g units with %d+ press AND lick trials (%gs or longer):\n' ...
    '\t%g (%.0f%%) are press-excited AND lick-excited;\n' ...
    '\t%g (%.0f%%) are press-inhibited AND lick-inhibited;\n' ...
    '\t%g (%.0f%%) are press-excited AND lick-inhibited;\n' ...
    '\t%g (%.0f%%) are press-inhibited AND lick-excited;\n'], ...
    nnz(c.hasPress & c.hasLick), p.minNumTrials, p.minTrialDuration, ...
    nnz(c.isPressUp & c.isLickUp), 100*nnz(c.isPressUp & c.isLickUp)/nTotal, ...
    nnz(c.isPressDown & c.isLickDown), 100*nnz(c.isPressDown & c.isLickDown)/nTotal, ...
    nnz(c.isPressUp & c.isLickDown), 100*nnz(c.isPressUp & c.isLickDown)/nTotal, ...
    nnz(c.isPressDown & c.isLickUp), 100*nnz(c.isPressDown & c.isLickUp)/nTotal)
sel = c.hasLick & c.hasPress;
fprintf('05 Calculate: Of %i: %i (%i%%) showed modulation for BOTH, %i (%i%%) showed modulation for lick only, %i (%i%%) showed modulation for reach only, %i (%i%%) for neither.', ...
    nnz(sel), ...
    nnz(sel & c.isPressResponsive & c.isLickResponsive), round(nnz(sel & c.isPressResponsive & c.isLickResponsive)/nnz(sel)*100), ...
    nnz(sel & c.isLickResponsive & ~c.isPressResponsive), round(nnz(sel & c.isLickResponsive & ~c.isPressResponsive)/nnz(sel)*100), ...
    nnz(sel & c.isPressResponsive & ~c.isLickResponsive), round(nnz(sel & c.isPressResponsive & ~c.isLickResponsive)/nnz(sel)*100), ...
    nnz(sel & ~c.isPressResponsive & ~c.isLickResponsive), round(nnz(sel & ~c.isPressResponsive & ~c.isLickResponsive)/nnz(sel)*100) ...
    )

%%
save('C:\SERVER\Units\meta_Lite_NonDuplicate_NonDrift_20250705.mat', 'p', 'c', 'eta', 'etaSmooth', 'euPos', 'meta', 'msr', 'onset', 'boot', 'ai')