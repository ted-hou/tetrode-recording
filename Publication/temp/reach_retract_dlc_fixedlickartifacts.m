% eu = EphysUnit.load('C:\SERVER\Units\PressVsLick_ArtifactsRemoved');
% load('C:\SERVER\Units\meta_PressVsLick_ArtifactsRemoved.mat');
% 
% %% Read Arduino events (servo/pressOff)
% clear ac
% [ac.ac, ac.euIndicesFirstInSession, ac.expIndices, ac.uniqueExpNames] = eu.loadArduinoConnection();
% 
% eu.alignTimestamps(["LEVER_RELEASED", "LEVER_PRESSED", "LEVER_RETRACT_START", "LEVER_RETRACTED", "LEVER_RETRACT_END", "TUBE_RETRACT_START", "TUBE_RETRACT_END"], ac={ac.ac, ac.euIndicesFirstInSession, ac.expIndices, ac.uniqueExpNames});
% 
% eu.loadDigitalEventsFromTetrodeRecording(["PressOn", "PressOff", "LickOn", "LickOff", "RewardOn", "RewardOff"]);
% 
% %% Find and correct bad lick labels
% clear nLicks
% nLicks.EUvReward = arrayfun(@(eu) [length(eu.EventTimes.LickOn), length(eu.EventTimes.RewardOn)], eu, UniformOutput=false);
% nLicks.EUvReward = cat(1, nLicks.EUvReward{:});
% 
% 
% lickIsReward = diff(nLicks.EUvReward, 1, 2) == 0; lickIsReward = lickIsReward(:)';
% 
% for iEu = 1:length(eu)
%     eu(iEu).EventTimes.LICK = [];
%     eu(iEu).EventTimes.LICK_OFF = [];
% end
% eu.alignTimestamps(["LICK", "LICK_OFF"], ac={ac.ac, ac.euIndicesFirstInSession, ac.expIndices, ac.uniqueExpNames});
% for iEu = find(lickIsReward)
%     eu(iEu).EventTimes.Bad_IsJustReward_Lick = eu(iEu).EventTimes.Lick;
%     eu(iEu).EventTimes.Bad_IsJustReward_LickOn = eu(iEu).EventTimes.LickOn;
%     eu(iEu).EventTimes.Bad_IsJustReward_LickOff = eu(iEu).EventTimes.LickOff;
%     eu(iEu).EventTimes.Lick = eu(iEu).EventTimes.LICK;
%     eu(iEu).EventTimes.LickOn = eu(iEu).EventTimes.LICK;
%     eu(iEu).EventTimes.LickOff = eu(iEu).EventTimes.LICK_OFF;
%     eu(iEu).Trials.Lick = Trial(eu(iEu).EventTimes.Cue, eu(iEu).EventTimes.LICK, 'first', eu(iEu).EventTimes.Press);
% end
% clear nLicks iEu
% 
% %%
% for iEu = 1:length(eu)
%     % try
%         % Discard false-lever-releases (minDuration=0.1)
%         releaseToPress = Trial([eu(iEu).EventTimes.PressOff, Inf], eu(iEu).EventTimes.PressOn, 'first');
%         releaseToPress = releaseToPress(releaseToPress.duration() > 5);
%         leverReleaseTimes = [releaseToPress.Start];
%         clear releaseToPress
% 
%         leverRetractTimes = eu(iEu).EventTimes.LEVER_RETRACT_START;
%         if isempty(leverRetractTimes)
%             leverRetractTimes = eu(iEu).EventTimes.LEVER_RETRACTED;
%             assert(~isempty(leverRetractTimes))
%         end
%         % Press trials
%         eu(iEu).Trials.Press = Trial(eu(iEu).EventTimes.Cue, eu(iEu).EventTimes.Press, 'first');
%         eu(iEu).Trials.PressValid = eu(iEu).Trials.Press(eu(iEu).Trials.Press.duration() >= 2);
%         eu(iEu).Trials.PressIncorrect = eu(iEu).Trials.Press(eu(iEu).Trials.Press.duration() < 4 & eu(iEu).Trials.Press.duration() >= 2);
%         eu(iEu).Trials.PressCorrect = eu(iEu).Trials.Press(eu(iEu).Trials.Press.duration() >= 4);
% 
%         % cue -> retract
%         eu(iEu).Trials.CueToLeverRetract = Trial([eu(iEu).Trials.PressValid.Start, Inf], leverRetractTimes, 'first');
%         eu(iEu).Trials.CueToLeverRetractIncorrect = Trial([eu(iEu).Trials.PressIncorrect.Start, Inf], leverRetractTimes, 'first');
%         eu(iEu).Trials.CueToLeverRetractCorrect = Trial([eu(iEu).Trials.PressCorrect.Start, Inf], leverRetractTimes, 'first');
%         % These assertions failed for 2 sessions ["daisy8_20210625", "desmond23_20220504"] b/c there were fewer retractions than cue (difference was only 1 or 12). WHY? 
%         % assert(length(eu(iEu).Trials.CueToLeverRetract) == length(eu(iEu).Trials.PressValid), "a: %i != %i", length(eu(iEu).Trials.CueToLeverRetract), length(eu(iEu).Trials.PressValid))
%         % assert(length(eu(iEu).Trials.CueToLeverRetractIncorrect) == length(eu(iEu).Trials.PressIncorrect), "b: %i != %i", length(eu(iEu).Trials.CueToLeverRetractIncorrect), length(eu(iEu).Trials.PressIncorrect))
%         % assert(length(eu(iEu).Trials.CueToLeverRetractCorrect) == length(eu(iEu).Trials.PressCorrect), "c: %i != %i", length(eu(iEu).Trials.CueToLeverRetractCorrect), length(eu(iEu).Trials.PressCorrect))
% 
%         % cue -> release
%         eu(iEu).Trials.CueToLeverRelease = Trial([eu(iEu).Trials.PressValid.Start, Inf], leverReleaseTimes, 'first');
%         eu(iEu).Trials.CueToLeverReleaseIncorrect = Trial([eu(iEu).Trials.PressIncorrect.Start, Inf], leverReleaseTimes, 'first'); % This should be almost identical to PressIncorrect/CueToLeverRetractIncorrect
%         eu(iEu).Trials.CueToLeverReleaseCorrect = Trial([eu(iEu).Trials.PressCorrect.Start, Inf], leverReleaseTimes, 'first');
% 
%         % correct trials: release before retract vs. release after retract
%         % cue = [eu(iEu).Trials.CueToLeverRetractCorrect.Start];
%         % retract = [eu(iEu).Trials.CueToLeverRetractCorrect.Stop];
%         % release = [eu(iEu).Trials.CueToLeverReleaseCorrect.Stop];
%         [isReleaseBeforeRetract, releaseBeforeRetract] = eu(iEu).Trials.CueToLeverRetractCorrect.inTrial(leverReleaseTimes);
%         releaseAfterRetract = leverReleaseTimes(~isReleaseBeforeRetract);
% 
%         eu(iEu).Trials.CueToLeverReleaseBeforeRetractCorrect = Trial([eu(iEu).Trials.CueToLeverRetractCorrect.Start, Inf], releaseBeforeRetract, 'first');
%         eu(iEu).Trials.RetractToLeverReleaseAfterRetractCorrect = Trial([eu(iEu).Trials.CueToLeverRetractCorrect.Stop, Inf], releaseAfterRetract, 'first', eu(iEu).EventTimes.PressOn);
%         clear isReleaseBeforeRetract releaseBeforeRetract releaseAfterRetract
% 
%         % [~, leverRetractTimesIncorrect] = eu(iEu).Trials.PressIncorrect.inTrial(leverRetractTimes, [0, 2], windowMode='stop');
%         % assert(nnz(leverRetractTimesIncorrect) > 0)
%         % [~, leverRetractTimesCorrect] = eu(iEu).Trials.PressCorrect.inTrial(leverRetractTimes, [-0.1, 8], windowMode='stop');
%         % eu(iEu).Trials.RetractReleaseIncorrect = Trial([leverRetractTimesIncorrect, Inf], leverReleaseTimes, 'first');
%         % eu(iEu).Trials.RetractReleaseCorrect = Trial([leverRetractTimesCorrect, Inf], leverReleaseTimes, 'first');
%         % eu(iEu).Trials.PressReleaseIncorrect = Trial([eu(iEu).Trials.PressIncorrect.Start, Inf], [eu(iEu).Trials.RetractReleaseIncorrect.Stop], 'first');
%         % eu(iEu).Trials.PressReleaseCorrect = Trial([eu(iEu).Trials.PressCorrect.Stop, Inf], [eu(iEu).Trials.RetractReleaseCorrect.Stop], 'first');
%         % eu(iEu).Trials.PressRetractIncorrect = Trial([eu(iEu).Trials.PressIncorrect.Stop], leverRetractTimesIncorrect, 'first');
%         % eu(iEu).Trials.PressRetractCorrect = Trial([eu(iEu).Trials.PressCorrect.Stop], leverRetractTimesCorrect, 'first');
%     % catch ME
%     %     warning('Counld not process iEu = %i, "%s"', iEu, eu(iEu).getName())
%     %     warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
%     % end
% end
% 
% nReleaseTrials = arrayfun(@(eu) struct(release=length(eu.Trials.CueToLeverRelease), releaseCorrect=length(eu.Trials.CueToLeverReleaseCorrect), ...
%     releaseBefore=length(eu.Trials.CueToLeverReleaseBeforeRetractCorrect), releaseAfter=length(eu.Trials.RetractToLeverReleaseAfterRetractCorrect)), eu);
% 
% clear ac iEu leverReleaseTimes leverRetractTimes
% 
% %% Save results
% eu.save('C:\SERVER\Units\PressVsLick_ArtifactsRemoved\FixedEventsAndTrials');

%% Load results
eu = EphysUnit.load('C:\SERVER\Units\PressVsLick_ArtifactsRemoved\FixedEventsAndTrials');
load('C:\SERVER\Units\meta_PressVsLick_ArtifactsRemoved.mat');

%% TODO: Load deeplabcut (Derin did all the older sessions)


%% Calculate ETA for correct reach vs. incorrect reach; correct vs. incorrect retract; correct vs. incorrect release; correct vs. incorrect lick
clear eta

clear artifactParams;
artifactParams(1) = struct(event='LickOn', length=10, lengthUnit='ms', direction='both');
artifactParams(2) = struct(event='LickOff', length=10, lengthUnit='ms', direction='both');
artifactParams(3) = struct(event='PressOn', length=10, lengthUnit='ms', direction='both');
artifactParams(4) = struct(event='PressOff', length=10, lengthUnit='ms', direction='both');

eta.artifactParams = artifactParams;
eta.baselineWindow = struct(press=[-4, -2], lick=[-4, -2], release=[-2, 0]);


eta.correctPress = eu.getETA('count', 'PressCorrect', [-4, 4], alignTo='stop', ...
    normalize='none', resolution=0.1, artifacts=artifactParams);
eta.incorrectPress = eu.getETA('count', 'PressIncorrect', [-4, 4], alignTo='stop', ...
    normalize='none', resolution=0.1, artifacts=artifactParams);

eta.correctRetract = eu.getETA('count', 'CueToLeverRetractCorrect', [-4, 4], alignTo='stop', ...
    normalize='none', resolution=0.1, artifacts=artifactParams);
eta.incorrectRetract = eu.getETA('count', 'CueToLeverRetractIncorrect', [-4, 4], alignTo='stop', ...
    normalize='none', resolution=0.1, artifacts=artifactParams);

eta.correctRelease = eu.getETA('count', 'CueToLeverReleaseCorrect', [-4, 4], alignTo='stop', ...
    normalize='none', resolution=0.1, artifacts=artifactParams);
eta.incorrectRelease = eu.getETA('count', 'CueToLeverReleaseIncorrect', [-4, 4], alignTo='stop', ...
    normalize='none', resolution=0.1, artifacts=artifactParams);

eta.correctReleaseBeforeRetract = eu.getETA('count', 'CueToLeverReleaseBeforeRetractCorrect', [-4, 4], alignTo='stop', ...
    normalize='none', resolution=0.1, artifacts=artifactParams);
eta.correctReleaseAfterRetract = eu.getETA('count', 'RetractToLeverReleaseAfterRetractCorrect', [-4, 4], alignTo='stop', ...
    normalize='none', resolution=0.1, artifacts=artifactParams);

eta.correctLick = eu.getETA('count', 'lick', [-4, 4], minTrialDuration=4, ...
    normalize='none', resolution=0.1, artifacts=artifactParams);
eta.incorrectLick = eu.getETA('count', 'lick', [-4, 4], minTrialDuration=2, maxTrialDuration=4, ...
    normalize='none', resolution=0.1, artifacts=artifactParams);

% Convert to sp/s
for field = ["correctPress", "incorrectPress", "correctLick", "incorrectLick", "correctRelease", "incorrectRelease", "correctReleaseBeforeRetract", "correctReleaseAfterRetract", "correctRetract", "incorrectRetract"]
    eta.(field).X = eta.(field).X ./ 0.1;
end

clear field

eta.pressNorm = eu.getETA('count', 'press', [-4, 4], alignTo='stop', minTrialDuration=2, normalize=eta.baselineWindow.press, resolution=0.1, artifacts=artifactParams);
eta.lickNorm = eu.getETA('count', 'lick', [-4, 4], alignTo='stop', minTrialDuration=2, normalize=eta.baselineWindow.lick, resolution=0.1, artifacts=artifactParams);
eta.lickNormToPress = eu.getETA('count', 'lick', [-4, 4], alignTo='stop', minTrialDuration=2, normalize=eta.pressNorm.stats, resolution=0.1, artifacts=artifactParams);
eta.correctReleaseNorm = eu.getETA('count', 'CueToLeverReleaseCorrect', [-4, 4], alignTo='stop', minTrialDuration=4, normalize=eta.baselineWindow.release, resolution=0.1, artifacts=artifactParams);
eta.correctReleaseNormToPress = eu.getETA('count', 'CueToLeverReleaseCorrect', [-4, 4], alignTo='stop', minTrialDuration=4, normalize=eta.pressNorm.stats, resolution=0.1, artifacts=artifactParams);

eta.correctPressNorm = eu.getETA('count', 'press', [-4, 4], alignTo='stop', minTrialDuration=4, normalize=eta.pressNorm.stats, resolution=0.1, artifacts=artifactParams);
eta.correctLickNormToPress = eu.getETA('count', 'lick', [-4, 4], alignTo='stop', minTrialDuration=4, normalize=eta.pressNorm.stats, resolution=0.1, artifacts=artifactParams);
eta.correctReleaseBeforeRetractNormToPress = eu.getETA('count', 'CueToLeverReleaseBeforeRetractCorrect', [-4, 4], alignTo='stop', minTrialDuration=4, normalize=eta.pressNorm.stats, resolution=0.1, artifacts=artifactParams);
eta.correctReleaseAfterRetractNormToPress = eu.getETA('count', 'RetractToLeverReleaseAfterRetractCorrect', [-4, 4], alignTo='stop', minTrialDuration=4, normalize=eta.pressNorm.stats, resolution=0.1, artifacts=artifactParams);

eta.incorrectPressNormToPress = eu.getETA('count', 'press', [-4, 4], alignTo='stop', minTrialDuration=2, maxTrialDuration=4, normalize=eta.pressNorm.stats, resolution=0.1, artifacts=artifactParams);
eta.incorrectLickNormToPress = eu.getETA('count', 'lick', [-4, 4], alignTo='stop', minTrialDuration=2, maxTrialDuration=4, normalize=eta.pressNorm.stats, resolution=0.1, artifacts=artifactParams);
%%
etaWithArti.pressNorm = eu.getETA('count', 'press', [-4, 4], alignTo='stop', minTrialDuration=2, normalize=eta.baselineWindow.press, resolution=0.1);
etaWithArti.lickNorm = eu.getETA('count', 'lick', [-4, 4], alignTo='stop', minTrialDuration=2, normalize=eta.baselineWindow.lick, resolution=0.1);
etaWithArti.lickNormToPress = eu.getETA('count', 'lick', [-4, 4], alignTo='stop', minTrialDuration=2, normalize=eta.pressNorm.stats, resolution=0.1);
etaWithArti.correctReleaseNorm = eu.getETA('count', 'CueToLeverReleaseCorrect', [-4, 4], alignTo='stop', minTrialDuration=4, normalize=eta.baselineWindow.release, resolution=0.1);
etaWithArti.correctReleaseNormToPress = eu.getETA('count', 'CueToLeverReleaseCorrect', [-4, 4], alignTo='stop', minTrialDuration=4, normalize=eta.pressNorm.stats, resolution=0.1);
etaWithArti.correctPressNorm = eu.getETA('count', 'press', [-4, 4], alignTo='stop', minTrialDuration=4, normalize=eta.pressNorm.stats, resolution=0.1);
etaWithArti.correctLickNormToPress = eu.getETA('count', 'lick', [-4, 4], alignTo='stop', minTrialDuration=4, normalize=eta.pressNorm.stats, resolution=0.1);
etaWithArti.correctReleaseBeforeRetractNormToPress = eu.getETA('count', 'CueToLeverReleaseBeforeRetractCorrect', [-4, 4], alignTo='stop', minTrialDuration=4, normalize=eta.pressNorm.stats, resolution=0.1);
etaWithArti.correctReleaseAfterRetractNormToPress = eu.getETA('count', 'RetractToLeverReleaseAfterRetractCorrect', [-4, 4], alignTo='stop', minTrialDuration=4, normalize=eta.pressNorm.stats, resolution=0.1);
etaWithArti.incorrectPressNormToPress = eu.getETA('count', 'press', [-4, 4], alignTo='stop', minTrialDuration=2, maxTrialDuration=4, normalize=eta.pressNorm.stats, resolution=0.1);
etaWithArti.incorrectLickNormToPress = eu.getETA('count', 'lick', [-4, 4], alignTo='stop', minTrialDuration=2, maxTrialDuration=4, normalize=eta.pressNorm.stats, resolution=0.1);

%% Calculate META
clear meta
meta.window.press = [-0.3, 0];
meta.window.lick = [-0.3, 0];
meta.window.correctRelease = [0, 0.3];

meta.pressNorm = mean(eta.pressNorm.X(:, eta.pressNorm.t > meta.window.press(1) & eta.pressNorm.t < meta.window.press(2)), 2, 'omitnan');
meta.lickNorm = mean(eta.lickNorm.X(:, eta.lickNorm.t > meta.window.lick(1) & eta.lickNorm.t < meta.window.lick(2)), 2, 'omitnan');
meta.lickNormToPress = mean(eta.lickNormToPress.X(:, eta.lickNormToPress.t > meta.window.lick(1) & eta.lickNormToPress.t < meta.window.lick(2)), 2, 'omitnan');
meta.correctReleaseNorm = mean(eta.correctReleaseNorm.X(:, eta.correctReleaseNorm.t > meta.window.correctRelease(1) & eta.correctReleaseNorm.t < meta.window.correctRelease(2)), 2, 'omitnan');
meta.correctReleaseNormToPress = mean(eta.correctReleaseNormToPress.X(:, eta.correctReleaseNormToPress.t > meta.window.correctRelease(1) & eta.correctReleaseNormToPress.t < meta.window.correctRelease(2)), 2, 'omitnan');

%% Bootstrap response directions
clear boot
boot.params = struct(...
    alpha = 0.01, ...
    nboot = 100000, ...
    responseWindow = meta.window, ...
    baselineWindow = struct(press=[-4, -2], lick='press', release='press'), ...
    artifacts = artifactParams ...
    );

[boot.press.h, boot.press.ci, boot.press.obs] = bootstrapMoveResponse(eu, ...
    'press', nboot=boot.params.nboot, ...
    baselineWindow=boot.params.baselineWindow.press, ...
    responseWindow=boot.params.responseWindow.press, ...
    alignTo='stop', allowedTrialDuration=[2, Inf], alpha=boot.params.alpha, ...
    artifacts=boot.params.artifacts);  

[boot.lick.h, boot.lick.ci, boot.lick.obs] = bootstrapMoveResponse(eu, ...
    struct(baseline='press', response='lick'), nboot=boot.params.nboot, ...
    baselineWindow=boot.params.baselineWindow.press, ...
    responseWindow=boot.params.responseWindow.lick, ...
    alignTo=struct(baseline='stop', response='stop'), ...
    allowedTrialDuration=struct(baseline=[2, Inf], response=[2, Inf]), ...
    alpha=boot.params.alpha, artifacts=boot.params.artifacts);  

[boot.correctRelease.h, boot.correctRelease.ci, boot.correctRelease.obs] = bootstrapMoveResponse(eu, ...
    struct(baseline='press', response='CueToLeverReleaseCorrect'), nboot=boot.params.nboot, ...
    baselineWindow=boot.params.baselineWindow.press, ...
    responseWindow=boot.params.responseWindow.correctRelease, ...
    alignTo=struct(baseline='stop', response='stop'), ...
    allowedTrialDuration=struct(baseline=[2, Inf], response=[4, Inf]), ...
    alpha=boot.params.alpha, artifacts=boot.params.artifacts);  

%% Categorize based on bootstrap results
cc.isPressUp = boot.press.h > 0 & ~isnan(boot.press.h);
cc.isPressDown = boot.press.h < 0 & ~isnan(boot.press.h);
cc.isLickUp = boot.lick.h > 0 & ~isnan(boot.lick.h);
cc.isLickDown = boot.lick.h < 0 & ~isnan(boot.lick.h);
cc.isCorrectReleaseUp = boot.correctRelease.h > 0 & ~isnan(boot.correctRelease.h);
cc.isCorrectReleaseDown = boot.correctRelease.h < 0 & ~isnan(boot.correctRelease.h);

cc.isPressResponsive = cc.isPressUp | cc.isPressDown;
cc.isLickResponsive = cc.isLickUp | cc.isLickDown;
cc.isCorrectReleaseResponsive = cc.isCorrectReleaseUp | cc.isCorrectReleaseDown;

fprintf('\n%i total units, %i (%.0f%%) reach modulated, %i (%.0f%%) lick modulated, %i (%.0f%%) correct-release modulated:\n', ...
    length(eu), ...
    nnz(cc.isPressResponsive), nnz(cc.isPressResponsive)/length(eu)*100, ...
    nnz(cc.isLickResponsive), nnz(cc.isLickResponsive)/length(eu)*100, ...
    nnz(cc.isCorrectReleaseResponsive), nnz(cc.isCorrectReleaseResponsive)/length(eu)*100)
fprintf('\t%i/%i units reach inc (%.0f%%);\n', nnz(cc.isPressUp), nnz(cc.isPressResponsive), nnz(cc.isPressUp)./nnz(cc.isPressResponsive)*100);
fprintf('\t%i/%i units reach dec (%.0f%%);\n', nnz(cc.isPressDown), nnz(cc.isPressResponsive), nnz(cc.isPressDown)./nnz(cc.isPressResponsive)*100);
fprintf('\t%i/%i units lick inc (%.0f%%);\n', nnz(cc.isLickUp), nnz(cc.isLickResponsive), nnz(cc.isLickUp)./nnz(cc.isLickResponsive)*100);
fprintf('\t%i/%i units lick dec (%.0f%%);\n', nnz(cc.isLickDown), nnz(cc.isLickResponsive), nnz(cc.isLickDown)./nnz(cc.isLickResponsive)*100);
fprintf('\t%i/%i units release inc (%.0f%%);\n', nnz(cc.isCorrectReleaseUp), nnz(cc.isCorrectReleaseResponsive), nnz(cc.isCorrectReleaseUp)./nnz(cc.isCorrectReleaseResponsive)*100);
fprintf('\t%i/%i units release dec (%.0f%%).\n', nnz(cc.isCorrectReleaseDown), nnz(cc.isCorrectReleaseResponsive), nnz(cc.isCorrectReleaseDown)./nnz(cc.isCorrectReleaseResponsive)*100);

sel = cc.isPressResponsive & cc.isLickResponsive & cc.isCorrectReleaseResponsive;
fprintf('\n%i/%i (%.0f%%) units modulated for reach, lick, and correct-release, of which:\n', nnz(sel), length(eu), nnz(sel)/length(eu)*100)
fprintf('\t%i/%i (%.0f%%) decrease for all three;\n', nnz(sel & cc.isPressDown & cc.isLickDown & cc.isCorrectReleaseDown), nnz(sel), nnz(sel & cc.isPressDown & cc.isLickDown & cc.isCorrectReleaseDown)./nnz(sel)*100)
fprintf('\t%i/%i (%.0f%%) increase for all three;\n\n', nnz(sel & cc.isPressUp & cc.isLickUp & cc.isCorrectReleaseUp), nnz(sel), nnz(sel & cc.isPressUp & cc.isLickUp & cc.isCorrectReleaseUp)./nnz(sel)*100)
fprintf('\t%i/%i (%.0f%%) decrease for lick and release;\n', nnz(sel & cc.isPressUp & cc.isLickDown & cc.isCorrectReleaseDown), nnz(sel), nnz(sel & cc.isPressUp & cc.isLickDown & cc.isCorrectReleaseDown)./nnz(sel)*100)
fprintf('\t%i/%i (%.0f%%) decrease for reach and release;\n', nnz(sel & cc.isPressDown & cc.isLickUp & cc.isCorrectReleaseDown), nnz(sel), nnz(sel & cc.isPressDown & cc.isLickUp & cc.isCorrectReleaseDown)./nnz(sel)*100)
fprintf('\t%i/%i (%.0f%%) decrease for reach and lick;\n\n', nnz(sel & cc.isPressDown & cc.isLickDown & cc.isCorrectReleaseUp), nnz(sel), nnz(sel & cc.isPressDown & cc.isLickDown & cc.isCorrectReleaseUp)./nnz(sel)*100)
fprintf('\t%i/%i (%.0f%%) decrease for just reach;\n', nnz(sel & cc.isPressDown & cc.isLickUp & cc.isCorrectReleaseUp), nnz(sel), nnz(sel & cc.isPressDown & cc.isLickUp & cc.isCorrectReleaseUp)./nnz(sel)*100)
fprintf('\t%i/%i (%.0f%%) decrease for just lick;\n', nnz(sel & cc.isPressUp & cc.isLickDown & cc.isCorrectReleaseUp), nnz(sel), nnz(sel & cc.isPressUp & cc.isLickDown & cc.isCorrectReleaseUp)./nnz(sel)*100)
fprintf('\t%i/%i (%.0f%%) decrease for just release;\n', nnz(sel & cc.isPressUp & cc.isLickUp & cc.isCorrectReleaseDown), nnz(sel), nnz(sel & cc.isPressUp & cc.isLickUp & cc.isCorrectReleaseDown)./nnz(sel)*100)

%% Save intermediate results because ETA with artifact blanking takes forever
% eu.save('C:\SERVER\Units\LickVsReach_FixedArtifacts')
save('C:\SERVER\Units\meta_LickVsReach_FixedArtifacts.mat', 'artifactParams', 'bootCLick', 'cc', 'circlick', 'eta', 'etaWithArti', 'oldSpikeTimes', 'meta', 'boot')

%% Load units, metadata and bootstrapping results.
% These units are good. They have lick vs reach trials, they have (Intan)
% lick artifacts filtered out and spikes redetected, units with drifting
% spike waveforms have been removed, Blackrock units were deemed good to
% keep without re-spike-sorting. I've fi    xed the missing digital events.
% I've blanked out [-10, 10]ms peri-lickOn/Off, peri-reachOn/Off for ETA
% calculation. Get these units with their fancy schmancy metadata now!
% eu = EphysUnit.load('C:\SERVER\Units\LickVsReach_FixedArtifacts');
% load('C:\SERVER\Units\meta_LickVsReach_FixedArtifacts.mat');
eu = EphysUnit.load('E:\DATA\Units\LickVsReach_FixedArtifacts');
load('E:\DATA\Units\meta_LickVsReach_FixedArtifacts.mat');

%% Scatter META vs META
% close all
% metaNames = ["pressNorm", "lickNorm", "lickNormToPress", "correctReleaseNorm", "correctReleaseNormToPress"];
% dispNames = {["press", "(norm to self)"], ["lick", "(norm to self)"], ["lick", "(norm to press)"], ["release", "(norm to self)"], ["release", "(norm to press)"]};
metaNames = ["pressNorm", "lickNormToPress", "correctReleaseNormToPress"];
bootResp = {cc.isPressResponsive, cc.isLickResponsive, cc.isCorrectReleaseResponsive};
dispNames = { ...
    ["Reach", sprintf("[%s]s", string(meta.window.press).join(', '))], ...
    ["Lick", sprintf("[%s]s", string(meta.window.lick).join(', '))], ...
    ["Release", sprintf("[%s]s", string(meta.window.correctRelease).join(', '))], ...
    };
n = length(metaNames);

clear tl
marginRatio = 2;
tl.parent = tiledlayout(figure(Units='inches', InnerPosition=[1 1 6 6]), marginRatio*n + 1, marginRatio*n + 1, TileSpacing='tight', TileIndexing='columnmajor', Padding='compact');
tl.yMargin = tiledlayout(tl.parent, n, 1, TileSpacing='none', Padding='none');
tl.yMargin.Layout.Tile = 1;
tl.yMargin.Layout.TileSpan = [marginRatio*n, 1];
tl.xMargin = tiledlayout(tl.parent, 1, n, TileSpacing='none', Padding='none');
tl.xMargin.Layout.Tile = 2*marginRatio*n + 2;
tl.xMargin.Layout.TileSpan = [1, marginRatio*n];
tl.center = tiledlayout(tl.parent, n, n, TileSpacing='none', Padding='none', TileIndexing='columnmajor');
tl.center.Layout.Tile = marginRatio*n + 2;
tl.center.Layout.TileSpan = [marginRatio*n, marginRatio*n];

for i = 1:n
    ii = n + 1 - i;
    y = meta.(metaNames(i));
    edges = -3:0.2:6;
    pdf = histcounts(y, edges, Normalization='pdf');
    centers = (edges(1:end-1) + edges(2:end)) / 2;
    ax = nexttile(tl.xMargin, i);
    bar(ax, centers, pdf, 1, EdgeColor='black', FaceColor='white');
    % set(ax, YDir='reverse', XAxisLocation='top')
    yticks(ax, [])
    xlim(ax, [-3, 6])
    
    ax = nexttile(tl.yMargin, ii);
    barh(ax, centers, pdf, 1, EdgeColor='black', FaceColor='white');
    xticks(ax, [])
    % set(ax, XDir='reverse', YAxisLocation='right')
    ylim(ax, [-3, 6])
    for j = 1:n
        ax = nexttile(tl.center, (j-1)*n + ii);
        hold(ax, 'on');
        x = meta.(metaNames(j));
        selSig = bootResp{i} & bootResp{j};
        scatter(ax, x(selSig), y(selSig), 5, 'black', 'filled');
        scatter(ax, x(~selSig), y(~selSig), 5, 'black');
        plot(ax, [-3, 6], [-3, 6], 'k--')
        xline(ax, 0, 'k--')
        yline(ax, 0, 'k--')
        % axis(ax, 'equal')
        xlim(ax, [-3, 6])
        ylim(ax, [-3, 6])
        hold(ax, 'off')
        xlabel(ax, dispNames{j}, Interpreter='tex');
        ylabel(ax, dispNames{i}, Interpreter='tex');

        if i > 1
            xticks(ax, [])
            xlabel(ax, '')
        end
        if j > 1
            yticks(ax, [])
            ylabel(ax, '')
        end
    end
end
fontsize(tl.parent, 9, 'points')
title(tl.parent, 'normalized to peri-reach [-4, -2]s')

clear metaNames i j ax n x y tl marginRatio selSig;

% Scatter META vs META
% close all
metaNames = ["pressNorm", "lickNorm", "correctReleaseNorm"];
dispNames = { ...
    ["Reach", sprintf("[%s]s", string(meta.window.press).join(', '))], ...
    ["Lick", sprintf("[%s]s", string(meta.window.lick).join(', '))], ...
    ["Release", sprintf("[%s]s", string(meta.window.correctRelease).join(', '))], ...
    };
n = length(metaNames);

clear tl
marginRatio = 2;
tl.parent = tiledlayout(figure(Units='inches', InnerPosition=[1 1 6 6]), marginRatio*n + 1, marginRatio*n + 1, TileSpacing='tight', TileIndexing='columnmajor', Padding='compact');
tl.yMargin = tiledlayout(tl.parent, n, 1, TileSpacing='none', Padding='none');
tl.yMargin.Layout.Tile = 1;
tl.yMargin.Layout.TileSpan = [marginRatio*n, 1];
tl.xMargin = tiledlayout(tl.parent, 1, n, TileSpacing='none', Padding='none');
tl.xMargin.Layout.Tile = 2*marginRatio*n + 2;
tl.xMargin.Layout.TileSpan = [1, marginRatio*n];
tl.center = tiledlayout(tl.parent, n, n, TileSpacing='none', Padding='none', TileIndexing='columnmajor');
tl.center.Layout.Tile = marginRatio*n + 2;
tl.center.Layout.TileSpan = [marginRatio*n, marginRatio*n];

for i = 1:n
    ii = n + 1 - i;
    y = meta.(metaNames(i));
    edges = -3:0.2:6;
    pdf = histcounts(y, edges, Normalization='pdf');
    centers = (edges(1:end-1) + edges(2:end)) / 2;
    ax = nexttile(tl.xMargin, i);
    bar(ax, centers, pdf, 1, EdgeColor='black', FaceColor='white');
    % set(ax, YDir='reverse', XAxisLocation='top')
    yticks(ax, [])
    xlim(ax, [-3, 6])
    
    ax = nexttile(tl.yMargin, ii);
    barh(ax, centers, pdf, 1, EdgeColor='black', FaceColor='white');
    xticks(ax, [])
    % set(ax, XDir='reverse', YAxisLocation='right')
    ylim(ax, [-3, 6])
    for j = 1:n
        ax = nexttile(tl.center, (j-1)*n + ii);
        hold(ax, 'on');
        x = meta.(metaNames(j));
        scatter(ax, x, y, 5, 'black');
        plot(ax, [-3, 6], [-3, 6], 'k--')
        xline(ax, 0, 'k--')
        yline(ax, 0, 'k--')
        % axis(ax, 'equal')
        xlim(ax, [-3, 6])
        ylim(ax, [-3, 6])
        hold(ax, 'off')
        xlabel(ax, dispNames{j}, Interpreter='tex');
        ylabel(ax, dispNames{i}, Interpreter='tex');

        if i > 1
            xticks(ax, [])
            xlabel(ax, '')
        end
        if j > 1
            yticks(ax, [])
            ylabel(ax, '')
        end
    end
end
fontsize(tl.parent, 9, 'points')
title(tl.parent, 'normalized independently')
clear metaNames i j ax n x y tl marginRatio;

%% Scatter3 animation
close all
fig = figure;
ax = axes(fig);
hold(ax, 'on')
plot3(ax, [-3, 6], [-3, 6], [-3, 6], 'k')

plot3(ax, [-3, 6], [0, 0], [0, 0], 'k--')
plot3(ax, [0, 0], [-3, 6], [0, 0], 'k--')
plot3(ax, [0, 0], [0, 0], [-3, 6], 'k--')

patch(ax, [-3, 6, 6, -3], [-3, -3, 6, 6], [0, 0, 0, 0], 'green', EdgeAlpha=1, LineStyle='--', FaceAlpha=0.6)
patch(ax, [-3, 6, 6, -3], [0, 0, 0, 0], [-3, -3, 6, 6], 'blue', EdgeAlpha=1, LineStyle='--', FaceAlpha=0.6)
patch(ax, [0, 0, 0, 0], [-3, 6, 6, -3], [-3, -3, 6, 6], 'red', EdgeAlpha=1, LineStyle='--', FaceAlpha=0.6)

selSig = cc.isPressResponsive & cc.isLickResponsive & cc.isCorrectReleaseResponsive;
scatter3(ax, meta.pressNorm(selSig), meta.lickNormToPress(selSig), meta.correctReleaseNormToPress(selSig), 10, 'black', 'filled', MarkerFaceAlpha=0.5);
scatter3(ax, meta.pressNorm(~selSig), meta.lickNormToPress(~selSig), meta.correctReleaseNormToPress(~selSig), 10, 'black', MarkerFaceAlpha=0.5);

axis(ax, 'equal')
xlabel(ax, 'Reach')
ylabel(ax, 'Lick')
zlabel(ax, 'Release')

filename = "D:\Downloads\New folder\scatter3.gif"; % Specify the output file name
if exist(filename, 'file')
    delete(filename)
end
angles = 45:405;
nImages = length(angles);
for idx = 1:nImages
    view(ax, [angles(idx), 30])
    drawnow();
    frame = getframe(fig);
    im{idx} = frame2im(frame);
end

for idx = 1:nImages
    [A, map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A, map, filename, "gif", "LoopCount", Inf, "DelayTime",1/30);
    else
        imwrite(A, map, filename, "gif", "WriteMode", "append", "DelayTime",1/30);
    end
end
clear filename fig ax angles idx fram im A map

%% Plot all units, PETH as heatmap: reach, release, lick
close all
clear etaOrder etaParams
etaParams.xlim = [-2, 2];
etaParams.clim = [-1.5, 1.5];

% Sorted independently, normed independently
tl = tiledlayout(figure(Units='inches', Position=[0 0 6 5]), 1, 3);
AX = gobjects(1, 3);
iAx = 1;

ax = nexttile(tl); AX(iAx) = ax; iAx = iAx + 1;
EphysUnit.plotETA(ax, eta.pressNorm, xlim=etaParams.xlim, clim=etaParams.clim, ...
    sortWindow=[-2, 0], signWindow=[-0.5, 0], sortThreshold=0.5, negativeSortThreshold=0.25, hideColorbar=true);
xline(ax, meta.window.press, 'k')
title(ax, 'Reach')
xlabel(ax, 'time to bar contact (s)')

ax = nexttile(tl); AX(iAx) = ax; iAx = iAx + 1;
EphysUnit.plotETA(ax, eta.lickNorm, xlim=etaParams.xlim, clim=etaParams.clim, ...
    sortWindow=[-2, 0], signWindow=[-0.5, 0], sortThreshold=0.5, negativeSortThreshold=0.25, hideColorbar=true);
xline(ax, meta.window.lick, 'k')
title(ax, 'Lick')
xlabel(ax, 'time to spout contact (s)')

ax = nexttile(tl); AX(iAx) = ax; iAx = iAx + 1;
EphysUnit.plotETA(ax, eta.correctReleaseNorm, xlim=etaParams.xlim, clim=etaParams.clim, ...
    sortWindow=[0, 1], signWindow=[0, 0.5], sortThreshold=0.5, negativeSortThreshold=0.25, hideColorbar=false);
xline(ax, meta.window.correctRelease, 'k')
title(ax, 'Release')
xlabel(ax, 'time to bar release (s)')

for iAx = 1:3
    applyCustomColormap(AX(iAx), [-1.5, 1.5], hlim=[0.375, 0, 0, -0.375], llim=[0.125, 0.5, 0.5, 0.25], hpwr=.5, lpwr=1, h0=0.33);
end

ax.Colorbar.Layout.Tile='east';
ylabel(AX, '')
title(tl, 'Ordered and normalized independently', FontWeight='bold')
fontsize(tl, 9, 'points')

% Sorted independently, normed to reach
tl = tiledlayout(figure(Units='inches', Position=[6 0 6 5]), 1, 3);
AX = gobjects(1, 3);
iAx = 1;

ax = nexttile(tl); AX(iAx) = ax; iAx = iAx + 1;
EphysUnit.plotETA(ax, eta.pressNorm, xlim=etaParams.xlim, clim=etaParams.clim, ...
    sortWindow=[-2, 0], signWindow=[-0.5, 0], sortThreshold=0.5, negativeSortThreshold=0.25, hideColorbar=true);
xline(ax, meta.window.press, 'k')
title(ax, 'Reach')
xlabel(ax, 'time to bar contact (s)')

ax = nexttile(tl); AX(iAx) = ax; iAx = iAx + 1;
EphysUnit.plotETA(ax, eta.lickNormToPress, xlim=etaParams.xlim, clim=etaParams.clim, ...
    sortWindow=[-2, 0], signWindow=[-0.5, 0], sortThreshold=0.5, negativeSortThreshold=0.25, hideColorbar=true);
xline(ax, meta.window.lick, 'k')
title(ax, 'Lick')
xlabel(ax, 'time to spout contact (s)')

ax = nexttile(tl); AX(iAx) = ax; iAx = iAx + 1;
EphysUnit.plotETA(ax, eta.correctReleaseNormToPress, xlim=etaParams.xlim, clim=etaParams.clim, ...
    sortWindow=[0, 1], signWindow=[0, 0.5], sortThreshold=0.5, negativeSortThreshold=0.25, hideColorbar=false);
xline(ax, meta.window.correctRelease, 'k')
title(ax, 'Release')
xlabel(ax, 'time to bar release (s)')

for iAx = 1:3
    applyCustomColormap(AX(iAx), [-1.5, 1.5], hlim=[0.375, 0, 0, -0.375], llim=[0.125, 0.5, 0.5, 0.25], hpwr=.5, lpwr=1, h0=0.33);
end

ax.Colorbar.Layout.Tile='east';
ylabel(AX, '')
title(tl, 'Ordered independently, normalized to peri-reach [-4, -2]s', FontWeight='bold')
fontsize(tl, 9, 'points')


% Sorted by reach, normed to reach
groupVar = zeros(size(eu));
groupVar(meta.pressNorm<0) = groupVar(meta.pressNorm<0) + 1;
groupVar(meta.lickNormToPress<0) = groupVar(meta.lickNormToPress<0) + 1;
groupVar(meta.correctReleaseNormToPress<0) = groupVar(meta.correctReleaseNormToPress<0) + 1;
groupVar(groupVar==0) = 4;

selJustOneDec = groupVar == 1;
groupVar(selJustOneDec & reshape(meta.pressNorm<0, size(groupVar))) = 1.0;
groupVar(selJustOneDec & reshape(meta.lickNormToPress<0, size(groupVar))) = 1.1;
groupVar(selJustOneDec & reshape(meta.correctReleaseNormToPress<0, size(groupVar))) = 1.2;

selJustTwoDec = groupVar == 2;
groupVar(selJustTwoDec & reshape(meta.pressNorm>0, size(groupVar))) = 2.0;
groupVar(selJustTwoDec & reshape(meta.lickNormToPress>0, size(groupVar))) = 2.1;
groupVar(selJustTwoDec & reshape(meta.correctReleaseNormToPress>0, size(groupVar))) = 2.2;


tl = tiledlayout(figure(Units='inches', Position=[0 6 6 5]), 1, 3);
AX = gobjects(1, 3);
iAx = 1;

ax = nexttile(tl); AX(iAx) = ax; iAx = iAx + 1;
[~, etaOrder.pressNormToPress] = EphysUnit.plotETA(ax, eta.pressNorm, xlim=etaParams.xlim, clim=etaParams.clim, ...
    sortGroup=groupVar, sortWindow=[-2, 0], signWindow=[-0.5, 0], sortThreshold=0.5, negativeSortThreshold=0.25, hideColorbar=true);
xline(ax, meta.window.press, 'k')
title(ax, 'Reach')
xlabel(ax, 'time to bar contact (s)')

ax = nexttile(tl); AX(iAx) = ax; iAx = iAx + 1;
EphysUnit.plotETA(ax, eta.lickNormToPress, xlim=etaParams.xlim, clim=etaParams.clim, ...
    order=etaOrder.pressNormToPress, hideColorbar=true);
xline(ax, meta.window.lick, 'k')
title(ax, 'Lick')
xlabel(ax, 'time to spout contact (s)')

ax = nexttile(tl); AX(iAx) = ax; iAx = iAx + 1;
EphysUnit.plotETA(ax, eta.correctReleaseNormToPress, xlim=etaParams.xlim, clim=etaParams.clim, ...
    order=etaOrder.pressNormToPress, hideColorbar=false);
xline(ax, meta.window.correctRelease, 'k')
title(ax, 'Release')
xlabel(ax, 'time to bar release (s)')

N = histcounts(groupVar, 1:5);
for iAx = 1:3
    applyCustomColormap(AX(iAx), [-1.5, 1.5], hlim=[0.375, 0, 0, -0.375], llim=[0.125, 0.5, 0.5, 0.25], hpwr=.5, lpwr=1, h0=0.33);
    yline(AX(iAx), cumsum(N(1:end-1)), '--', LineWidth=2, Color='magenta');
    yticks(AX(iAx), cumsum(N));
end

ax.Colorbar.Layout.Tile='east';
ylabel(AX, '')
title(tl, 'Same order, normalized to peri-reach [-4, -2]s', FontWeight='bold')
fontsize(tl, 9, 'points')


% Sorted by reach, normed independently

groupVar = zeros(size(eu));
groupVar(meta.pressNorm<0) = groupVar(meta.pressNorm<0) + 1;
groupVar(meta.lickNorm<0) = groupVar(meta.lickNorm<0) + 1;
groupVar(meta.correctReleaseNorm<0) = groupVar(meta.correctReleaseNorm<0) + 1;
groupVar(groupVar==0) = 4;

selJustOneDec = groupVar == 1;
groupVar(selJustOneDec & reshape(meta.pressNorm<0, size(groupVar))) = 1.0;
groupVar(selJustOneDec & reshape(meta.lickNorm<0, size(groupVar))) = 1.1;
groupVar(selJustOneDec & reshape(meta.correctReleaseNorm<0, size(groupVar))) = 1.2;

selJustTwoDec = groupVar == 2;
groupVar(selJustTwoDec & reshape(meta.pressNorm>0, size(groupVar))) = 2.0;
groupVar(selJustTwoDec & reshape(meta.lickNorm>0, size(groupVar))) = 2.1;
groupVar(selJustTwoDec & reshape(meta.correctReleaseNorm>0, size(groupVar))) = 2.2;

tl = tiledlayout(figure(Units='inches', Position=[6 6 6 5]), 1, 3);
AX = gobjects(1, 3);
iAx = 1;

ax = nexttile(tl); AX(iAx) = ax; iAx = iAx + 1;
[~, etaOrder.pressNormToSelf] = EphysUnit.plotETA(ax, eta.pressNorm, xlim=etaParams.xlim, clim=etaParams.clim, ...
    sortGroup=groupVar, sortWindow=[-2, 0], signWindow=[-0.5, 0], sortThreshold=0.5, negativeSortThreshold=0.25, hideColorbar=true);
xline(ax, meta.window.press, 'k')
title(ax, 'Reach')
xlabel(ax, 'time to bar contact (s)')

ax = nexttile(tl); AX(iAx) = ax; iAx = iAx + 1;
EphysUnit.plotETA(ax, eta.lickNorm, xlim=etaParams.xlim, clim=etaParams.clim, ...
    order=etaOrder.pressNormToSelf, hideColorbar=true);
xline(ax, meta.window.lick       , 'k')
title(ax, 'Lick')
xlabel(ax, 'time to spout contact (s)')

ax = nexttile(tl); AX(iAx) = ax; iAx = iAx + 1;
EphysUnit.plotETA(ax, eta.correctReleaseNorm, xlim=etaParams.xlim, clim=etaParams.clim, ...
    order=etaOrder.pressNormToSelf, hideColorbar=false);
xline(ax, meta.window.correctRelease, 'k')
title(ax, 'Release')
xlabel(ax, 'time to bar release (s)')

N = histcounts(groupVar, 1:5);
for iAx = 1:3
    applyCustomColormap(AX(iAx), [-1.5, 1.5], hlim=[0.375, 0, 0, -0.375], llim=[0.125, 0.5, 0.5, 0.25], hpwr=.5, lpwr=1, h0=0.33);
    yline(AX(iAx), cumsum(N(1:end-1)) + 1, '--', LineWidth=2, Color='magenta');
end

ax.Colorbar.Layout.Tile='east';
ylabel(AX, '')
title(tl, 'Same order, normalized independently', FontWeight='bold')
fontsize(tl, 9, 'points')



clear ax iAx AX tl N groupVar selJustOneDec selJustTwoDec

%% Heatmap with fancy sort

% Sorted by reach, normed to reach
groupVar = zeros(size(eu));
groupVar(cc.isPressDown) = groupVar(cc.isPressDown) + 1;
groupVar(cc.isLickDown) = groupVar(cc.isLickDown) + 1;
groupVar(cc.isCorrectReleaseDown) = groupVar(cc.isCorrectReleaseDown) + 1;
groupVar(groupVar==0) = 4;

selJustOneDec = groupVar == 1;
groupVar(selJustOneDec & reshape(cc.isPressDown, size(groupVar))) = 1.0;
groupVar(selJustOneDec & reshape(cc.isLickDown, size(groupVar))) = 1.1;
groupVar(selJustOneDec & reshape(cc.isCorrectReleaseDown, size(groupVar))) = 1.2;

selJustTwoDec = groupVar == 2;
groupVar(selJustTwoDec & reshape(~cc.isCorrectReleaseDown, size(groupVar))) = 2.0;
groupVar(selJustTwoDec & reshape(~cc.isLickDown, size(groupVar))) = 2.1;
groupVar(selJustTwoDec & reshape(~cc.isPressDown, size(groupVar))) = 2.2;

tl = tiledlayout(figure(Units='inches', Position=[0 1 6 5]), 1, 3);
AX = gobjects(1, 3);
iAx = 1;

ax = nexttile(tl); AX(iAx) = ax; iAx = iAx + 1;
[~, etaOrder.pressNormToPress] = EphysUnit.plotETA(ax, eta.pressNorm, xlim=etaParams.xlim, clim=etaParams.clim, ...
    sortGroup=groupVar, sortWindow=[-2, 0], signWindow=[-0.5, 0], sortThreshold=0.5, negativeSortThreshold=0.25, hideColorbar=true);
xline(ax, meta.window.press, 'k')
title(ax, 'Reach')
xlabel(ax, 'time to bar contact (s)')

ax = nexttile(tl); AX(iAx) = ax; iAx = iAx + 1;
EphysUnit.plotETA(ax, eta.lickNormToPress, xlim=etaParams.xlim, clim=etaParams.clim, ...
    order=etaOrder.pressNormToPress, hideColorbar=true);
xline(ax, meta.window.lick, 'k')
title(ax, 'Lick')
xlabel(ax, 'time to spout contact (s)')

ax = nexttile(tl); AX(iAx) = ax; iAx = iAx + 1;
EphysUnit.plotETA(ax, eta.correctReleaseNormToPress, xlim=etaParams.xlim, clim=etaParams.clim, ...
    order=etaOrder.pressNormToPress, hideColorbar=false);
xline(ax, meta.window.correctRelease, 'k')
title(ax, 'Release')
xlabel(ax, 'time to bar release (s)')

N = histcounts(groupVar, 1:5);
for iAx = 1:3
    applyCustomColormap(AX(iAx), [-1.5, 1.5], hlim=[0.375, 0, 0, -0.375], llim=[0.125, 0.5, 0.5, 0.25], hpwr=.5, lpwr=1, h0=0.33);
    yline(AX(iAx), cumsum(N(1:end-1)), '--', LineWidth=2, Color='magenta');
    yticks(AX(iAx), cumsum(N));
end

ax.Colorbar.Layout.Tile='east';
ylabel(AX, '')
title(tl, 'Same order, normalized to peri-reach [-4, -2]s', FontWeight='bold')
fontsize(tl, 9, 'points')



%% Plot individial units, PETH as trace: lick, reach, release
close all
path = 'E:\Figures\noArtifacts_pressAndLickBlancking_10ms\';
if ~exist(path, 'dir')
    mkdir(path)
end
% ETA = {eta.correctLick, eta.correctPress, eta.correctReleaseAfterRetract, eta.correctRelease; 
%     eta.incorrectLick, eta.incorrectPress, eta.correctReleaseBeforeRetract, eta.incorrectRelease};
% NAME = ["self-timed lick", "self-timed reach", "release (correct)", "release"];
% OUTCOME = ["correct", "correct", "after bar retract", "correct"; ...
%     "incorrect", "incorrect", "before bar retract", "incorrect"];
% XLIM = {[-2, 2], [-2, 2], [-1, 2], [-1, 2]};

ETA = {eta.correctLick, eta.correctPress, eta.correctReleaseAfterRetract; 
    eta.incorrectLick, eta.incorrectPress, eta.correctReleaseBeforeRetract};
NAME = ["self-timed lick", "self-timed reach", "release (correct)"];
OUTCOME = ["correct", "correct", "on retract"; ...
    "incorrect", "incorrect", "early"];
XLIM = {[-2, 2], [-2, 2], [-1, 2]};
ZEROLABEL = ["lick", "touch", "release"];

fig = figure(Units='inches', OuterPosition=[1, 1, 8, 3]);
tl = tiledlayout(fig, 1, sum(cellfun(@diff, XLIM)), TileSpacing='compact', Padding='tight');
ax = gobjects(1, length(XLIM));
h = gobjects(2, 1);
for i = 1:length(XLIM)
    ax(i) = nexttile(tl, [1, diff(XLIM{i})]);
end
for iEu = 0:length(eu)
    yl = [Inf, -Inf];
    for i = 1:length(XLIM)
        cla(ax(i));
        hold(ax(i), 'on')
        t = ETA{1, i}.t;
        sel = t > -2 & t < 1.9;

        if iEu == 0
            h(1) = plot(ax(i), ETA{1, i}.t(sel), mean(ETA{1, i}.X(:, sel), 1, 'omitnan'), Color=[0.2 0.6 0.2], DisplayName=sprintf("%s", OUTCOME(1, i)), LineWidth=1.5);
            h(2) = plot(ax(i), ETA{2, i}.t(sel), mean(ETA{2, i}.X(:, sel), 1, 'omitnan'), Color=[0.2 0.2 0.2], DisplayName=sprintf("%s", OUTCOME(2, i)), LineWidth=1.5);
            yl(1) = 20;
            yl(2) = 80;
        else
            h(1) = plot(ax(i), ETA{1, i}.t(sel), ETA{1, i}.X(iEu, sel), Color=[0.2 0.6 0.2], DisplayName=sprintf("%s (%i trials)", OUTCOME(1, i), ETA{1, i}.N(iEu)), LineWidth=1.5);
            h(2) = plot(ax(i), ETA{2, i}.t(sel), ETA{2, i}.X(iEu, sel), Color=[0.2 0.2 0.2], DisplayName=sprintf("%s (%i trials)", OUTCOME(2, i), ETA{2, i}.N(iEu)), LineWidth=1.5);
            yl(1) = min([yl(1), ETA{1, i}.X(iEu, sel), ETA{2, i}.X(iEu, sel)], [], 'omitnan');
            yl(2) = max([yl(2), ETA{1, i}.X(iEu, sel), ETA{2, i}.X(iEu, sel)], [], 'omitnan');
        end
        hold(ax(i), 'off')
        lgd = legend(ax(i), h, AutoUpdate=false, Location='northoutside');
        title(lgd, NAME(i))
        xlim(ax(i), XLIM{i})
        xline(ax(i), 0, '-')
        if iEu == 0
            yline(ax(i), mean([eta.pressNorm.stats.mean]./0.1), 'k:')
        else
            yline(ax(i), eta.pressNorm.stats(iEu).mean./0.1, 'k:')
        end
        xticks(ax(i), -4:4)
        xticklabels(ax(i), [string(-4:-1), ZEROLABEL(i), string(1:4)])
        xtickangle(ax(i), 0)
    end

    ylim(ax, yl + 0.05*diff(yl)*[-1, 1]);
    yticks(ax(2:end), [])

    xlabel(tl, 'Time to event (s)')
    ylabel(tl, 'sp/s')

    fontsize(tl, 9, 'points')

    if iEu == 0
        title(tl, sprintf("%i unit average", length(eu)), Interpreter='none', FontWeight='bold')
        print(fig, sprintf('%s\\grand_average_%i_units.png', path, length(eu)), '-dpng', '-r0')
    else
        title(tl, eu(iEu).getName(), Interpreter='none', FontWeight='bold')
        print(fig, sprintf('%s\\%s.png', path, eu(iEu).getName()), '-dpng', '-r0')
    end

end

clear fig tl ax ETA NAME OUTCOME iEu i lgd yl path h