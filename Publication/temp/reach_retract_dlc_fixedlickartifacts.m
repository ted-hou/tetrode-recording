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

%% Find (last lick after cue correct, tube retract correct, tube retract incorrect)
for iEu = 1:length(eu)
    eu(iEu).Trials.LickValid = eu(iEu).Trials.Lick(eu(iEu).Trials.Lick.duration() >= 2);
    eu(iEu).Trials.LickIncorrect = eu(iEu).Trials.Lick(eu(iEu).Trials.Lick.duration() < 4 & eu(iEu).Trials.Lick.duration() >= 2);
    eu(iEu).Trials.LickCorrect = eu(iEu).Trials.Lick(eu(iEu).Trials.Lick.duration() >= 4);

    % For each self-timed Lick, find the next cue (of any kind)
    firstLickToCueValid = Trial([eu(iEu).Trials.LickValid.Stop, Inf], eu(iEu).EventTimes.Cue, 'first');
    firstLickToCueIncorrect = Trial([eu(iEu).Trials.LickIncorrect.Stop, Inf], eu(iEu).EventTimes.Cue, 'first');
    firstLickToCueCorrect = Trial([eu(iEu).Trials.LickCorrect.Stop, Inf], eu(iEu).EventTimes.Cue, 'first');
    % Find the last lickOff between firstLick and nextCue
    [~, lickOff, trialIndices] = firstLickToCueValid.inTrial(eu(iEu).EventTimes.LickOff);
    [~, lastLickOff, ~] = unique(trialIndices, 'last');
    lastLickOff = lickOff(lastLickOff);
    eu(iEu).Trials.CueToLastLickOff = Trial([eu(iEu).Trials.LickValid.Start, Inf], lastLickOff(:)', 'last');

    % CueToLastLickOffIncorrect
    [~, lickOff, trialIndices] = firstLickToCueIncorrect.inTrial(eu(iEu).EventTimes.LickOff);
    [~, lastLickOff, ~] = unique(trialIndices, 'last');
    lastLickOff = lickOff(lastLickOff);
    eu(iEu).Trials.CueToLastLickOffIncorrect = Trial([eu(iEu).Trials.LickIncorrect.Start, Inf], lastLickOff(:)', 'last');

    % CueToLastLickOffCorrect
    [~, lickOff, trialIndices] = firstLickToCueCorrect.inTrial(eu(iEu).EventTimes.LickOff);
    [~, lastLickOff, ~] = unique(trialIndices, 'last');
    lastLickOff = lickOff(lastLickOff);
    eu(iEu).Trials.CueToLastLickOffCorrect = Trial([eu(iEu).Trials.LickCorrect.Start, Inf], lastLickOff(:)', 'last');

    % LickToLastLickOffIncorrect
    [~, lickOff, trialIndices] = firstLickToCueIncorrect.inTrial(eu(iEu).EventTimes.LickOff);
    [~, lastLickOff, ~] = unique(trialIndices, 'last');
    lastLickOff = lickOff(lastLickOff);
    eu(iEu).Trials.LickToLastLickOffIncorrect = Trial([eu(iEu).Trials.LickIncorrect.Stop, Inf], lastLickOff(:)', 'last');

    % LickToLastLickOffCorrect
    [~, lickOff, trialIndices] = firstLickToCueCorrect.inTrial(eu(iEu).EventTimes.LickOff);
    [~, lastLickOff, ~] = unique(trialIndices, 'last');
    lastLickOff = lickOff(lastLickOff);
    eu(iEu).Trials.LickToLastLickOffCorrect = Trial([eu(iEu).Trials.LickCorrect.Stop, Inf], lastLickOff(:)', 'last');


    % CorrectPressToFirstRewardLick
    % First lick after rewarded self-timed reach, exclude trials with cues in between (i.e. rewarded press, but did not lick before next cue)
    eu(iEu).Trials.CorrectPressToFirstRewardLick = Trial([eu(iEu).Trials.PressCorrect.Stop, Inf], eu(iEu).EventTimes.LickOn, 'first', eu(iEu).EventTimes.Cue);

    % For each self-timed reach, find the next cue (of any kind)
    correctPressToCue = Trial([eu(iEu).Trials.CorrectPressToFirstRewardLick.Start, Inf], eu(iEu).EventTimes.Cue, 'first');
    % CorrectPressToLastLickOff
    [~, lickOff, trialIndices] = correctPressToCue.inTrial(eu(iEu).EventTimes.LickOff);
    [~, lastLickOff, ~] = unique(trialIndices, 'last');
    lastLickOff = lickOff(lastLickOff);
    eu(iEu).Trials.CorrectPressToLastLickOff = Trial([eu(iEu).Trials.CorrectPressToFirstRewardLick.Start, Inf], lastLickOff(:)', 'last');

end
clear iEu firstLickToCueValid firstLickToCueIncorrect firstLickToCueCorrect lickOff trialIndices lastLickOff correctPressToCue

% 
% nReleaseTrials = arrayfun(@(eu) struct(release=length(eu.Trials.CueToLeverRelease), releaseCorrect=length(eu.Trials.CueToLeverReleaseCorrect), ...
%     releaseBefore=length(eu.Trials.CueToLeverReleaseBeforeRetractCorrect), releaseAfter=length(eu.Trials.RetractToLeverReleaseAfterRetractCorrect)), eu);
% 
% clear ac iEu leverReleaseTimes leverRetractTimes
% 
% %% Save results
% eu.save('C:\SERVER\Units\PressVsLick_ArtifactsRemoved\FixedEventsAndTrials');

%% Load results
eu = EphysUnit.load('C:\SERVER\Units\LickVsReach_FixedArtifacts');
load('C:\SERVER\Units\meta_LickVsReach_FixedArtifacts.mat');

%% TODO: Load deeplabcut (Derin did all the older sessions)


%% Source of spikes switcheroo
% 
% for iEu = 1:length(eu)
%     newSpikeTimes{iEu} = eu(iEu).SpikeTimes;
% end

for iEu = 1:length(eu)
    eu(iEu).SpikeTimes = oldSpikeTimes(iEu).data;
end

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

eta.correctReleaseBeforeRetract = eu.getETA('count', 'CueToLeverReleaseBeforeRetractCorrect', [-4, 4], alignTo='stop', ...
    normalize='none', resolution=0.1, artifacts=artifactParams);
eta.correctReleaseAfterRetract = eu.getETA('count', 'RetractToLeverReleaseAfterRetractCorrect', [-4, 4], alignTo='stop', ...
    normalize='none', resolution=0.1, artifacts=artifactParams);

eta.correctLick = eu.getETA('count', 'lick', [-4, 4], minTrialDuration=4, ...
    normalize='none', resolution=0.1, artifacts=artifactParams);
eta.incorrectLick = eu.getETA('count', 'lick', [-4, 4], minTrialDuration=2, maxTrialDuration=4, ...
    normalize='none', resolution=0.1, artifacts=artifactParams);

eta.correctLickLastLickOff = eu.getETA('count', 'CueToLastLickOffCorrect', [-4, 4], alignTo='stop', ...
    normalize='none', resolution=0.1, artifacts=artifactParams);

eta.correctPressFirstLick = eu.getETA('count', 'CorrectPressToFirstRewardLick', [-4, 4], alignTo='stop', normalize='none', resolution=0.1, artifacts=artifactParams);
eta.correctPressLastLickOff = eu.getETA('count', 'CorrectPressToLastLickOff', [-4, 4], alignTo='stop', normalize='none', resolution=0.1, artifacts=artifactParams);

% Convert to sp/s
for field = ["correctPress", "incorrectPress", "correctLick", "incorrectLick", "correctRelease", "correctReleaseBeforeRetract", "correctReleaseAfterRetract", "correctRetract", "incorrectRetract", ...
        "correctLickLastLickOff", "correctPressFirstLick", "correctPressLastLickOff"]
    eta.(field).X = eta.(field).X ./ 0.1;
end

clear field

eta.pressNorm = eu.getETA('count', 'press', [-4, 4], alignTo='stop', minTrialDuration=2, normalize=eta.baselineWindow.press, resolution=0.1, artifacts=artifactParams);
eta.lickNorm = eu.getETA('count', 'lick', [-4, 4], alignTo='stop', minTrialDuration=2, normalize=eta.pressNorm.stats, resolution=0.1, artifacts=artifactParams);
eta.correctReleaseNorm = eu.getETA('count', 'CueToLeverReleaseCorrect', [-4, 4], alignTo='stop', normalize=eta.pressNorm.stats, resolution=0.1, artifacts=artifactParams);

eta.lickNormToSelf = eu.getETA('count', 'lick', [-4, 4], alignTo='stop', minTrialDuration=2, normalize=eta.baselineWindow.lick, resolution=0.1, artifacts=artifactParams);
eta.correctReleaseNormToSelf = eu.getETA('count', 'CueToLeverReleaseCorrect', [-4, 4], alignTo='stop', normalize=eta.baselineWindow.release, resolution=0.1, artifacts=artifactParams);

eta.correctPressNorm = eu.getETA('count', 'press', [-4, 4], alignTo='stop', minTrialDuration=4, normalize=eta.pressNorm.stats, resolution=0.1, artifacts=artifactParams);
eta.correctLickNorm = eu.getETA('count', 'lick', [-4, 4], alignTo='stop', minTrialDuration=4, normalize=eta.pressNorm.stats, resolution=0.1, artifacts=artifactParams);
eta.correctReleaseBeforeRetractNorm = eu.getETA('count', 'CueToLeverReleaseBeforeRetractCorrect', [-4, 4], alignTo='stop', normalize=eta.pressNorm.stats, resolution=0.1, artifacts=artifactParams);
eta.correctReleaseAfterRetractNorm = eu.getETA('count', 'RetractToLeverReleaseAfterRetractCorrect', [-4, 4], alignTo='stop', normalize=eta.pressNorm.stats, resolution=0.1, artifacts=artifactParams);

eta.incorrectPressNorm = eu.getETA('count', 'press', [-4, 4], alignTo='stop', minTrialDuration=2, maxTrialDuration=4, normalize=eta.pressNorm.stats, resolution=0.1, artifacts=artifactParams);
eta.incorrectLickNorm = eu.getETA('count', 'lick', [-4, 4], alignTo='stop', minTrialDuration=2, maxTrialDuration=4, normalize=eta.pressNorm.stats, resolution=0.1, artifacts=artifactParams);

% Lick trials, last lickOff before next cue
eta.correctLickLastLickOffNorm = eu.getETA('count', 'CueToLastLickOffCorrect', [-4, 4], alignTo='stop', normalize=eta.pressNorm.stats, resolution=0.1, artifacts=artifactParams);
% etaArtiFree.correctLickToLastLickOffNorm = eu.getETA('count', 'LickToLastLickOffCorrect', [-4, 20], alignTo='start', normalize=etaArtiFree.pressNorm.stats, resolution=0.1, artifacts=artifactParams, includeInvalid=false);
eta.incorrectLickToLastLickOffNorm = eu.getETA('count', 'LickToLastLickOffIncorrect', [-4, 20], alignTo='start', normalize=eta.pressNorm.stats, resolution=0.1, artifacts=artifactParams, includeInvalid=false);

% Reach trials, correct, reward licking
eta.correctPressFirstLickNorm = eu.getETA('count', 'CorrectPressToFirstRewardLick', [-4, 4], alignTo='stop', normalize=eta.pressNorm.stats, resolution=0.1, artifacts=artifactParams);
eta.correctPressLastLickOffNorm = eu.getETA('count', 'CorrectPressToLastLickOff', [-4, 4], alignTo='stop', normalize=eta.pressNorm.stats, resolution=0.1, artifacts=artifactParams);

%%
clear etaWithArti
etaWithArti.artifactParams = [];
etaWithArti.baselineWindow = struct(press=[-4, -2], lick=[-4, -2], release=[-2, 0]);


etaWithArti.correctPress = eu.getETA('count', 'PressCorrect', [-4, 4], alignTo='stop', ...
    normalize='none', resolution=0.1, artifacts=[]);
etaWithArti.incorrectPress = eu.getETA('count', 'PressIncorrect', [-4, 4], alignTo='stop', ...
    normalize='none', resolution=0.1, artifacts=[]);

etaWithArti.correctRetract = eu.getETA('count', 'CueToLeverRetractCorrect', [-4, 4], alignTo='stop', ...
    normalize='none', resolution=0.1, artifacts=[]);
etaWithArti.incorrectRetract = eu.getETA('count', 'CueToLeverRetractIncorrect', [-4, 4], alignTo='stop', ...
    normalize='none', resolution=0.1, artifacts=[]);

etaWithArti.correctRelease = eu.getETA('count', 'CueToLeverReleaseCorrect', [-4, 4], alignTo='stop', ...
    normalize='none', resolution=0.1, artifacts=[]);

etaWithArti.correctReleaseBeforeRetract = eu.getETA('count', 'CueToLeverReleaseBeforeRetractCorrect', [-4, 4], alignTo='stop', ...
    normalize='none', resolution=0.1, artifacts=[]);
etaWithArti.correctReleaseAfterRetract = eu.getETA('count', 'RetractToLeverReleaseAfterRetractCorrect', [-4, 4], alignTo='stop', ...
    normalize='none', resolution=0.1, artifacts=[]);

etaWithArti.correctLick = eu.getETA('count', 'lick', [-4, 4], minTrialDuration=4, ...
    normalize='none', resolution=0.1, artifacts=[]);
etaWithArti.incorrectLick = eu.getETA('count', 'lick', [-4, 4], minTrialDuration=2, maxTrialDuration=4, ...
    normalize='none', resolution=0.1, artifacts=[]);

etaWithArti.correctLickLastLickOff = eu.getETA('count', 'CueToLastLickOffCorrect', [-4, 4], alignTo='stop', ...
    normalize='none', resolution=0.1, artifacts=[]);

etaWithArti.correctPressFirstLick = eu.getETA('count', 'CorrectPressToFirstRewardLick', [-4, 4], alignTo='stop', normalize='none', resolution=0.1, artifacts=[]);
etaWithArti.correctPressLastLickOff = eu.getETA('count', 'CorrectPressToLastLickOff', [-4, 4], alignTo='stop', normalize='none', resolution=0.1, artifacts=[]);

% Convert to sp/s
for field = ["correctPress", "incorrectPress", "correctLick", "incorrectLick", "correctRelease", "correctReleaseBeforeRetract", "correctReleaseAfterRetract", "correctRetract", "incorrectRetract", ...
        "correctLickLastLickOff", "correctPressFirstLick", "correctPressLastLickOff"]
    etaWithArti.(field).X = etaWithArti.(field).X ./ 0.1;
end

clear field

etaWithArti.pressNorm = eu.getETA('count', 'press', [-4, 4], alignTo='stop', minTrialDuration=2, normalize=etaWithArti.baselineWindow.press, resolution=0.1, artifacts=[]);
etaWithArti.lickNorm = eu.getETA('count', 'lick', [-4, 4], alignTo='stop', minTrialDuration=2, normalize=etaWithArti.pressNorm.stats, resolution=0.1, artifacts=[]);
etaWithArti.correctReleaseNorm = eu.getETA('count', 'CueToLeverReleaseCorrect', [-4, 4], alignTo='stop', normalize=etaWithArti.pressNorm.stats, resolution=0.1, artifacts=[]);

etaWithArti.lickNormToSelf = eu.getETA('count', 'lick', [-4, 4], alignTo='stop', minTrialDuration=2, normalize=etaWithArti.baselineWindow.lick, resolution=0.1, artifacts=[]);
etaWithArti.correctReleaseNormToSelf = eu.getETA('count', 'CueToLeverReleaseCorrect', [-4, 4], alignTo='stop', normalize=etaWithArti.baselineWindow.release, resolution=0.1, artifacts=[]);

etaWithArti.correctPressNorm = eu.getETA('count', 'press', [-4, 4], alignTo='stop', minTrialDuration=4, normalize=etaWithArti.pressNorm.stats, resolution=0.1, artifacts=[]);
etaWithArti.correctLickNorm = eu.getETA('count', 'lick', [-4, 4], alignTo='stop', minTrialDuration=4, normalize=etaWithArti.pressNorm.stats, resolution=0.1, artifacts=[]);
etaWithArti.correctReleaseBeforeRetractNorm = eu.getETA('count', 'CueToLeverReleaseBeforeRetractCorrect', [-4, 4], alignTo='stop', normalize=etaWithArti.pressNorm.stats, resolution=0.1, artifacts=[]);
etaWithArti.correctReleaseAfterRetractNorm = eu.getETA('count', 'RetractToLeverReleaseAfterRetractCorrect', [-4, 4], alignTo='stop', normalize=etaWithArti.pressNorm.stats, resolution=0.1, artifacts=[]);

etaWithArti.incorrectPressNorm = eu.getETA('count', 'press', [-4, 4], alignTo='stop', minTrialDuration=2, maxTrialDuration=4, normalize=etaWithArti.pressNorm.stats, resolution=0.1, artifacts=[]);
etaWithArti.incorrectLickNorm = eu.getETA('count', 'lick', [-4, 4], alignTo='stop', minTrialDuration=2, maxTrialDuration=4, normalize=etaWithArti.pressNorm.stats, resolution=0.1, artifacts=[]);

% Lick trials, last lickOff before next cue
etaWithArti.correctLickLastLickOffNorm = eu.getETA('count', 'CueToLastLickOffCorrect', [-4, 4], alignTo='stop', normalize=etaWithArti.pressNorm.stats, resolution=0.1, artifacts=[]);
% eta.correctLickToLastLickOffNorm = eu.getETA('count', 'LickToLastLickOffCorrect', [-4, 20], alignTo='start', normalize=eta.pressNorm.stats, resolution=0.1, artifacts=[], includeInvalid=false);
etaWithArti.incorrectLickToLastLickOffNorm = eu.getETA('count', 'LickToLastLickOffIncorrect', [-4, 20], alignTo='start', normalize=etaWithArti.pressNorm.stats, resolution=0.1, artifacts=[], includeInvalid=false);

% Reach trials, correct, reward licking
etaWithArti.correctPressFirstLickNorm = eu.getETA('count', 'CorrectPressToFirstRewardLick', [-4, 4], alignTo='stop', normalize=etaWithArti.pressNorm.stats, resolution=0.1, artifacts=[]);
etaWithArti.correctPressLastLickOffNorm = eu.getETA('count', 'CorrectPressToLastLickOff', [-4, 4], alignTo='stop', normalize=etaWithArti.pressNorm.stats, resolution=0.1, artifacts=[]);


%% Calculate META
clear meta
meta.window.press = [-0.3, 0.3];
meta.window.lick = [-0.3, 0.3];
meta.window.release = [-0.3, 0.3];
meta.window.lastLickOff = [-0.3, 0.3];
meta.window.correctPressFirstLick = [-0.3, 0.3];

meta.localBaselineWindow.press = [-4, 2];
meta.localBaselineWindow.lick = [-4, 2];
meta.localBaselineWindow.release = [-1, -0.3];
meta.localBaselineWindow.lastLickOff = [-1, -0.3];
meta.localBaselineWindow.correctPressFirstLick = [-1, -0.3];

meta.pressNorm = mean(eta.pressNorm.X(:, eta.pressNorm.t > meta.window.press(1) & eta.pressNorm.t < meta.window.press(2)), 2, 'omitnan');
meta.lickNorm = mean(eta.lickNorm.X(:, eta.lickNorm.t > meta.window.lick(1) & eta.lickNorm.t < meta.window.lick(2)), 2, 'omitnan');
meta.incorrectLickNorm = mean(eta.incorrectLickNorm.X(:, eta.incorrectLickNorm.t > meta.window.lick(1) & eta.incorrectLickNorm.t < meta.window.lick(2)), 2, 'omitnan');
meta.correctLickNorm = mean(eta.correctLickNorm.X(:, eta.correctLickNorm.t > meta.window.lick(1) & eta.correctLickNorm.t < meta.window.lick(2)), 2, 'omitnan');
meta.correctReleaseNorm = mean(eta.correctReleaseNorm.X(:, eta.correctReleaseNorm.t > meta.window.release(1) & eta.correctReleaseNorm.t < meta.window.release(2)), 2, 'omitnan');
meta.correctLickLastLickOffNorm = mean(eta.correctLickLastLickOffNorm.X(:, eta.correctLickLastLickOffNorm.t > meta.window.lastLickOff(1) & eta.correctLickLastLickOffNorm.t < meta.window.lastLickOff(2)), 2, 'omitnan');
meta.correctPressLastLickOffNorm = mean(eta.correctPressLastLickOffNorm.X(:, eta.correctPressLastLickOffNorm.t > meta.window.lastLickOff(1) & eta.correctPressLastLickOffNorm.t < meta.window.lastLickOff(2)), 2, 'omitnan');
meta.correctPressFirstLickNorm = mean(eta.correctPressFirstLickNorm.X(:, eta.correctPressFirstLickNorm.t > meta.window.correctPressFirstLick(1) & eta.correctPressFirstLickNorm.t < meta.window.correctPressFirstLick(2)), 2, 'omitnan');

selBaseline = eta.pressNorm.t < meta.localBaselineWindow.press(2) & eta.pressNorm.t > meta.localBaselineWindow.press(1); meta.localNorm.pressNorm = (meta.pressNorm - mean(eta.pressNorm.X(:, selBaseline), 2, 'omitnan')) ./ std(eta.pressNorm.X(:, selBaseline), 1, 2, 'omitnan');
selBaseline = eta.lickNorm.t < meta.localBaselineWindow.lick(2) & eta.lickNorm.t > meta.localBaselineWindow.lick(1); meta.localNorm.lickNorm = (meta.lickNorm - mean(eta.lickNorm.X(:, selBaseline), 2, 'omitnan')) ./ std(eta.lickNorm.X(:, selBaseline), 1, 2, 'omitnan');
selBaseline = eta.incorrectLickNorm.t < meta.localBaselineWindow.lick(2) & eta.incorrectLickNorm.t > meta.localBaselineWindow.lick(1); meta.localNorm.incorrectLickNorm = (meta.incorrectLickNorm - mean(eta.incorrectLickNorm.X(:, selBaseline), 2, 'omitnan')) ./ std(eta.incorrectLickNorm.X(:, selBaseline), 1, 2, 'omitnan');
selBaseline = eta.correctLickNorm.t < meta.localBaselineWindow.lick(2) & eta.correctLickNorm.t > meta.localBaselineWindow.lick(1); meta.localNorm.correctLickNorm = (meta.correctLickNorm - mean(eta.correctLickNorm.X(:, selBaseline), 2, 'omitnan')) ./ std(eta.correctLickNorm.X(:, selBaseline), 1, 2, 'omitnan');
selBaseline = eta.correctReleaseNorm.t < meta.localBaselineWindow.release(2) & eta.correctReleaseNorm.t > meta.localBaselineWindow.release(1); meta.localNorm.correctReleaseNorm = (meta.correctReleaseNorm - mean(eta.correctReleaseNorm.X(:, selBaseline), 2, 'omitnan')) ./ std(eta.correctReleaseNorm.X(:, selBaseline), 1, 2, 'omitnan');
selBaseline = eta.correctLickLastLickOffNorm.t < meta.localBaselineWindow.lastLickOff(2) & eta.correctLickLastLickOffNorm.t > meta.localBaselineWindow.lastLickOff(1); meta.localNorm.correctLickLastLickOffNorm = (meta.correctLickLastLickOffNorm - mean(eta.correctLickLastLickOffNorm.X(:, selBaseline), 2, 'omitnan')) ./ std(eta.correctLickLastLickOffNorm.X(:, selBaseline), 1, 2, 'omitnan');
selBaseline = eta.correctPressLastLickOffNorm.t < meta.localBaselineWindow.lastLickOff(2) & eta.correctPressLastLickOffNorm.t > meta.localBaselineWindow.lastLickOff(1); meta.localNorm.correctPressLastLickOffNorm = (meta.correctPressLastLickOffNorm - mean(eta.correctPressLastLickOffNorm.X(:, selBaseline), 2, 'omitnan')) ./ std(eta.correctPressLastLickOffNorm.X(:, selBaseline), 1, 2, 'omitnan');
selBaseline = eta.correctPressFirstLickNorm.t < meta.localBaselineWindow.correctPressFirstLick(2) & eta.correctPressFirstLickNorm.t > meta.localBaselineWindow.correctPressFirstLick(1); meta.localNorm.correctPressFirstLickNorm = (meta.correctPressFirstLickNorm - mean(eta.correctPressFirstLickNorm.X(:, selBaseline), 2, 'omitnan')) ./ std(eta.correctPressFirstLickNorm.X(:, selBaseline), 1, 2, 'omitnan');
clear selBaseline


%% Bootstrap response directions
clear boot
boot.params = struct(...
    alpha = 0.01, ...
    nboot = 100000, ...
    responseWindow = meta.window, ...
    baselineWindow = struct(press=[-4, -2], lick='press', release='press', lastLickOff='press'), ...
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
    responseWindow=boot.params.responseWindow.release, ...
    alignTo=struct(baseline='stop', response='stop'), ...
    allowedTrialDuration=struct(baseline=[2, Inf], response=[4, Inf]), ...
    alpha=boot.params.alpha, artifacts=boot.params.artifacts);  

[boot.correctLickLastLickOff.h, boot.correctLickLastLickOff.ci, boot.correctLickLastLickOff.obs] = bootstrapMoveResponse(eu, ...
    struct(baseline='press', response='CueToLastLickOffCorrect'), nboot=boot.params.nboot, ...
    baselineWindow=boot.params.baselineWindow.press, ...
    responseWindow=boot.params.responseWindow.lastLickOff, ...
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
cc.isCorrectLastLickOffUp = boot.correctLickLastLickOff.h > 0 & ~isnan(boot.correctLickLastLickOff.h);
cc.isCorrectLastLickOffDown = boot.correctLickLastLickOff.h < 0 & ~isnan(boot.correctLickLastLickOff.h);

cc.isPressResponsive = cc.isPressUp | cc.isPressDown;
cc.isLickResponsive = cc.isLickUp | cc.isLickDown;
cc.isCorrectReleaseResponsive = cc.isCorrectReleaseUp | cc.isCorrectReleaseDown;
cc.isCorrectLastLickOffResponsive = cc.isCorrectLastLickOffUp | cc.isCorrectLastLickOffDown;

fprintf('\n%i total units, %i (%.0f%%) reach modulated, %i (%.0f%%) lick modulated, %i (%.0f%%) correct-release modulated, %i (%.0f%%) correct-lastLickOff modulated:\n', ...
    length(eu), ...
    nnz(cc.isPressResponsive), nnz(cc.isPressResponsive)/length(eu)*100, ...
    nnz(cc.isLickResponsive), nnz(cc.isLickResponsive)/length(eu)*100, ...
    nnz(cc.isCorrectReleaseResponsive), nnz(cc.isCorrectReleaseResponsive)/length(eu)*100, ...
    nnz(cc.isCorrectLastLickOffResponsive), nnz(cc.isCorrectLastLickOffResponsive)/length(eu)*100)
fprintf('\t%i/%i units reach inc (%.0f%%);\n', nnz(cc.isPressUp), nnz(cc.isPressResponsive), nnz(cc.isPressUp)./nnz(cc.isPressResponsive)*100);
fprintf('\t%i/%i units reach dec (%.0f%%);\n', nnz(cc.isPressDown), nnz(cc.isPressResponsive), nnz(cc.isPressDown)./nnz(cc.isPressResponsive)*100);
fprintf('\t%i/%i units lick inc (%.0f%%);\n', nnz(cc.isLickUp), nnz(cc.isLickResponsive), nnz(cc.isLickUp)./nnz(cc.isLickResponsive)*100);
fprintf('\t%i/%i units lick dec (%.0f%%);\n', nnz(cc.isLickDown), nnz(cc.isLickResponsive), nnz(cc.isLickDown)./nnz(cc.isLickResponsive)*100);
fprintf('\t%i/%i units release inc (%.0f%%);\n', nnz(cc.isCorrectReleaseUp), nnz(cc.isCorrectReleaseResponsive), nnz(cc.isCorrectReleaseUp)./nnz(cc.isCorrectReleaseResponsive)*100);
fprintf('\t%i/%i units release dec (%.0f%%).\n', nnz(cc.isCorrectReleaseDown), nnz(cc.isCorrectReleaseResponsive), nnz(cc.isCorrectReleaseDown)./nnz(cc.isCorrectReleaseResponsive)*100);
fprintf('\t%i/%i units lastLickOff inc (%.0f%%);\n', nnz(cc.isCorrectLastLickOffUp), nnz(cc.isCorrectLastLickOffResponsive), nnz(cc.isCorrectLastLickOffUp)./nnz(cc.isCorrectLastLickOffResponsive)*100);
fprintf('\t%i/%i units lastLickOff dec (%.0f%%).\n', nnz(cc.isCorrectLastLickOffDown), nnz(cc.isCorrectLastLickOffResponsive), nnz(cc.isCorrectLastLickOffDown)./nnz(cc.isCorrectLastLickOffResponsive)*100);

sel = cc.isPressResponsive & cc.isLickResponsive & cc.isCorrectReleaseResponsive & cc.isCorrectLastLickOffResponsive;
fprintf('\n%i/%i (%.0f%%) units modulated for reach, lick, and correct-release, of which:\n', nnz(sel), length(eu), nnz(sel)/length(eu)*100)
fprintf('\t%i/%i (%.0f%%) decrease for all four;\n', nnz(sel & cc.isPressDown & cc.isLickDown & cc.isCorrectReleaseDown & cc.isCorrectLastLickOffDown), nnz(sel), nnz(sel & cc.isPressDown & cc.isLickDown & cc.isCorrectReleaseDown & cc.isCorrectLastLickOffDown)./nnz(sel)*100)
fprintf('\t%i/%i (%.0f%%) increase for all four;\n\n', nnz(sel & cc.isPressUp & cc.isLickUp & cc.isCorrectReleaseUp & cc.isCorrectLastLickOffUp), nnz(sel), nnz(sel & cc.isPressUp & cc.isLickUp & cc.isCorrectReleaseUp & cc.isCorrectLastLickOffUp)./nnz(sel)*100)
fprintf('\t%i/%i (%.0f%%) decrease for just reach;\n', nnz(sel & cc.isPressDown & cc.isLickUp & cc.isCorrectReleaseUp & cc.isCorrectLastLickOffUp), nnz(sel), nnz(sel & cc.isPressDown & cc.isLickUp & cc.isCorrectReleaseUp & cc.isCorrectLastLickOffUp)./nnz(sel)*100)
fprintf('\t%i/%i (%.0f%%) decrease for just lick;\n', nnz(sel & cc.isPressUp & cc.isLickDown & cc.isCorrectReleaseUp & cc.isCorrectLastLickOffUp), nnz(sel), nnz(sel & cc.isPressUp & cc.isLickDown & cc.isCorrectReleaseUp & cc.isCorrectLastLickOffUp)./nnz(sel)*100)
fprintf('\t%i/%i (%.0f%%) decrease for just release;\n', nnz(sel & cc.isPressUp & cc.isLickUp & cc.isCorrectReleaseDown & cc.isCorrectLastLickOffUp), nnz(sel), nnz(sel & cc.isPressUp & cc.isLickUp & cc.isCorrectReleaseDown & cc.isCorrectLastLickOffUp)./nnz(sel)*100)
fprintf('\t%i/%i (%.0f%%) decrease for just lastLickOff;\n', nnz(sel & cc.isPressUp & cc.isLickUp & cc.isCorrectReleaseUp & cc.isCorrectLastLickOffDown), nnz(sel), nnz(sel & cc.isPressUp & cc.isLickUp & cc.isCorrectReleaseUp & cc.isCorrectLastLickOffDown)./nnz(sel)*100)

%% Save intermediate results because ETA with artifact blanking takes forever
% eu.save('C:\SERVER\Units\LickVsReach_FixedArtifacts')
save('C:\SERVER\Units\meta_LickVsReach_FixedArtifacts.mat', 'artifactParams', 'bootCLick', 'cc', 'circlick', 'eta', 'etaWithArti', 'oldSpikeTimes', 'newSpikeTimes', 'meta', 'boot')

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

%% Spike-sorting approach to categorize PETH shapes
ETANAME = ["correctPressNorm", "incorrectPressNorm", "correctLickNorm", "incorrectLickNorm", "correctPressFirstLickNorm", "correctLickLastLickOffNorm", "correctPressLastLickOffNorm" ,"correctReleaseNorm"];
ETA = {eta.correctPressNorm, eta.incorrectPressNorm, eta.correctLickNorm, eta.incorrectLickNorm, eta.correctPressFirstLickNorm, eta.correctLickLastLickOffNorm, eta.correctPressLastLickOffNorm, eta.correctReleaseNorm};
NAME = ["Reach(correct)", "Reach(incorect)", "Lick(correct)", "Lick(incorect)", "FirstLick(Reach)", "LastLick(Lick)", "LastLick(Reach)", "Release(Reach)"];
ZEROLABEL = ["touch", "touch", "lick", "lick", "lick", "lick-off", "lick-off", "release"];
WINDOW = {[-0.4, 0.4], [-0.4, 0.4], [-0.4, 0.4], [-0.4, 0.4], [-0.3, 0.3], [-0.3, 0.3], [-0.3, 0.3], [-0.3, 0.3]};
XLIM = {[-2, 2], [-2, 2], [-2, 2], [-2, 2], [-1, 1], [-1, 1], [-1, 1], [-1, 1]};
USEFORGROUPVAR = [false, true, false, true, false, false, false, true];
% plotMode = 'mean+sd'; % mean, raw
plotMode = 'mean+raw'; % mean, raw
nDims = 5;
K = [6, 6, 6, 6, 6, 6, 6, 6];

close all
clear cluster clusterOrder
for iFig = 1:length(ETA)
    k = K(iFig);
    etaTemp = ETA{iFig};
    X = etaTemp.X(:, isin(etaTemp.t, WINDOW{iFig}));
    t = etaTemp.t(isin(etaTemp.t, WINDOW{iFig}));
    % interp over nans
    for i = 1:size(X, 1)
        selnan = isnan(X(i, :));
        X(i, selnan) = interp1(t(~selnan), X(i, ~selnan), t(selnan), 'linear', 'extrap');
        assert(all(~isnan(X(i, :))))
    end
    [coeff, score, ~, ~, explained] = pca(X);
    [clusterTemp, C] = kmeans(score(:, 1:nDims), k);
    % [cluster, C] = kmeans(X, k);
    % Sort clusters by size
    clusterSize = zeros(1, k);
    for iCluster = 1:k
        clusterSize(iCluster) = nnz(clusterTemp==iCluster);
    end
    [~, clusterOrderTemp] = sort(clusterSize, 'descend');
    cluster.(ETANAME(iFig)) = clusterTemp;
    clusterOrder.(ETANAME(iFig)) = clusterOrderTemp;
    
    
    
    w = [2, 3, 5];
    cw = cumsum([1, w]);
    fig = figure(Units='inches', Position=[1, 1, 15, 5]);
    tl = tiledlayout(fig, 1, sum(w));
    ax = gobjects(1, 2);
    ax(1) = nexttile(tl, cw(1), [1, w(1)]);
    ax(2) = nexttile(tl, cw(2), [1, w(2)]);
    ax(3) = nexttile(tl, cw(3), [1, w(3)]);
    hold(ax, 'on')
    
    h = gobjects(1, nDims);
    for iDim = 1:nDims
        h(iDim) = plot(ax(1), t, coeff(:, iDim), LineWidth=1.5, Color=[getColor(iDim, nDims, 0.7), 0.7], DisplayName=sprintf("pc%i (%i%%)", iDim, round(explained(iDim))));
    end
    xline(ax(1), 0, 'k--')
    yline(ax(1), 0, 'k--')
    legend(ax(1), h, Location='northoutside')
    
    h = gobjects(k, 1);
    for iCluster = 1:k
        color = getColor(iCluster, k, 0.7);
        selCluster = clusterTemp==clusterOrderTemp(iCluster);
        scatter3(ax(2), score(selCluster, 1), score(selCluster, 2), score(selCluster, 3), 15, color, 'filled');
    
        x = etaTemp.X(selCluster, :);
        tt = etaTemp.t;
        switch plotMode
            case 'mean'
                mu = mean(x, 1, 'omitnan');
                h(iCluster) = plot(ax(3), tt, mu, LineWidth=1.5, Color=color, DisplayName=sprintf('cluster%i (%i%%, n=%i)', iCluster, round(100*nnz(selCluster)/length(eu)), nnz(selCluster)));
            case 'mean+sd'
                mu = mean(x, 1, 'omitnan');
                sd = std(x, 1, 'omitnan');
                h(iCluster) = plot(ax(3), tt, mu, LineWidth=1.5, Color=color, DisplayName=sprintf('cluster%i (%i%%, n=%i)', iCluster, round(100*nnz(selCluster)/length(eu)), nnz(selCluster)));
                patch(ax(3), [tt, flip(tt)], [mu+sd, flip(mu-sd)], color, FaceAlpha=0.1, EdgeColor='none')
            case 'mean+qt'
                mu = mean(x, 1, 'omitnan');
                qt = quantile(x, [0.05, 0.95], 1);
                h(iCluster) = plot(ax(3), tt, mu, LineWidth=1.5, Color=color, DisplayName=sprintf('cluster%i (%i%%, n=%i)', iCluster, round(100*nnz(selCluster)/length(eu)), nnz(selCluster)));
                patch(ax(3), [tt, flip(tt)], [mu+qt(1, :), flip(mu-qt(2, :))], color, FaceAlpha=0.1, EdgeColor='none')
            case 'raw'
                hh = plot(ax(3), tt, x, LineWidth=0.5, Color=[color, 0.1], DisplayName=sprintf('cluster%i (%i%%, n=%i)', iCluster, round(100*nnz(selCluster)/length(eu)), nnz(selCluster)));
                h(iCluster) = hh(1);
            case 'mean+raw'
                mu = mean(x, 1, 'omitnan');
                h(iCluster) = plot(ax(3), tt, mu, LineWidth=1.5, Color=color, DisplayName=sprintf('cluster%i (%i%%, n=%i)', iCluster, round(100*nnz(selCluster)/length(eu)), nnz(selCluster)));
                plot(ax(3), tt, x, LineWidth=0.5, Color=[color, 0.1], DisplayName=sprintf('cluster%i (%i%%, n=%i)', iCluster, round(100*nnz(selCluster)/length(eu)), nnz(selCluster)));
        end
        xline(ax(3), 0, 'k--')
        yline(ax(3), 0, 'k--')
    end
    hold(ax, 'off')
    axis(ax(2), 'equal')
    xlabel(ax(2), 'pc1')
    ylabel(ax(2), 'pc2')
    zlabel(ax(2), 'pc3')
    xlabel(ax(3), 'time (s)')
    ylabel(ax(3), 'normalized spike rate (a.u.)')
    xlim(ax(3), XLIM{iFig})
    xticks(ax(3), -1:1)
    xticklabels(ax(3), ["-1", ZEROLABEL(iFig), "1"])
    ylim(ax(3), [-3, 6])
    ax(2).XAxis.Color = getColor(1, nDims, 0.7);
    ax(2).YAxis.Color = getColor(2, nDims, 0.7);
    ax(2).ZAxis.Color = getColor(3, nDims, 0.7);
    legend(ax(3), h, Location='eastoutside')
    title(tl, NAME(iFig), FontSize=9, FontWeight='bold')
end

% Convert to groupvar
ii = 0;
base = 10;
assert(all(K<base));
groupVar = zeros(length(eu), 1);

for i = find(USEFORGROUPVAR)
    groupVar = groupVar + 10^ii*(cluster.(ETANAME(i)) - 1);
    ii = ii + 1;
end


%% Scatter META vs META
% close all
metaNames = ["pressNorm", "lickNorm", "correctReleaseNorm"];
bootResp = {cc.isPressResponsive, cc.isLickResponsive, cc.isCorrectReleaseResponsive};
dispNames = { ...
    ["Reach", sprintf("[%s]s", string(meta.window.press).join(', '))], ...
    ["Lick", sprintf("[%s]s", string(meta.window.lick).join(', '))], ...
    ["Release", sprintf("[%s]s", string(meta.window.release).join(', '))], ...
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
scatter3(ax, meta.pressNorm(selSig), meta.lickNorm(selSig), meta.correctReleaseNorm(selSig), 10, 'black', 'filled', MarkerFaceAlpha=0.5);
scatter3(ax, meta.pressNorm(~selSig), meta.lickNorm(~selSig), meta.correctReleaseNorm(~selSig), 10, 'black', MarkerFaceAlpha=0.5);

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

%% Expanded heatmap (correct/incorrect reach/lick/barrelease/tuberelease with fancy sort

path = 'E:\Figures\reach_lick_firstLickReach_lastLickLick_lastLickReach_release_fancygrouping';
useLocalNorm = false;
if useLocalNorm
    path = sprintf("%s\\useLocalNorm", path);
end
if ~exist(path, 'dir')
    mkdir(path)
end

close all
% Whole timecourse
% ETA = ["incorrectPressNormToPress", "correctPressNorm", "correctReleaseNormToPress", "incorrectLickNormToPress", "correctLickNormToPress", "correctLastLickOffNormToPress"];
% WINDOW = ["press", "press", "release", "lick", "lick", "lastLickOff"];
% TITLE = ["Reach (incorrect)", "Reach (correct)", "Release (correct)", "Lick (incorrect)", "Lick (correct)", "LastLickOff (correct)"];
% EVENT = ["bar-contact", "bar-contact", "bar-release", "spout-contact", "spout-contact", "spout-disengage"];
% XLIM = {[-2, 1], [-2, 1], [-1, 2], [-2, 1], [-2, 1], [-1, 2]};

% Just meta window
ETA = [ ...
    "correctPressNorm",     "correctLickNorm",  "correctPressFirstLickNorm", "correctLickLastLickOffNorm", "correctPressLastLickOffNorm", "correctReleaseNorm"; ...
    "incorrectPressNorm",   "incorrectLickNorm", "", "", "", ""];
WINDOW = [ ...
    "press", "lick", "correctPressFirstLick", "lastLickOff", "lastLickOff", "release"; ...
    "press", "lick", "correctPressFirstLick", "lastLickOff", "lastLickOff", "release" ...
    ];
TITLE = [ ...
    "Reach (correct)", "Lick (correct)", "FirstLick(Reach)", "LastLick(Lick)", "LastLick(Reach)", "Release(Reach)"; ...
    "Reach (incorrect)", "Lick (incorrect)", "", "", "", "" ...
    ];
EVENT = [ ...
    "bar-contact", "lick", "lick", "lick-off", "lick-off", "release"; ...
    "bar-contact", "lick", "", "", "", ""];
XLIM = {[-2, 2], [-2, 2], [-1, 1], [-1, 1], [-1, 1], [-1, 1]};
nCols = size(ETA, 2);
nRows = size(ETA, 1);

% % Sorted by reach, normed to reach
% groupVar = zeros(size(eu));
% groupVar(meta.pressNorm<0) = groupVar(meta.pressNorm<0) + 1;
% groupVar(meta.correctReleaseNormToPress<0) = groupVar(meta.correctReleaseNormToPress<0) + 1;
% groupVar(meta.lickNormToPress<0) = groupVar(meta.lickNormToPress<0) + 1;
% groupVar(meta.correctLastLickOffNormToPress<0) = groupVar(meta.correctLastLickOffNormToPress<0) + 1;
% groupVar(groupVar==0) = 5;
% 
% selJustOneDec = groupVar == 1;
% groupVar(selJustOneDec & reshape(meta.pressNorm<0, size(groupVar))) = 1.1;
% groupVar(selJustOneDec & reshape(meta.correctReleaseNormToPress<0, size(groupVar))) = 1.2;
% groupVar(selJustOneDec & reshape(meta.lickNormToPress<0, size(groupVar))) = 1.3;
% groupVar(selJustOneDec & reshape(meta.correctLastLickOffNormToPress<0, size(groupVar))) = 1.4;
% 
% N = histcounts(groupVar, 1:6);
% N1 = histcounts(groupVar, 1:0.1:1.5);

% Alt sorting groups
if useLocalNorm
    groupVar = zeros(size(eu));
    groupVar(meta.localNorm.pressNorm<0) = groupVar(meta.localNorm.pressNorm<0) + 1;
    groupVar(meta.localNorm.lickNorm<0) = groupVar(meta.localNorm.lickNorm<0) + 1;
    groupVar(meta.localNorm.correctLickLastLickOffNorm<0) = groupVar(meta.localNorm.correctLickLastLickOffNorm<0) + 1;
    groupVar(meta.localNorm.correctReleaseNorm<0) = groupVar(meta.localNorm.correctReleaseNorm<0) + 1;
    groupVar(groupVar==0) = 5;
    
    selJustOneDec = groupVar == 1;
    groupVar(selJustOneDec & reshape(meta.localNorm.pressNorm<0, size(groupVar))) = 1.1;
    groupVar(selJustOneDec & reshape(meta.localNorm.lickNorm<0, size(groupVar))) = 1.2;
    groupVar(selJustOneDec & reshape(meta.localNorm.correctLickLastLickOffNorm<0, size(groupVar))) = 1.3;
    groupVar(selJustOneDec & reshape(meta.localNorm.correctReleaseNorm<0, size(groupVar))) = 1.4;
else
    groupVar = zeros(size(eu));
    groupVar(meta.pressNorm<0) = groupVar(meta.pressNorm<0) + 1;
    groupVar(meta.lickNorm<0) = groupVar(meta.lickNorm<0) + 1;
    groupVar(meta.correctLickLastLickOffNorm<0) = groupVar(meta.correctLickLastLickOffNorm<0) + 1;
    groupVar(meta.correctReleaseNorm<0) = groupVar(meta.correctReleaseNorm<0) + 1;
    groupVar(groupVar==0) = 5;
    
    selJustOneDec = groupVar == 1;
    groupVar(selJustOneDec & reshape(meta.pressNorm<0, size(groupVar))) = 1.1;
    groupVar(selJustOneDec & reshape(meta.lickNorm<0, size(groupVar))) = 1.2;
    groupVar(selJustOneDec & reshape(meta.correctLickLastLickOffNorm<0, size(groupVar))) = 1.3;
    groupVar(selJustOneDec & reshape(meta.correctReleaseNorm<0, size(groupVar))) = 1.4;
end

groupNames = dictionary([], string([]));
groupNames(1.1) = "reach-dec";
groupNames(1.2) = "lick-dec";
groupNames(1.3) = "last-lick-dec";
groupNames(1.4) = "release-dec";
groupNames(2) = "2xdec";
groupNames(3) = "3xdec";
groupNames(4) = "4xdec";
groupNames(5) = "all-inc";

N = histcounts(groupVar, 1:6);
N1 = histcounts(groupVar, 1:0.1:1.5);


w = round(10*(cellfun(@diff, XLIM)));
cw = cumsum([0, w]) + 1;
tl = tiledlayout(figure(Units='inches', Position=[0 1 12 10]), nRows, sum(w), TileSpacing='tight', Padding='tight');
AX = gobjects(nRows, nCols);

for iRow = 1:nRows
    for iCol = 1:nCols
        ax = nexttile(tl, (iRow-1)*sum(w) + cw(iCol), [1, w(iCol)]); AX(iRow, iCol) = ax;
        if ETA(iRow, iCol) == ""
            ax.Visible = false;
            continue
        end
        etaTemp = eta.(ETA(iRow, iCol));
        if iCol == 1 && iRow == 1
            [~, tempOrder] = EphysUnit.plotETA(ax, etaTemp, xlim=XLIM{iCol}, clim=etaParams.clim, ...
                sortGroup=groupVar, sortWindow=[-2, 0], signWindow=[-0.5, 0], sortThreshold=0.5, negativeSortThreshold=0.25, hideColorbar=true);
        elseif iCol == nCols && iRow == 1
            EphysUnit.plotETA(ax, etaTemp, xlim=XLIM{iCol}, clim=etaParams.clim, order=tempOrder, hideColorbar=false);
            ax.Colorbar.Layout.Tile = 'east';
        else
            EphysUnit.plotETA(ax, etaTemp, xlim=XLIM{iCol}, clim=etaParams.clim, order=tempOrder, hideColorbar=true);
        end
        xline(ax, meta.window.(WINDOW(iRow, iCol)), 'k')
        yline(ax, cumsum(N(1:end-1)), '-', LineWidth=1.5, Color=hsl2rgb([0.4, 0.8, 0.2]));
        yline(ax, cumsum(N1(1:end-1)), '--', LineWidth=1.5, Color=hsl2rgb([0.4, 0.8, 0.2]));
        yticks(ax, cumsum(N));
        % applyCustomColormap(ax, [-1.5, 1.5], hlim=[0.375, 0, 0, -0.375], llim=[0.125, 0.5, 0.5, 0.25], hpwr=.5, lpwr=1, h0=0.33);
        applyCustomColormap(ax, [-1.5, 1.5], hlim=[0.375, 0, 0, -0.375], llim=[0.125, 1, 1, 0.25], hpwr=.5, lpwr=0.33, h0=0.33);
    
        title(ax, TITLE(iRow, iCol))
        xlabel(ax, sprintf("time to %s (s)", EVENT(iRow, iCol)))
    end
end

ylabel(AX, '')
yticks(AX(:, 2:end), [])
title(tl, 'Same order, normalized to peri-reach [-4, -2]s', FontWeight='bold')
fontsize(tl, 9, 'points')

print(tl.Parent, sprintf('%s\\PETH_Heatmap.png', path), '-dpng', '-r0')

% Mean PETH
TITLE = ["Reach", "Lick", "FirstLick(Reach)", "LastLick(Lick)", "LastLick(Reach)", "Release(Reach)"];
LEGEND = [ ...
    "correct", "correct", "correct", "correct", "correct", "correct"; ...
    "incorrect", "incorrect", "incorrect", "incorrect", "incorrect", "incorrect"; ...
    ];
nCols = length(ETA);

uniqueGroupVars = unique(groupVar);
nGroups = length(uniqueGroupVars);

tl = tiledlayout(figure(Units='inches', Position=[12 1 12 5]), nGroups, sum(w), TileSpacing='tight', Padding='tight');
AX = gobjects(nGroups, 3);

h = gobjects(2, 1);
colors = [0.2, 0.8, 0.2; 0.2, 0.2, 0.2];
for iGroup = 1:nGroups
    sel = groupVar == uniqueGroupVars(iGroup);
    for iCol = 1:nCols
        ax = nexttile(tl, (iGroup-1)*sum(w) + cw(iCol), [1, w(iCol)]); AX(iGroup, iCol) = ax;
        hold(ax, 'on')
        for iRow = 1:nRows
            if ETA(iRow, iCol) == ""
                continue
            end
            etaTemp = eta.(ETA(iRow, iCol));
            selT = etaTemp.t>XLIM{iCol}(1) & etaTemp.t<XLIM{iCol}(2);
            if iCol == 1 && iGroup == 1
                h(iRow) = plot(ax, etaTemp.t(selT), mean(etaTemp.X(sel, selT), 1, 'omitnan'), Color=colors(iRow, :), LineWidth=1.5, DisplayName=LEGEND(iRow, iCol));
            else
                plot(ax, etaTemp.t(selT), mean(etaTemp.X(sel, selT), 1, 'omitnan'), Color=colors(iRow, :), LineWidth=1.5, DisplayName=LEGEND(iRow, iCol));
            end
        end
        hold(ax, 'off')
        if iCol == 1
            ylabel(ax, [sprintf("%s", groupNames(uniqueGroupVars(iGroup))), sprintf("(n=%i)", nnz(sel))], Rotation=0)
        end
        xline(ax, 0, 'k:')
        yline(ax, 0, 'k:')
    end
    % xlim(AX(iGroup, :), [-2, 2])
    drawnow();
    yl = vertcat(AX(iGroup, :).YLim);
    ylim(AX(iGroup, :), [min(yl(:)), max(yl(:))]);
end

for iCol = 1:nCols
    title(AX(1, iCol), TITLE(iCol))
    xlabel(AX(end, iCol), sprintf('time to %s (s)', EVENT(iCol)))
end
lgd = legend(h, Orientation='horizontal');
lgd.Layout.Tile = 'north';
xticks(AX(1:end-1, :), [])
yticks(AX(:, 2:end), [])

fontsize(AX, 9, 'points')

ylabel(tl, "Spike rate (sp/s)", FontSize=9);
print(tl.Parent, sprintf('%s\\per_group_average.png', path), '-dpng', '-r0')

% clear tempOrder etaTemp iAx ax tl AX groupVar ETA WINDOW nPlots
% clear iGroup uniqueGrroupVars groupVars selJustOneDec selJustTwoDec tl AX iAx ax N iAx


%% Plot individial units, PETH as trace: lick, reach, release

grpNames = groupNames.values;
grpNames = grpNames(:)';
for grpName = grpNames
    if ~exist(sprintf("%s\\%s", path, grpName), 'dir')
        mkdir(sprintf("%s\\%s", path, grpName))
    end
end

% ETA = {eta.correctLick, eta.correctPress, eta.correctReleaseAfterRetract, eta.correctRelease; 
%     eta.incorrectLick, eta.incorrectPress, eta.correctReleaseBeforeRetract, eta.incorrectRelease};
% NAME = ["self-timed lick", "self-timed reach", "release (correct)", "release"];
% OUTCOME = ["correct", "correct", "after bar retract", "correct"; ...
%     "incorrect", "incorrect", "before bar retract", "incorrect"];
% XLIM = {[-2, 2], [-2, 2], [-1, 2], [-1, 2]};

ETA = {eta.correctPress, eta.correctLick, eta.correctPressFirstLick, eta.correctLickLastLickOff, eta.correctPressLastLickOff, eta.correctRelease; ...
        eta.incorrectPress, eta.incorrectLick, [], [], [], []};
NAME = ["Reach", "Lick", "FirstLick(Reach)", "LastLick(Lick)", "LastLick(Reach)", "Release(Reach)"];
ZEROLABEL = ["touch", "lick", "lick", "lick-off", "lick-off", "release"];
META = ["pressNorm", "lickNorm", "", "correctLickLastLickOffNorm", "", "correctReleaseNorm"];
WINDOW = ["press", "lick", "", "lastLickOff", "", "release"];
OUTCOME = ["cor"; "inc"];
XLIM = {[-2, 2], [-2, 2], [-1, 1], [-1, 1], [-1, 1], [-1, 1]};

fig = figure(Units='inches', OuterPosition=[1, 1, 12, 3]);
tl = tiledlayout(fig, 1, sum(cellfun(@diff, XLIM)), TileSpacing='compact', Padding='tight');
ax = gobjects(1, length(XLIM));
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

        h = gobjects(2, 1);
        if iEu == 0
            h(1) = plot(ax(i), ETA{1, i}.t(sel), mean(ETA{1, i}.X(:, sel), 1, 'omitnan'), Color=[0.2 0.6 0.2], DisplayName=sprintf("%s", OUTCOME(1)), LineWidth=1.5);
            if ~isempty(ETA{2, i})
                h(2) = plot(ax(i), ETA{2, i}.t(sel), mean(ETA{2, i}.X(:, sel), 1, 'omitnan'), Color=[0.2 0.2 0.2], DisplayName=sprintf("%s", OUTCOME(2)), LineWidth=1.5);
            else
                h = h(1);
            end
            yl(1) = 20;
            yl(2) = 80;
        else
            h(1) = plot(ax(i), ETA{1, i}.t(sel), ETA{1, i}.X(iEu, sel), Color=[0.2 0.6 0.2], DisplayName=sprintf("%s(n=%i)", OUTCOME(1), ETA{1, i}.N(iEu)), LineWidth=1.5);
            if ~isempty(ETA{2, i})
                h(2) = plot(ax(i), ETA{2, i}.t(sel), ETA{2, i}.X(iEu, sel), Color=[0.2 0.2 0.2], DisplayName=sprintf("%s(n=%i)", OUTCOME(2), ETA{2, i}.N(iEu)), LineWidth=1.5);
                yl(1) = min([yl(1), ETA{1, i}.X(iEu, sel), ETA{2, i}.X(iEu, sel)], [], 'omitnan');
                yl(2) = max([yl(2), ETA{1, i}.X(iEu, sel), ETA{2, i}.X(iEu, sel)], [], 'omitnan');
            else
                yl(1) = min([yl(1), ETA{1, i}.X(iEu, sel)], [], 'omitnan');
                yl(2) = max([yl(2), ETA{1, i}.X(iEu, sel)], [], 'omitnan');       
                h = h(1);         
            end
            if META(i) ~= ""
                window = meta.window.(WINDOW(i));
                if useLocalNorm
                    metaTemp = meta.localNorm.(META(i));
                else
                    metaTemp = meta.(META(i));
                end
                if metaTemp(iEu) < 0
                    color = 'blue';
                else
                    color = 'red';
                end
                patch(ax(i), window([1 2 2 1]), yl([1 1 2 2]), color, EdgeColor=color, FaceAlpha=0.33)
            end
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
        grpName = groupNames(groupVar(iEu));
        title(tl, sprintf("%s (%s)", eu(iEu).getName(), grpName), Interpreter='none', FontWeight='bold')
        print(fig, sprintf('%s\\%s\\%s.png', path, grpName, eu(iEu).getName()), '-dpng', '-r0')
    end
end

clear fig tl ax ETA NAME OUTCOME iEu i lgd yl path h grpName