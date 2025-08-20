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

clear artifactParams;
artifactParams(1) = struct(event='LickOn', length=10, lengthUnit='ms', direction='both');
artifactParams(2) = struct(event='LickOff', length=10, lengthUnit='ms', direction='both');
artifactParams(3) = struct(event='PressOn', length=10, lengthUnit='ms', direction='both');
artifactParams(4) = struct(event='PressOff', length=10, lengthUnit='ms', direction='both');


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

eta.pressNorm = eu.getETA('count', 'press', [-4, 4], alignTo='stop', minTrialDuration=2, normalize=[-4, -2], resolution=0.1, artifacts=artifactParams);
eta.lickNorm = eu.getETA('count', 'lick', [-4, 4], alignTo='stop', minTrialDuration=2, normalize=eta.pressNorm.stats, resolution=0.1, artifacts=artifactParams);
eta.correctReleaseNorm = eu.getETA('count', 'CueToLeverReleaseCorrect', [-4, 4], alignTo='stop', minTrialDuration=4, normalize=eta.pressNorm.stats, resolution=0.1, artifacts=artifactParams);

eta.correctPressNorm = eu.getETA('count', 'press', [-4, 4], alignTo='stop', minTrialDuration=4, normalize=eta.pressNorm.stats, resolution=0.1, artifacts=artifactParams);
eta.correctLickNorm = eu.getETA('count', 'lick', [-4, 4], alignTo='stop', minTrialDuration=4, normalize=eta.pressNorm.stats, resolution=0.1, artifacts=artifactParams);
eta.correctReleaseBeforeRetractNorm = eu.getETA('count', 'CueToLeverReleaseBeforeRetractCorrect', [-4, 4], alignTo='stop', minTrialDuration=4, normalize=eta.pressNorm.stats, resolution=0.1, artifacts=artifactParams);
eta.correctReleaseAfterRetractNorm = eu.getETA('count', 'RetractToLeverReleaseAfterRetractCorrect', [-4, 4], alignTo='stop', minTrialDuration=4, normalize=eta.pressNorm.stats, resolution=0.1, artifacts=artifactParams);

eta.incorrectPressNorm = eu.getETA('count', 'press', [-4, 4], alignTo='stop', minTrialDuration=2, maxTrialDuration=4, normalize=eta.pressNorm.stats, resolution=0.1, artifacts=artifactParams);
eta.incorrectLickNorm = eu.getETA('count', 'lick', [-4, 4], alignTo='stop', minTrialDuration=2, maxTrialDuration=4, normalize=eta.pressNorm.stats, resolution=0.1, artifacts=artifactParams);


%% Plot all units, PETH as heatmap: lick, reach, release
tl = tiledlayout(figure(Units='inches', Position=[1 1 10 5]), 1, 3);
clear etaOrder ax

% PETH, lick
ax = nexttile(tl);
[~, etaOrder.lick] = EphysUnit.plotETA(ax, eta.lickNorm);

ax = nexttile(tl);
[~, etaOrder.lick] = EphysUnit.plotETA(ax, eta.correctLick);


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