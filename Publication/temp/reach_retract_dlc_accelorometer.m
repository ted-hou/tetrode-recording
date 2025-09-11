%% Load ephys units
eu = EphysUnit.load('\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\Units\TwoColor_SNr_SCRetro\ReverseInjection\SingleUnit_NonDuplicate_NonDrift_SNr', waveforms=false, spikecounts=false, spikerates=false);

%%
sessions = [...
    % "daisy27_20250624", ... SNr, traditional
    % "daisy27_20250626", ... SNr, Accel Lick
    ... "daisy27_20250707", ... SC, Accel?
    ... "daisy27_20250715", ... SC, Accel?
    "daisy27_20250717", ... SNr, Accel Lick and Reach
    "daisy27_20250721", ... SNr, Accel Lick and Reach
    "daisy27_20250724", ... SNr, Accel Lick and Reach
    "daisy28_20250701", ... SNr, Accel Lick and Reach
    "daisy28_20250702", ... SNr, Accel Lick and Reach
    ... "daisy28_20250714", ... SC
    "daisy28_20250716", ... SNr, Accel
    "daisy28_20250718", ... SNr, Accel
    "daisy28_20250723", ... SNr, Accel
    "daisy28_20250725", ... SNr, Accel
    ... "daisy28_20250728", ... SNr, Accel (Camera 2 bad for first 14 min)
    "daisy28_20250729" ... SNr, Accel
];

selUnits = ismember(string({eu.ExpName}), sessions);
eu = eu(selUnits);

%% Load deeplabcut
dlcPath = "\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\DeepLabCut\Results";
exp = CompleteExperiment3(eu, cameras='lr', deeplabcutPath=dlcPath);
exp.alignTimestamps(refEventNameArduino={'REWARD_ON', 'OPTO1_ON', 'OPTO2_ON', 'LEVER_PRESSED'}, refEventNameEphys={'RewardOn', 'LaserModBlueOn', 'LaserModRedOn', 'Press'}, trialDurationTolerance=1);

%% Both were recording in left hemisphere, and in the DLC labeling, HandContra/HandIpsi refers to camera-side/far-side, we need to fix this for both the left and right cameras. Left/right cameras are placed to the left/right side of mouse.
% lcam: HandContra -> handContra, HandIpsi -> handIpsi
% rcam: HandContra -> handIpsi, HandIpsi -> handContra

for iExp = 1:length(exp)
    exp(iExp).vtdL = renamevars(exp(iExp).vtdL, ...
        [ ...
            "HandContra_X", "HandContra_Y", "HandContra_Likelihood", ...
            "HandIpsi_X", "HandIpsi_Y", "HandIpsi_Likelihood", ...
            "Tongue_X", "Tongue_Y", "Tongue_Likelihood", ...
            "Jaw_X", "Jaw_Y", "Jaw_Likelihood", ...
        ], ...
        [ ...
            "handContra_X", "handContra_Y", "handContra_Likelihood", ...
            "handIpsi_X", "handIpsi_Y", "handIpsi_Likelihood", ...
            "tongue_X", "tongue_Y", "tongue_Likelihood", ...
            "jaw_X", "jaw_Y", "jaw_Likelihood", ...
        ]);
    exp(iExp).vtdR = renamevars(exp(iExp).vtdR, ...
        [ ...
            "HandContra_X", "HandContra_Y", "HandContra_Likelihood", ...
            "HandIpsi_X", "HandIpsi_Y", "HandIpsi_Likelihood", ...
            "Tongue_X", "Tongue_Y", "Tongue_Likelihood", ...
            "Jaw_X", "Jaw_Y", "Jaw_Likelihood", ...
        ], ...
        [ ...
            "handIpsi_X", "handIpsi_Y", "handIpsi_Likelihood", ...
            "handContra_X", "handContra_Y", "handContra_Likelihood", ...
            "tongue_X", "tongue_Y", "tongue_Likelihood", ...
            "jaw_X", "jaw_Y", "jaw_Likelihood", ...
        ]);
end

%% Make lick and reach trials, but first cleanup the bad accel artifacts
pulseWidthThreshold = 1.5e-3;

fig = figure(Units='inches', Position=[1, 1, 7, 7]);
tl = tiledlayout(fig, length(exp), 1, TileSpacing='tight', Padding='tight');
for iExp = 1:length(exp)
    ax = nexttile(tl); hold(ax, 'on')
    histogram(ax, 1e3*(exp(iExp).eu(1).EventTimes.PressOff - exp(iExp).eu(1).EventTimes.PressOn), [0:0.1:10, Inf], EdgeColor='red', FaceColor='red', FaceAlpha=0.2, Normalization='pdf', DisplayName='reach')
    histogram(ax, 1e3*(exp(iExp).eu(1).EventTimes.LickOff - exp(iExp).eu(1).EventTimes.LickOn), [0:0.1:10, Inf], EdgeColor='blue', FaceColor='blue', FaceAlpha=0.2, Normalization='pdf', DisplayName='lick')
    xline(ax, 1e3*pulseWidthThreshold, 'k:', LineWidth=1.5, DisplayName='threshold')
end
lgd = legend(ax);
lgd.Layout.Tile = 'north';
xlabel(tl, 'pulse width (ms)')
ylabel(tl, 'prob')
for iEu = 1:length(eu)
    lick = eu(iEu).EventTimes.Lick;
    press = eu(iEu).EventTimes.Press;
    selLick = eu(iEu).EventTimes.LickOff - eu(iEu).EventTimes.LickOn > 1.5e-3;
    selPress = eu(iEu).EventTimes.PressOff - eu(iEu).EventTimes.PressOn > 1.5e-3;
    eu(iEu).Trials.Press = Trial(eu(iEu).EventTimes.TIMEOUT_START, eu(iEu).EventTimes.Press(selPress), stopMode='first', exclude=eu(iEu).EventTimes.Lick(selLick));
    eu(iEu).Trials.Lick = Trial(eu(iEu).EventTimes.TIMEOUT_START, eu(iEu).EventTimes.Lick(selLick), stopMode='first', exclude=eu(iEu).EventTimes.Press(selPress));
end
clear iEu fig tl iExp ax lick press selLick selPress

fprintf("Using a threshold of %gms, we have the following number of trials for each session:\n(A press trial is determined as the first PressOn after TIMEOUT_START, without any LickOn in between)\n(A lick trial is determined as the first LickOn after TIMEOUT_START, without any PressOn in between)\n", 1e3*pulseWidthThreshold)
disp(table(arrayfun(@(exp) length(exp.eu(1).Trials.Press), exp(:)), arrayfun(@(exp) length(exp.eu(1).Trials.Lick), exp(:)), VariableNames=["Press", "Lick"]))

% %%
% close all
% for iExp = 1:length(exp)
%     ax = axes(figure);
%     hold(ax, 'on')
%     [X, Y, L, t] = exp(iExp).getTrajectoryByTrial('r', 'handContra', trialType='press', window=[-1, 0.5], likelihoodThreshold=0.4);
%     X = X - X(:, 1); Y = Y - Y(:, 1);
%     plot(ax, median(X, 1, 'omitnan'), median(Y, 1, 'omitnan'), 'r-', Marker='x', MarkerSize=25, LineWidth=1.5, DisplayName='reach')
%     [X, Y, L, t] = exp(iExp).getTrajectoryByTrial('r', 'handContra', trialType='lick', window=[-1, 0.5], likelihoodThreshold=0.4);
%     X = X - X(:, 1); Y = Y - Y(:, 1);
%     plot(ax, median(X, 1, 'omitnan'), median(Y, 1, 'omitnan'), 'b-', Marker='x', MarkerSize=25, LineWidth=1.5, DisplayName='lick')
%     xline(ax, 0, 'k--')
%     yline(ax, 0, 'k--')
%     title(ax, exp(iExp).name, Interpreter='none')
%     axis(ax, 'image')
%     axis(ax, 'equal')
% end

% Calculate arm traversal distance during reach vs. lick trials.
clear traj
traj(length(exp)) = struct(window=-[], minLikelihood=[], press=[], lick=[]);
badExpIndices = [];
for iExp = 1:length(exp)
    try
        traj(iExp).window = [-1, 0.5];
        traj(iExp).minLikelihood = 0.4;
        [traj(iExp).press.handContra.X, traj(iExp).press.handContra.Y, traj(iExp).press.handContra.L, traj(iExp).press.handContra.t] = exp(iExp).getTrajectoryByTrial('r', 'handContra', trialType='press', window=traj(iExp).window, likelihoodThreshold=traj(iExp).minLikelihood);
        [traj(iExp).lick.handContra.X, traj(iExp).lick.handContra.Y, traj(iExp).lick.handContra.L, traj(iExp).lick.handContra.t] = exp(iExp).getTrajectoryByTrial('r', 'handContra', trialType='lick', window=traj(iExp).window, likelihoodThreshold=traj(iExp).minLikelihood);
        [traj(iExp).press.handIpsi.X, traj(iExp).press.handIpsi.Y, traj(iExp).press.handIpsi.L, traj(iExp).press.handIpsi.t] = exp(iExp).getTrajectoryByTrial('l', 'handIpsi', trialType='press', window=traj(iExp).window, likelihoodThreshold=traj(iExp).minLikelihood);
        [traj(iExp).lick.handIpsi.X, traj(iExp).lick.handIpsi.Y, traj(iExp).lick.handIpsi.L, traj(iExp).lick.handIpsi.t] = exp(iExp).getTrajectoryByTrial('l', 'handIpsi', trialType='lick', window=traj(iExp).window, likelihoodThreshold=traj(iExp).minLikelihood);
        [traj(iExp).press.tongue.X, traj(iExp).press.tongue.Y, traj(iExp).press.tongue.L, traj(iExp).press.tongue.t] = exp(iExp).getTrajectoryByTrial('r', 'tongue', trialType='press', window=traj(iExp).window, likelihoodThreshold=traj(iExp).minLikelihood);
        [traj(iExp).lick.tongue.X, traj(iExp).lick.tongue.Y, traj(iExp).lick.tongue.L, traj(iExp).lick.tongue.t] = exp(iExp).getTrajectoryByTrial('r', 'tongue', trialType='lick', window=traj(iExp).window, likelihoodThreshold=traj(iExp).minLikelihood);
    
        % sum the traversal distance 
        for trialType = ["press", "lick"]
            traj(iExp).(trialType).handContra.traversal = sum(sqrt(diff(traj(iExp).(trialType).handContra.X, 1, 2).^2 + diff(traj(iExp).(trialType).handContra.Y, 1, 2).^2), 2, 'omitnan');
            traj(iExp).(trialType).handIpsi.traversal = sum(sqrt(diff(traj(iExp).(trialType).handIpsi.X, 1, 2).^2 + diff(traj(iExp).(trialType).handIpsi.Y, 1, 2).^2), 2, 'omitnan');
        end
    catch
        badExpIndices = [badExpIndices, iExp];
        warning('Could not calculate trajectory for session %i', iExp)
    end
end
badExpIndices = unique(badExpIndices);

isGoodUnit = true(size(eu));
for iExp = badExpIndices
    isGoodUnit(ismember(eu, exp(iExp).eu)) = false;
end

fprintf('%i/%i sessions (%i/%i units) failed to generate trajectories and will not be included.\n', length(badExpIndices), length(exp), nnz(~isGoodUnit), length(eu));

% %% Plot them trajectories
% close all
% fig = figure(Units='inches', Position=[1, 1, 7, 7]);
% tl = tiledlayout(fig, length(exp), 1, TileSpacing='tight');
% for iExp = 1:length(exp)
%     try
%         ax = nexttile(tl); hold(ax, 'on');
%         colors = struct(press='red', lick='blue');
%         for trialType = ["press", "lick"]
%             histogram(ax, traj(iExp).(trialType).handContra.traversal, 0:1:100, DisplayName=trialType, FaceColor='none', EdgeColor=colors.(trialType), Normalization='pdf');
%         end
%         xline(ax, quantile(traj(iExp).press.handContra.traversal, [0.05, 0.25, 0.5]), colors.press, LineStyle=':')
%     end
% end
% 
% clear fig tl iExp trialType ax colors
% 
% %% Pick 1 random reach trials with arm traversal > 50prct of reach
% close all
% for iExp = 1:length(exp)
%     try
%         iTrial = find(isin(traj(iExp).press.handContra.traversal, quantile(traj(iExp).press.handContra.traversal, [0.5, 1]), true));
%         iTrial = iTrial(randi([1, length(iTrial)]));
%         [clip, t] = exp(iExp).getVideoClip(exp(iExp).eu(1).Trials.Press(iTrial).Stop, side='r', numFramesBefore=30, numFramesAfter=30, bodyParts={'handIpsi', 'handContra', 'jaw', 'tongue'}, ...
%             minLikelihood=traj(iExp).minLikelihood);
%         implay(clip, 30)
%     catch
%         warning('could not process exp %i: %s', iExp, exp(iExp).name)
%     end
% end
% 
% %% Pick 1 random lick trials with arm traversal < 25prct of reach
% close all
% for iExp = 1:length(exp)
%     try
%         iTrial = find(isin(traj(iExp).lick.handContra.traversal, quantile(traj(iExp).press.handContra.traversal, [0, 0.25]), true));
%         iTrial = iTrial(randi([1, length(iTrial)]));
%         [clip, t] = exp(iExp).getVideoClip(exp(iExp).eu(1).Trials.Lick(iTrial).Stop, side='r', numFramesBefore=30, numFramesAfter=30, bodyParts={'handIpsi', 'handContra', 'jaw', 'tongue'}, ...
%             minLikelihood=traj(iExp).minLikelihood);
%         implay(clip, 30)
%     catch
%         warning('could not process exp %i: %s', iExp, exp(iExp).name)
%     end
% end

minNumTrials = 4;
% Restrict lick trials and reach trials by arm traversal
for iEu = 1:length(eu)
    iExp = find(arrayfun(@(exp) ismember(eu(iEu), exp.eu), exp));
    if ismember(iExp, badExpIndices)
        continue
    end

    clear selTrials
    selTrials.press = isin(traj(iExp).press.handContra.traversal, quantile(traj(iExp).press.handContra.traversal, [0.5, 1]), true) | isin(traj(iExp).press.handIpsi.traversal, quantile(traj(iExp).press.handIpsi.traversal, [0.5, 1]), true);
    selTrials.lick = isin(traj(iExp).lick.handContra.traversal, quantile(traj(iExp).press.handContra.traversal, [0, 0.5]), true) & isin(traj(iExp).lick.handIpsi.traversal, quantile(traj(iExp).press.handIpsi.traversal, [0, 0.5]), true);
    
    if eu(iEu) == exp(iExp).eu(1)
        fprintf('Exp%i %s: kept %i/%i reach trials, kept %i/%i lick trials.\n', iExp, exp(iExp).name, nnz(selTrials.press), length(selTrials.press), nnz(selTrials.lick), length(selTrials.lick))
    end

    if nnz(selTrials.press) < minNumTrials || nnz(selTrials.lick) < minNumTrials
        isGoodUnit(iEu) = false;
        badExpIndices = [badExpIndices, iExp];
    end

    eu(iEu).Trials.Press = eu(iEu).Trials.Press(selTrials.press);
    eu(iEu).Trials.Lick = eu(iEu).Trials.Lick(selTrials.lick);
end
badExpIndices = unique(badExpIndices);

%% Now that we've cleanup'd press/lick trials (no confounding movement, no false-positives), let's get the following
%%% Bar retract, Rewarded/unrewarded, First/last lick
[ac.ac, ac.euIndicesFirstInSession, ac.expIndices, ac.uniqueExpNames] = eu.loadArduinoConnection();
eu.alignTimestamps(["LEVER_RETRACT_START", "TUBE_RETRACT_START"], ac={ac.ac, ac.euIndicesFirstInSession, ac.expIndices, ac.uniqueExpNames}, acRefEventName="REWARD_ON", euRefEventName="RewardOn");


for iExp = 1:length(exp)
    if ismember(iExp, badExpIndices)
        continue
    end

    % Rewarded/unrewarded trials
    pressTrials = exp(iExp).eu(1).Trials.Press;
    lickTrials = exp(iExp).eu(1).Trials.Lick;
    rewardTimes = exp(iExp).eu(1).EventTimes.RewardOn;
    pressWaitfortouchToReward = Trial([pressTrials.Start, Inf], rewardTimes, 'first', exclude=[lickTrials.Start]);
    lickWaitfortouchToReward = Trial([lickTrials.Start, Inf], rewardTimes, 'first', exclude=[pressTrials.Start]);
    isRewardedPress = pressWaitfortouchToReward.inTrial([pressTrials.Stop], [-0.5, 0.5], windowMode='stop');
    isRewardedLick = lickWaitfortouchToReward.inTrial([lickTrials.Stop], [-0.5, 0.5], windowMode='stop');

    % CorrectPressRetract (retract=bar retract, release=arm release)
    leverRetractTimes = exp(iExp).eu(1).EventTimes.LEVER_RETRACT_START;
    correctPressToLeverRetract = Trial([pressTrials(isRewardedPress).Stop, Inf], leverRetractTimes, 'first', exclude=[[pressTrials.Start], [lickTrials.Start]]);

    % CorrectPressToFirstLick
    lickTimes = exp(iExp).eu(1).EventTimes.LickOn;
    correctPressToFirstLick = Trial([pressTrials(isRewardedPress).Stop, Inf], lickTimes, 'first', exclude=[[pressTrials.Start], [lickTrials.Start]]);

    % CorrectPressToLastLick
    TODO!
        % % For each self-timed Lick, find the next cue (of any kind)
        % firstLickToCueValid = Trial([eu(iEu).Trials.LickValid.Stop, Inf], eu(iEu).EventTimes.Cue, 'first');
        % firstLickToCueIncorrect = Trial([eu(iEu).Trials.LickIncorrect.Stop, Inf], eu(iEu).EventTimes.Cue, 'first');
        % firstLickToCueCorrect = Trial([eu(iEu).Trials.LickCorrect.Stop, Inf], eu(iEu).EventTimes.Cue, 'first');
        % % Find the last lickOff between firstLick and nextCue
        % [~, lickOff, trialIndices] = firstLickToCueValid.inTrial(eu(iEu).EventTimes.LickOff);
        % [~, lastLickOff, ~] = unique(trialIndices, 'last');
        % lastLickOff = lickOff(lastLickOff);
        % eu(iEu).Trials.CueToLastLickOff = Trial([eu(iEu).Trials.LickValid.Start, Inf], lastLickOff(:)', 'last');
    
    % CorrectLickToLastLick

    for iEu = 1:length(exp(iExp).eu)
        exp(iExp).eu(iEu).Trials.CorrectPress = pressTrials(isRewardedPress);
        exp(iExp).eu(iEu).Trials.IncorrectPress = pressTrials(~isRewardedPress);
        exp(iExp).eu(iEu).Trials.CorrectLick = lickTrials(isRewardedLick);
        exp(iExp).eu(iEu).Trials.IncorrectLick = lickTrials(~isRewardedLick);
        exp(iExp).eu(iEu).Trials.CorrectPressToLeverRetract = correctPressToLeverRetract;
        exp(iExp).eu(iEu).Trials.CorrectPressToFirstLick = correctPressToFirstLick;
    end

    fprintf('iExp=%i, %s, reach (%i rewarded, %i unrewarded), lick (%i rewarded, %i unrewarded).\n', iExp, exp(iExp).name, nnz(isRewardedPress), nnz(~isRewardedPress), nnz(isRewardedLick), nnz(~isRewardedLick));
end
clear pressTrials lickTrials rewardTimes pressWaitfortouchToReward lickWaitfortouchToReward isRewardedPress isRewardedLick iExp iEu

%% Make/Plot ETA Heatmap

normWindow=[-4, -2];
xlDisp = [-2, 2];

clear eta
selUnits = isGoodUnit;
eta.pressNorm = eu.getETA('count', 'press', [-4, 4], selUnits=selUnits, resolution=0.1, alignTo='stop', includeInvalid=true, normalize=normWindow, minTrialDuration=0, maxTrialDuration=Inf);
eta.lickNorm = eu.getETA('count', 'lick', [-4, 4], selUnits=selUnits, resolution=0.1, alignTo='stop', includeInvalid=true, normalize=eta.pressNorm.stats, minTrialDuration=0, maxTrialDuration=Inf);
eta.correctPressNorm = eu.getETA('count', 'CorrectPress', [-4, 4], selUnits=selUnits, resolution=0.1, alignTo='stop', includeInvalid=true, normalize=eta.pressNorm.stats, minTrialDuration=0, maxTrialDuration=Inf);
eta.correctLickNorm = eu.getETA('count', 'CorrectLick', [-4, 4], selUnits=selUnits, resolution=0.1, alignTo='stop', includeInvalid=true, normalize=eta.pressNorm.stats, minTrialDuration=0, maxTrialDuration=Inf);
eta.incorrectPressNorm = eu.getETA('count', 'IncorrectPress', [-4, 4], selUnits=selUnits, resolution=0.1, alignTo='stop', includeInvalid=true, normalize=eta.pressNorm.stats, minTrialDuration=0, maxTrialDuration=Inf);
eta.incorrectLickNorm = eu.getETA('count', 'IncorrectLick', [-4, 4], selUnits=selUnits, resolution=0.1, alignTo='stop', includeInvalid=true, normalize=eta.pressNorm.stats, minTrialDuration=0, maxTrialDuration=Inf);

%% For sorting, make templates to dot-product with
clear template
template(length(ETASORT)) = struct(t=[], x=[]);
for iETA = 1:length(ETASORT)
    template(iETA).t = etaCombined.t;
    template(iETA).x = zeros(1, length(etaCombined.t));
    template(iETA).x(1, isin(etaCombined.t, SORTWINDOW{iETA}) & etaCombined.epoch==iETA) = 1;
end

%%
close all
% EphysUnit.plotDoubleETA(eta.lick, eta.press, selUnits, 'lick', 'reach', clim=[-1, 1], xlim=xlDisp, sortWindow=[-3, 1], sortThreshold=0.25)

tl = tiledlayout(figure(Units='normalized', OuterPosition=[0.5, 0, 0.5, 1]), 1, 4);
ax = nexttile(tl);
[~, orderA] = EphysUnit.plotETA(ax, eta.lick, isGoodUnit, event='lick', clim=[-1.5, 1.5], xlim=xlDisp, sortWindow=[-3, 1], signWindow=[-.1, 0], sortThreshold=0.25, hideColorbar=true);
xline(ax, 0)
title(ax, 'Lick (sort A)')

ax = nexttile(tl);
[~, ~] = EphysUnit.plotETA(ax, eta.press, isGoodUnit, event='reach', order=orderA, clim=[-1.5, 1.5], xlim=xlDisp, sortWindow=[-3, 1], signWindow=[-.1, 0], sortThreshold=0.25, hideColorbar=true);
xline(ax, 0)
title(ax, 'Reach (sort A)')

ax = nexttile(tl);
[~, orderB] = EphysUnit.plotETA(ax, eta.press, isGoodUnit, event='reach', clim=[-1.5, 1.5], xlim=xlDisp, sortWindow=[-3, 1], signWindow=[-0.5, 0], sortThreshold=0.25, hideColorbar=false);
ax.Colorbar.Layout.Tile = 'east';
xline(ax, 0)
title(ax, 'Reach (sort B)')

ax = nexttile(tl);
[~, ~] = EphysUnit.plotETA(ax, eta.lick, isGoodUnit, event='lick', order=orderB, clim=[-1.5, 1.5], xlim=xlDisp, sortWindow=[-3, 1], signWindow=[-.1, 0], sortThreshold=0.25, hideColorbar=true);
xline(ax, 0)
title(ax, 'Lick (sort B)')

clear normWindow xlDisp selUnits tl ax orderA orderB