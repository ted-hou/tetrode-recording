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

[ac.ac, ac.euIndicesFirstInSession, ac.expIndices, ac.uniqueExpNames] = eu.loadArduinoConnection();
eu.alignTimestamps(["LEVER_RETRACT_START", "TUBE_RETRACT_START"], ac={ac.ac, ac.euIndicesFirstInSession, ac.expIndices, ac.uniqueExpNames}, acRefEventName="REWARD_ON", euRefEventName="RewardOn");

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
pulseWidthThreshold = struct(press = 1.5e-3, lick=0);
minTrialLength = 2;

fig = figure(Units='inches', Position=[1, 1, 7, 7]);
tl = tiledlayout(fig, length(exp), 1, TileSpacing='tight', Padding='tight');
for iExp = 1:length(exp)
    ax = nexttile(tl); hold(ax, 'on')
    histogram(ax, 1e3*(exp(iExp).eu(1).EventTimes.PressOff - exp(iExp).eu(1).EventTimes.PressOn), [0:0.1:10, Inf], EdgeColor='red', FaceColor='red', FaceAlpha=0.2, Normalization='pdf', DisplayName='reach')
    histogram(ax, 1e3*(exp(iExp).eu(1).EventTimes.LickOff - exp(iExp).eu(1).EventTimes.LickOn), [0:0.1:10, Inf], EdgeColor='blue', FaceColor='blue', FaceAlpha=0.2, Normalization='pdf', DisplayName='lick')
    xline(ax, 1e3*pulseWidthThreshold.press, 'k:', LineWidth=1.5, DisplayName='threshold')
end
lgd = legend(ax);
lgd.Layout.Tile = 'north';
xlabel(tl, 'pulse width (ms)')
ylabel(tl, 'prob')
for iEu = 1:length(eu)
    selPress = eu(iEu).EventTimes.PressOff - eu(iEu).EventTimes.PressOn > pulseWidthThreshold.press;
    selLick = eu(iEu).EventTimes.LickOff - eu(iEu).EventTimes.LickOn > pulseWidthThreshold.lick;
    eu(iEu).Trials.Press = Trial(eu(iEu).EventTimes.TIMEOUT_START, eu(iEu).EventTimes.Press(selPress), stopMode='first', exclude=eu(iEu).EventTimes.Lick(selLick));
    eu(iEu).Trials.Lick = Trial(eu(iEu).EventTimes.TIMEOUT_START, eu(iEu).EventTimes.Lick(selLick), stopMode='first', exclude=eu(iEu).EventTimes.Press(selPress));
    eu(iEu).Trials.Press = eu(iEu).Trials.Press(eu(iEu).Trials.Press.duration() >= minTrialLength);
    eu(iEu).Trials.Lick = eu(iEu).Trials.Lick(eu(iEu).Trials.Lick.duration() >= minTrialLength);
    eu(iEu).EventTimes.ValidPress = eu(iEu).EventTimes.Press(selPress);
    eu(iEu).EventTimes.ValidLick = eu(iEu).EventTimes.Lick(selLick);
end
clear fig tl iExp lgd iEu ax selLick selPress

fprintf("Using a threshold of %gms for press and %gms for lick, we have the following number of trials for each session:\n(A press trial is determined as the first PressOn after TIMEOUT_START, without any LickOn in between)\n(A lick trial is determined as the first LickOn after TIMEOUT_START, without any PressOn in between)\n", 1e3*pulseWidthThreshold.press, 1e3*pulseWidthThreshold.lick)
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
clear iExp trialType

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

minNumTrials = 8;
% Restrict lick trials and reach trials by arm traversal
for iEu = 1:length(eu)
    iExp = find(arrayfun(@(exp) ismember(eu(iEu), exp.eu), exp));
    if ismember(iExp, badExpIndices)
        continue
    end

    clear selTrials
    selTrials.press = isin(traj(iExp).press.handContra.traversal, quantile(traj(iExp).press.handContra.traversal, [0.5, 1]), true) | isin(traj(iExp).press.handIpsi.traversal, quantile(traj(iExp).press.handIpsi.traversal, [0.5, 1]), true);
    selTrials.lick = isin(traj(iExp).lick.handContra.traversal, quantile(traj(iExp).press.handContra.traversal, [0, 0.3]), true) & isin(traj(iExp).lick.handIpsi.traversal, quantile(traj(iExp).press.handIpsi.traversal, [0, 0.3]), true);
    
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

clear minNumTrials iEu ixp selTrials

%% Now that we've cleanup'd press/lick trials (no confounding movement, no false-positives), let's get the following trial types:
% CorrectPress
% IncorrectPress
% CorrectLick
% IncorrectLick
% CorrectPressToLeverRetract
% CorrectPressToFirstLick
% CorrectPressToLastLick
% CorrectLickToLastLick

clear nTrials
nTrials(length(exp)) = struct(iExp=[], name="", nUnits=[], correctPress=[], incorrectPress=[], correctLick=[], incorrectLick=[], correctPressToLeverRetract=[], correctPressToFirstLick=[], correctPressToLastLick=[], correctLickToLastLick=[]);

for iExp = 1:length(exp)
    if ismember(iExp, badExpIndices)
        continue
    end

    % Rewarded/unrewarded trials
    % wait = WAITFORTOUCH (guaranteed reward) % timeout = TIMEOUT_START (either rewareded or not)
    timeoutToPress = exp(iExp).eu(1).Trials.Press;
    timeoutToLick = exp(iExp).eu(1).Trials.Lick;

    assert(issorted([timeoutToPress.Start]))
    assert(issorted([timeoutToLick.Start]))
    [~, ~, selCorrectPress] = timeoutToPress.inTrial(exp(iExp).eu(1).EventTimes.WAITFORTOUCH);
    [~, ~, selCorrectLick] = timeoutToLick.inTrial(exp(iExp).eu(1).EventTimes.WAITFORTOUCH);
    isCorrectPress = ismember(1:length(timeoutToPress), selCorrectPress);
    isCorrectLick = ismember(1:length(timeoutToLick), selCorrectLick);

    correctTimeoutToPress = timeoutToPress(isCorrectPress);
    incorrectTimeoutToPress = timeoutToPress(~isCorrectPress);
    correctTimeoutToLick = timeoutToLick(isCorrectLick);
    incorrectTimeoutToLick = timeoutToLick(~isCorrectLick);

    % CorrectPressRetract (retract=bar retract, release=arm release)
    correctPressToLeverRetract = Trial([correctTimeoutToPress.Stop, Inf], exp(iExp).eu(1).EventTimes.LEVER_RETRACT_START, 'first', exclude=exp(iExp).eu(1).EventTimes.TIMEOUT_START);

    % CorrectPressToFirstLick
    correctPressToFirstLick = Trial([correctTimeoutToPress.Stop, Inf], exp(iExp).eu(1).EventTimes.ValidLick, 'first', exclude=exp(iExp).eu(1).EventTimes.TIMEOUT_START);

    % CorrectPressToLastLick
    % For each correctPress, find the next timeout (of any kind)
    correctPressToTimeout = Trial([correctTimeoutToPress.Stop, Inf], exp(iExp).eu(1).EventTimes.TIMEOUT_START, 'first');
    % Find the last lick between correctPress and timeoutstart
    [~, lick, trialIndices] = correctPressToTimeout.inTrial(exp(iExp).eu(1).EventTimes.ValidLick);
    [~, lastLick, ~] = unique(trialIndices, 'last');
    lastLick = lick(lastLick);
    correctPressToLastLick = Trial([correctTimeoutToPress.Stop, Inf], lastLick(:)', 'last');

    % CorrectLickToLastLick
    % For each correctLick, find the next timeout (of any kind)
    correctLickToTimeout = Trial([correctTimeoutToLick.Stop, Inf], exp(iExp).eu(1).EventTimes.TIMEOUT_START, 'first');
    % Find the last lick between correctLick and timeoutstart
    [~, lick, trialIndices] = correctLickToTimeout.inTrial(exp(iExp).eu(1).EventTimes.ValidLick);
    [~, lastLick, ~] = unique(trialIndices, 'last');
    lastLick = lick(lastLick);
    correctLickToLastLick = Trial([correctTimeoutToLick.Stop, Inf], lastLick(:)', 'last');

    for iEu = 1:length(exp(iExp).eu)
        exp(iExp).eu(iEu).Trials.CorrectPress = correctTimeoutToPress;
        exp(iExp).eu(iEu).Trials.IncorrectPress = incorrectTimeoutToPress;
        exp(iExp).eu(iEu).Trials.CorrectLick = correctTimeoutToLick;
        exp(iExp).eu(iEu).Trials.IncorrectLick = incorrectTimeoutToLick;
        exp(iExp).eu(iEu).Trials.CorrectPressToLeverRetract = correctPressToLeverRetract;
        exp(iExp).eu(iEu).Trials.CorrectPressToFirstLick = correctPressToFirstLick;
        exp(iExp).eu(iEu).Trials.CorrectPressToLastLick = correctPressToLastLick;
        exp(iExp).eu(iEu).Trials.CorrectLickToLastLick = correctLickToLastLick;
    end

    nTrials(iExp).iExp = iExp;
    nTrials(iExp).nUnits = length(exp(iExp).eu);
    nTrials(iExp).name = string(exp(iExp).name);
    nTrials(iExp).correctPress = length(correctTimeoutToPress);
    nTrials(iExp).incorrectPress = length(incorrectTimeoutToPress);
    nTrials(iExp).correctLick = length(correctTimeoutToLick);
    nTrials(iExp).incorrectLick = length(incorrectTimeoutToLick);
    nTrials(iExp).correctPressToLeverRetract = length(correctPressToLeverRetract);
    nTrials(iExp).correctPressToFirstLick = length(correctPressToFirstLick);
    nTrials(iExp).correctPressToLastLick = length(correctPressToLastLick);
    nTrials(iExp).correctLickToLastLick = length(correctLickToLastLick);
end

disp(struct2table(nTrials(~ismember(1:length(exp), badExpIndices))))

clear iExp timeoutToPress timeoutToLick selCorrectPress selCorrectLick isCorrectPress isCorrectLick correctTimeoutToPress incorrectTimeoutToLick incorrectTimeoutToPress correctTimeoutToLick correctPressToLeverRetract correctPressToFirstLick correctPressToTimeout lick trialIndices lastLick correctPressToLastLick correctLickToLastLick iEu

%% Make/Plot ETA Heatmap


clear eta
eta.normWindow = [-3, -1.5];
eta.resolution = 0.05;
selUnits = isGoodUnit;
eta.pressNorm = eu.getETA('count', 'press', [-4, 4], selUnits=selUnits, resolution=eta.resolution, alignTo='stop', includeInvalid=false, normalize=eta.normWindow, minTrialDuration=0, maxTrialDuration=Inf);
eta.lickNorm = eu.getETA('count', 'lick', [-4, 4], selUnits=selUnits, resolution=eta.resolution, alignTo='stop', includeInvalid=false, normalize=eta.normWindow, minTrialDuration=0, maxTrialDuration=Inf);
eta.correctPressNorm = eu.getETA('count', 'CorrectPress', [-4, 4], selUnits=selUnits, resolution=eta.resolution, alignTo='stop', includeInvalid=false, normalize=eta.pressNorm.stats, minTrialDuration=0, maxTrialDuration=Inf);
eta.incorrectPressNorm = eu.getETA('count', 'IncorrectPress', [-4, 4], selUnits=selUnits, resolution=eta.resolution, alignTo='stop', includeInvalid=false, normalize=eta.pressNorm.stats, minTrialDuration=0, maxTrialDuration=Inf);
eta.correctLickNorm = eu.getETA('count', 'CorrectLick', [-4, 4], selUnits=selUnits, resolution=eta.resolution, alignTo='stop', includeInvalid=false, normalize=eta.pressNorm.stats, minTrialDuration=0, maxTrialDuration=Inf);
eta.incorrectLickNorm = eu.getETA('count', 'IncorrectLick', [-4, 4], selUnits=selUnits, resolution=eta.resolution, alignTo='stop', includeInvalid=false, normalize=eta.pressNorm.stats, minTrialDuration=0, maxTrialDuration=Inf);
eta.correctPressToLeverRetractNorm = eu.getETA('count', 'CorrectPressToLeverRetract', [-4, 4], selUnits=selUnits, resolution=eta.resolution, alignTo='stop', includeInvalid=false, normalize=eta.pressNorm.stats, minTrialDuration=0, maxTrialDuration=Inf);
eta.correctPressToFirstLickNorm = eu.getETA('count', 'CorrectPressToFirstLick', [-4, 4], selUnits=selUnits, resolution=eta.resolution, alignTo='stop', includeInvalid=false, normalize=eta.pressNorm.stats, minTrialDuration=0, maxTrialDuration=Inf);
eta.correctPressToLastLickNorm = eu.getETA('count', 'CorrectPressToLastLick', [-4, 4], selUnits=selUnits, resolution=eta.resolution, alignTo='stop', includeInvalid=false, normalize=eta.pressNorm.stats, minTrialDuration=0, maxTrialDuration=Inf);
eta.correctLickToLastLickNorm = eu.getETA('count', 'CorrectLickToLastLick', [-4, 4], selUnits=selUnits, resolution=eta.resolution, alignTo='stop', includeInvalid=false, normalize=eta.pressNorm.stats, minTrialDuration=0, maxTrialDuration=Inf);
% 
% eta.pressNorm = eu.getETA('count', 'press', [-4, 4], selUnits=selUnits, resolution=eta.resolution, alignTo='stop', includeInvalid=false, normalize=eta.normWindow, minTrialDuration=0, maxTrialDuration=Inf);
% eta.lickNorm = eu.getETA('count', 'lick', [-4, 4], selUnits=selUnits, resolution=eta.resolution, alignTo='stop', includeInvalid=false, normalize=eta.normWindow, minTrialDuration=0, maxTrialDuration=Inf);
% eta.correctPressNorm = eu.getETA('count', 'CorrectPress', [-4, 4], selUnits=selUnits, resolution=eta.resolution, alignTo='stop', includeInvalid=false, normalize=[-1, -0.3], minTrialDuration=0, maxTrialDuration=Inf);
% eta.incorrectPressNorm = eu.getETA('count', 'IncorrectPress', [-4, 4], selUnits=selUnits, resolution=eta.resolution, alignTo='stop', includeInvalid=false, normalize=[-1, -0.3], minTrialDuration=0, maxTrialDuration=Inf);
% eta.correctLickNorm = eu.getETA('count', 'CorrectLick', [-4, 4], selUnits=selUnits, resolution=eta.resolution, alignTo='stop', includeInvalid=false, normalize=[-1, -0.3], minTrialDuration=0, maxTrialDuration=Inf);
% eta.incorrectLickNorm = eu.getETA('count', 'IncorrectLick', [-4, 4], selUnits=selUnits, resolution=eta.resolution, alignTo='stop', includeInvalid=false, normalize=[-1, -0.3], minTrialDuration=0, maxTrialDuration=Inf);
% eta.correctPressToLeverRetractNorm = eu.getETA('count', 'CorrectPressToLeverRetract', [-4, 4], selUnits=selUnits, resolution=eta.resolution, alignTo='stop', includeInvalid=false, normalize=[-1, -0.3], minTrialDuration=0, maxTrialDuration=Inf);
% eta.correctPressToFirstLickNorm = eu.getETA('count', 'CorrectPressToFirstLick', [-4, 4], selUnits=selUnits, resolution=eta.resolution, alignTo='stop', includeInvalid=false, normalize=[-1, -0.3], minTrialDuration=0, maxTrialDuration=Inf);
% eta.correctPressToLastLickNorm = eu.getETA('count', 'CorrectPressToLastLick', [-4, 4], selUnits=selUnits, resolution=eta.resolution, alignTo='stop', includeInvalid=false, normalize=[-1, -0.3], minTrialDuration=0, maxTrialDuration=Inf);
% eta.correctLickToLastLickNorm = eu.getETA('count', 'CorrectLickToLastLick', [-4, 4], selUnits=selUnits, resolution=eta.resolution, alignTo='stop', includeInvalid=false, normalize=[-1, -0.3], minTrialDuration=0, maxTrialDuration=Inf);

% Lick bouts (norm to pre-press [-4, -2])
eta.lickBoutNaive = eu.getETA('count', 'lickbout_naive', window=[0, 2*pi*4], resolution=2*pi/30, normalize='none',  minInterval=0.05, maxInterval=0.20, ...
    minBoutCycles=2, maxBoutCycles=4, artifacts=artifactParams);
eta.lickBoutNaiveNorm = eta.lickBoutNaive;
eta.lickBoutNaiveNorm.X = (eta.lickBoutNaiveNorm.X - vertcat(eta.pressNorm.stats.mean)/0.1) ./ (vertcat(eta.pressNorm.stats.sd)/0.1);

%% Plottings
% close all
ETASORT = {eta.pressNorm, eta.lickNorm, eta.correctPressToFirstLickNorm, eta.correctPressToLastLickNorm, eta.correctPressToLeverRetractNorm};
SORTWINDOW = {[-0.3, 0.3], [-0.3, 0.3], [-0.1, 0.3], [-0.1, 0.3], [-0.1, 0.3]};
ETA = {eta.pressNorm, eta.lickNorm, eta.correctPressToFirstLickNorm, eta.correctPressToLastLickNorm, eta.correctLickToLastLickNorm, eta.correctPressToLeverRetractNorm};
NAME = ["Reach", "Lick", "FirstLick\n(correct-reach)", "LastLick\n(correct-reach)", "LastLick\n(correct-lick)", "BarRetract\n(correct-reach)"];
ZEROLABEL = ["touch", "lick", "lick", "last-lick", "last-lick", "bar-retract"];
XLIM = {[-1, 0.5], [-1, 0.5], [-0.3, 0.3], [-0.3, 0.3], [-0.3, 0.3], [-0.3, 0.3]};
w = cellfun(@(xl) round(10*diff(xl)), XLIM);
cw = [0, cumsum(w)];

% Combine ETA, PCA, and sort along 1st dimension
etaCombined = struct(X=[], t=[]);
etaCombined.X = cellfun(@(eta) eta.X, ETASORT, UniformOutput=false);
etaCombined.X = cat(2, etaCombined.X{:});
etaCombined.t = cellfun(@(eta) eta.t, ETASORT, UniformOutput=false);
etaCombined.t = cat(2, etaCombined.t{:});
etaCombined.epoch = arrayfun(@(i) i*ones(1, length(ETASORT{i}.t)), 1:length(ETASORT), UniformOutput=false);
etaCombined.epoch = cat(2, etaCombined.epoch{:});
etaCombined.X(etaCombined.X>3) = 3;
etaCombined.X(etaCombined.X<-1.5) = -1.5;

etaCombined.X = etaCombined.X(selUnits, :);

% For sorting, make templates to dot-product with
clear template
template(length(ETASORT)) = struct(t=[], x=[]);
for iETA = 1:length(ETASORT)
    template(iETA).t = etaCombined.t;
    template(iETA).x = zeros(1, length(etaCombined.t));
    template(iETA).x(1, isin(etaCombined.t, SORTWINDOW{iETA}) & etaCombined.epoch==iETA) = 1;
end


score = zeros(size(etaCombined.X, 1), length(ETASORT));
etaCombined.X(isnan(etaCombined.X)) = 0;
for iETA = 1:length(ETASORT)
    score(:, iETA) = etaCombined.X * template(iETA).x';
end
groupVar = arrayfun(@(i) bitshift(int16(score(:, i)>0), length(ETASORT)-i), 1:size(score, 2), UniformOutput=false);
groupVar = sum(horzcat(groupVar{:}), 2);

% First, sort by number of negative modulations
numNeg = sum(score<0, 2);
[uniqueGroupVars, ia] = unique(groupVar);
[~, I] = sort(numNeg(ia), 'ascend');
groupVar = changem(groupVar, 0:length(uniqueGroupVars)-1, uniqueGroupVars(I));

% Then, put all small groups (excluding single neg ones) at the bottom
[uniqueGroupVars, ia] = unique(groupVar);
assert(length(uniqueGroupVars) == max(groupVar)+1);
groupSize = histcounts(groupVar, 0:length(uniqueGroupVars));

numUnitsInSameGroup = arrayfun(@(gv) nnz(groupVar==gv), groupVar);
isRare = numUnitsInSameGroup < 3;
isSingleNeg = numNeg==1;
groupVar(isRare & ~isSingleNeg) = max(groupVar)+1;
% Tighten up the groupvars
uniqueGroupVars = unique(groupVar);
groupVar = changem(groupVar, 0:length(uniqueGroupVars)-1, uniqueGroupVars);
groupSize = histcounts(groupVar, 0:length(uniqueGroupVars));
groupSizeCum = cumsum(groupSize);
[~, sortOrder] = sort(double(groupVar)*10 + score(:, 1)./max(abs(score(:, 1))), 'ascend');


fig = figure(Units='inches', Position=[1, 1, 8, 5]);
tl = tiledlayout(fig, 1, sum(w), TileSpacing='loose', Padding='compact');
ax = gobjects(1, length(ETA));
for iAx = 1:length(ETA)
    hidecb = iAx < length(ETA);
    ax(iAx) = nexttile(tl, 1 + cw(iAx), [1, w(iAx)]);
    EphysUnit.plotETA(ax(iAx), ETA{iAx}, selUnits, xlim=XLIM{iAx}, clim=[-1.5, 1.5], order=sortOrder, hidecolorbar=hidecb);
    % applyCustomColormap(ax(iAx), [-1.5, 1.5], hlim=[0.375, 0, 0, -0.375], llim=[0.125, 0.5, 0.5, 0.25], hpwr=.5, lpwr=1, h0=0.33);
    applyCustomColormap(ax(iAx), [-1.5, 3], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);
    if ~hidecb
        ax(iAx).Colorbar.Layout.Tile = 'east';
    end
    if iAx > 1
        yticks(ax(iAx), [])
    else
        yticks(ax(iAx), groupSizeCum(1:end)+0.5)
        yticklabels(ax(iAx), string(groupSizeCum(1:end)))
    end
    title(ax(iAx), strsplit(NAME(iAx), "\\n"))
    xlabel(ax(iAx), "")
    ylabel(ax(iAx), "")
    xticks(ax(iAx), [XLIM{iAx}(1), 0, XLIM{iAx}(2)])
    xticklabels(ax(iAx), [string(1000*XLIM{iAx}(1)), ZEROLABEL(iAx), string(1000*XLIM{iAx}(2))])
    xtickangle(ax(iAx), 0)
    xline(ax(iAx), 0, 'k-')
    yline(ax(iAx), groupSizeCum(1:end-1)+0.5, 'k:', LineWidth=1.5)
    yline(ax(iAx), groupSizeCum([1, 1 + length(ETASORT)])+0.5, 'k-', LineWidth=1.5)
    ax(iAx).YAxis.TickLength = [0, 0];
end
xlabel(tl, "Time (ms)")
ylabel(tl, "Unit")
fontsize(fig, 9, 'points')



% %% Old heatmaps
% close all
% % EphysUnit.plotDoubleETA(eta.lick, eta.press, selUnits, 'lick', 'reach', clim=[-1, 1], xlim=xlDisp, sortWindow=[-3, 1], sortThreshold=0.25)
% 
% tl = tiledlayout(figure(Units='normalized', OuterPosition=[0.5, 0, 0.5, 1]), 1, 4);
% ax = nexttile(tl);
% [~, orderA] = EphysUnit.plotETA(ax, eta.lickNorm, isGoodUnit, event='lick', clim=[-1.5, 1.5], xlim=[-2, 2], sortWindow=[-3, 1], signWindow=[-.1, 0], sortThreshold=0.25, hideColorbar=true);
% xline(ax, 0)
% title(ax, 'Lick (sort A)')
% 
% ax = nexttile(tl);
% [~, ~] = EphysUnit.plotETA(ax, eta.pressNorm, isGoodUnit, event='reach', order=orderA, clim=[-1.5, 1.5], xlim=[-2, 2], sortWindow=[-3, 1], signWindow=[-.1, 0], sortThreshold=0.25, hideColorbar=true);
% xline(ax, 0)
% title(ax, 'Reach (sort A)')
% 
% ax = nexttile(tl);
% [~, orderB] = EphysUnit.plotETA(ax, eta.pressNorm, isGoodUnit, event='reach', clim=[-1.5, 1.5], xlim=[-2, 2], sortWindow=[-3, 1], signWindow=[-0.5, 0], sortThreshold=0.25, hideColorbar=false);
% ax.Colorbar.Layout.Tile = 'east';
% xline(ax, 0)
% title(ax, 'Reach (sort B)')
% 
% ax = nexttile(tl);
% [~, ~] = EphysUnit.plotETA(ax, eta.lickNorm, isGoodUnit, event='lick', order=orderB, clim=[-1.5, 1.5], xlim=[-2, 2], sortWindow=[-3, 1], signWindow=[-.1, 0], sortThreshold=0.25, hideColorbar=true);
% xline(ax, 0)
% title(ax, 'Lick (sort B)')
% 
% clear normWindow tl ax orderA orderB

%% Optotaggins

% ETA Stim
p.isiWindow = [-0.7, 0.3];
p.isiRes = 1e-3;
p.isiBaselineWindow = [-0.1, 0];
p.stimBluePowers = [2000, 8000, 16000]*1e-6; 
p.stimRedPowers = [2000, 8000, 16000]*1e-6;
p.stimBlueDurations = [20]*1e-3;
p.stimRedDurations = [20]*1e-3;

p.xlim.stim = [-0.1, 0.3];
p.xlim.move = [-4, 2];
p.rasterSzStim = 1;
p.rasterSzMove = 1;

close all
XBlue = cell(length(eu), 1);
XRed = cell(length(eu), 1);
for iEu = 1:length(eu)
    groupsBlue = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'duration'}, selectBy=struct(power=p.stimBluePowers, duration=p.stimBlueDurations, location=[], wavelength=[470, 473]));
    groupsRed = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'duration'}, selectBy=struct(power=p.stimRedPowers, duration=p.stimRedDurations, location=[], wavelength=635));

    % rd = eu(iEu).getRasterData('stimtwocolor', p.isiWindow, trials=[groupsRed.trials], alignTo='start', shutterDelay=0, sort=false, photoelectricBlankDuration=0.5e-3);
    % EphysUnit.plotRaster(rd)

    if ~isempty(groupsBlue)
        [isi, t] = eu(iEu).getMeanPEISI('stimtwocolor', [groupsBlue.trials], window=p.isiWindow, resolution=p.isiRes, photoelectricBlankDuration=1.5e-3);
        selBaseline = t>p.isiBaselineWindow(1) & t<p.isiBaselineWindow(2);
        normSR = (1./isi - mean(1./isi(:, selBaseline), 'omitnan')) ./ std(1./isi(:, selBaseline), 0, 2, 'omitnan');
        XBlue{iEu} = normSR;
    else
        XBlue{iEu} = NaN(size(t));
    end

    if ~isempty(groupsRed)
        [isi, t] = eu(iEu).getMeanPEISI('stimtwocolor', [groupsRed.trials], window=p.isiWindow, resolution=p.isiRes, ...
            photoelectricBlankDuration=1.5e-3);
        selBaseline = t>p.isiBaselineWindow(1) & t<p.isiBaselineWindow(2);
        normSR = (1./isi - mean(1./isi(:, selBaseline), 'omitnan')) ./ std(1./isi(:, selBaseline), 0, 2, 'omitnan');
        XRed{iEu} = normSR;
    else
        XRed{iEu} = NaN(size(t));
    end

    % fprintf('nan=%i, nan=%i\n', nnz(isnan(XBlue{iEu})), nnz(isnan(XRed{iEu})))
end

eta.stimBlue = struct(X=cat(1, XBlue{:}), t=t, N=[], D=[], stats=[]);
eta.stimRed = struct(X=cat(1, XRed{:}), t=t, N=[], D=[], stats=[]);

clear XBlue XRed iEu groupsBlue groupsRed isi t selBaseline normSR

%% Calculate META
clear meta
p.metaWindowStim = [0.005, 0.050];
p.posRespThresholdStim = 2;
p.negRespThresholdStim = -1.5;

t = eta.stimBlue.t;
meta.stimBlue = mean(eta.stimBlue.X(:, t>=p.metaWindowStim(1) & t<=p.metaWindowStim(2)), 2, 'omitnan');
t = eta.stimRed.t;
meta.stimRed = mean(eta.stimRed.X(:, t>=p.metaWindowStim(1) & t<=p.metaWindowStim(2)), 2, 'omitnan');
clear t

c.isStimBlueUp = meta.stimBlue >= p.posRespThresholdStim;
c.isStimBlueDown = meta.stimBlue <= p.negRespThresholdStim;
c.isStimRedUp = meta.stimRed >= p.posRespThresholdStim;
c.isStimRedDown = meta.stimRed <= p.negRespThresholdStim;

c.isStimBlueUpRedUpThereforeChrimsonMaybe = c.isStimBlueUp & c.isStimRedUp;
c.isStimBlueUpRedNotUpThereforeCoChrMaybe = c.isStimBlueUp & ~c.isStimRedUp;
c.isStimBlueNotUpRedUpThereforeChrimsonMaybe = ~c.isStimBlueUp & c.isStimRedUp;
c.isStimBlueNotUpRedNotUp = ~c.isStimBlueUp & ~c.isStimRedUp;

fprintf('ChrimsonR %i, CoChR %i\n', nnz(c.isStimRedUp), nnz(c.isStimBlueUpRedNotUpThereforeCoChrMaybe))
fprintf('Red %i, Blue %i\n', nnz(c.isStimRedUp), nnz(c.isStimBlueUp))


%% Plot stim heatmap
selUnitsStim = (c.isStimBlueUp | c.isStimRedUp) & isGoodUnit(:);
close all
ETASORT = {eta.stimRed, eta.stimBlue};
SORTWINDOW = {[5, 50]*1e-3, [5, 50]*1e-3};
ETA = {eta.stimRed, eta.stimBlue, eta.pressNorm, eta.lickNorm};
NAME = ["Red", "Blue", "Reach", "Lick"];
ZEROLABEL = ["stim", "stim", "touch", "lick"];
XLIM = {[-100, 100]*1e-3, [-100, 100]*1e-3, [-1000, 500]*1e-3, [-1000, 500]*1e-3};
CLIM = {[-10, 10], [-10, 10], [-1.5, 3], [-1.5, 3]};
CBTILE = {'west', '', 'east', ''};
w = cellfun(@(xl) round(10*diff(xl)), XLIM);
w(1:2) = w(1:2)*3;
cw = [0, cumsum(w)];

% Combine ETA, PCA, and sort along 1st dimension
etaCombined = struct(X=[], t=[]);
etaCombined.X = cellfun(@(eta) eta.X, ETASORT, UniformOutput=false);
etaCombined.X = cat(2, etaCombined.X{:});
etaCombined.t = cellfun(@(eta) eta.t, ETASORT, UniformOutput=false);
etaCombined.t = cat(2, etaCombined.t{:});
etaCombined.epoch = arrayfun(@(i) i*ones(1, length(ETASORT{i}.t)), 1:length(ETASORT), UniformOutput=false);
etaCombined.epoch = cat(2, etaCombined.epoch{:});
etaCombined.X(etaCombined.X>3) = 3;
etaCombined.X(etaCombined.X<-1.5) = -1.5;

etaCombined.X = etaCombined.X(selUnitsStim, :);

% For sorting, make templates to dot-product with
clear template
template(length(ETASORT)) = struct(t=[], x=[]);
for iETA = 1:length(ETASORT)
    template(iETA).t = etaCombined.t;
    template(iETA).x = zeros(1, length(etaCombined.t));
    template(iETA).x(1, isin(etaCombined.t, SORTWINDOW{iETA}) & etaCombined.epoch==iETA) = 1;
end


score = zeros(size(etaCombined.X, 1), length(ETASORT));
etaCombined.X(isnan(etaCombined.X)) = 0;
for iETA = 1:length(ETASORT)
    score(:, iETA) = etaCombined.X * template(iETA).x';
end
groupVar = arrayfun(@(i) bitshift(int16(score(:, i)>40), length(ETASORT)-i), 1:size(score, 2), UniformOutput=false);
groupVar = sum(horzcat(groupVar{:}), 2);

% First, sort by number of negative modulations
numNeg = sum(score>40, 2);
[uniqueGroupVars, ia] = unique(groupVar);
[~, I] = sort(numNeg(ia), 'ascend');
groupVar = changem(groupVar, 0:length(uniqueGroupVars)-1, uniqueGroupVars(I));

% Then, put all small groups (excluding single neg ones) at the bottom
[uniqueGroupVars, ia] = unique(groupVar);
assert(length(uniqueGroupVars) == max(groupVar)+1);
groupSize = histcounts(groupVar, 0:length(uniqueGroupVars));

numUnitsInSameGroup = arrayfun(@(gv) nnz(groupVar==gv), groupVar);
isRare = numUnitsInSameGroup < 3;
isSingleNeg = numNeg==1;
groupVar(isRare & ~isSingleNeg) = max(groupVar)+1;
% Tighten up the groupvars
uniqueGroupVars = unique(groupVar);
groupVar = changem(groupVar, 0:length(uniqueGroupVars)-1, uniqueGroupVars);
groupSize = histcounts(groupVar, 0:length(uniqueGroupVars));
groupSizeCum = cumsum(groupSize);
[~, sortOrder] = sort(double(groupVar)*10 + score(:, 1)./max(abs(score(:, 1))), 'ascend');


fig = figure(Units='inches', Position=[1, 1, 8, 5]);
tl = tiledlayout(fig, 1, sum(w), TileSpacing='loose', Padding='compact');
ax = gobjects(1, length(ETA));
for iAx = 1:length(ETA)
    hidecb = isempty(CBTILE{iAx});
    ax(iAx) = nexttile(tl, 1 + cw(iAx), [1, w(iAx)]);
    EphysUnit.plotETA(ax(iAx), ETA{iAx}, selUnitsStim, xlim=XLIM{iAx}, clim=[-1.5, 1.5], order=sortOrder, hidecolorbar=hidecb);
    % applyCustomColormap(ax(iAx), [-1.5, 1.5], hlim=[0.375, 0, 0, -0.375], llim=[0.125, 0.5, 0.5, 0.25], hpwr=.5, lpwr=1, h0=0.33);
    applyCustomColormap(ax(iAx), CLIM{iAx}, hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);
    if ~hidecb
        ax(iAx).Colorbar.Layout.Tile = CBTILE{iAx};
    end
    if iAx > 1
        yticks(ax(iAx), [])
    else
        yticks(ax(iAx), groupSizeCum(1:end)+0.5)
        yticklabels(ax(iAx), string(groupSizeCum(1:end)))
    end
    title(ax(iAx), strsplit(NAME(iAx), "\\n"))
    xlabel(ax(iAx), "")
    ylabel(ax(iAx), "")
    xticks(ax(iAx), [XLIM{iAx}(1), 0, XLIM{iAx}(2)])
    xticklabels(ax(iAx), [string(1000*XLIM{iAx}(1)), ZEROLABEL(iAx), string(1000*XLIM{iAx}(2))])
    xtickangle(ax(iAx), 0)
    xline(ax(iAx), 0, 'k-')
    if iAx <= 2
        xline(ax(iAx), 0.02, 'k-')
    end
    yline(ax(iAx), groupSizeCum(1:end-1)+0.5, 'k:', LineWidth=1.5)
    ax(iAx).YAxis.TickLength = [0, 0];
end
xlabel(tl, "Time (ms)")
ylabel(tl, "Unit")
fontsize(fig, 9, 'points')


%%
fig = figure(Units="inches", Position=[0 0 12.5, 7.5]);
ax = arrayfun(@(i) subplot(1, 2, i), 1:2);

[~, order] = EphysUnit.plotETA(ax(1), eta.stimBlue, sortWindow=[0, 50], signWindow=[0, 25], ...
    sortThreshold=0.25, negativeSortThreshold=0.25, ...
    clim=[-10, 10], xlim=[-50, 100], hidecolorbar=true, timeUnit='ms');
EphysUnit.plotETA(ax(2), eta.stimRed, order=order, ...
    clim=[-10, 10], xlim=[-50, 100], hidecolorbar=true, timeUnit='ms');

xline(ax(1), [0, 20])
xline(ax(2), [0, 20])