%% Load ephys units
eu = EphysUnit.load('\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\Units\TwoColor_SNr_SCRetro\ReverseInjection\SingleUnit_NonDuplicate_NonDrift_SNr', waveforms=false, spikecounts=false, spikerates=false);

%%
sessions = [...
    % "daisy27_20250624", ... SNr, traditional
    "daisy27_20250626", ... SNr, Accel Lick
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
exp.alignTimestamps(refEventNameArduino='LICK', refEventNameEphys='LickOn', trialDurationTolerance=0.1);

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

%% Verify aligment by looking at clips
iExp = 1;
iReward = 10;
[clip, t] = exp(iExp).getVideoClip(exp(iExp).eu(1).EventTimes.RewardOn(iReward), side='r', numFramesBefore=60, numFramesAfter=60, bodyParts={'handIpsi', 'handContra', 'jaw', 'tongue'}, ...
    minLikelihood=0.6);
implay(clip, 30)

%%
close all
for iExp = 1:length(exp)
    ax = axes(figure);
    hold(ax, 'on')
    [X, Y, L, t] = exp(iExp).getTrajectoryByTrial('r', 'handContra', trialType='press', window=[-1, 0], likelihoodThreshold=0.4);
    X = X - X(:, 1); Y = Y - Y(:, 1);
    plot(ax, median(X(:, end), 1, 'omitnan'), median(Y(:, end), 1, 'omitnan'), 'r-', Marker='x', MarkerSize=25, LineWidth=1.5, DisplayName='reach')
    scatter(ax, X(:, end), Y(:, end), 10, [1, 0, 0], DisplayName='reach')
    [X, Y, L, t] = exp(iExp).getTrajectoryByTrial('r', 'handContra', trialType='lick', window=[-1, 0], likelihoodThreshold=0.4);
    X = X - X(:, 1); Y = Y - Y(:, 1);
    plot(ax, median(X(:, end), 1, 'omitnan'), median(Y(:, end), 1, 'omitnan'), 'b-', Marker='x', MarkerSize=25, LineWidth=1.5, DisplayName='lick')
    scatter(ax, X(:, end), Y(:, end), 10, [0, 0, 1], DisplayName='lick')
    xline(ax, 0, 'k--')
    yline(ax, 0, 'k--')
    title(ax, exp(iExp).name, Interpreter='none')
    axis(ax, 'image')
    axis(ax, 'equal')
end