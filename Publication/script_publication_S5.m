
folders = [ ...
    ... "C:\SERVER\daisy26\daisy26_20250424", ...
    ... "C:\SERVER\daisy26\daisy26_20250425", ...
    ... "C:\SERVER\desmond38\desmond38_20250401", ...
    ... "C:\SERVER\desmond38\desmond38_20250402", ...
    ... "C:\SERVER\desmond38\desmond38_20250407", ...
    ... "C:\SERVER\desmond38\desmond38_20250417", ...
    ... "C:\SERVER\desmond39\desmond39_20250416", ...
    ... "C:\SERVER\desmond39\desmond39_20250423", ...
    "C:\SERVER\daisy27\daisy27_20250624", ...
    "C:\SERVER\daisy27\daisy27_20250626", ...
    "C:\SERVER\daisy27\daisy27_20250707", ...
    "C:\SERVER\daisy27\daisy27_20250715", ...
    "C:\SERVER\daisy27\daisy27_20250717", ...
    "C:\SERVER\daisy27\daisy27_20250721", ...
    "C:\SERVER\daisy27\daisy27_20250724", ...
    "C:\SERVER\daisy28\daisy28_20250701", ...
    "C:\SERVER\daisy28\daisy28_20250702", ...
    "C:\SERVER\daisy28\daisy28_20250716", ...
    "C:\SERVER\daisy28\daisy28_20250718", ...
    "C:\SERVER\daisy28\daisy28_20250723", ...
    "C:\SERVER\daisy28\daisy28_20250725", ...
    "C:\SERVER\daisy28\daisy28_20250729", ...
    "C:\SERVER\daisy28\daisy28_20250714", ...
    "C:\SERVER\daisy29\daisy29_20251023", ...
    "C:\SERVER\daisy29\daisy29_20251024", ...
    "C:\SERVER\daisy29\daisy29_20251025", ...
    "C:\SERVER\daisy29\daisy29_20251027", ...
    "C:\SERVER\daisy29\daisy29_20251028", ...
    "C:\SERVER\daisy29\daisy29_20251031", ...
    "C:\SERVER\daisy29\daisy29_20251118", ...
    "C:\SERVER\daisy29\daisy29_20251120", ...
    "C:\SERVER\daisy29\daisy29_20251103", ...
    "C:\SERVER\daisy30\daisy30_20251029", ...
    "C:\SERVER\daisy30\daisy30_20251030", ...
    "C:\SERVER\daisy30\daisy30_20251119", ...
    "C:\SERVER\daisy30\daisy30_20251121", ...
    "C:\SERVER\desmond41\desmond41_20251028", ...
    "C:\SERVER\desmond41\desmond41_20251029", ...
    "C:\SERVER\desmond42\desmond42_20251030", ...
    "C:\SERVER\desmond42\desmond42_20251031", ...
    "C:\SERVER\desmond41\desmond41_20251104", ...
    "C:\SERVER\desmond41\desmond41_20251117", ...
    "C:\SERVER\desmond42\desmond42_20251121", ...
    "C:\SERVER\desmond41\desmond41_20251124", ...
    ];

euFolders = [ ...
    % "C:\SERVER\Units\TwoColor_SNr_SCRetro\SingleUnit_NonDuplicate_NonDrift_SNr", ... daisy26, desmond38, desmond39
    "C:\SERVER\Units\TwoColor_SNr_SCRetro\ReverseInjection\SingleUnit_NonDuplicate_NonDrift_SNr_withTrials", ... daisy27, 28
    "C:\SERVER\Units\TwoColor_Striatonigral\SingleUnit_NonDuplicate_NonDrift_SNr", ... daisy29, 30, desmond41, 42
];

dlcResultsPath = 'C:\SERVER\DeepLabCut\Results\FourPawsTongueJawSpine';

excludeAnimals = {'daisy26', 'desmond38', 'desmond39'}; % Lick artifacts, this is pre-accelerometer, we should remove anything before 20250530 (first lick/reach accelerometer test)

%% Load ephysunits
eu = cell(length(euFolders), 1);
for iDir = 1:length(euFolders)
    eu{iDir} = EphysUnit.load(euFolders(iDir));
end

eu = cat(2, eu{:});


%% Make CompleteExperiment3

exp = CompleteExperiment3(eu, cameras='lr', deeplabcutPath=dlcResultsPath);

exp.alignTimestamps(refEventNameArduino={'OPTO1_ON', 'OPTO2_ON', 'REWARD_ON'}, refEventNameEphys={'StimOn', 'RewardOn'}, trialDurationTolerance=2);

clear results
results(length(exp)) = struct(name=[], varsL=[], hasTimestampsL=[], varsR=[], hasTimestampsR=[], isValid=[]);
for iExp = 1:length(exp)
    results(iExp).name = exp(iExp).name;
    results(iExp).varsR = string(exp(iExp).vtdR.Properties.VariableNames)';
    results(iExp).hasTimestampsR = ismember('Timestamp', exp(iExp).vtdR.Properties.VariableNames);

    if ~isempty(exp(iExp).vtdL)
        results(iExp).varsL = string(exp(iExp).vtdL.Properties.VariableNames)';
        results(iExp).hasTimestampsL = ismember('Timestamp', exp(iExp).vtdL.Properties.VariableNames);
    else
        results(iExp).hasTimestampsL = false;
    end

    results(iExp).isValid = length(results(iExp).varsL) == 23 && length(results(iExp).varsR) == 23;
end
clear iExp

% Cull bad experiments (some are still pending DLC)
assert(isequal({exp.name}, {results.name}));
exp = exp([results.isValid]);
results = results([results.isValid]);
eu = [exp.eu];

% Rename variables
from = ["HandIpsiCam", "HandContraCam", "LegIpsiCam", "LegContraCam", "Jaw", "Tongue", "Spine"];
toL = ["HandL", "HandR", "FootL", "FootR", "Jaw", "Tongue", "Spine"];
toR = ["HandR", "HandL", "FootR", "FootL", "Jaw", "Tongue", "Spine"];
for iExp = 1:length(exp)
    if results(iExp).hasTimestampsL
        vtd = exp(iExp).vtdL;
        for i = 1:length(toL)
            if ismember(sprintf("%s_X", from(i)), results(iExp).varsL)
                vtd = renamevars(vtd, ...
                    [sprintf("%s_X", from(i)), sprintf("%s_Y", from(i)), sprintf("%s_Likelihood", from(i))], ...
                    [sprintf("%s_X", toL(i)), sprintf("%s_Y", toL(i)), sprintf("%s_Likelihood", toL(i))]);
            end
        end
        exp(iExp).vtdL = vtd;
    end
    if results(iExp).hasTimestampsR
        vtd = exp(iExp).vtdR;
        for i = 1:length(toR)
            if ismember(sprintf("%s_X", from(i)), results(iExp).varsR)
                vtd = renamevars(vtd, ...
                    [sprintf("%s_X", from(i)), sprintf("%s_Y", from(i)), sprintf("%s_Likelihood", from(i))], ...
                    [sprintf("%s_X", toR(i)), sprintf("%s_Y", toR(i)), sprintf("%s_Likelihood", toR(i))]);
            end
        end
        exp(iExp).vtdR = vtd;
    end
end

%% Reset results
clearvars -except eu exp results

expIndices = cellfun(@(name) find(strcmpi(name, {exp.name}), 1, 'first'), {eu.ExpName});

%%

[ac.ac, ac.euIndicesFirstInSession, ac.expIndices, ac.uniqueExpNames] = eu.loadArduinoConnection();
eu.alignTimestamps(["WAITFORTOUCH", "TIMEOUT_START"], ac={ac.ac, ac.euIndicesFirstInSession, ac.expIndices, ac.uniqueExpNames}, acRefEventName="REWARD_ON", euRefEventName="RewardOn");


%% Both were recording in left hemisphere, and in the DLC labeling, HandContra/HandIpsi refers to camera-side/far-side, we need to fix this for both the left and right cameras. Left/right cameras are placed to the left/right side of mouse.
% lcam: HandContra -> handContra, HandIpsi -> handIpsi
% rcam: HandContra -> handIpsi, HandIpsi -> handContra

for iExp = 1:length(exp)
    switch exp(iExp).animalName
        case {'daisy27', 'daisy28'} % Left hemi implant
            exp(iExp).vtdL = renamevars(exp(iExp).vtdL, ...
                [ ...
                    "HandR_X", "HandR_Y", "HandR_Likelihood", ...
                    "HandL_X", "HandL_Y", "HandL_Likelihood", ...
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
                    "HandR_X", "HandR_Y", "HandR_Likelihood", ...
                    "HandL_X", "HandL_Y", "HandL_Likelihood", ...
                    "Tongue_X", "Tongue_Y", "Tongue_Likelihood", ...
                    "Jaw_X", "Jaw_Y", "Jaw_Likelihood", ...
                ], ...
                [ ...
                    "handContra_X", "handContra_Y", "handContra_Likelihood", ...
                    "handIpsi_X", "handIpsi_Y", "handIpsi_Likelihood", ...
                    "tongue_X", "tongue_Y", "tongue_Likelihood", ...
                    "jaw_X", "jaw_Y", "jaw_Likelihood", ...
                ]);            
        case {'daisy29', 'daisy30', 'daisy33', 'daisy34', 'desmond41', 'desmond42', 'desmond43'} % Right hemi implant
            exp(iExp).vtdL = renamevars(exp(iExp).vtdL, ...
                [ ...
                    "HandL_X", "HandL_Y", "HandL_Likelihood", ...
                    "HandR_X", "HandR_Y", "HandR_Likelihood", ...
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
                    "HandL_X", "HandL_Y", "HandL_Likelihood", ...
                    "HandR_X", "HandR_Y", "HandR_Likelihood", ...
                    "Tongue_X", "Tongue_Y", "Tongue_Likelihood", ...
                    "Jaw_X", "Jaw_Y", "Jaw_Likelihood", ...
                ], ...
                [ ...
                    "handContra_X", "handContra_Y", "handContra_Likelihood", ...
                    "handIpsi_X", "handIpsi_Y", "handIpsi_Likelihood", ...
                    "tongue_X", "tongue_Y", "tongue_Likelihood", ...
                    "jaw_X", "jaw_Y", "jaw_Likelihood", ...
                ]);           
        otherwise
            error('Unknown animal %s', exp(iExp).animalName)
    end
end

%% Do Lick detection by video
% %% Make BennyHill-style videos of lick frames, grouped by tongue likelihoods
% llhThresholds = 0.4:0.2:1;
% iExp = 2;
% 
% side = 'l';
% clear frames
% frames(length(llhThresholds)-1) = struct(frames=[], t=[], llhWindow=[]);
% for i = 1:length(llhThresholds)-1
%     frames(i).llhWindow = llhThresholds(i:i+1);
%     selFrames = isin(exp(iExp).vtdL.tongue_Likelihood, frames(i).llhWindow);
%     switch side
%         case 'l'
%             [frames(i).frames, frames(i).t] = exp(iExp).getVideoFrames(exp(iExp).vtdL.Timestamp(selFrames), 'l');    
%         case 'r'
%             [frames(i).frames, frames(i).t] = exp(iExp).getVideoFrames(exp(iExp).vtdR.Timestamp(selFrames), 'r');
%     end
% end
% 
% clear llhThresholds iExp i selFrames

% Let's choose a tongue llh>0.5 as lick
tongueLikelihoodThreshold = 0.5;
for iExp = 1:length(exp)
    % try
        l1 = exp(iExp).vtdL.tongue_Likelihood;
        l2 = interp1(exp(iExp).vtdR.Timestamp, exp(iExp).vtdR.tongue_Likelihood, exp(iExp).vtdL.Timestamp, 'linear');

        isTongueVisible = l1' > tongueLikelihoodThreshold | l2' > tongueLikelihoodThreshold;
        lickOn = exp(iExp).vtdL.Timestamp(strfind(isTongueVisible, [false, true]) + 1);
        lickOff = exp(iExp).vtdL.Timestamp(strfind(isTongueVisible, [true, false]) + 1);

        for iEu = 1:length(exp(iExp).eu)
            exp(iExp).eu(iEu).EventTimes.Lick = lickOn;
            exp(iExp).eu(iEu).EventTimes.LickOn = lickOn;
            exp(iExp).eu(iEu).EventTimes.LickOff = lickOff;
        end

        fprintf('%i - %s, %i(%i) lickOn, %i(%i) lickOff by video(accelorometer).\n', iExp, exp(iExp).name, length(lickOn), length(exp(iExp).eu(iEu).EventTimes.LickOn), length(lickOff), length(exp(iExp).eu(iEu).EventTimes.LickOff));
    % catch
    %     warning('exp %i failed', iExp)
    % end
end

%% Make lick and reach trials, but first cleanup the bad accel artifacts
pulseWidthThreshold = struct(press = 1.5e-3, lick=0);
minTrialLength = 1;

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
    try
        if any(isfield(eu(iEu).EventTimes, {'ValidPress', 'ValidLick'}))
            continue
        end
        selPress = eu(iEu).EventTimes.PressOff - eu(iEu).EventTimes.PressOn > pulseWidthThreshold.press;
        eu(iEu).Trials.Press = Trial(eu(iEu).EventTimes.TIMEOUT_START, eu(iEu).EventTimes.Press(selPress), stopMode='first', exclude=eu(iEu).EventTimes.Lick);
        eu(iEu).Trials.Press = eu(iEu).Trials.Press(eu(iEu).Trials.Press.duration() >= minTrialLength);
        eu(iEu).EventTimes.ValidPress = eu(iEu).EventTimes.Press(selPress);
    
        eu(iEu).EventTimes.PressOn = eu(iEu).EventTimes.PressOn(selPress);
        eu(iEu).EventTimes.PressOff = eu(iEu).EventTimes.PressOff(selPress);
    
        % selLick = eu(iEu).EventTimes.LickOff - eu(iEu).EventTimes.LickOn > pulseWidthThreshold.lick;
        eu(iEu).Trials.Lick = Trial(eu(iEu).EventTimes.TIMEOUT_START, eu(iEu).EventTimes.Lick, stopMode='first', exclude=eu(iEu).EventTimes.ValidPress);
        eu(iEu).Trials.Lick = eu(iEu).Trials.Lick(eu(iEu).Trials.Lick.duration() >= minTrialLength);
        eu(iEu).EventTimes.ValidLick = eu(iEu).EventTimes.Lick;  
    catch
        warning('Could not process unit %i', iEu)
    end
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

%% Calculate arm traversal distance during reach vs. lick trials.
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

c.isGoodUnit = true(size(eu));
for iExp = badExpIndices
    if isempty(iExp)
        continue
    end
    c.isGoodUnit(ismember(eu, exp(iExp).eu)) = false;
end

fprintf('%i/%i sessions (%i/%i units) failed to generate trajectories and will not be included.\n', length(badExpIndices), length(exp), nnz(~c.isGoodUnit), length(eu));
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
    selTrials.lick = isin(traj(iExp).lick.handContra.traversal, quantile(traj(iExp).press.handContra.traversal, [0, 0.8]), true) & isin(traj(iExp).lick.handIpsi.traversal, quantile(traj(iExp).press.handIpsi.traversal, [0, 0.8]), true);

    if eu(iEu) == exp(iExp).eu(1)
        fprintf('Exp%i %s: kept %i/%i reach trials, kept %i/%i lick trials.\n', iExp, exp(iExp).name, nnz(selTrials.press), length(selTrials.press), nnz(selTrials.lick), length(selTrials.lick))
    end

    if nnz(selTrials.press) < minNumTrials || nnz(selTrials.lick) < minNumTrials
        c.isGoodUnit(iEu) = false;
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
    % correctPressToLeverRetract = Trial([correctTimeoutToPress.Stop, Inf], exp(iExp).eu(1).EventTimes.LEVER_RETRACT_START, 'first', exclude=exp(iExp).eu(1).EventTimes.TIMEOUT_START);

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

%%
nTrialsPress = arrayfun(@(eu) length(eu.Trials.Press), eu);
nTrialsLick = arrayfun(@(eu) length(eu.Trials.Lick), eu);
c.hasPress = nTrialsPress >= 15;
c.hasLick = nTrialsLick >= 15;

clear eta
eta.normWindow = [-3, -1.5];
eta.resolution = 0.1;
selUnits = c.hasPress & c.hasLick;
eta.pressNorm = eu.getETA('count', 'press', [-4, 4], selUnits=selUnits, resolution=eta.resolution, alignTo='stop', includeInvalid=false, normalize=eta.normWindow, minTrialDuration=0, maxTrialDuration=Inf);
eta.lickNorm = eu.getETA('count', 'lick', [-4, 4], selUnits=selUnits, resolution=eta.resolution, alignTo='stop', includeInvalid=false, normalize=eta.normWindow, minTrialDuration=0, maxTrialDuration=Inf);

meta.press = mean(eta.pressNorm.X(:, isin(eta.pressNorm.t, p.responseWindowPress)), 2, 'omitnan');
meta.lick = mean(eta.lickNorm.X(:, isin(eta.lickNorm.t, p.responseWindowLick)), 2, 'omitnan');

%%
close all
fig = figure(Units='inches', Position=[1, 1, 4, 4]);
tl = tiledlayout(fig, 1, 2);

nnz(sel)

ax = gobjects(1, 2);
ax(1) = nexttile(tl);
ax(2) = nexttile(tl);

sel = c.hasPress & c.hasLick;

groupVar = NaN(length(sel), 1);
groupVar(meta.press<0 & meta.lick>0) = 0;
groupVar(meta.press>0 & meta.lick<0) = 1;
groupVar(meta.press<0 & meta.lick<0) = 2;
groupVar(meta.press>0 & meta.lick>0) = 3;
groupVar = groupVar(sel);
[~, order] = EphysUnit.plotETA(ax(1), eta.pressNorm, sel, sortGroup=groupVar, ...
    clim=[-1.5, 1.5], xlim=[-2.5, 0.5], sortWindow=[-3, 0], signWindow=[-0.3, 0], ...
    sortThreshold=0.25, negativeSortThreshold=0.25, hidecolorbar=true);
EphysUnit.plotETA(ax(2), eta.lickNorm, sel, order=order, sortGroup=groupVar, ...
    clim=[-1.5, 1.5], xlim=[-2.5, 0.5], hidecolorbar=true);

for iAx = 1:2
    applyCustomColormap(ax(iAx), [-1.5, 1.5], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);
    % applyCustomColormap(ax(iAx), [-1.5, 1.5], hlim=[0.375, 0, -0.375, -0.375], llim=[0.1, 1, 1, 0.3], hpwr=.3, lpwr=0.25, h0=0.33);
    % applyCustomColormap(ax(iAx), [-1.5, 1.5], hlim=[0.375, 0, 0, -0.375], llim=[0.125, 0.5, 0.5, 0.25], hpwr=.5, lpwr=1, h0=0.33);
end

N = histcounts(groupVar, -0.5:2:3.5);
yline(ax(1), cumsum(N(1:end-1)) + 1, 'k--');
yline(ax(2), cumsum(N(1:end-1)) + 1, 'k--');

xline(ax(1), 0, 'k--')
xline(ax(2), 0, 'k--')
ylim(ax(1:2), [0, nnz(sel)+1])
yt = 0:100:nnz(sel);
yt(1) = 1;
if round(yt(end)./100) == round(nnz(sel)./100)
    yt(end) = nnz(sel);
else
    yt(end + 1) = nnz(sel);
end
yt = unique(yt);
yticks(ax(1), [1, cumsum(N(1:end-1)) + 1, nnz(sel)])
yticks(ax(2), [])

title(ax(1), 'reach')
title(ax(2), 'lick')
ylabel(tl, 'unit', FontSize=p.fontSize)
ylabel(ax, '')
xlabel(tl, 'time to bar/spout contact (s)', FontSize=p.fontSize)
xlabel(ax, '')
fontsize(ax, p.fontSize, 'points')
fontname(ax, 'Arial')
axc = ax(2);

cb = colorbar(ax(1));
cb.Layout.Tile = 'east';
cb.Label.String = 'normalized spike rate (a.u.)';
cb.Label.Position(1) = 0;
cb.Label.VerticalAlignment = 'bottom';

% hLetter = text(ax(1), 0, 0, 'b', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
% ax(1).Units = 'inches';
% hLetter.HorizontalAlignment = 'right';
% hLetter.VerticalAlignment = 'top';
% hLetter.Position = [-0.3, ax(1).Position(4) + 0.1, 0];

copygraphics(fig, ContentType='vector', BackgroundColor='none')



%% Boot movement responses
p.bootAlpha = 0.01;
p.nboot = 100000;
p.responseWindowPress = [-0.3, 0];
p.responseWindowLick = [-0.3, 0];
assert(isequal(p.responseWindowPress, [-0.3, 0]))
assert(isequal(p.responseWindowLick, [-0.3, 0]))
boot.press = struct('h', NaN(length(eu), 1), 'muDiffCI', NaN(length(eu), 2), 'muDiffObs', NaN(length(eu), 1));
boot.lick = struct('h', NaN(length(eu), 1), 'muDiffCI', NaN(length(eu), 2), 'muDiffObs', NaN(length(eu), 1));
[boot.press.h(c.isGoodUnit), boot.press.muDiffCI(c.isGoodUnit, :), boot.press.muDiffObs(c.isGoodUnit)] = bootstrapMoveResponse( ...
    eu(c.isGoodUnit), 'press', nboot=p.nboot, alpha=p.bootAlpha, withReplacement=false, oneSided=false, ...
    responseWindow=p.responseWindowPress);
[boot.lick.h(c.isGoodUnit), boot.lick.muDiffCI(c.isGoodUnit, :), boot.lick.muDiffObs(c.isGoodUnit)] = bootstrapMoveResponse( ...
    eu(c.isGoodUnit), 'lick', nboot=p.nboot, alpha=p.bootAlpha, withReplacement=false, oneSided=false, ...
    responseWindow=p.responseWindowLick);
fprintf(1, '\nAll done\n')

%% Report bootstraped movement response direction


% figure, histogram(boot.press.h)
c.isPressUp = boot.press.h' == 1 & c.isGoodUnit & c.hasPress;
c.isPressDown = boot.press.h' == -1 & c.isGoodUnit & c.hasPress;
c.isPressResponsive = (c.isPressUp | c.isPressDown);

% figure, histogram(boot.lick.h)
c.isLickUp = boot.lick.h' == 1 & c.isGoodUnit & c.hasLick;
c.isLickDown = boot.lick.h' == -1 & c.isGoodUnit & c.hasLick;
c.isLickResponsive = (c.isLickUp | c.isLickDown);

fprintf(1, ['%g good units, %g modulated for reach:\n' ...
    '\t%g (%.0f%%) are excited (p<%g);\n' ...
    '\t%g (%.0f%%) are inhibited (p<%g).\n'], ...
    nnz(c.hasPress), nnz(c.isPressResponsive), ...
    nnz(c.isPressUp), 100*nnz(c.isPressUp)/nnz(c.isPressResponsive), p.bootAlpha, ...
    nnz(c.isPressDown), 100*nnz(c.isPressDown)/nnz(c.isPressResponsive), p.bootAlpha);

fprintf(1, ['%g good units, %g modulated for lick:\n' ...
    '\t%g (%.0f%%) are excited (p<%g);\n' ...
    '\t%g (%.0f%%) are inhibited (p<%g).\n'], ...
    nnz(c.hasLick), nnz(c.isLickResponsive), ...
    nnz(c.isLickUp), 100*nnz(c.isLickUp)/nnz(c.isLickResponsive), p.bootAlpha, ...
    nnz(c.isLickDown), 100*nnz(c.isLickDown)/nnz(c.isLickResponsive), p.bootAlpha);

nTotal = nnz(c.isPressResponsive & c.isLickResponsive);
fprintf(1, ['%g good units, %g modulated for reach AND lick:\n' ...
    '\t%g (%.0f%%) are press-excited AND lick-excited;\n' ...
    '\t%g (%.0f%%) are press-inhibited AND lick-inhibited;\n' ...
    '\t%g (%.0f%%) are press-excited AND lick-inhibited;\n' ...
    '\t%g (%.0f%%) are press-inhibited AND lick-excited;\n'], ...
    nnz(c.hasPress & c.hasLick), nTotal, ...
    nnz(c.isPressUp & c.isLickUp), 100*nnz(c.isPressUp & c.isLickUp)/nTotal, ...
    nnz(c.isPressDown & c.isLickDown), 100*nnz(c.isPressDown & c.isLickDown)/nTotal, ...
    nnz(c.isPressUp & c.isLickDown), 100*nnz(c.isPressUp & c.isLickDown)/nTotal, ...
    nnz(c.isPressDown & c.isLickUp), 100*nnz(c.isPressDown & c.isLickUp)/nTotal)   

sel = c.hasPress & c.hasLick;
fprintf('05 Calculate: Of %i: %i (%i%%) showed modulation for BOTH, %i (%i%%) showed modulation for lick only, %i (%i%%) showed modulation for reach only, %i (%i%%) for neither.', ...
    nnz(sel), ...
    nnz(sel & c.isPressResponsive & c.isLickResponsive), round(nnz(sel & c.isPressResponsive & c.isLickResponsive)/nnz(sel)*100), ...
    nnz(sel & c.isLickResponsive & ~c.isPressResponsive), round(nnz(sel & c.isLickResponsive & ~c.isPressResponsive)/nnz(sel)*100), ...
    nnz(sel & c.isPressResponsive & ~c.isLickResponsive), round(nnz(sel & c.isPressResponsive & ~c.isLickResponsive)/nnz(sel)*100), ...
    nnz(sel & ~c.isPressResponsive & ~c.isLickResponsive), round(nnz(sel & ~c.isPressResponsive & ~c.isLickResponsive)/nnz(sel)*100) ...
    )

clear nTotal sel

save('C:\SERVER\Units\boot_1225_reach_vs_lick.mat', 'boot', 'c', 'eta')



% 627 good units, 207 modulated for reach:
% 	164 (79%) are excited (p<0.01);
% 	43 (21%) are inhibited (p<0.01).
% 1084 good units, 216 modulated for lick:
% 	164 (76%) are excited (p<0.01);
% 	52 (24%) are inhibited (p<0.01).
% 600 good units, 48 modulated for reach AND lick:
% 	35 (73%) are press-excited AND lick-excited;
% 	9 (19%) are press-inhibited AND lick-inhibited;
% 	2 (4%) are press-excited AND lick-inhibited;
% 	2 (4%) are press-inhibited AND lick-excited;
% 05 Calculate: Of 600: 48 (8%) showed modulation for BOTH, 53 (9%) showed modulation for lick only, 151 (25%) showed modulation for reach only, 348 (58%) for neither.>> 

%% Save results
eu.save('C:\SERVER\Units\ReachVsLick_1225\')
save('C:\SERVER\Units\meta_ReachVsLick_1225_20260728.mat', 'boot', 'c', 'eta')