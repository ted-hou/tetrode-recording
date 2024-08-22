zThreshold = 2;
expName = 'daisy24_20240627';
%% Load data
load(sprintf('C:\\SERVER\\PawAnalysis\\EMG_%s.mat', expName))
files = dir(sprintf('C:\\SERVER\\Units\\acute_3cam_reach_direction_2tgts\\SingleUnits_NonDuplicate\\%s_Channel*_Unit*.mat', expName));
assert(~isempty(files))
files = arrayfun(@(f) sprintf('%s\\%s', f.folder, f.name), files, UniformOutput=false);

eu = EphysUnit.load(files);
trials = eu(1).makeTrials('press_spontaneous2');
goodTrials = trials(trials.duration() > 4);% & arrayfun(@(t) all(t <= eu.EventTimes.RewardTimes | t >= eu.EventTimes.RewardTimes + 4), [trials.Stop])) ;
eu(1).Trials.PressSpontaneous = goodTrials;
pa = Pawnalyzer2(eu, refEvent='press');
pa.getClips(noImage=true, nFramesBefore=15, nFramesAfter=0, keepData=false, trials='PressSpontaneous');
pa.load('auto')

%% Smooth EMG
ai.SmoothData = smoothdata(ai.Data, 'movmean', 1000); % 1000 sample (33.3ms) moving average window
ai.SmoothNormData = (ai.SmoothData - median(ai.SmoothData))./mad(ai.SmoothData)*0.6745;

%% %Calculate paw trajectories
clear traj
traj = struct(contra=[], ipsi=[], target=[], t=[], pca=[], kmeans=[]);
nt = 10;
nFrames = 16;
[traj.contra, traj.ipsi, traj.target, traj.t] = pa.getTrajectories(1, zero=nFrames - nt + 1);
if ~ismember(pa.exp.animalName, {'desmond29', 'daisy25'})
    traj.contra.x = -traj.contra.x;
    traj.ipsi.x = -traj.ipsi.x;
end

%% Calculate video-movement onset times
nTrials = length(traj.target);
trueStartTime = NaN(nTrials, 1);
t = traj.t;
for iTrial = 1:nTrials
    data = [ ...
        traj.contra.x(iTrial, :); ...
        traj.contra.y(iTrial, :); ...
        traj.contra.z(iTrial, :); ...
        traj.ipsi.x(iTrial, :); ...
        traj.ipsi.y(iTrial, :); ...
        traj.ipsi.z(iTrial, :); ...
        ];
    dist = sqrt(sum(data.^2, 1));
    normDist = normalize(dist, 'zscore', 'robust')*0.6745;
    if any(isnan(normDist))
        continue
    end
    lastIndexAboveThreshold = find(normDist > zThreshold, 1, 'last'); % We don't do abs since dist is already positive, although z-scoring will generate negatives, true movement should be positive.
    if isempty(lastIndexAboveThreshold)
        continue
    end
    trueStartIndex = find((normDist < zThreshold) & ((1:nFrames) <= lastIndexAboveThreshold), 1, 'last');
    if ~isempty(trueStartIndex)
        trueStartTime(iTrial) = t(trueStartIndex) - t(end);
    end
end
figure
histogram(trueStartTime, -0.5:0.05:0, Normalization='probability')
xlabel('Time to touch (s)')
ylabel('Probability')
title(sprintf('Video detected movement onset time (when z(dist) > %g)', zThreshold))

%% Calculate touch aligned EMG
baseWindow = [-4, -2];
clear emg
tLocal = -4:1/30000:0;
tEvent = reshape([eu(1).Trials.PressSpontaneous.Stop], [], 1);
tGlobal = tLocal + tEvent;
nTrials = length(tEvent);
X = zeros(nTrials, length(tLocal));
for iTrial = 1:nTrials
    X(iTrial, :) = interp1(ai.Timestamps, ai.SmoothData, tGlobal(iTrial, :), 'linear');
end

selBase = tLocal >= baseWindow(1) & tLocal <= baseWindow(2);
baseMeanByTrial = mean(X(:, selBase), 2, 'omitnan');
baseSdByTrial = std(X(:, selBase), 1, 2, 'omitnan');
normX = (X - baseMeanByTrial)./baseSdByTrial;

emg.touchAligned.t = tLocal;
emg.touchAligned.T = tGlobal;
emg.touchAligned.X = X;
emg.touchAligned.normX = normX;
clear tLocal tEvent tGlobal nTrials X iTrial

%% Calculate video-movement-onset aligned EMG
tLocal = -4:1/30000:0.5;
tEvent = reshape([eu(1).Trials.PressSpontaneous.Stop], [], 1) + reshape(trueStartTime, [], 1);
tGlobal = tLocal + tEvent;
nTrials = length(tEvent);
X = zeros(nTrials, length(tLocal));
for iTrial = 1:nTrials
    X(iTrial, :) = interp1(ai.Timestamps, ai.SmoothData, tGlobal(iTrial, :), 'linear');
end

selBase = tLocal >= baseWindow(1) & tLocal <= baseWindow(2);
baseMeanByTrial = mean(X(:, selBase), 2, 'omitnan');
baseSdByTrial = std(X(:, selBase), 1, 2, 'omitnan');
normX = (X - baseMeanByTrial)./baseSdByTrial;

emg.vidOnsetAligned.t = tLocal;
emg.vidOnsetAligned.T = tGlobal;
emg.vidOnsetAligned.X = X;
emg.vidOnsetAligned.normX = normX;
clear tLocal tEvent tGlobal nTrials X iTrial

%% Trial-average EMG aligned to touch
figure
plot(emg.touchAligned.t, mean(emg.touchAligned.normX, 1, 'omitnan'))
xlabel('Time to touch (s)')
ylabel('Smoothed normalzied EMG')
title('Trial-average EMG aligned to touch')

%% Detect movement onset time based on z(EMG)>2.5
nTrials = length(traj.target);
trueStartTimeEMG = NaN(nTrials, 1);
t = emg.touchAligned.t;
for iTrial = 1:nTrials
    x = emg.touchAligned.normX(iTrial, :);
    lastIndexAboveThreshold = find(x > zThreshold, 1, 'last'); % We don't do abs since dist is already positive, although z-scoring will generate negatives, true movement should be positive.
    if isempty(lastIndexAboveThreshold)
        continue
    end
    trueStartIndex = find((x < zThreshold) & ((1:length(t)) <= lastIndexAboveThreshold), 1, 'last');
    if ~isempty(trueStartIndex)
        trueStartTimeEMG(iTrial) = t(trueStartIndex) - t(end);
    end
end
figure
histogram(trueStartTimeEMG, -0.5:0.05:0, Normalization='probability')
xlabel('Time to touch (s)')
ylabel('Probability')
title(sprintf('EMG detected movement onset time (when z(EMG) > %g)', zThreshold))

%% Calculate EMG-onset aligned EMG
tLocal = -4:1/30000:0.5;
tEvent = reshape([eu(1).Trials.PressSpontaneous.Stop], [], 1) + reshape(trueStartTimeEMG, [], 1);
tGlobal = tLocal + tEvent;
nTrials = length(tEvent);
X = zeros(nTrials, length(tLocal));
for iTrial = 1:nTrials
    X(iTrial, :) = interp1(ai.Timestamps, ai.SmoothData, tGlobal(iTrial, :), 'linear');
end

selBase = tLocal >= baseWindow(1) & tLocal <= baseWindow(2);
baseMeanByTrial = mean(X(:, selBase), 2, 'omitnan');
baseSdByTrial = std(X(:, selBase), 1, 2, 'omitnan');
normX = (X - baseMeanByTrial)./baseSdByTrial;

emg.emgOnsetAligned.t = tLocal;
emg.emgOnsetAligned.T = tGlobal;
emg.emgOnsetAligned.X = X;
emg.emgOnsetAligned.normX = normX;
clear tLocal tEvent tGlobal nTrials X iTrial
%% Plot EMG-Video latency difference
figure
df = trueStartTimeEMG - trueStartTime;
histogram(df, -0.5:0.05:0.5, Normalization='probability')
xlabel('tOnset(EMG) - tOnset(Video) (s)')
ylabel('Prob')
title(sprintf('EMG onset preceeds video-onset by (\\mu=%.1f, \\sigma=%.1f) ms', mean(df, 'omitnan')*1000, std(df, 'omitnan')*1000))


%% Trial-average EMG aligned to vid onset
ax = axes(figure);
hold(ax, 'on')
plot(ax, emg.vidOnsetAligned.t, emg.vidOnsetAligned.normX(1:10, :)')
plot(ax, emg.vidOnsetAligned.t, mean(emg.vidOnsetAligned.normX, 1, 'omitnan'), 'k', LineWidth=3)
hold on
ylim(ax, 'auto')
xlim(ax, 'auto')
plot([0, 0], ax.YLim, 'k--')
plot(ax.XLim, [2.5, 2.5], 'k--')
xlabel(ax, 'Time to touch (s)')
ylabel(ax, 'Smoothed normalzied EMG')
title(ax, 'Trial-average EMG aligned to vid-onset')

%% Do some cherry picking: find trials without spurious pre-reach movements, compare to those with
close all
hThreshold = 0.04;
nTrials = size(emg.vidOnsetAligned.X, 1);
tThreshold = -0.5;
hMax = max(emg.vidOnsetAligned.X(:, emg.vidOnsetAligned.t <= tThreshold), [], 2, 'omitnan');
selTrials = hMax <= hThreshold;

ax = axes(figure); hold(ax, 'on');
plot(ax, emg.vidOnsetAligned.t, mean(emg.vidOnsetAligned.X(selTrials, :), 1, 'omitnan'), DisplayName=sprintf('Good trials (n=%i)', nnz(selTrials)))
plot(ax, emg.vidOnsetAligned.t, mean(emg.vidOnsetAligned.X(~selTrials, :), 1, 'omitnan'), DisplayName=sprintf('Naughty trials (n=%i)', nnz(~selTrials)))
hold(ax, 'off')
legend(ax)

ax = axes(figure); hold(ax, 'on')
histogram(ax, trueStartTimeEMG(selTrials) - trueStartTime(selTrials), -5:0.1:0)
histogram(ax, trueStartTimeEMG(selTrials), -5:0.1:0)
histogram(ax, trueStartTimeEMG(~selTrials), -5:0.1:0)
title(ax, sprintf('mean=%g s, %g s, %g s', mean(trueStartTimeEMG(selTrials) - trueStartTime(selTrials)), mean(trueStartTimeEMG(selTrials)), mean(trueStartTimeEMG(~selTrials))))

ax = axes(figure); hold(ax, 'on');
plot(ax, emg.emgOnsetAligned.t, mean(emg.emgOnsetAligned.X(selTrials, :), 1, 'omitnan'), DisplayName=sprintf('Good trials (n=%i)', nnz(selTrials)))
plot(ax, emg.emgOnsetAligned.t, mean(emg.emgOnsetAligned.X(~selTrials, :), 1, 'omitnan'), DisplayName=sprintf('Naughty trials (n=%i)', nnz(~selTrials)))
hold(ax, 'off')
legend(ax)

%% Load eu and make PETH
% expName = 'daisy24_20240626';
% files = dir(sprintf('C:\\SERVER\\Units\\acute_3cam_reach_direction_2tgts\\SingleUnits_NonDuplicate\\%s*.mat', expName));
% files = arrayfun(@(f) [f.folder, '\', f.name], files, 'UniformOutput', false);
% eu = EphysUnit.load(files);
%%
eta.good = eu.getETA('count', 'press_spontaneous', window=[-4, 0.5], resolution=0.1, normalize='iti', alignTo='stop', includeInvalid=true, ...
            trials=goodTrials(selTrials), correction=trueStartTimeEMG(selTrials));
eta.bad = eu.getETA('count', 'press_spontaneous', window=[-4, 0.5], resolution=0.1, normalize='iti', alignTo='stop', includeInvalid=true, ...
            trials=goodTrials(~selTrials), correction=trueStartTimeEMG(~selTrials));
EphysUnit.plotDoubleETA(eta.good, eta.bad, clim=[-1.5, 1.5], xlim=[-4, 0], sortWindow=[-4, 0])

