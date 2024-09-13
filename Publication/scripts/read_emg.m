function [exp, onset, isGoodTrial, etaMerged, onsetThreshold, spuriousMovementThreshold] = read_emg()

onsetThreshold = 0.25;
expName = {'daisy24_20240626', 'daisy24_20240627', 'daisy24_20240628'};

% Load data
exp(length(expName)) = struct(ai=[], pa=[], eu=[], goodTrials=[]);
for iExp = 1:length(exp)
    ai = load(sprintf('C:\\SERVER\\PawAnalysis\\EMG_%s.mat', expName{iExp}));
    exp(iExp).ai = ai.ai; 
    clear ai;
    files = dir(sprintf('C:\\SERVER\\Units\\acute_3cam_reach_direction_2tgts\\SingleUnits_NonDuplicate\\%s_Channel*_Unit*.mat', expName{iExp}));
    assert(~isempty(files))
    files = arrayfun(@(f) sprintf('%s\\%s', f.folder, f.name), files, UniformOutput=false);
    
    exp(iExp).eu = EphysUnit.load(files);
    trials = exp(iExp).eu(1).makeTrials('press_spontaneous2');
    exp(iExp).goodTrials = trials(trials.duration() > 4);% & arrayfun(@(t) all(t <= eu.EventTimes.RewardTimes | t >= eu.EventTimes.RewardTimes + 4), [trials.Stop])) ;
    exp(iExp).eu(1).Trials.PressSpontaneous = exp(iExp).goodTrials;
    exp(iExp).pa = Pawnalyzer2(exp(iExp).eu, refEvent='press');
    exp(iExp).pa.getClips(noImage=true, nFramesBefore=15, nFramesAfter=0, keepData=false, trials='PressSpontaneous');
    exp(iExp).pa.load('auto')

    % Smooth EMG
    exp(iExp).ai.SmoothData = smoothdata(exp(iExp).ai.Data, 'movmean', 1000); % 1000 sample (33.3ms) moving average window
    exp(iExp).ai.SmoothNormData = (exp(iExp).ai.SmoothData - median(exp(iExp).ai.SmoothData))./mad(exp(iExp).ai.SmoothData)*0.6745;
end
clear iExp files trials

%% Calculate paw trajectories
nt = 10;
nFrames = 16;
for iExp = 1:length(exp)
    clear traj
    traj = struct(contra=[], ipsi=[], target=[], t=[], pca=[], kmeans=[]);
    [traj.contra, traj.ipsi, traj.target, traj.t] = exp(iExp).pa.getTrajectories(1, zero=nFrames - nt + 1);
    if ~ismember(exp(iExp).pa.exp.animalName, {'desmond29', 'daisy25'})
        traj.contra.x = -traj.contra.x;
        traj.ipsi.x = -traj.ipsi.x;
    end
    exp(iExp).traj = traj;
end
clear iExp traj
%% Calculate video-movement onset times
for iExp = 1:length(exp)
    nTrials = length(exp(iExp).traj.target);
    exp(iExp).trueStartTime = NaN(nTrials, 1);
    t = exp(iExp).traj.t;
    dist = zeros(nTrials, length(t));
    for iTrial = 1:nTrials
        pos = [ ...
            exp(iExp).traj.contra.x(iTrial, :); ...
            exp(iExp).traj.contra.y(iTrial, :); ...
            exp(iExp).traj.contra.z(iTrial, :); ...            
            ];
        dist(iTrial, :) = sqrt(sum(pos.^2, 1));
    end
    distNorm = (dist - median(dist(:), 'omitnan'))./(mad(dist(:), 1)./0.6745);
    for iTrial = 1:nTrials
        if any(isnan(distNorm(iTrial, :)))
            continue
        end
        isAbove = distNorm(iTrial, :) >= onsetThreshold;
        iLastAbove = find(isAbove, 1, 'last'); % We don't do abs since dist is already positive, although z-scoring will generate negatives, true movement should be positive.
        if isempty(iLastAbove)
            continue
        end
        iOnset = strfind(isAbove(1:iLastAbove), [0 1 1]) + 1;
        if isempty(iOnset)
            continue
        end
        iOnset = iOnset(end);
        exp(iExp).trueStartTime(iTrial) = t(iOnset) - t(end);
    end
end

for iExp = 1:length(exp)
    nTrials = length(exp(iExp).traj.target);
    exp(iExp).trueStartTimeLeft = NaN(nTrials, 1);
    t = exp(iExp).traj.t;
    dist = zeros(nTrials, length(t));
    for iTrial = 1:nTrials
        pos = [ ...
            exp(iExp).traj.contra.x(iTrial, :); ...
            exp(iExp).traj.contra.y(iTrial, :); ...
            exp(iExp).traj.contra.z(iTrial, :); ...            
            ];
        dist(iTrial, :) = sqrt(sum(pos.^2, 1));
    end
    distNorm = (dist - median(dist(:), 'omitnan'))./(mad(dist(:), 1)./0.6745);
    for iTrial = 1:nTrials
        if any(isnan(distNorm(iTrial, :)))
            continue
        end
        isAbove = distNorm(iTrial, :) >= onsetThreshold;
        iOnset = strfind(isAbove, [0 1 1]) + 1;
        if isempty(iOnset)
            continue
        end
        iOnset = iOnset(1);
        exp(iExp).trueStartTimeLeft(iTrial) = t(iOnset) - t(end);
    end
end



figure
histogram(vertcat(exp.trueStartTime), -0.5:0.05:0, Normalization='probability')
hold on

histogram(vertcat(exp.trueStartTimeLeft), -0.5:0.05:0, Normalization='probability')
xlabel('Time to touch (s)')
ylabel('Probability')
title(sprintf('Video detected movement onset time (when z(dist) > %g for 2 conseq bins)', onsetThreshold))
%% Calculate touch aligned EMG
baseWindow = [-4, -2];
for iExp = 1:length(exp)
    tLocal = -4:1/30000:0;
    tEvent = reshape([exp(iExp).eu(1).Trials.PressSpontaneous.Stop], [], 1);
    tGlobal = tLocal + tEvent;
    nTrials = length(tEvent);
    X = zeros(nTrials, length(tLocal));
    for iTrial = 1:nTrials
        X(iTrial, :) = interp1(exp(iExp).ai.Timestamps, exp(iExp).ai.SmoothData, tGlobal(iTrial, :), 'linear');
    end
    
    selBase = tLocal >= baseWindow(1) & tLocal <= baseWindow(2);
    baseMeanByTrial = mean(X(:, selBase), 2, 'omitnan');
    baseSdByTrial = std(X(:, selBase), 1, 2, 'omitnan');
    normX = (X - baseMeanByTrial)./baseSdByTrial;
    
    exp(iExp).emg.touchAligned.t = tLocal;
    exp(iExp).emg.touchAligned.T = tGlobal;
    exp(iExp).emg.touchAligned.X = X;
    exp(iExp).emg.touchAligned.normX = normX;
end
clear tLocal tEvent tGlobal nTrials X iTrial emg iExp selBase baseMeanByTrial baseSdByTrial normX

%% Calculate video-movement-onset aligned EMG
for iExp = 1:length(exp)
    tLocal = -4:1/30000:0.5;
    tEvent = reshape([exp(iExp).eu(1).Trials.PressSpontaneous.Stop], [], 1) + reshape(exp(iExp).trueStartTime, [], 1);
    tGlobal = tLocal + tEvent;
    nTrials = length(tEvent);
    X = zeros(nTrials, length(tLocal));
    for iTrial = 1:nTrials
        X(iTrial, :) = interp1(exp(iExp).ai.Timestamps, exp(iExp).ai.SmoothData, tGlobal(iTrial, :), 'linear');
    end
    
    selBase = tLocal >= baseWindow(1) & tLocal <= baseWindow(2);
    baseMeanByTrial = mean(X(:, selBase), 2, 'omitnan');
    baseSdByTrial = std(X(:, selBase), 1, 2, 'omitnan');
    normX = (X - baseMeanByTrial)./baseSdByTrial;
    
    exp(iExp).emg.vidOnsetAligned.t = tLocal;
    exp(iExp).emg.vidOnsetAligned.T = tGlobal;
    exp(iExp).emg.vidOnsetAligned.X = X;
    exp(iExp).emg.vidOnsetAligned.normX = normX;
end
clear tLocal tEvent tGlobal nTrials X iTrial iExp selBase baseMeanByTrial baseSdByTrial normX
%% PLOT Trial-average EMG aligned to touch
% figure
% normX = arrayfun(@(exp) exp.emg.touchAligned.normX, expEMG, UniformOutput=false);
% normX = cat(1, normX{:});
% plot(expEMG(1).emg.touchAligned.t, mean(normX, 1, 'omitnan'))
% xlabel('Time to touch (s)')
% ylabel('Smoothed normalzied EMG')
% title(sprintf('Trial-average EMG aligned to touch (when z(dist) > %g for 2 conseq bins)', size(normX, 1), length(expEMG)))

%% Detect movement onset time based on z(EMG)>0.25 for 2 conseq bins (from right)
for iExp = 1:length(exp)
    nTrials = length(exp(iExp).traj.target);
    trueStartTimeEMG = NaN(nTrials, 1);
    t = -2:1/30:0;
    for iTrial = 1:nTrials
        x = interp1(exp(iExp).emg.touchAligned.t, exp(iExp).emg.touchAligned.normX(iTrial, :), t);
        isAbove = x >= onsetThreshold;
        iLastAbove = find(isAbove, 1, 'last'); % We don't do abs since dist is already positive, although z-scoring will generate negatives, true movement should be positive.
        if isempty(iLastAbove)
            continue
        end
        iOnset = strfind(isAbove(1:iLastAbove), [0 1 1]) + 1;
        if isempty(iOnset)
            continue
        end
        iOnset = iOnset(end);
        trueStartTimeEMG(iTrial) = t(iOnset) - t(end);
    end
    exp(iExp).trueStartTimeEMG = trueStartTimeEMG;
end

% Detect video and EMG onset from the left
for iExp = 1:length(exp)
    nTrials = length(exp(iExp).traj.target);
    exp(iExp).trueStartTimeEMGLeft = NaN(nTrials, 1);
    t = -2:1/30:0;
    for iTrial = 1:nTrials
        x = interp1(exp(iExp).emg.touchAligned.t, exp(iExp).emg.touchAligned.normX(iTrial, :), t);
        isAbove = x >= 1;
        iOnset = strfind(isAbove, [0 1 1]) + 1;
        if isempty(iOnset)
            continue
        end
        iOnset = iOnset(1);
        exp(iExp).trueStartTimeEMGLeft(iTrial) = t(iOnset) - t(end);
    end
end

% figure
% histogram(vertcat(expEMG.trueStartTimeEMG), -2:0.05:0, Normalization='probability')
% hold on
% histogram(vertcat(expEMG.trueStartTimeEMGLeft), -2:0.05:0, Normalization='probability')
% hold on
% t = vertcat(expEMG.trueStartTimeEMGLeft);
% histogram(t(selTrials), -2:0.05:0, Normalization='probability')
% xlabel('Time to touch (s)')
% ylabel('Probability')
% title(sprintf('EMG detected movement onset time (when z(EMG) > %g)', onsetThreshold))


%% Plot EMG-Video latency difference
% figure
% df = trueStartTimeEMG - trueStartTime;
% histogram(df, -0.5:0.05:0.5, Normalization='probability')
% xlabel('tOnset(EMG) - tOnset(Video) (s)')
% ylabel('Prob')
% title(sprintf('EMG onset preceeds video-onset by (\\mu=%.1f, \\sigma=%.1f) ms', mean(df, 'omitnan')*1000, std(df, 'omitnan')*1000))


%% PLOT Trial-average EMG aligned to vid onset

% normX = arrayfun(@(exp) exp.emg.vidOnsetAligned.normX, expEMG, UniformOutput=false);
% normX = cat(1, normX{:});
% 
% ax = axes(figure);
% hold(ax, 'on')
% % plot(ax, exp(1).emg.vidOnsetAligned.t, normX(1:10, :)')
% plot(ax, expEMG(1).emg.vidOnsetAligned.t, mean(normX, 1, 'omitnan'), 'k', LineWidth=3)
% hold on
% ylim(ax, 'auto')
% xlim(ax, 'auto')
% plot([0, 0], ax.YLim, 'k--')
% plot(ax.XLim, [0.25, 0.25], 'k--')
% xlabel(ax, 'Time to touch (s)')
% ylabel(ax, 'Smoothed normalzied EMG')
% title(ax, 'Trial-average EMG aligned to vid-onset')

%% Do some cherry picking: find trials without spurious pre-reach movements, compare to those with
% close all
spuriousMovementThreshold = 5;

normX = arrayfun(@(exp) exp.emg.touchAligned.normX, exp, UniformOutput=false);
normX = cat(1, normX{:});
t = exp(1).emg.touchAligned.t;

trueStartTimeEMG = vertcat(exp.trueStartTimeEMG);
trueStartTime = vertcat(exp.trueStartTime);


nTrials = size(normX, 1);
tThreshold = -0.5;
hMaxPre = max(normX(:, t <= tThreshold), [], 2, 'omitnan');
hMaxPost = max(normX(:, t > tThreshold), [], 2, 'omitnan');
isGoodTrial = hMaxPre <= spuriousMovementThreshold & hMaxPost >= spuriousMovementThreshold;

ax = axes(figure); hold(ax, 'on');
plot(ax, t, mean(normX(isGoodTrial, :), 1, 'omitnan'), DisplayName=sprintf('Good trials (n=%i)', nnz(isGoodTrial)))
plot(ax, t, mean(normX(~isGoodTrial, :), 1, 'omitnan'), DisplayName=sprintf('Naughty trials (n=%i)', nnz(~isGoodTrial)))
plot(ax, [-4, 0], [0.25, 0.25], 'k--', DisplayName=sprintf('Threshold=%g', onsetThreshold))
hold(ax, 'off')
legend(ax)
xlabel(ax, 'Time to touch (s)')
ylabel(ax, 'EMG')
ylim(ax, [-1, 15])


normX = arrayfun(@(exp) exp.emg.vidOnsetAligned.normX, exp, UniformOutput=false);
normX = cat(1, normX{:});
t = exp(1).emg.vidOnsetAligned.t;

ax = axes(figure); hold(ax, 'on');
plot(ax, t, mean(normX(isGoodTrial, :), 1, 'omitnan'), DisplayName=sprintf('Good trials (n=%i)', nnz(isGoodTrial)))
plot(ax, t, mean(normX(~isGoodTrial, :), 1, 'omitnan'), DisplayName=sprintf('Naughty trials (n=%i)', nnz(~isGoodTrial)))
plot(ax, [-4, 0], [0.25, 0.25], 'k--', DisplayName=sprintf('Threshold=%g', onsetThreshold))
hold(ax, 'off')
legend(ax)
xlabel(ax, 'Time to video onset (s)')
ylabel(ax, 'EMG')
xlim(ax, [-4, 0])
ylim(ax, [0-1, 15])

ax = axes(figure); hold(ax, 'on')
% histogram(ax, trueStartTimeEMG(selTrials) - trueStartTime(selTrials), -5:0.1:0)
histogram(ax, trueStartTimeEMG(isGoodTrial), -0.5:0.05:0)
% histogram(ax, trueStartTimeEMG(~selTrials), -5:0.1:0)
title(ax, sprintf('mean=%g s, %g s, %g s', mean(trueStartTimeEMG(isGoodTrial) - trueStartTime(isGoodTrial), 'omitnan'), mean(trueStartTimeEMG(isGoodTrial), 'omitnan'), mean(trueStartTimeEMG(~isGoodTrial), 'omitnan')))

%% PLOT ETAs side by side (good vs bad trials), aligned to EMG onset
% Good trials should still have early activation
eta(length(exp)) = struct(good=[], bad=[]);
csnTrials = cumsum([0, arrayfun(@(exp) length(exp.goodTrials), exp)]);
for iExp = 1:length(exp)
    selTrialInSession = csnTrials(iExp) + 1:csnTrials(iExp + 1);
    eta(iExp).good = exp(iExp).eu.getETA('count', 'press_spontaneous', window=[-4, 0.5], resolution=0.1, normalize=[-4, -2], alignTo='stop', includeInvalid=true, ...
                trials=exp(iExp).goodTrials(isGoodTrial(selTrialInSession)), correction=trueStartTimeEMG(isGoodTrial(selTrialInSession)));
    eta(iExp).bad = exp(iExp).eu.getETA('count', 'press_spontaneous', window=[-4, 0.5], resolution=0.1, normalize=[-4, -2], alignTo='stop', includeInvalid=true, ...
                trials=exp(iExp).goodTrials(~isGoodTrial(selTrialInSession)), correction=trueStartTimeEMG(~isGoodTrial(selTrialInSession)));
end

etaGood = [eta.good];
etaBad = [eta.bad];
etaMerged.good = struct(X=cat(1, etaGood.X), t=etaGood(1).t, N=cat(1, etaGood.N), D=cat(1, etaGood.D));
etaMerged.bad = struct(X=cat(1, etaBad.X), t=etaBad(1).t, N=cat(1, etaBad.N), D=cat(1, etaBad.D));
clear etaGood etaBad csnTrials iExp selTrialInSession

[ax, order] = EphysUnit.plotDoubleETA(etaMerged.good, etaMerged.bad, clim=[-1.5, 1.5], xlim=[-4, 0], sortWindow=[-4, 0], signWindow=[-0.3, 0], sortThreshold=0.25, negativeSortThreshold=0.25);
title(ax(1), sprintf('%i good trials', nnz(isGoodTrial)))
title(ax(2), sprintf('%i bad trials', nnz(~isGoodTrial)))
xlabel(ax, 'Time to EMG onset (s)')

onset.video = trueStartTime;
onset.emg = trueStartTimeEMG;


