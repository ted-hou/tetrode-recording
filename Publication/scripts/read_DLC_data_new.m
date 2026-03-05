if ~exist('euArtiFree', 'var')
    if exist('E:\Data\Units\PressVsLick_ArtifactsRemoved_Full\FixedEventsAndTrials', 'dir')
        euArtiFree = EphysUnit.load('E:\Data\Units\PressVsLick_ArtifactsRemoved_Full\FixedEventsAndTrials');
        load('E:\Data\Units\meta_PressVsLick_ArtifactsRemoved_Full_20260107.mat');
        etaArtiFree = metaArtiFree.eta;
    else
        euArtiFree = EphysUnit.load('C:\SERVER\Units\PressVsLick_ArtifactsRemoved_Full\FixedEventsAndTrials');
        load('C:\SERVER\Units\meta_PressVsLick_ArtifactsRemoved_Full_20260107.mat');
        etaArtiFree = metaArtiFree.eta;
    end
end
%% Daisy 2, 3, 8, 9, 10, 13, 14, 15, desmond10, 11, 22, 23, 24, 25, 26, 27
% Load DLC from 
% paths = [ ...
%     "\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\DeepLabCut\Results", ...
%     "\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\DeepLabCut\Projects\July14", ...
%     "\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\DeepLabCut\Projects\July3-DT-2025-07-03" ...
%     ];

paths = [ ...
    "\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\DeepLabCut\Results" ...
    % "\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\DeepLabCut\Projects\July14", ...
    % "\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\DeepLabCut\Projects\July3-DT-2025-07-03" ...
    ];

clear sessions
sessionNames = unique(string({euArtiFree.ExpName}));
nFound = 0;
for iSession = 1:length(sessionNames)
    sessions(iSession).name = sessionNames(iSession);
    for iPath = 1:length(paths)
        files = dir(sprintf("%s\\%s*.csv", paths(iPath), sessionNames(iSession)));
        if ~isempty(files)
            sessions(iSession).files = files;
            sessions(iSession).path = paths(iPath);
            sessions(iSession).eu = euArtiFree(ismember(string({euArtiFree.ExpName}), sessionNames(iSession)));
            break
        end
    end
end
clear iSession iPath files

%% Make CompleteExperiment3 objects
clc
exp = CompleteExperiment3([sessions.eu], cameras='lr', deeplabcutPath='\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\DeepLabCut\Results');

exp.alignTimestamps(refEventNameArduino={'CUE_ON'}, refEventNameEphys={'Cue'}, trialDurationTolerance=2);


%%
clear results
results(length(exp)) = struct(name=[], varsL=[], hasTimestampsL=[], varsR=[], hasTimestampsR=[]);
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
end
clear iExp


% Rename variables
from = ["HandCameraSide", "handIpsi", "handCont", "footIpsi", "footCont", "tongue"];
toL = ["HandL", "HandL", "HandR", "FootL", "FootR", "Tongue"];
toR = ["HandR", "HandR", "HandL", "FootR", "FootL", "Tongue"];
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

clear results
results(length(exp)) = struct(name=[], varsL=[], hasTimestampsL=[], varsR=[], hasTimestampsR=[]);
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
end
clear iExp toL toR vtd i from 

%% Remove bad sessions
sel = [results.hasTimestampsL] | [results.hasTimestampsR];
exp = exp(sel);
results = results(sel);
eu = [exp.eu];

clear sel

%% Get onset of reaching, offset of retraction
clear sequence trials
for iExp = 1:length(exp)
    trials(iExp).CueToLeverReleaseCorrect = exp(iExp).eu(1).Trials.CueToLeverReleaseCorrect;
    trials(iExp).PressCorrect = exp(iExp).eu(1).Trials.PressCorrect;
    % theseTrials.CorrectPressToFirstRewardLick
    % theseTrials.CorrectPressToLastLickOff
    
    [lia, locb] = ismember([trials(iExp).PressCorrect.Start], [trials(iExp).CueToLeverReleaseCorrect.Start]);
    if ~all(lia)
        trials(iExp).PressCorrect = trials(iExp).PressCorrect(lia);
    end
    sequence(iExp).Cue = [trials(iExp).PressCorrect.Start];
    sequence(iExp).ReachEnd = [trials(iExp).PressCorrect.Stop];
    sequence(iExp).RetractStart = [trials(iExp).CueToLeverReleaseCorrect.Stop];
end
clear iExp lia locb

%% Get reach start
for iExp = 1:length(exp)
    [X, Y, L, t] = exp(iExp).getTrajectoryByTrial('r', 'HandR', trialType='press', trials=trials(iExp).PressCorrect, alignTo='stop', window=[-4, 0], includeInvalid=true, likelihoodThreshold=0.5);
    selBaseline = t <= -2 & t >= -4;
    X = (X - mean(X(:, selBaseline), 2, 'omitnan')) ./ std(X(:, selBaseline), 0, 2, 'omitnan');
    Y = (Y - mean(Y(:, selBaseline), 2, 'omitnan')) ./ std(Y(:, selBaseline), 0, 2, 'omitnan');
    theta = 3;
    B = abs(X) >= theta | abs(Y) >= theta;
    tOnset = NaN(1, size(B, 1));
    for iTrial = 1:size(B, 1)
        iOnset = strfind(B(iTrial, :), [0, 1, 1]) + 1;
        if isempty(iOnset)
            tOnset(1, iTrial) = NaN;
        else
            iOnset = iOnset(end);
            tOnset(1, iTrial) = t(iOnset);
        end
    end
    sequence(iExp).ReachStart = tOnset + sequence(iExp).ReachEnd;
end
clear iExp X Y L t selBaseline theta B tOnset iTrial iOnset

% Get retract end
for iExp = 1:length(exp)
    [X, Y, L, t] = exp(iExp).getTrajectoryByTrial('r', 'HandR', trialType='press', trials=trials(iExp).CueToLeverReleaseCorrect, alignTo='stop', window=[0, 4], includeInvalid=true, likelihoodThreshold=0.5);
    selBaseline = t <= 4 & t >= 2;
    X = (X - mean(X(:, selBaseline), 2, 'omitnan')) ./ std(X(:, selBaseline), 0, 2, 'omitnan');
    Y = (Y - mean(Y(:, selBaseline), 2, 'omitnan')) ./ std(Y(:, selBaseline), 0, 2, 'omitnan');
    theta = 3;
    B = abs(X) >= theta | abs(Y) >= theta;
    tOffset = NaN(1, size(B, 1));
    for iTrial = 1:size(B, 1)
        iOffset = strfind(B(iTrial, :), [1, 0, 0]) + 1;
        if isempty(iOffset)
            tOffset(1, iTrial) = NaN;
        else
            iOffset = iOffset(end);
            tOffset(1, iTrial) = t(iOffset);
        end
    end
    sequence(iExp).RetractEnd = tOffset + sequence(iExp).RetractStart;
end
clear iExp X Y L t selBaseline theta B tOffset iTrial iOffset

%% Reshape and select common elements
for iExp = 1:length(exp)
    sel = true(length(sequence(iExp).Cue), 1);
    for fn = ["Cue", "ReachStart", "ReachEnd", "RetractStart", "RetractEnd"]
        sel = sel & ~isnan(sequence(iExp).(fn)');
    end
    sequence(iExp).isValid = sel;
    for fn = ["Cue", "ReachStart", "ReachEnd", "RetractStart", "RetractEnd"]
        sequence(iExp).(fn) = sequence(iExp).(fn)(sel)';
    end
end
clear iExp sel fn

%% Filter out invalid sessions
minNumTrials = 5;
sel = cellfun(@(b) nnz(b) >= minNumTrials, {sequence.isValid});
sequence = sequence(sel);
results = results(sel);
exp = exp(sel);
eu = [exp.eu];

%% Get binned spike counts and trajectories
[lia, expIndex] = ismember({eu.ExpName}, {exp.name});
assert(all(lia));
clear lia sc sr msr

nSamples = 10;
nSamplesBaseline = 20;
durationPreMove = 2;

% shuffle = [-100, 100];
shuffle = [0, 0];

tLocal = unique([linspace(-2, -1, nSamples), linspace(-1, 0, nSamples), linspace(0, 1, nSamples), linspace(1, 2, nSamples), linspace(2, 3, nSamples)]);
msr = NaN(length(eu), length(tLocal)-1);
for iEu = 1:length(eu)
    iExp = expIndex(iEu);
    nTrials = length(sequence(iExp).Cue);
    sr = NaN(nTrials, length(tLocal)-1);
    for iTrial = 1:nTrials
        shuffledTimeShift = shuffle(1) + (rand(1) * diff(shuffle));
        edges = unique([linspace(sequence(iExp).ReachStart(iTrial)-durationPreMove, sequence(iExp).ReachStart(iTrial), nSamples), linspace(sequence(iExp).ReachStart(iTrial), sequence(iExp).ReachEnd(iTrial), nSamples), linspace(sequence(iExp).ReachEnd(iTrial), sequence(iExp).RetractStart(iTrial), nSamples), linspace(sequence(iExp).RetractStart(iTrial), sequence(iExp).RetractEnd(iTrial), nSamples), linspace(sequence(iExp).RetractEnd(iTrial), sequence(iExp).RetractEnd(iTrial)+2, nSamples)]) + shuffledTimeShift;
        [sc, t] = eu(iEu).getSpikeCounts(edges);

        sc = double(sc)./diff(edges);
        baselineEdges = linspace(sequence(iExp).ReachStart(iTrial)-4, sequence(iExp).ReachStart(iTrial)-2, nSamplesBaseline);
        [bsc, ~] = eu(iEu).getSpikeCounts(baselineEdges);
        bsc = double(bsc)./diff(baselineEdges);

        sr(iTrial, :) = (sc - mean(bsc, 'all')) ./ std(bsc, 0, 'all');
    end
    msr(iEu, :) = mean(sr, 1);
end
clear iEu iExp nTrials sr iTrial edges sc t baselineEdges bsc


% Get trajectories per experiment
clear traj
traj(length(exp)) = struct(X=[], Y=[], L=[], t=[], T=[]);
for iExp = 1:length(exp)
    nTrials = length(sequence(iExp).Cue);
    X = zeros(nTrials, length(tLocal));
    Y = X;
    L = X;
    T = X;
    XBase = zeros(nTrials, nSamplesBaseline);
    YBase = XBase;

    for iTrial = 1:nTrials
        t = unique([linspace(sequence(iExp).ReachStart(iTrial)-durationPreMove, sequence(iExp).ReachStart(iTrial), nSamples), linspace(sequence(iExp).ReachStart(iTrial), sequence(iExp).ReachEnd(iTrial), nSamples), linspace(sequence(iExp).ReachEnd(iTrial), sequence(iExp).RetractStart(iTrial), nSamples), linspace(sequence(iExp).RetractStart(iTrial), sequence(iExp).RetractEnd(iTrial), nSamples), linspace(sequence(iExp).RetractEnd(iTrial), sequence(iExp).RetractEnd(iTrial)+2, nSamples)]);
        [X(iTrial, :), Y(iTrial, :), L(iTrial, :)] = exp(iExp).getTrajectory(t, 'r', 'HandR', likelihoodThreshold=0.5);
        T(iTrial, :) = t;

        % Baseline [-4, -2]
        tBaseline = linspace(sequence(iExp).ReachStart(iTrial)-4, sequence(iExp).ReachStart(iTrial)-2, nSamplesBaseline);
        [XBase(iTrial, :), YBase(iTrial, :), ~] =  exp(iExp).getTrajectory(tBaseline, 'r', 'HandR', likelihoodThreshold=0.5);
    end
    traj(iExp).X = (X - mean(XBase, 'all', 'omitnan')) ./ std(XBase, 0, 'all', 'omitnan');
    traj(iExp).Y = (Y - mean(YBase, 'all', 'omitnan')) ./ std(YBase, 0, 'all', 'omitnan');
    traj(iExp).L = L;
    traj(iExp).t = tLocal;
    traj(iExp).T = T;
end
clear iExp nTrials X Y L T XBase YBase iTrial t tBaseline

%% Plot PETH
% close all
fig = figure;
layout.h = [1, 3];
tl = tiledlayout(fig, sum(layout.h), 1, TileSpacing='compact', Padding='loose');
AX = gobjects(2, 1);
ax = nexttile(tl, [layout.h(1), 1]);
AX(1) = ax;
hold(ax, 'on')
for iExp = 1:length(exp)
    x = mean(traj(iExp).X, 1, 'omitnan');
    y = mean(traj(iExp).Y, 1, 'omitnan');
    plot(ax, traj(iExp).t, sqrt(x.^2 + y.^2), Color=[0.15, 0.15, 0.15, 0.3])
end

ax = nexttile(tl, layout.h(1) + 1, [layout.h(2), 1]);
AX(2) = ax;
hold(ax, 'on')
for iEu = 1:length(eu)
    % plot(ax, (tLocal(1:end-1) + tLocal(2:end))*0.5, (msr(iEu, :) - eu(iEu).SpikeRateStats.median)./eu(iEu).SpikeRateStats.madITI, Color=[0.15, 0.15, 0.15, 0.33])
    plot(ax, (tLocal(1:end-1) + tLocal(2:end))*0.5, msr(iEu, :))
end

for iAx = 1:2
    xline(AX(iAx), -1, 'k-', 'reach start')
    xline(AX(iAx), 0, 'k-', 'reach end')
    xline(AX(iAx), 1, 'k-', 'retract start')
    xline(AX(iAx), 2, 'k-', 'retract end')
end
ylabel(ax, 'normalized spike rate (a.u.)')
xlabel(ax, '"phase"')
xticks(ax, [])
title(ax, sprintf('n=%i SNr units, rewarded reach trials', length(eu)))
clear fig tl AX ax iExp x y iEu iAx


% Plot heatmap
fig = figure;
ax = axes(fig);
eta = struct(X=msr, t=(tLocal(1:end-1) + tLocal(2:end))*0.5, N=[]);

% EphysUnit.plotETA(ax, eta, signWindow=[-0.5, 0], sortWindow=[-2, 0], sortThreshold=0.25, negativeSortThreshold=0.25);

% [~, IValley] = min(eta.X(:, eta.t < -1), [], 2, 'omitnan');
[~, IValley] = min(eta.X, [], 2, 'omitnan');
[~, order] = sort(IValley);
EphysUnit.plotETA(ax, eta, order=order);

applyCustomColormap(ax, [-3, 3], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);    

xline(ax, -1, 'k-', 'reach start')
xline(ax, 0, 'k-', 'reach end')
xline(ax, 1, 'k-', 'retract start')
xline(ax, 2, 'k-', 'retract end')



%% Bootstrap the blue diagonal
% H0: the min spike rate for each neuron/PETH trace can be generated by
% chance (if we shuffle the timeseries)

clear bootMinSR
nBoot = 100;
bootMinSR = struct(obs=[], boot=NaN(length(eu), nBoot), p=[], h=[], ci95=[], ci99=[]);
tLocal = unique([linspace(-2, -1, nSamples), linspace(-1, 0, nSamples), linspace(0, 1, nSamples), linspace(1, 2, nSamples), linspace(2, 3, nSamples)]);


% 1. Do shuffled minSR

msr = NaN(length(eu), length(tLocal)-1);
shuffle = [-100, 100];
lineLength = 0;
ticTotal = tic();
for iBoot = 1:nBoot
    ticBoot = tic();
    for iEu = 1:length(eu)
        fprintf(repmat('\b', [1, lineLength]));
        lineLength = fprintf("iBoot=%i/%i, iEu=%i/%i (%.2fs, %.2fs total)\n", iBoot, nBoot, iEu, length(eu), toc(ticBoot), toc(ticTotal));
        iExp = expIndex(iEu);
        nTrials = length(sequence(iExp).Cue);
        sr = NaN(nTrials, length(tLocal)-1);
        for iTrial = 1:nTrials
            shuffledTimeShift = shuffle(1) + (rand(1) * diff(shuffle));
            edges = unique([linspace(sequence(iExp).ReachStart(iTrial)-durationPreMove, sequence(iExp).ReachStart(iTrial), nSamples), linspace(sequence(iExp).ReachStart(iTrial), sequence(iExp).ReachEnd(iTrial), nSamples), linspace(sequence(iExp).ReachEnd(iTrial), sequence(iExp).RetractStart(iTrial), nSamples), linspace(sequence(iExp).RetractStart(iTrial), sequence(iExp).RetractEnd(iTrial), nSamples), linspace(sequence(iExp).RetractEnd(iTrial), sequence(iExp).RetractEnd(iTrial)+2, nSamples)]) + shuffledTimeShift;
            [sc, ~] = eu(iEu).getSpikeCounts(edges);
    
            sc = double(sc)./diff(edges);
            baselineEdges = linspace(sequence(iExp).ReachStart(iTrial)-4, sequence(iExp).ReachStart(iTrial)-2, nSamplesBaseline);
            [bsc, ~] = eu(iEu).getSpikeCounts(baselineEdges);
            bsc = double(bsc)./diff(baselineEdges);
    
            sr(iTrial, :) = (sc - mean(bsc, 'all')) ./ std(bsc, 0, 'all');
        end
        msr(iEu, :) = mean(sr, 1);
    end
    clear iEu iExp nTrials sr iTrial edges sc t baselineEdges bsc
    [bootMinSR.boot(:, iBoot), ~] = min(msr, [], 2, 'omitnan');
end

% 2. Do observed minSR
shuffle = [0, 0];

msr = NaN(length(eu), length(tLocal)-1);
for iEu = 1:length(eu)
    iExp = expIndex(iEu);
    nTrials = length(sequence(iExp).Cue);
    sr = NaN(nTrials, length(tLocal)-1);
    for iTrial = 1:nTrials
        shuffledTimeShift = shuffle(1) + (rand(1) * diff(shuffle));
        edges = unique([linspace(sequence(iExp).ReachStart(iTrial)-durationPreMove, sequence(iExp).ReachStart(iTrial), nSamples), linspace(sequence(iExp).ReachStart(iTrial), sequence(iExp).ReachEnd(iTrial), nSamples), linspace(sequence(iExp).ReachEnd(iTrial), sequence(iExp).RetractStart(iTrial), nSamples), linspace(sequence(iExp).RetractStart(iTrial), sequence(iExp).RetractEnd(iTrial), nSamples), linspace(sequence(iExp).RetractEnd(iTrial), sequence(iExp).RetractEnd(iTrial)+2, nSamples)]) + shuffledTimeShift;
        [sc, t] = eu(iEu).getSpikeCounts(edges);

        sc = double(sc)./diff(edges);
        baselineEdges = linspace(sequence(iExp).ReachStart(iTrial)-4, sequence(iExp).ReachStart(iTrial)-2, nSamplesBaseline);
        [bsc, ~] = eu(iEu).getSpikeCounts(baselineEdges);
        bsc = double(bsc)./diff(baselineEdges);

        sr(iTrial, :) = (sc - mean(bsc, 'all')) ./ std(bsc, 0, 'all');
    end
    msr(iEu, :) = mean(sr, 1);
end
clear iEu iExp nTrials sr iTrial edges sc t baselineEdges bsc
[bootMinSR.obs, ~] = min(msr, [], 2, 'omitnan');
clear iBoot ticBoot ticTotal lineLength shuffle

%%
bootMinSR.ci95=quantile(bootMinSR.boot, [0.05, 1], 2);
bootMinSR.ci99=quantile(bootMinSR.boot, [0.01, 1], 2);

for iEu = 1:length(eu)
    bootMinSR.h95(iEu, 1) = bootMinSR.obs(iEu) < bootMinSR.ci95(1);
    bootMinSR.h99(iEu, 1) = bootMinSR.obs(iEu) < bootMinSR.ci99(1);
end
clear iEu


% Plot heatmap
fig = figure;
ax = axes(fig);
eta = struct(X=msr(bootMinSR.h95, :), t=(tLocal(1:end-1) + tLocal(2:end))*0.5, N=[]);

% EphysUnit.plotETA(ax, eta, signWindow=[-0.5, 0], sortWindow=[-2, 0], sortThreshold=0.25, negativeSortThreshold=0.25);

[~, IValley] = min(eta.X, [], 2, 'omitnan');
[~, order] = sort(IValley);
EphysUnit.plotETA(ax, eta, order=order);

applyCustomColormap(ax, [-3, 3], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);    

xline(ax, -1, 'k-', 'reach start')
xline(ax, 0, 'k-', 'reach end')
xline(ax, 1, 'k-', 'retract start')
xline(ax, 2, 'k-', 'retract end')
