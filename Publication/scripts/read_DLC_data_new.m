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

%% Reset
clearvars -except euArtiFree paths sessionNames sessions exp results

%% Collect sequence of events for rewarded reach trials
clear sequence trials
sequence(length(exp)) = struct(CorrectReach=[], CorrectLick=[]);
for iExp = 1:length(exp)
    trials(iExp).CueToLeverReleaseCorrect = exp(iExp).eu(1).Trials.CueToLeverReleaseCorrect;
    trials(iExp).PressCorrect = exp(iExp).eu(1).Trials.PressCorrect;
    trials(iExp).CorrectPressToFirstRewardLick = exp(iExp).eu(1).Trials.CorrectPressToFirstRewardLick;
    trials(iExp).LickCorrect = exp(iExp).eu(1).Trials.LickCorrect;
    trials(iExp).CueToLastLickOffCorrect = exp(iExp).eu(1).Trials.CueToLastLickOffCorrect;
    
    [lia, locb] = ismember([trials(iExp).PressCorrect.Start], [trials(iExp).CueToLeverReleaseCorrect.Start]);
    if ~all(lia)
        trials(iExp).PressCorrect = trials(iExp).PressCorrect(lia);
    end
    sequence(iExp).CorrectReach.Cue = [trials(iExp).PressCorrect.Start];
    sequence(iExp).CorrectReach.ReachEnd = [trials(iExp).PressCorrect.Stop];
    sequence(iExp).CorrectReach.RetractStart = [trials(iExp).CueToLeverReleaseCorrect.Stop];

    % First lick after bar-contact, before arm-retract
    lickOn = [trials(iExp).CorrectPressToFirstRewardLick.Stop];
    [lia, locb] = ismember([trials(iExp).CorrectPressToFirstRewardLick.Start], [trials(iExp).PressCorrect.Stop]);
    if ~all(lia)
        trials(iExp).CorrectPressToFirstRewardLick = trials(iExp).CorrectPressToFirstRewardLick(lia);
        lickOn = lickOn(lia);
        locb = locb(lia);
    end
    sequence(iExp).CorrectReach.FirstLick = NaN(1, length(trials(iExp).PressCorrect));
    sequence(iExp).CorrectReach.FirstLick(locb) = lickOn;

    % CorrectLick
    sequence(iExp).CorrectLick.Cue = [trials(iExp).LickCorrect.Start];
    sequence(iExp).CorrectLick.Lick = [trials(iExp).LickCorrect.Stop];
    sequence(iExp).CorrectLick.LastLickOff = [trials(iExp).CueToLastLickOffCorrect.Stop];
end
clear iExp lia locb trialsBarHeld B I lickOn

%% Get the onset of arm reach/offset of arm retraction via video
% Get ReachStart
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
    sequence(iExp).CorrectReach.ReachStart = tOnset + sequence(iExp).CorrectReach.ReachEnd;
end
clear iExp X Y L t selBaseline theta B tOnset iTrial iOnset

% Get RetractEnd
for iExp = 1:length(exp)
    [X, Y, L, t] = exp(iExp).getTrajectoryByTrial('r', 'HandR', trialType='press', trials=trials(iExp).CueToLeverReleaseCorrect, alignTo='stop', window=[0, 4], includeInvalid=true, likelihoodThreshold=0.5);
    selBaseline = t <= 4 & t >= 2;
    X = (X - mean(X(:, selBaseline), 2, 'omitnan')) ./ std(X(:, selBaseline), 0, 2, 'omitnan');
    Y = (Y - mean(Y(:, selBaseline), 2, 'omitnan')) ./ std(Y(:, selBaseline), 0, 2, 'omitnan');
    theta = 0.5;
    B = abs(X) <= theta & abs(Y) <= theta;
    tOffset = NaN(1, size(B, 1));
    for iTrial = 1:size(B, 1)
        iOffset = strfind(B(iTrial, :), [0, 1, 1]) + 1;
        if isempty(iOffset)
            tOffset(1, iTrial) = NaN;
        else
            iOffset = iOffset(1);
            tOffset(1, iTrial) = t(iOffset);
        end
    end
    sequence(iExp).CorrectReach.RetractEnd = tOffset + sequence(iExp).CorrectReach.RetractStart;
end
clear iExp X Y L t selBaseline theta B tOffset iTrial iOffset

%% Reshape and select common elements
for iExp = 1:length(exp)
    sel = true(length(sequence(iExp).CorrectReach.Cue), 1);
    for fn = ["Cue", "ReachStart", "ReachEnd", "FirstLick", "RetractStart", "RetractEnd"]
        sel = sel & ~isnan(sequence(iExp).CorrectReach.(fn)');
    end
    sequence(iExp).CorrectReach.isValid = sel;
    for fn = ["Cue", "ReachStart", "ReachEnd", "FirstLick", "RetractStart", "RetractEnd"]
        sequence(iExp).CorrectReach.(fn) = sequence(iExp).CorrectReach.(fn)(sel)';
    end
end
clear iExp sel fn

%% Filter out invalid sessions
minNumTrials = 5;
thisSeq = [sequence.CorrectReach];
sel = cellfun(@(b) nnz(b) >= minNumTrials, {thisSeq.isValid});
sequence = sequence(sel);
results = results(sel);
exp = exp(sel);
eu = [exp.eu];
clear sel

%% Get binned spike counts and trajectories
[lia, expIndex] = ismember({eu.ExpName}, {exp.name});
assert(all(lia));
clear lia sc sr msr

nSamples = 8;
nSamplesBaseline = 20;
durationPreMove = 2;

% shuffle = [-100, 100];
shuffle = [0, 0];

tLocal.CorrectReach = unique([ ...
    linspace(-2, -1, nSamples), ... ReachStart-2 -> ReachStart
    linspace(-1, 0, nSamples), ... ReachStart -> ReachEnd
    linspace(0, 1, nSamples), ... ReachEnd -> FirstLick
    linspace(1, 2, nSamples), ... FirstLick -> FirstLick+1
    linspace(2, 3, nSamples), ... FirstLick+1 -> RetractStart
    linspace(3, 4, nSamples), ... RetractStart -> RetractEnd
    linspace(4, 5, nSamples) ... RetractEnd -> RetractEnd+2
    ]);
tLocal.CorrectLick = unique([ ...
    linspace(-1, 0, nSamples), ... Lick-2 -> Lick
    linspace(0, 1, nSamples), ... Lick -> Lick+1
    linspace(1, 2, nSamples), ... Lick+1 -> LastLickOff
    linspace(2, 3, nSamples), ... LastLickOff -> LastLickOff+2
    ]);
msr.CorrectReach = NaN(length(eu), length(tLocal.CorrectReach)-1);
msr.CorrectLick = NaN(length(eu), length(tLocal.CorrectLick)-1);
for iEu = 1:length(eu)
    iExp = expIndex(iEu);

    % CorrectReach
    nTrials = length(sequence(iExp).CorrectReach.Cue);
    sr = NaN(nTrials, length(tLocal.CorrectReach)-1);
    for iTrial = 1:nTrials
        try
            shuffledTimeShift = shuffle(1) + (rand(1) * diff(shuffle));
            edges = [ ...
                linspace(sequence(iExp).CorrectReach.ReachStart(iTrial)-durationPreMove, sequence(iExp).CorrectReach.ReachStart(iTrial), nSamples), ...ReachStart-2 -> ReachStart
                linspace(sequence(iExp).CorrectReach.ReachStart(iTrial), sequence(iExp).CorrectReach.ReachEnd(iTrial), nSamples), ... ReachStart -> ReachEnd
                linspace(sequence(iExp).CorrectReach.ReachEnd(iTrial), sequence(iExp).CorrectReach.FirstLick(iTrial), nSamples), ... ReachEnd -> FirstLick
                linspace(sequence(iExp).CorrectReach.FirstLick(iTrial), min(sequence(iExp).CorrectReach.FirstLick(iTrial)+1, sequence(iExp).CorrectReach.RetractStart(iTrial)), nSamples), ... FirstLick -> FirstLick+1
                linspace(min(sequence(iExp).CorrectReach.FirstLick(iTrial)+1, sequence(iExp).CorrectReach.RetractStart(iTrial)), sequence(iExp).CorrectReach.RetractStart(iTrial), nSamples), ... FirstLick+1 -> RetractStart
                linspace(sequence(iExp).CorrectReach.RetractStart(iTrial), sequence(iExp).CorrectReach.RetractEnd(iTrial), nSamples), ... RetractStart -> RetractEnd
                linspace(sequence(iExp).CorrectReach.RetractEnd(iTrial), sequence(iExp).CorrectReach.RetractEnd(iTrial)+2, nSamples) ... RetractEnd -> RetractEnd+2
                ];
            if ~all(diff(edges) >= 0)
                error()
            end
            edges = unique(edges) + shuffledTimeShift;
            [sc, t] = eu(iEu).getSpikeCounts(edges);
    
            sc = double(sc)./diff(edges);
            baselineEdges = linspace(sequence(iExp).CorrectReach.ReachStart(iTrial)-4, sequence(iExp).CorrectReach.ReachStart(iTrial)-2, nSamplesBaseline);
            [bsc, ~] = eu(iEu).getSpikeCounts(baselineEdges);
            bsc = double(bsc)./diff(baselineEdges);
    
            sr(iTrial, :) = (sc - mean(bsc, 'all', 'omitnan')) ./ std(bsc, 0, 'all', 'omitnan');
        catch
            warning('Cannot process CorrectReach PETH for iExp=%i, iEu=%i, iTrial=%i, %s', iExp, iEu, iTrial, eu(iEu).getName)
        end
    end
    msr.CorrectReach(iEu, :) = mean(sr, 1, 'omitnan');

    % CorrectLick
    nTrials = length(sequence(iExp).CorrectLick.Cue);
    sr = NaN(nTrials, length(tLocal.CorrectLick)-1);
    for iTrial = 1:nTrials
        try
            shuffledTimeShift = shuffle(1) + (rand(1) * diff(shuffle));
            edges = [ ...
                    linspace(sequence(iExp).CorrectLick.Lick(iTrial)-durationPreMove, sequence(iExp).CorrectLick.Lick(iTrial), nSamples), ... Lick-2 -> Lick
                    linspace(sequence(iExp).CorrectLick.Lick(iTrial), min(sequence(iExp).CorrectLick.Lick(iTrial)+1, sequence(iExp).CorrectLick.LastLickOff(iTrial)), nSamples), ... Lick -> Lick+1
                    linspace(min(sequence(iExp).CorrectLick.Lick(iTrial)+1, sequence(iExp).CorrectLick.LastLickOff(iTrial)), sequence(iExp).CorrectLick.LastLickOff(iTrial), nSamples), ... Lick+1 -> LastLickOff
                    linspace(sequence(iExp).CorrectLick.LastLickOff(iTrial), sequence(iExp).CorrectLick.LastLickOff(iTrial)+durationPreMove, nSamples), ... LastLickOff -> LastLickOff+2
                ];
            if ~all(diff(edges) >= 0)
                error()
            end
            edges = unique(edges) + shuffledTimeShift;
            [sc, t] = eu(iEu).getSpikeCounts(edges);
    
            sc = double(sc)./diff(edges);
            baselineEdges = linspace(sequence(iExp).CorrectLick.Lick(iTrial)-4, sequence(iExp).CorrectLick.Lick(iTrial)-2, nSamplesBaseline);
            [bsc, ~] = eu(iEu).getSpikeCounts(baselineEdges);
            bsc = double(bsc)./diff(baselineEdges);
    
            sr(iTrial, :) = (sc - mean(bsc, 'all', 'omitnan')) ./ std(bsc, 0, 'all', 'omitnan');
        catch
            warning('Cannot process CorrectLick PETH for iExp=%i, iEu=%i, iTrial=%i, %s', iExp, iEu, iTrial, eu(iEu).getName)
        end
    end
    msr.CorrectLick(iEu, :) = mean(sr, 1, 'omitnan');


end
clear iEu iExp nTrials sr iTrial edges sc t baselineEdges bsc



%% Get trajectories per experiment
clear traj
traj(length(exp)) = struct(CorrectReach=[], CorrectLick=[]);
for iExp = 1:length(exp)
    nTrials = length(sequence(iExp).CorrectReach.Cue);

    for bodypart = ["HandR", "Tongue"]
        X = zeros(nTrials, length(tLocal.CorrectReach));
        Y = X;
        L = X;
        T = X;
        XBase = zeros(nTrials, nSamplesBaseline);
        YBase = XBase;
    
        for iTrial = 1:nTrials
            try
                t = [ ...
                    linspace(sequence(iExp).CorrectReach.ReachStart(iTrial)-durationPreMove, sequence(iExp).CorrectReach.ReachStart(iTrial), nSamples), ...ReachStart-2 -> ReachStart
                    linspace(sequence(iExp).CorrectReach.ReachStart(iTrial), sequence(iExp).CorrectReach.ReachEnd(iTrial), nSamples), ... ReachStart -> ReachEnd
                    linspace(sequence(iExp).CorrectReach.ReachEnd(iTrial), sequence(iExp).CorrectReach.FirstLick(iTrial), nSamples), ... ReachEnd -> FirstLick
                    linspace(sequence(iExp).CorrectReach.FirstLick(iTrial), min(sequence(iExp).CorrectReach.FirstLick(iTrial)+1, sequence(iExp).CorrectReach.RetractStart(iTrial)), nSamples), ... FirstLick -> FirstLick+1
                    linspace(min(sequence(iExp).CorrectReach.FirstLick(iTrial)+1, sequence(iExp).CorrectReach.RetractStart(iTrial)), sequence(iExp).CorrectReach.RetractStart(iTrial), nSamples), ... FirstLick+1 -> RetractStart
                    linspace(sequence(iExp).CorrectReach.RetractStart(iTrial), sequence(iExp).CorrectReach.RetractEnd(iTrial), nSamples), ... RetractStart -> RetractEnd
                    linspace(sequence(iExp).CorrectReach.RetractEnd(iTrial), sequence(iExp).CorrectReach.RetractEnd(iTrial)+2, nSamples) ... RetractEnd -> RetractEnd+2
                    ];
                if ~all(diff(t) >= 0)
                    error()
                end
                t = unique(t);
                [X(iTrial, :), Y(iTrial, :), ~] = exp(iExp).getTrajectory(t, 'r', char(bodypart), likelihoodThreshold=0.5);
                [~, ~, L(iTrial, :)] = exp(iExp).getTrajectory(t, 'r', char(bodypart), likelihoodThreshold=0);
                T(iTrial, :) = t;
        
                % Baseline [-4, -2]
                tBaseline = linspace(sequence(iExp).CorrectReach.ReachStart(iTrial)-4, sequence(iExp).CorrectReach.ReachStart(iTrial)-2, nSamplesBaseline);
                [XBase(iTrial, :), YBase(iTrial, :), ~] =  exp(iExp).getTrajectory(tBaseline, 'r', char(bodypart), likelihoodThreshold=0.5);
            catch
                warning('Cannot process "%s" trajectory for iExp=%i, iTrial=%i, %s', bodypart, iExp, iTrial, exp(iExp).name)
            end
        end
        traj(iExp).CorrectReach.(bodypart).X = (X - mean(XBase, 'all', 'omitnan')) ./ std(XBase, 0, 'all', 'omitnan');
        traj(iExp).CorrectReach.(bodypart).Y = (Y - mean(YBase, 'all', 'omitnan')) ./ std(YBase, 0, 'all', 'omitnan');
        traj(iExp).CorrectReach.(bodypart).L = L;
        traj(iExp).CorrectReach.(bodypart).t = tLocal.CorrectReach;
        traj(iExp).CorrectReach.(bodypart).T = T;
    end


    nTrials = length(sequence(iExp).CorrectLick.Cue);
    for bodypart = ["HandR", "Tongue"]
        X = zeros(nTrials, length(tLocal.CorrectLick));
        Y = X;
        L = X;
        T = X;
        XBase = zeros(nTrials, nSamplesBaseline);
        YBase = XBase;
    
        for iTrial = 1:nTrials
            try
                t = [ ...
                    linspace(sequence(iExp).CorrectLick.Lick(iTrial)-durationPreMove, sequence(iExp).CorrectLick.Lick(iTrial), nSamples), ... Lick-2 -> Lick
                    linspace(sequence(iExp).CorrectLick.Lick(iTrial), min(sequence(iExp).CorrectLick.Lick(iTrial)+1, sequence(iExp).CorrectLick.LastLickOff(iTrial)), nSamples), ... Lick -> Lick+1
                    linspace(min(sequence(iExp).CorrectLick.Lick(iTrial)+1, sequence(iExp).CorrectLick.LastLickOff(iTrial)), sequence(iExp).CorrectLick.LastLickOff(iTrial), nSamples), ... Lick+1 -> LastLickOff
                    linspace(sequence(iExp).CorrectLick.LastLickOff(iTrial), sequence(iExp).CorrectLick.LastLickOff(iTrial)+durationPreMove, nSamples), ... LastLickOff -> LastLickOff+2
                    ];
                if ~all(diff(t) >= 0)
                    error()
                end
                t = unique(t);
                [X(iTrial, :), Y(iTrial, :), ~] = exp(iExp).getTrajectory(t, 'r', char(bodypart), likelihoodThreshold=0.5);
                [~, ~, L(iTrial, :)] = exp(iExp).getTrajectory(t, 'r', char(bodypart), likelihoodThreshold=0);
                T(iTrial, :) = t;
        
                % Baseline [-4, -2]
                tBaseline = linspace(sequence(iExp).CorrectLick.Lick(iTrial)-4, sequence(iExp).CorrectLick.Lick(iTrial)-2, nSamplesBaseline);
                [XBase(iTrial, :), YBase(iTrial, :), ~] =  exp(iExp).getTrajectory(tBaseline, 'r', char(bodypart), likelihoodThreshold=0.5);
            catch
                warning('Cannot process "%s" trajectory for iExp=%i, iTrial=%i, %s', bodypart, iExp, iTrial, exp(iExp).name)
            end
        end
        traj(iExp).CorrectLick.(bodypart).X = (X - mean(XBase, 'all', 'omitnan')) ./ std(XBase, 0, 'all', 'omitnan');
        traj(iExp).CorrectLick.(bodypart).Y = (Y - mean(YBase, 'all', 'omitnan')) ./ std(YBase, 0, 'all', 'omitnan');
        traj(iExp).CorrectLick.(bodypart).L = L;
        traj(iExp).CorrectLick.(bodypart).t = tLocal.CorrectLick;
        traj(iExp).CorrectLick.(bodypart).T = T;        
    end
end
clear iExp nTrials X Y L T XBase YBase iTrial t tBaseline bodypart

%% Plot PETH
close all
layout.w = [4, 7];
layout.h = [1, 3];
fig = figure(Units='inches', Position=[0.5, 0.5, 6, 4]);
tl = tiledlayout(fig, sum(layout.h), sum(layout.w), TileSpacing='compact', Padding='loose', TileIndexing='columnmajor');
AX = gobjects(2, 2);

% CorrectLick trajectories
ax = nexttile(tl, 1, [layout.h(1), layout.w(1)]);
AX(1, 1) = ax;
hold(ax, 'on')

% CorrectLick PETH
ax = nexttile(tl, layout.h(1) + 1, [layout.h(2), layout.w(1)]);
AX(2, 1) = ax;
hold(ax, 'on')
for iEu = 1:length(eu)
    plot(ax, (tLocal.CorrectLick(1:end-1) + tLocal.CorrectLick(2:end))*0.5, msr.CorrectLick(iEu, :))
end

for iAx = 1:2
    xline(AX(iAx, 1), 0, 'k-', 'first lick')
    xline(AX(iAx, 1), 1, 'k-', 'first lick+1')
    xline(AX(iAx, 1), 2, 'k-', 'last lick')
end
ylabel(ax, 'normalized spike rate (a.u.)')
xlabel(ax, '"phase"')
xticks(ax, [])
title(ax, sprintf('n=%i SNr units, rewarded lick trials', length(eu)))

% CorrectReach trajectories
ax = nexttile(tl, layout.w(1)*sum(layout.h)+1, [layout.h(1), layout.w(2)]);
AX(1, 2) = ax;
hold(ax, 'on')
trajTemp = [traj.CorrectReach];
HandR = [trajTemp.HandR];
Tongue = [trajTemp.Tongue];
x = mean(cat(1, HandR.X), 1, 'omitnan');
y = mean(cat(1, HandR.Y), 1, 'omitnan');
l = cat(1, Tongue.L);
l(l>0.5) = 1;
l(l<0.5) = 0;
l = mean(l, 1, 'omitnan');
l = l./max(l);
colororder(ax, [0.75, 0.15, 0.15; 0.15, 0.15, 0.75])
yyaxis(ax, 'left')
plot(ax, HandR(1).t, sqrt(x.^2 + y.^2), Color=[0.75, 0.15, 0.15], DisplayName='HandR')
ylabel(ax, 'HandR')
yyaxis(ax, 'right')
plot(ax, Tongue(1).t, l, Color=[0.15, 0.15, 0.75], DisplayName='Tongue')
ylim(ax, [0, 1])
ylabel(ax, 'Tongue')

% CorrectReach PETH
ax = nexttile(tl, layout.w(1)*sum(layout.h) + layout.h(1) + 1, [layout.h(2), layout.w(2)]);
AX(2, 2) = ax;
hold(ax, 'on')
for iEu = 1:length(eu)
    % plot(ax, (tLocal.CorrectReach(1:end-1) + tLocal.CorrectReach(2:end))*0.5, (msr.CorrectReach(iEu, :) - eu(iEu).SpikeRateStats.median)./eu(iEu).SpikeRateStats.madITI, Color=[0.15, 0.15, 0.15, 0.33])
    plot(ax, (tLocal.CorrectReach(1:end-1) + tLocal.CorrectReach(2:end))*0.5, msr.CorrectReach(iEu, :))
end

for iAx = 1:2
    xline(AX(iAx, 2), -1, 'k-', 'reach start')
    xline(AX(iAx, 2), 0, 'k-', 'reach end')
    xline(AX(iAx, 2), 1, 'k-', 'lick start')
    xline(AX(iAx, 2), 2, 'k-', 'lick start+1')
    xline(AX(iAx, 2), 3, 'k-', 'retract start')
    xline(AX(iAx, 2), 4, 'k-', 'retract end')
end
ylabel(ax, 'normalized spike rate (a.u.)')
xlabel(ax, '"phase"')
xticks(ax, [])
title(ax, sprintf('n=%i SNr units, rewarded reach trials', length(eu)))

clear fig tl AX ax iExp iEu iAx


%% Plot heatmap
layout.w = [7, 4];
layout.h = [1, 5];
fig = figure(Units='inches', Position=[6.5, 0.5, 6, 9.5]);
tl = tiledlayout(fig, sum(layout.h), sum(layout.w), TileSpacing='compact', Padding='loose', TileIndexing='columnmajor');
AX = gobjects(2, 2);

[~, IValley] = min(msr.CorrectLick, [], 2, 'omitnan');
nValley = sum(msr.CorrectLick < -0.5, 2);
[~, order] = sort(nValley*max(IValley)*10 + IValley, 'ascend');

meta.lick = mean(msr.CorrectLick(:, tLocal.CorrectLick>0 & tLocal.CorrectLick<1), 2, 'omitnan');
meta.reachLick = mean(msr.CorrectReach(:, tLocal.CorrectReach>1 & tLocal.CorrectReach<2), 2, 'omitnan');

metaSign.lick = meta.lick >= 0;
metaSign.reachLick = meta.reachLick >= 0;

% groupVar = zeros(size(meta.lick));
% groupVar(~metaSign.lick & metaSign.reachLick) = 0;
% groupVar(metaSign.lick & ~metaSign.reachLick) = 1;
% groupVar(~metaSign.lick & ~metaSign.reachLick) = 2;
% groupVar(metaSign.lick & metaSign.reachLick) = 3;
groupVar = (metaSign.lick)*2 + metaSign.reachLick;
sortVar = groupVar*10*max(abs(meta.lick)) + meta.lick;
[~, order] = sort(sortVar);

% CorrectReach trajectories
ax = nexttile(tl, 1, [layout.h(1), layout.w(1)]);
AX(1, 1) = ax;
hold(ax, 'on')
trajTemp = [traj.CorrectReach];
HandR = [trajTemp.HandR];
Tongue = [trajTemp.Tongue];
x = mean(cat(1, HandR.X), 1, 'omitnan');
y = mean(cat(1, HandR.Y), 1, 'omitnan');
l = cat(1, Tongue.L);
l(l>0.5) = 1;
l(l<0.5) = 0;
l = mean(l, 1, 'omitnan');
l = l./max(l);
colororder(ax, [0.75, 0.15, 0.15; 0.15, 0.15, 0.75])
yyaxis(ax, 'left')
plot(ax, HandR(1).t, sqrt(x.^2 + y.^2), Color=[0.75, 0.15, 0.15], DisplayName='HandR')
ylabel(ax, 'HandR')
yyaxis(ax, 'right')
plot(ax, Tongue(1).t, l, Color=[0.15, 0.15, 0.75], DisplayName='Tongue')
ylim(ax, [0, 1])
ylabel(ax, 'Tongue')

% CorrectReach PETH
ax = nexttile(tl, layout.h(1) + 1, [layout.h(2), layout.w(1)]);
AX(2, 1) = ax;

eta = struct(X=msr.CorrectReach, t=(tLocal.CorrectReach(1:end-1) + tLocal.CorrectReach(2:end))*0.5, N=[]);
% [~, IValley] = min(eta.X, [], 2, 'omitnan');
% nValley = sum(eta.X < -0.5, 2);
% [~, order] = sort(nValley*max(IValley)*10 + IValley, 'ascend');
EphysUnit.plotETA(ax, eta, order=order);
applyCustomColormap(ax, [-3, 3], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);    

for iAx = 1:2
    xline(AX(iAx, 1), -1, 'k-', 'reach start')
    xline(AX(iAx, 1), 0, 'k-', 'reach end')
    xline(AX(iAx, 1), 1, 'k-', 'lick start')
    xline(AX(iAx, 1), 2, 'k-', 'lick start+1')
    xline(AX(iAx, 1), 3, 'k-', 'retract start')
    xline(AX(iAx, 1), 4, 'k-', 'retract end')
end

% CorrectLick trajectories
ax = nexttile(tl, layout.w(1)*sum(layout.h)+1, [layout.h(1), layout.w(2)]);
AX(1, 2) = ax;
hold(ax, 'on')
trajTemp = [traj.CorrectLick];
HandR = [trajTemp.HandR];
Tongue = [trajTemp.Tongue];
x = mean(cat(1, HandR.X), 1, 'omitnan');
y = mean(cat(1, HandR.Y), 1, 'omitnan');
l = cat(1, Tongue.L);
l(l>0.5) = 1;
l(l<0.5) = 0;
l = mean(l, 1, 'omitnan');
l = l./max(l);
colororder(ax, [0.75, 0.15, 0.15; 0.15, 0.15, 0.75])
yyaxis(ax, 'left')
plot(ax, HandR(1).t, sqrt(x.^2 + y.^2), Color=[0.75, 0.15, 0.15], DisplayName='HandR')
ylabel(ax, 'HandR')
yyaxis(ax, 'right')
plot(ax, Tongue(1).t, l, Color=[0.15, 0.15, 0.75], DisplayName='Tongue')
ylim(ax, [0, 1])
ylabel(ax, 'Tongue')
clear HandR Tongue trajTemp x y l

% CorrectLick PETH
ax = nexttile(tl, layout.w(1)*sum(layout.h) + layout.h(1) + 1, [layout.h(2), layout.w(2)]);
AX(2, 2) = ax;

eta = struct(X=msr.CorrectLick, t=(tLocal.CorrectLick(1:end-1) + tLocal.CorrectLick(2:end))*0.5, N=[]);
% [~, IValley] = min(eta.X, [], 2, 'omitnan');
% nValley = sum(eta.X < -0.25, 2);
% [~, order] = sort(nValley*max(IValley)*10 + IValley, 'ascend');
EphysUnit.plotETA(ax, eta, order=order);
applyCustomColormap(ax, [-3, 3], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);    

for iAx = 1:2
    xline(AX(iAx, 2), 0, 'k-', 'first lick')
    xline(AX(iAx, 2), 1, 'k-', 'first lick+1')
    xline(AX(iAx, 2), 2, 'k-', 'last lick')
end

title(AX, '')
title(AX(:, 1), 'Reach (rewarded)')
title(AX(:, 2), 'Lick (rewarded)')

xlabel(AX(2, :), '')
xlim(AX(:, 1), [-2, 5])
xlim(AX(:, 2), [-1, 3])
% xticks(AX, [])

yl = [min(arrayfun(@(ax) ax.YAxis(1).Limits(1), AX(1, :))), max(arrayfun(@(ax) ax.YAxis(1).Limits(2), AX(1, :)))];
AX(1, 1).YAxis(1).Limits = yl;
AX(1, 2).YAxis(1).Limits = yl;

xticks(AX, [])
xticks(AX(2, 1), -1:4)
xticklabels(AX(2, 1), ...
    [ ...
        "reach start", ...
        "reach end", ...
        "first lick", ...
        "first lick + 1", ...
        "retract start", ...
        "retract end", ...
    ])
xticks(AX(2, 2), 0:2)
xticklabels(AX(2, 2), ...
    [ ...
        "first lick", ...
        "first lick + 1", ...
        "last lick", ...
    ])
xtickangle(AX, 90)
fontsize(fig, 9, 'points')

clear ax AX eta iAx fig tl IValley nValley order meta


%% Do some neuron categorization by response
[~, expIndex] = ismember({eu.ExpName}, {exp.name});

meta.lick = mean(msr.CorrectLick(:, tLocal.CorrectLick>=0 & tLocal.CorrectLick<=1), 2, 'omitnan');
meta.reachLick = mean(msr.CorrectReach(:, tLocal.CorrectReach>=1 & tLocal.CorrectReach<=2), 2, 'omitnan');
meta.reach = mean(msr.CorrectReach(:, tLocal.CorrectReach>=-1 & tLocal.CorrectReach<=0), 2, 'omitnan');

tl = tiledlayout(figure(Units='inches', Position=[1, 1, 10, 3]), 1, 3);
epochs = ["lick", "reachLick", "reach"];
for i = 1:3
    ax = nexttile(tl);
    histogram(ax, meta.(epochs(i)), -30:1:30)
    title(ax, epochs(i))
end

% theta.pos = 1;
% theta.neg = -0.5;

theta.pos = 0;
theta.neg = 0;

metaSign = struct(lick=[], reachLick=[], reach=[]);
for i = 1:3
    x = meta.(epochs(i));
    y = ones(size(x));
    y(x<=theta.neg) = 0;
    y(x>=theta.pos) = 2;
    metaSign.(epochs(i)) = y;
end

base = 3;
metaHash = base^2*metaSign.reach + base^1*metaSign.reachLick + base^0*metaSign.lick;
descElements = [
    "reach-dec", "reachLick-dec", "lick-dec"; ...
    "reach-flat", "reachLick-flat", "lick-flat"; ...
    "reach-inc", "reachLick-inc", "lick-inc"; ...
    ];
hashDesc = configureDictionary("double", "string");
for i = 1:3
    for j = 1:3
        for k = 1:3
            hash = base^2*(i-1) + base^1*(j-1) + base^0*(k-1);
            desc = sprintf("%s %s %s", descElements(i, 1), descElements(j, 2), descElements(k, 3));
            hashDesc(hash) = desc;
        end
    end
end

exportPath = 'E:\Figures\reach_vs_lick_continuous\heatmaps_by_hash';
uniqueMetaHash = reshape(unique(metaHash), 1, []);
if exist(exportPath)
    rmdir(exportPath, 's')
end

layout.w = [7, 4];
close all
iFig = 0;
wFig = 4.5;

hashGroups = {[18, 8], [24, 20, 6, 2], [0, 26]};
H0 = [0.5, 0.5, 0.5];
hashGroupNames = ["action_selection", "weird", "all_congruent"];

for hash = uniqueMetaHash
    sel = metaHash == hash;

    layout.h = [1, nnz(sel)];

    iGroup = find(cellfun(@(hg) ismember(hash, hg), hashGroups));

    fig = figure(Units='inches', InnerPosition=[0+(iGroup-1)*wFig+0.25*(iGroup-1), H0(iGroup), wFig, 1+nnz(sel)*0.1], MenuBar='none'); 
    H0(iGroup) = H0(iGroup) + 1.5+nnz(sel)*0.1;
    iFig = iFig + 1;
    tl = tiledlayout(fig, 1, sum(layout.w), TileSpacing='tight', Padding='tight', TileIndexing='columnmajor');
    ax = gobjects(1, 2);
    ax(1) = nexttile(tl, 1, [1, layout.w(1)]);
    ax(2) = nexttile(tl, layout.w(1) + 1, [1, layout.w(2)]);

    tasks = ["CorrectReach", "CorrectLick"];
    for iTask = 1:2
        task = tasks(iTask);
        eta = struct(X=msr.(task), t=(tLocal.(task)(1:end-1) + tLocal.(task)(2:end))*0.5, N=[]);
        if iTask == 1
            % [~, IValley] = min(eta.X(sel, :), [], 2, 'omitnan');
            % nValley = sum(eta.X(sel, :) < 0, 2);
            tValley = zeros(size(eta.X, 1), 1);
            for iEu = 1:size(eta.X, 1)
                selT = eta.X(iEu, :)<-0.25;
                if nnz(selT) == 0
                    tValley(iEu) = Inf;
                else
                    w = abs(eta.X(iEu, selT));
                    w = w./sum(w);
                    tValley(iEu) = eta.t(selT)*w';
                end
            end
            tValley = tValley(sel);
            [~, order] = sort(tValley, 'ascend');
        end
        EphysUnit.plotETA(ax(iTask), eta, sel, order=order, hidecolorbar=iTask==1);
        applyCustomColormap(ax(iTask), [-3, 3], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);    
        
        if iTask == 1
            xline(ax(iTask), -1, 'k--')
            xline(ax(iTask), 0, 'k--')
            xline(ax(iTask), 1, 'k-', LineWidth=2)
            xline(ax(iTask), 2, 'k-', LineWidth=2)
            xline(ax(iTask), 3, 'k--')
            xline(ax(iTask), 4, 'k--')
            xticks(ax(1), -1:4)
            xticklabels(ax(1), ...
                [ ...
                    "reach start", ...
                    "reach end", ...
                    "first lick", ...
                    "first lick + 1", ...
                    "retract start", ...
                    "retract end", ...
                ])
            title(ax(iTask), 'Reach (rewarded)')
            xlim(ax(iTask), [-2, 5])
        else
            xline(ax(iTask), 0, 'k-', LineWidth=2)
            xline(ax(iTask), 1, 'k-', LineWidth=2)
            xline(ax(iTask), 2, 'k--')
            title(ax(iTask), 'Lick (rewarded)')
            xlim(ax(iTask), [-1, 3])   
            xticks(ax(iTask), 0:2)
            xticklabels(ax(iTask), ...
                [ ...
                    "first lick", ...
                    "first lick + 1", ...
                    "last lick", ...
                ])
        end
    end

    title(ax, '')
    xlabel(ax, '')
    ylabel(ax(2), '')
    yticks(ax(2), []) 
    xtickangle(ax, 90)
    fontsize(fig, 9, 'points')
        

    thisPath = sprintf("%s\\%s", exportPath, hashGroupNames(iGroup));
    if ~exist(thisPath, 'dir')
        mkdir(thisPath);
    end

    title(tl, sprintf("%i (n=%i) - %s", hash, nnz(sel), hashDesc(hash)), FontWeight='bold', FontSize=11)
    print(fig, sprintf("%s\\%i (n=%i) - %s.png", thisPath, hash, nnz(sel), hashDesc(hash)), '-dpng', '-r0')

    clear fig ax tl
    % Make the PETHs and movement traces for single neurons
    for iEu = reshape(find(sel), 1, [])
        iExp = expIndex(iEu);
        layout.h = [2, 4];
        fig = figure(Units='inches', InnerPosition=[2, 2, 6, 6]);
        tl = tiledlayout(fig, sum(layout.h), sum(layout.w), TileSpacing='tight', Padding='tight', TileIndexing='columnmajor');
        ax = gobjects(2, 2);
        for iTask = 1:2
            ax(1, iTask) = nexttile(tl, 1 + (iTask-1)*sum(layout.h)*layout.w(1), [layout.h(1), layout.w(iTask)]);
            ax(2, iTask) = nexttile(tl, 1 + (iTask-1)*sum(layout.h)*layout.w(1) + layout.h(1), [layout.h(2), layout.w(iTask)]);
        end
        
        for iTask = 1:2
            task = tasks(iTask);

            etaTemp = struct(X=msr.(task)(iEu, :), t=(tLocal.(task)(1:end-1) + tLocal.(task)(2:end))*0.5, N=[]);

            trajTemp = traj(iExp).(task);
            HandR = [trajTemp.HandR];
            Tongue = [trajTemp.Tongue];
            x = mean(cat(1, HandR.X), 1, 'omitnan');
            y = mean(cat(1, HandR.Y), 1, 'omitnan');
            l = cat(1, Tongue.L);
            l(l>0.5) = 1;
            l(l<0.5) = 0;
            l = mean(l, 1, 'omitnan');
            l = l./max(l);

            % Plot traj
            colororder(ax(1, iTask), [0.75, 0.15, 0.15; 0.15, 0.15, 0.75])
            yyaxis(ax(1, iTask), 'left')
            plot(ax(1, iTask), HandR(1).t, sqrt(x.^2 + y.^2), Color=[0.75, 0.15, 0.15], DisplayName='HandR')
            ylabel(ax(1, iTask), 'HandR')
            yyaxis(ax(1, iTask), 'right')
            plot(ax(1, iTask), Tongue(1).t, l, Color=[0.15, 0.15, 0.75], DisplayName='Tongue')
            ylim(ax(1, iTask), [0, 1])
            ylabel(ax(1, iTask), 'Tongue')
    
            % Plot PETH
            plot(ax(2, iTask), etaTemp.t, etaTemp.X, Color='k', LineWidth=1.5)
            yline(ax(2, iTask), 0, 'k--', Alpha=0.2)

            if iTask == 1
                for iRow = 1:2
                    xline(ax(iRow, iTask), -1, 'k--', Alpha=0.2)
                    xline(ax(iRow, iTask), 0, 'k--', Alpha=0.2)
                    xline(ax(iRow, iTask), 1, 'k-', Alpha=0.2, LineWidth=1)
                    xline(ax(iRow, iTask), 2, 'k-', Alpha=0.2, LineWidth=1)
                    xline(ax(iRow, iTask), 3, 'k--', Alpha=0.2)
                    xline(ax(iRow, iTask), 4, 'k--', Alpha=0.2)
                end
                xticks(ax(2, 1), -1:4)
                xticklabels(ax(2, 1), ...
                    [ ...
                        "reach start", ...
                        "reach end", ...
                        "first lick", ...
                        "first lick + 1", ...
                        "retract start", ...
                        "retract end", ...
                    ])
                title(ax(1, iTask), 'Displacement (rewarded reach)')
                title(ax(2, iTask), 'PETH (rewarded reach)')
                xlim(ax(:, iTask), [-2, 5])
            else
                for iRow = 1:2
                    xline(ax(iRow, iTask), 0, 'k-', Alpha=0.2, LineWidth=1)
                    xline(ax(iRow, iTask), 1, 'k-', Alpha=0.2, LineWidth=1)
                    xline(ax(iRow, iTask), 2, 'k--', Alpha=0.2)
                end
                title(ax(1, iTask), 'Displacement (rewarded lick)')
                title(ax(2, iTask), 'PETH (rewarded lick)')
                xlim(ax(:, iTask), [-1, 3])   
                xticks(ax(2, iTask), 0:2)
                xticklabels(ax(2, iTask), ...
                    [ ...
                        "first lick", ...
                        "first lick + 1", ...
                        "last lick", ...
                    ])
            end
            xticks(ax(1, :), [])

            clear etaTemp task trajTemp HandR Tongue x y l
        end
        xtickangle(ax, 90)
        yl = [min(arrayfun(@(ax) ax.YAxis(1).Limits(1), ax(1, :))), max(arrayfun(@(ax) ax.YAxis(1).Limits(2), ax(1, :)))];
        ax(1, 1).YAxis(1).Limits = yl;
        ax(1, 2).YAxis(1).Limits = yl;
        yl = [min(arrayfun(@(ax) ax.YLim(1), ax(2, :))), max(arrayfun(@(ax) ax.YLim(2), ax(2, :)))];
        yl = [min(-3, yl(1)), max(3, yl(2))];
        ylim(ax(2, :), yl)
        fontsize(fig, 9, 'points')
        title(tl, eu(iEu).getName(), FontSize=11, FontWeight='bold', Interpreter='none')

        thisPath = sprintf("%s\\%s\\%i", exportPath, hashGroupNames(iGroup), hash);
        if ~exist(thisPath, 'dir')
            mkdir(thisPath);
        end
        print(fig, sprintf("%s\\%i - %s.png", thisPath, hash, eu(iEu).getName()), '-dpng', '-r0')

        close(fig)
    end
end

clear tl epochs i ax theta metaSign meta i x y base metaHash descElements hashDesc i j k exportPath uniqueMetaHash layout iFig wFig hashGroups H0 hashGroupNames hash sel iGroup fig  ax AX eta IValley nValley tValley iEu selT w order desc