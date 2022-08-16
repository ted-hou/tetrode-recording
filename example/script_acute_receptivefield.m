%% Load
clear
% Load data if necessary
if ~exist('exp_D1', 'var')
    exp_D1.ar = AcuteRecording.load('C:\SERVER\Experiment_Galvo_D1Cre;DlxFlp;Ai80\AcuteRecording');
end
if ~exist('exp_A2A', 'var')
    exp_A2A.ar = AcuteRecording.load('C:\SERVER\Experiment_Galvo_A2ACre\AcuteRecording');
end

%% Plot move vs stim, color by ml
clear ax
close all
selection = AcuteRecording.makeSelection('Light', [0.4, 0.5], 'Duration', 0.01);
[~, ax{1}] = exp_D1.ar.plotStimVsMoveResponse('Lick', 'Select', selection, 'StimThreshold', 0.5, 'MoveThreshold', 0.5, 'Highlight', 'move', 'MergeGroups', 'max', 'Hue', 'ml');
[~, ax{2}] = exp_A2A.ar.plotStimVsMoveResponse('Lick', 'Select', selection, 'StimThreshold', 0.5, 'MoveThreshold', 0.5, 'Highlight', 'move', 'MergeGroups', 'max', 'Hue', 'ml');
AcuteRecording.unifyAxesLims(ax{1})
AcuteRecording.unifyAxesLims(ax{2})
AcuteRecording.drawLines(ax{1}, true, true)
AcuteRecording.drawLines(ax{2}, true, true)
clear ax


%% Temp
obj = exp_D1.ar;
critLight = [0.4, 0.5];
critDuration = 0.01;
[statsPadded, pooledGroupsPadded, stats, groups, pooledGroups] = getStuff(obj, critLight, critDuration);

function [statsPadded, pooledGroupsPadded, stats, groups, pooledGroups] = getStuff(obj, critLight, critDuration)
    % Extract groups and grouped stats for each experiment
    groups = cell(length(obj), 1);
    stats = cell(length(obj), 1);
    for iExp = 1:length(obj)
        [~, I] = obj(iExp).selectStimResponse('Light', critLight, 'Duration', critDuration);
        groups{iExp} = obj(iExp).groupByStimCondition(I, {'light', 'duration', 'ml', 'dv'});
        stats{iExp} = obj(iExp).summarizeStimResponse(groups{iExp}, 'peak');
    end
    
    % Find all unique group hashes across experiments and prepare for data merging. 
    % If a condition is not tested in certain experiments, missing data will be padded with NaNs.
    h = cellfun(@(g) [g.groupHash], groups, 'UniformOutput', false);
    uniqueHashes = unique(cat(2, h{:}));
    nGroups = length(uniqueHashes);
    statsPadded = cell(length(obj), 1);
    for iExp = 1:length(obj)
        nUnits = size(stats{iExp}, 1);
        statsPadded{iExp} = NaN(nUnits, nGroups);
        for iGrp = 1:nGroups
            sel = [groups{iExp}.groupHash] == uniqueHashes(iGrp);
            if any(sel)
                assert(nnz(sel) == 1)
                statsPadded{iExp}(:, iGrp) = stats{iExp}(:, sel);
            end
        end
    end
    
    pooledGroups = AcuteRecording.poolGroups(groups);
    for iGrp = 1:nGroups
        sel = find([pooledGroups.groupHash] == uniqueHashes(iGrp));
        if ~isempty(sel)
            pooledGroupsPadded(iGrp) = pooledGroups(sel);
        end
    end

    statsPadded = cat(1, statsPadded{:});
end