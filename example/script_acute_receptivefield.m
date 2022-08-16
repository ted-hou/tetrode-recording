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
[~, ax{1}] = exp_D1.ar.plotStimVsMoveResponse('Lick', 'Light', [0.4, 0.5], 'Duration', 0.01, 'StimThreshold', 0.5, 'MoveThreshold', 0.5, 'Highlight', 'move', 'MergeGroups', 'max', 'Hue', 'ml');
[~, ax{2}] = exp_A2A.ar.plotStimVsMoveResponse('Lick', 'Light', [0.4, 0.5], 'Duration', 0.01, 'StimThreshold', 0.5, 'MoveThreshold', 0.5, 'Highlight', 'move', 'MergeGroups', 'max', 'Hue', 'ml');
AcuteRecording.unifyAxesLims(ax{1})
AcuteRecording.unifyAxesLims(ax{2})
AcuteRecording.drawLines(ax{1}, true, true)
AcuteRecording.drawLines(ax{2}, true, true)
clear ax


%% Temp
obj = exp_D1.ar;
[stimResponse, stimGroups] = obj.getStimResponse(2, 0.01);
pressResponse = obj.getMoveResponse('Press');
lickResponse = obj.getMoveResponse('Lick');
