%% Load
clear
% Load data if necessary
if ~exist('exp_D1', 'var')
    exp_D1.ar = AcuteRecording.load('C:\SERVER\Experiment_Galvo_D1Cre;DlxFlp;Ai80\AcuteRecording');
end
if ~exist('exp_A2A', 'var')
    exp_A2A.ar = AcuteRecording.load('C:\SERVER\Experiment_Galvo_A2ACre\AcuteRecording');
end

%% Extract data
[exp_D1.sr, exp_D1.groups] = exp_D1.ar.getStimResponse([0.4 0.5], 0.01);
exp_D1.pr = exp_D1.ar.getMoveResponse('Press');
exp_D1.lr = exp_D1.ar.getMoveResponse('Lick');
exp_D1.pos = exp_D1.ar.getProbeCoords();

[exp_A2A.sr, exp_A2A.groups] = exp_A2A.ar.getStimResponse([0.4 0.5], 0.01);
exp_A2A.pr = exp_A2A.ar.getMoveResponse('Press');
exp_A2A.lr = exp_A2A.ar.getMoveResponse('Lick');
exp_A2A.pos = exp_A2A.ar.getProbeCoords();

ar = [exp_A2A.ar, exp_D1.ar];
[~, ax(1)] = ar.plotPressVsLickResponse('Hue', 'ml', 'PressThreshold', 0.67, 'LickThreshold', 0.67);
[~, ax(2)] = ar.plotPressVsLickResponse('Hue', 'dv', 'PressThreshold', 0.67, 'LickThreshold', 0.67);
AcuteRecording.unifyAxesLims(ax)
AcuteRecording.drawLines(ax, true, true)
clear ar ax

%%
sel = max(abs(stimResp), [], 2, 'omitnan') > 0.5;
stimResp = stimResp(sel, :);
coords = coords(sel, :);
normStimResp = stimResp ./ max(abs(stimResp), [], 2, 'omitnan');

bar(mean(stimResp, 'omitnan'));
hold on
errorbar(1:8, mean(stimResp, 'omitnan'), std(stimResp, 'omitnan'))
xticks(1:8)
xticklabels({stimGroups.label})

strCoords = [1.3, -4.15; 1.3, -3.48; 1.3, -2.81; 1.3, -2.15; 3.4, -4.15; 3.4, -3.48; 3.4, -2.81; 3.4, -2.15];

[dr, d] = rfscore(normStimResp, strCoords);

%%
chunkSize = 10;
for i = 1:chunkSize:length(normStimResp)
    figure()
    plot(normStimResp(i:min(i+chunkSize-1, length(normStimResp)), :)');
end

%%
data = [0 0 0 0 0 0 0 0;
    0 0 0 0 2 2 2 2;
    0 0 0 0 0 2 0 0;
    0 0 2 0 0 0 0 0;
    0 2 0 0 0 0 0 0;];
pos = (1:8)';

[score, deltaR, dist] = rfscore(data, pos);

tbl = table(data, deltaR, deltaR./dist, score);

%%
function [score, deltaR, dist] = rfscore(data, pos)
    p = 0;
    for i = 1:size(data, 2)-1
        for j = i+1:size(data, 2)
            p = p + 1;
            dist(p) = sqrt(sum((pos(j, :) - pos(i, :)).^2));
            deltaR(:, p) = abs(data(:, j) - data(:, i));
        end
    end
    score = sum(deltaR ./ dist, 2, 'omitnan');
end