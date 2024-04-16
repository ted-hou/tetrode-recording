% 1. Read bar touch times (from ephys)
% 2. Read move init times (already done in VTD)
% 3. Get video frames between init and touch
% 4. Play video
% 5. Make start frame


%% Load EU
euReachDir = EphysUnit.load('C:\SERVER\Units\acute_3cam_reach_direction');

% Make CompleteExperiment3 from EU
expReachDir = CompleteExperiment3(euReachDir);
expReachDir.alignTimestamps();
posOrder = zeros(length(expReachDir), 4);
posNames = {'contra-out', 'contra-front', 'contra-in', 'ipsi-front'};
expIndex = zeros(length(euReachDir), 1);

i = 1;
% Correct feature names
for iExp = 1:length(expReachDir)
    expIndex(i:i+length(expReachDir(iExp).eu)-1) = iExp;
    i = i + length(expReachDir(iExp).eu);
    switch expReachDir(iExp).animalName
        case {'desmond28', 'desmond30'}
            expReachDir(iExp).vtdL = renamevars(expReachDir(iExp).vtdL, ...
                {'handIpsiCam_X', 'handIpsiCam_Y', 'handIpsiCam_Likelihood', 'footIpsiCam_X', 'footIpsiCam_Y', 'footIpsiCam_Likelihood'}, ...
                {'handContra_X', 'handContra_Y', 'handContra_Likelihood', 'footContra_X', 'footContra_Y', 'footContra_Likelihood'});
            expReachDir(iExp).vtdR = renamevars(expReachDir(iExp).vtdR, ...
                {'handIpsiCam_X', 'handIpsiCam_Y', 'handIpsiCam_Likelihood', 'footIpsiCam_X', 'footIpsiCam_Y', 'footIpsiCam_Likelihood'}, ...
                {'handIpsi_X', 'handIpsi_Y', 'handIpsi_Likelihood', 'footIpsi_X', 'footIpsi_Y', 'footIpsi_Likelihood'});
            posOrder(iExp, :) = [4, 3, 2, 1];

        case 'desmond29'
            expReachDir(iExp).vtdF = renamevars(expReachDir(iExp).vtdF, ...
                {'handIpsi_X', 'handIpsi_Y', 'handIpsi_Likelihood', 'footIpsi_X', 'footIpsi_Y', 'footIpsi_Likelihood', ...
                'handContra_X', 'handContra_Y', 'handContra_Likelihood', 'footContra_X', 'footContra_Y', 'footContra_Likelihood'}, ...
                {'handContra_X', 'handContra_Y', 'handContra_Likelihood', 'footContra_X', 'footContra_Y', 'footContra_Likelihood', ...
                'handIpsi_X', 'handIpsi_Y', 'handIpsi_Likelihood', 'footIpsi_X', 'footIpsi_Y', 'footIpsi_Likelihood'});
            expReachDir(iExp).vtdL = renamevars(expReachDir(iExp).vtdL, ...
                {'handIpsiCam_X', 'handIpsiCam_Y', 'handIpsiCam_Likelihood', 'footIpsiCam_X', 'footIpsiCam_Y', 'footIpsiCam_Likelihood'}, ...
                {'handIpsi_X', 'handIpsi_Y', 'handIpsi_Likelihood', 'footIpsi_X', 'footIpsi_Y', 'footIpsi_Likelihood'});
            expReachDir(iExp).vtdR = renamevars(expReachDir(iExp).vtdR, ...
                {'handIpsiCam_X', 'handIpsiCam_Y', 'handIpsiCam_Likelihood', 'footIpsiCam_X', 'footIpsiCam_Y', 'footIpsiCam_Likelihood'}, ...
                {'handContra_X', 'handContra_Y', 'handContra_Likelihood', 'footContra_X', 'footContra_Y', 'footContra_Likelihood'});
            posOrder(iExp, :) = [1, 2, 3, 4];            
    end
end
clear iExp i

%% Get some video clips at lever touch
clear clips
trials = cell(length(exp), 1);
clips = cell(length(exp));
t = cell(length(exp), 1);
nFramseBefore = 10;
nFramesAfter = 5;
for iExp = 1:length(exp)
    trials{iExp} = VideoDirReachTrial(exp(iExp));
    trials{iExp} = trials{iExp}(1:5);
    timestamps = [trials{iExp}.Stop];
    [clipsF, t{iExp}.front] = exp(iExp).getVideoClip(timestamps, side='f', bodyparts={'handContra'}, ...
        numFramesBefore=nFramseBefore, numFramesAfter=nFramesAfter);

    switch exp(iExp).animalName
        case {'desmond28', 'desmond30'}
            [clipsIpsi, t{iExp}.ipsi] = exp(iExp).getVideoClip(timestamps, side='r', bodyparts={'handIpsi'}, numFramesBefore=nFramseBefore, numFramesAfter=nFramesAfter);
            [clipsContra, t{iExp}.contra] = exp(iExp).getVideoClip(timestamps, side='l', bodyparts={'handContra'}, numFramesBefore=nFramseBefore, numFramesAfter=nFramesAfter);
        case 'desmond29'
            [clipsIpsi, t{iExp}.ipsi] = exp(iExp).getVideoClip(timestamps, side='l', bodyparts={'handIpsi'}, numFramesBefore=nFramseBefore, numFramesAfter=nFramesAfter);
            [clipsContra, t{iExp}.contra] = exp(iExp).getVideoClip(timestamps, side='r', bodyparts={'handContra'}, numFramesBefore=nFramseBefore, numFramesAfter=nFramesAfter);
    end
    clips{iExp} = cell(length(trials{iExp}), 1);
    for iTrial = 1:length(timestamps)
        clips{iExp}{iTrial} = horzcat(clipsF{iTrial}, clipsIpsi{iTrial}, clipsContra{iTrial});
    end
end
clear iExp timestamps clipsF clipsIpsi clipsContra

%% Make video annoatator dialog thingy
nFrames = nFramesAfter + nFramseBefore + 1;

fig = uifigure();
fig.UserData.iExp = 1;
fig.UserData.iTrial = 1;
fig.UserData.iFrame = 1;
fig.UserData.nFrames = nFramesAfter + nFramseBefore + 1;
ax = uiaxes(fig);
fig.UserData.ax = ax;
imagesc(ax, clips{fig.UserData.iExp}{fig.UserData.iTrial}(:, :, :, fig.UserData.iFrame))
axis(ax, 'image')
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';

fig.WindowKeyPressFcn = @onKeyPressed;
fig.WindowScrollWheelFcn = @onScrollWheeled;
fig.UserData.clips = clips;

function onKeyPressed(fig, data)
    
end

function onScrollWheeled(fig, data)
    if data.VerticalScrollCount == 1
        showFrame(fig, fig.UserData.iExp, fig.UserData.iTrial, min(fig.UserData.nFrames, fig.UserData.iFrame + 1))
    elseif data.VerticalScrollCount == -1
        showFrame(fig, fig.UserData.iExp, fig.UserData.iTrial, max(1, fig.UserData.iFrame - 1))
    end
end

function showFrame(fig, iExp, iTrial, iFrame)
    imagesc(fig.UserData.ax, fig.UserData.clips{iExp}{iTrial}(:, :, :, iFrame))
    fig.UserData.iFrame = iFrame;
end