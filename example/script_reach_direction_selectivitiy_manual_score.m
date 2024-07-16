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
trials = cell(length(expReachDir), 1);
clips = cell(length(expReachDir), 1);
t = cell(length(expReachDir), 1);
nFramseBefore = 10;
nFramesAfter = 5;
for iExp = 1:length(expReachDir)
    trials{iExp} = VideoDirReachTrial(expReachDir(iExp));
    trials{iExp} = trials{iExp}(1:5);
    timestamps = [trials{iExp}.Stop];
    [clipsF, t{iExp}.front] = expReachDir(iExp).getVideoClip(timestamps, side='f', bodyparts={'handContra'}, ...
        numFramesBefore=nFramseBefore, numFramesAfter=nFramesAfter);

    switch expReachDir(iExp).animalName
        case {'desmond28', 'desmond30'}
            [clipsIpsi, t{iExp}.ipsi] = expReachDir(iExp).getVideoClip(timestamps, side='r', bodyparts={'handIpsi'}, numFramesBefore=nFramseBefore, numFramesAfter=nFramesAfter);
            [clipsContra, t{iExp}.contra] = expReachDir(iExp).getVideoClip(timestamps, side='l', bodyparts={'handContra'}, numFramesBefore=nFramseBefore, numFramesAfter=nFramesAfter);
        case 'desmond29'
            [clipsIpsi, t{iExp}.ipsi] = expReachDir(iExp).getVideoClip(timestamps, side='l', bodyparts={'handIpsi'}, numFramesBefore=nFramseBefore, numFramesAfter=nFramesAfter);
            [clipsContra, t{iExp}.contra] = expReachDir(iExp).getVideoClip(timestamps, side='r', bodyparts={'handContra'}, numFramesBefore=nFramseBefore, numFramesAfter=nFramesAfter);
    end
    clips{iExp} = cell(length(trials{iExp}), 1);
    for iTrial = 1:length(timestamps)
        clips{iExp}{iTrial} = horzcat(clipsF{iTrial}, clipsContra{iTrial}, clipsIpsi{iTrial});
    end
end
clear iExp timestamps clipsF clipsIpsi clipsContra

%% Make video annoatator dialog thingy
% Mouse scroll wheel to change frame
% Left click to label contra paw
% Right click to label ipsi paw
% Drag to edit
nFrames = nFramesAfter + nFramseBefore + 1;

fig = uifigure();
fig.UserData.iExp = 1;
fig.UserData.iTrial = 1;
fig.UserData.iFrame = 1;
fig.UserData.nFrames = nFramesAfter + nFramseBefore + 1;
fig.UserData.nTrials = cellfun(@length, clips);
ax = uiaxes(fig);
fig.UserData.ax = ax;
imagesc(ax, clips{fig.UserData.iExp}{fig.UserData.iTrial}(:, :, :, fig.UserData.iFrame))
axis(ax, 'image')
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
hold(ax, 'on')
fig.UserData.clips = clips;
fig.UserData.currentMode = 'none';
fig.UserData.ax = ax;
fig.UserData.hLabel = gobjects(2, 3);
fig.UserData.isMouseDown = false;
fig.UserData.precisionMode = struct(enabled=false, x0=NaN, y0=NaN);
for iExp = 1:length(clips)
    for iTrial = 1:length(clips{iExp})
        fig.UserData.output{iExp}{iTrial} = NaN(nFrames, 8); 
    end
end

fig.WindowButtonMotionFcn = @onMouseMove;
fig.WindowButtonDownFcn = @onMouseDown;
fig.WindowButtonUpFcn = @onMouseUp;
fig.KeyPressFcn = @onKeyDown;
fig.KeyReleaseFcn = @onKeyUp;
fig.WindowScrollWheelFcn = @onScrollWheeled;

function label(fig, x, y, side)
    if isnan(x) || isnan(y)
        return;
    end
    w = size(fig.UserData.clips{fig.UserData.iExp}{fig.UserData.iTrial}, 2) ./ 3;
    if x <= w
        iCam = 1;
    elseif x <= 2*w
        iCam = 2;
    else
        iCam = 3;
    end

    switch side
        case 'contra'
            iSide = 1;
            style = 'ro';
            if iCam == 3
                return;
            end
        case 'ipsi'
            iSide = 2;
            style = 'yo';
            if iCam == 2
                return;
            end
    end
    if ~isvalid(fig.UserData.hLabel(iSide, iCam)) || ~isgraphics(fig.UserData.hLabel(iSide, iCam))
        fig.UserData.hLabel(iSide, iCam) = plot(fig.UserData.ax, x, y, style, MarkerSize=25);
    else
        set(fig.UserData.hLabel(iSide, iCam), XData=x, YData=y);
    end
end

function onMouseDown(fig, data)
    fig.UserData.isMouseDown = true;
    x = fig.UserData.ax.CurrentPoint(1, 1);
    y = fig.UserData.ax.CurrentPoint(1, 2);
    if fig.UserData.precisionMode.enabled
        fig.UserData.precisionMode.x0 = x;
        fig.UserData.precisionMode.y0 = y;
    else
        fig.UserData.precisionMode.x0 = NaN;
        fig.UserData.precisionMode.y0 = NaN;
    end
    switch fig.UserData.currentMode
        case 'contra'
            label(fig, x, y, 'contra')
        case 'ipsi'
            label(fig, x, y, 'ipsi')
    end
end

function onMouseUp(fig, data)
    fig.UserData.isMouseDown = false;

    % Store data
    writeData(fig, fig.UserData.iExp, fig.UserData.iTrial, fig.UserData.iFrame, fig.UserData.hLabel, 'contra', 1);
    writeData(fig, fig.UserData.iExp, fig.UserData.iTrial, fig.UserData.iFrame, fig.UserData.hLabel, 'contra', 3);
    writeData(fig, fig.UserData.iExp, fig.UserData.iTrial, fig.UserData.iFrame, fig.UserData.hLabel, 'ipsi', 1);
    writeData(fig, fig.UserData.iExp, fig.UserData.iTrial, fig.UserData.iFrame, fig.UserData.hLabel, 'ipsi', 2);
end

% (iFrame, iFeat)
% iFeat 1-8 are contraX_front, contraY_front, contraX_side, contraY_side,
% ipsiX_front, ipsiY_front, ipsiX_side, ipsiY_side,
function writeData(fig, iExp, iTrial, iFrame, hLabel, side, iCam)
    switch side
        case 'contra'
            iSide = 1;
            switch iCam
                case 1
                    iX = 1;
                    iY = 2;
                case 2
                    iX = 3;
                    iY = 4;
            end
        case 'ipsi'
            iSide = 2;
            switch iCam
                case 1
                    iX = 5;
                    iY = 6;
                case 3
                    iX = 7;
                    iY = 8;
            end
    end

    if isgraphics(hLabel(iSide, iCam)) && isvalid(hLabel(iSide, iCam))
        fig.UserData.output{iExp}{iTrial}(iFrame, iX) = hLabel(iSide, iCam).XData;
        fig.UserData.output{iExp}{iTrial}(iFrame, iY) = hLabel(iSide, iCam).YData;
    end
end

function [x, y] = readData(fig, iExp, iTrial, iFrame, side, iCam)
    switch side
        case 'contra'
            iSide = 1;
            switch iCam
                case 1
                    iX = 1;
                    iY = 2;
                case 2
                    iX = 3;
                    iY = 4;
            end
        case 'ipsi'
            iSide = 2;
            switch iCam
                case 1
                    iX = 5;
                    iY = 6;
                case 3
                    iX = 7;
                    iY = 8;
            end
    end

    x = fig.UserData.output{iExp}{iTrial}(iFrame, iX);
    y = fig.UserData.output{iExp}{iTrial}(iFrame, iY);
end

function onMouseMove(fig, data)
    if ~fig.UserData.isMouseDown
        return
    end
    x = fig.UserData.ax.CurrentPoint(1, 1);
    y = fig.UserData.ax.CurrentPoint(1, 2);
    if fig.UserData.precisionMode.enabled
        x = (x - fig.UserData.precisionMode.x0) ./ 10 + fig.UserData.precisionMode.x0;
        y = (y - fig.UserData.precisionMode.y0) ./ 10 + fig.UserData.precisionMode.y0;
    end
    switch fig.UserData.currentMode
        case 'contra'
            label(fig, x, y, 'contra')
        case 'ipsi'
            label(fig, x, y, 'ipsi')
    end
end

function onKeyDown(fig, data)
    switch data.Key
        case {'rightarrow', 'period', 'd'}
            showFrame(fig, fig.UserData.iExp, fig.UserData.iTrial, min(fig.UserData.nFrames, fig.UserData.iFrame + 1))
        case {'leftarrow', 'comma', 'a'}
            showFrame(fig, fig.UserData.iExp, fig.UserData.iTrial, max(1, fig.UserData.iFrame - 1))
        case {'downarrow', 'pagedown', 's'}
            if ismember('shift', data.Modifier)
                showFrame(fig, min(length(fig.UserData.nTrials), fig.UserData.iExp + 1), 1, 1)
            else
                showFrame(fig, fig.UserData.iExp, min(fig.UserData.nTrials(fig.UserData.iExp), fig.UserData.iTrial + 1), 1)
            end
        case {'uparrow', 'pageup', 'w'}
            if ismember('shift', data.Modifier)
                showFrame(fig, max(1, fig.UserData.iExp - 1), 1, 1)
            else
                showFrame(fig, fig.UserData.iExp, max(1, fig.UserData.iTrial - 1), 1)
            end
        case 'z'
            fig.UserData.currentMode = 'contra';
        case 'c'
            fig.UserData.currentMode = 'ipsi';
        case 'alt'
            if ~fig.UserData.precisionMode.enabled
                fig.UserData.precisionMode.enabled = true;
                fig.UserData.precisionMode.x0 = fig.UserData.ax.CurrentPoint(1, 1);
                fig.UserData.precisionMode.y0 = fig.UserData.ax.CurrentPoint(1, 2);
            end
    end
end

function onKeyUp(fig, data)
    switch data.Key
        case 'alt'
            fig.UserData.precisionMode.enabled = false;
            fig.UserData.precisionMode.x0 = NaN;
            fig.UserData.precisionMode.y0 = NaN;
    end
end

function onScrollWheeled(fig, data)
    if data.VerticalScrollCount == 1
        showFrame(fig, fig.UserData.iExp, fig.UserData.iTrial, min(fig.UserData.nFrames, fig.UserData.iFrame + 1))
    elseif data.VerticalScrollCount == -1
        showFrame(fig, fig.UserData.iExp, fig.UserData.iTrial, max(1, fig.UserData.iFrame - 1))
    end
end

function showFrame(fig, iExp, iTrial, iFrame)
    cla(fig.UserData.ax)
    imagesc(fig.UserData.ax, fig.UserData.clips{iExp}{iTrial}(:, :, :, iFrame))
    fig.UserData.iExp = iExp;
    fig.UserData.iTrial = iTrial;
    fig.UserData.iFrame = iFrame;

    side = 'contra';
    for iCam = [1, 2]
        [x, y] = readData(fig, iExp, iTrial, iFrame, side, iCam);
        label(fig, x, y, side)
    end
    side = 'ipsi';
    for iCam = [1, 3]
        [x, y] = readData(fig, iExp, iTrial, iFrame, side, iCam);
        label(fig, x, y, side)
    end
end