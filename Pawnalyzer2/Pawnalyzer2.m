% Used for manually labeling reach trajectories for the directional reach
% task (4 targets, 3 cameras)
classdef Pawnalyzer2 < handle
    properties
        t
        data
        clipParams = struct(nFramesBefore=NaN, nFramesAfter=NaN, trials=[])
        path
    end

    properties (Transient)
        eu = EphysUnit.empty;
        exp = CompleteExperiment3.empty;
        gui
        clips
        dataEMG = struct([])
    end

    methods
        %% Data
        function obj = Pawnalyzer2(varargin)
            if nargin == 0
                obj = Pawnalyzer2(uigetdir());
                return
            end

            p = inputParser();
            if ischar(varargin{1})
                p.addRequired('path', @ischar);
                p.addParameter('refEvent', 'cue', @(x) ismember(x, {'cue', 'press', 'reward'}));
                p.parse(varargin{:});
                path = p.Results.path;
                eu = EphysUnit.load(path);
                exp = CompleteExperiment3(eu);
                switch p.Results.refEvent
                    case 'cue'
                        exp.alignTimestamps(refEventNameArduino='CUE_ON', refEventNameEphys='Cue', trialDurationTolerance=0.1);
                    case 'press'
                        exp.alignTimestamps(refEventNameArduino='LEVER_PRESSED', refEventNameEphys='Press', trialDurationTolerance=1);
                    case 'reward'
                        exp.alignTimestamps(refEventNameArduino='REWARD_ON', refEventNameEphys='RewardTimes', trialDurationTolerance=1);
                    otherwise
                        error();
                end
            else
                if isa(varargin{1}, 'EphysUnit')
                    p.addRequired('eu', @(x) isa(x, 'EphysUnit'))
                    p.addOptional('exp', [], @(x) isa(x, 'CompleteExperiment3'))
                    p.addParameter('refEvent', 'cue', @(x) ismember(x, {'cue', 'press', 'reward'}));
                    p.parse(varargin{:})
                    eu = p.Results.eu;
                    exp = p.Results.exp;
                    if isempty(exp)
                        exp = CompleteExperiment3(eu);
                        switch p.Results.refEvent
                            case 'cue'
                                exp.alignTimestamps(refEventNameArduino='CUE_ON', refEventNameEphys='Cue', trialDurationTolerance=0.1);
                            case 'press'
                                exp.alignTimestamps(refEventNameArduino='LEVER_PRESSED', refEventNameEphys='Press', trialDurationTolerance=1.5);
                            case 'reward'
                                exp.alignTimestamps(refEventNameArduino='REWARD_ON', refEventNameEphys='RewardTimes', trialDurationTolerance=1.5);
                            otherwise
                                error();
                        end
                    end
                elseif isa(varargin{1}, 'CompleteExperiment3')
                    p.addRequired('exp', @(x) isa(x, 'CompleteExperiment3'))
                    p.parse(varargin{:})
                    exp = p.Results.exp;
                    eu = [exp.eu];
                end
            end

            
            % Correct feature names
%             for iExp = 1:length(exp)
%                 switch exp(iExp).animalName
%                     case {'desmond28', 'desmond30'}
%                         exp(iExp).vtdL = renamevars(exp(iExp).vtdL, ...
%                             {'handIpsiCam_X', 'handIpsiCam_Y', 'handIpsiCam_Likelihood', 'footIpsiCam_X', 'footIpsiCam_Y', 'footIpsiCam_Likelihood'}, ...
%                             {'handContra_X', 'handContra_Y', 'handContra_Likelihood', 'footContra_X', 'footContra_Y', 'footContra_Likelihood'});
%                         exp(iExp).vtdR = renamevars(exp(iExp).vtdR, ...
%                             {'handIpsiCam_X', 'handIpsiCam_Y', 'handIpsiCam_Likelihood', 'footIpsiCam_X', 'footIpsiCam_Y', 'footIpsiCam_Likelihood'}, ...
%                             {'handIpsi_X', 'handIpsi_Y', 'handIpsi_Likelihood', 'footIpsi_X', 'footIpsi_Y', 'footIpsi_Likelihood'});
%             
%                     case 'desmond29'
%                         exp(iExp).vtdF = renamevars(exp(iExp).vtdF, ...
%                             {'handIpsi_X', 'handIpsi_Y', 'handIpsi_Likelihood', 'footIpsi_X', 'footIpsi_Y', 'footIpsi_Likelihood', ...
%                             'handContra_X', 'handContra_Y', 'handContra_Likelihood', 'footContra_X', 'footContra_Y', 'footContra_Likelihood'}, ...
%                             {'handContra_X', 'handContra_Y', 'handContra_Likelihood', 'footContra_X', 'footContra_Y', 'footContra_Likelihood', ...
%                             'handIpsi_X', 'handIpsi_Y', 'handIpsi_Likelihood', 'footIpsi_X', 'footIpsi_Y', 'footIpsi_Likelihood'});
%                         exp(iExp).vtdL = renamevars(exp(iExp).vtdL, ...
%                             {'handIpsiCam_X', 'handIpsiCam_Y', 'handIpsiCam_Likelihood', 'footIpsiCam_X', 'footIpsiCam_Y', 'footIpsiCam_Likelihood'}, ...
%                             {'handIpsi_X', 'handIpsi_Y', 'handIpsi_Likelihood', 'footIpsi_X', 'footIpsi_Y', 'footIpsi_Likelihood'});
%                         exp(iExp).vtdR = renamevars(exp(iExp).vtdR, ...
%                             {'handIpsiCam_X', 'handIpsiCam_Y', 'handIpsiCam_Likelihood', 'footIpsiCam_X', 'footIpsiCam_Y', 'footIpsiCam_Likelihood'}, ...
%                             {'handContra_X', 'handContra_Y', 'handContra_Likelihood', 'footContra_X', 'footContra_Y', 'footContra_Likelihood'});
%                     otherwise
%                         error('Unknown animal %s, please manually specify camera side vs. contra/ipsi.', exp(iExp).animalName)
%                 end
%             end

            obj.eu = eu;
            obj.exp = exp;
        end

        function [clips, t] = getClips(obj, varargin)
            p = inputParser();
            p.addParameter('nFramesBefore', 15, @isnumeric);
            p.addParameter('nFramesAfter', 0, @isnumeric);
            p.addParameter('trials', [], @(x) isnumeric(x) || ismember(x, {'PressSpontaneous', 'Press', 'PressSpontaneousMedial', 'PressSpontaneousLateral'}));
            p.addParameter('keepData', true, @islogical);
            p.addParameter('noImage', false, @islogical);
            p.parse(varargin{:})
            nFramesBefore = p.Results.nFramesBefore;
            nFramesAfter = p.Results.nFramesAfter;
            selTrials = p.Results.trials;
            keepData = p.Results.keepData;
            noImage = p.Results.noImage;

            trials = cell(length(obj.exp), 1);
            clips = cell(length(obj.exp), 1);
            t = cell(length(obj.exp), 1);
            for iExp = 1:length(obj.exp)
                if ischar(selTrials)
                    trials{iExp} = obj.exp(iExp).eu(1).Trials.(selTrials);
                else
                    trials{iExp} = VideoDirReachTrial(obj.exp(iExp));
                    if ~isempty(selTrials)
                        trials{iExp} = trials{iExp}(selTrials);
                    end
                end
                timestamps = [trials{iExp}.Stop];
                if iExp == 3
                    disp(3)
                end
                [clipsF, t{iExp}.front] = obj.exp(iExp).getVideoClip(timestamps, side='f', bodyparts={}, ...
                    numFramesBefore=nFramesBefore, numFramesAfter=nFramesAfter, noImage=noImage);
            
                switch obj.exp(iExp).animalName
                    case {'desmond28', 'desmond30', 'daisy23', 'daisy24'}
                        [clipsIpsi, t{iExp}.ipsi] = obj.exp(iExp).getVideoClip(timestamps, side='r', bodyparts={}, numFramesBefore=nFramesBefore, numFramesAfter=nFramesAfter, noImage=noImage);
                        [clipsContra, t{iExp}.contra] = obj.exp(iExp).getVideoClip(timestamps, side='l', bodyparts={}, numFramesBefore=nFramesBefore, numFramesAfter=nFramesAfter, noImage=noImage);
                    case {'desmond29', 'daisy25'}
                        [clipsIpsi, t{iExp}.ipsi] = obj.exp(iExp).getVideoClip(timestamps, side='l', bodyparts={}, numFramesBefore=nFramesBefore, numFramesAfter=nFramesAfter, noImage=noImage);
                        [clipsContra, t{iExp}.contra] = obj.exp(iExp).getVideoClip(timestamps, side='r', bodyparts={}, numFramesBefore=nFramesBefore, numFramesAfter=nFramesAfter, noImage=noImage);
                end
                clips{iExp} = cell(length(trials{iExp}), 1);
                for iTrial = 1:length(timestamps)
                    clips{iExp}{iTrial} = horzcat(clipsF{iTrial}, clipsContra{iTrial}, clipsIpsi{iTrial});
                end
            end
            obj.clips = clips;
            obj.t = t;
            if keepData
                assert(obj.clipParams.nFramesBefore == nFramesBefore)
                assert(obj.clipParams.nFramesAfter == nFramesAfter)
                assert(obj.clipParams.trials == selTrials)
            end

            if ~keepData
                obj.clipParams = struct(nFramesBefore=nFramesBefore, nFramesAfter=nFramesAfter, trials=selTrials);
                for iExp = 1:obj.getLength('exp')
                    for iTrial = 1:obj.getLength('trial', exp=iExp)
                        obj.data{iExp}{iTrial} = NaN(obj.getLength('frame', exp=iExp, trial=iTrial), 8);
                    end
                end
            end
        end

        function save(obj, path)
            if nargin < 2 || isempty(path)
                if isempty(obj.path) || isequal(obj.path, [0, 0])
                    if length(obj.exp) == 1
                        [file, path] = uiputfile(sprintf('C:\\SERVER\\%s\\%s\\pawnalyzer2_%s.mat', obj.exp.animalName, obj.exp.name, obj.exp.name));
                        path = [path, file];
                        if isequal(path, [0, 0])
                            return
                        end
                    else
                        [file, path] = uiputfile('C:\SERVER\pawnalyzer2_data.mat');
                        path = [path, file];
                        if isequal(path, [0, 0])
                            return
                        end
                    end
                else
                    path = obj.path;
                end
            end
            obj.path = path;
            save(path, 'obj', '-v7.3')
            if isfield(obj.gui, 'ax') || isvalid(obj.gui.ax)
                text(obj.gui.ax, 20, 20, 'Saved!', BackgroundColor='black', Color='white', FontName='Arial', FontSize=24, VerticalAlignment='bottom', HorizontalAlignment='left')
            end
        end

        function load(obj, path)
            if nargin < 2
                [file, path] = uigetfile();
                path = [path, file];
            end
            if strcmpi(path, 'auto')
                for iExp = 1:length(obj.exp)
                    expName = obj.exp(iExp).name;
                    animalName = obj.exp(iExp).eu(1).getAnimalName();
                    dataPath = sprintf('C:\\SERVER\\%s\\%s\\pawnalyzer2_%s*.mat', animalName, expName, expName);
                    files = dir(dataPath);
                    if isempty(files)
                        continue
                    end
                    filePath = sprintf('%s\\%s', files.folder, files.name);
                    S = load(filePath);
                    assert(isa(S.obj, 'Pawnalyzer2'))
                    obj.data{iExp} = S.obj.data{1};
                end
                return
            end

            S = load(path);
            assert(isa(S.obj, 'Pawnalyzer2'))
            assert(obj.clipParams.nFramesBefore == S.obj.clipParams.nFramesBefore)
            assert(obj.clipParams.nFramesAfter == S.obj.clipParams.nFramesAfter)
%             assert(isequal(obj.clipParams.trials, S.obj.clipParams.trials))
            obj.data = S.obj.data;
        end

        %% GUI
        function start(obj, varargin)
            p = inputParser();
            p.addParameter('index', struct(exp=[], trial=[], frame=[]), @(x) isstruct(x))
            p.addParameter('forceRestart', false, @islogical)
            p.parse(varargin{:})
            index = p.Results.index;
            forceRestart = p.Results.forceRestart;

            if isfield(obj.gui, 'fig') && isgraphics(obj.gui.fig) && isvalid(obj.gui.fig)
                if forceRestart
                    close(obj.gui.fig);
                else
                    movegui(obj.gui.fig);
                    return
                end
            end

            % Use previous indices if possible
            if isempty(index.exp) || isempty(index.trial) || isempty(index.frame)
                if isfield(obj.gui, 'index')
                    index = obj.gui.index;
                else
                    index = struct(exp=1, trial=1, frame=1);
                end
            end

            obj.gui = [];

            obj.gui.fig = uifigure(Pointer='crosshair', Units='normalized', Position=[0, 0, 1, 1], WindowState='maximized');
            obj.gui.ax = uiaxes(obj.gui.fig, Units='normalized', Position=[0, 0, 1, 1]);
            axis(obj.gui.ax, 'image');
            hold(obj.gui.ax, 'on');
            obj.gui.label = gobjects(2, 3);
            obj.gui.text.title = title(obj.gui.ax, '', FontName='Arial', FontSize=30, FontWeight='bold', Interpreter='none');
            obj.gui.text.side = gobjects(1, 3);
            obj.gui.ax.XAxis.Visible = 'off';
            obj.gui.ax.YAxis.Visible = 'off';
            obj.gui.ax.YAxis.Direction = 'reverse';

            obj.gui.precisionMode = struct(enabled=false, x0=NaN, y0=NaN);
            obj.gui.trackingMode = struct(enabled=false, side='');
            obj.gui.index = struct(exp=index.exp, trial=index.trial, frame=index.frame);
            obj.gui.isMouseDown = false;
            obj.gui.isKeyDown = struct(control=false, shift=false, alt=false, capslock=false);

            obj.gui.fig.WindowButtonDownFcn = @obj.onMouseDown;
            obj.gui.fig.WindowButtonUpFcn = @obj.onMouseUp;
            obj.gui.fig.WindowButtonMotionFcn = @obj.onMouseMove;
            obj.gui.fig.WindowScrollWheelFcn = @obj.onMouseScroll;
            obj.gui.fig.WindowKeyPressFcn = @obj.onKeyDown;
            obj.gui.fig.WindowKeyReleaseFcn = @obj.onKeyUp;

            if obj.hasEMG()
                obj.gui.axEMG = uiaxes(obj.gui.fig, Units='normalized', Position=[0.33, 0, 0.33, 0.25]);
            end

            obj.showFrame(exp=index.exp, trial=index.trial, frame=index.frame);
        end

        function b = hasEMG(obj, iExp)
            if nargin < 2
                b = ~isempty(obj.dataEMG);
                return
            end
            if isempty(obj.dataEMG)
                b = false;
                return
            end
            if length(obj.dataEMG) < iExp
                b = false;
                return
            end
            b = ~isempty(obj.dataEMG(iExp).X);
        end

        function [x, y, side] = saveCursorPos(obj)
            x = obj.gui.ax.CurrentPoint(1, 1);
            y = obj.gui.ax.CurrentPoint(1, 2);
            if obj.gui.precisionMode.enabled
                x = (x - obj.gui.precisionMode.x0)./5 + obj.gui.precisionMode.x0;
                y = (y - obj.gui.precisionMode.y0)./5 + obj.gui.precisionMode.y0;
            end

            if obj.gui.trackingMode.enabled
                side = obj.gui.trackingMode.side;
            else
                switch obj.gui.fig.SelectionType
                    % Left click
                    case 'normal'
                        side = 'contra';
                    % Ctrl-click or right-click
                    case 'alt'
                        side = 'ipsi';
                end
            end

            [iExp, iTrial, iFrame, ~] = obj.drawLabel(x, y, side);
            obj.setData(x, y, iExp, iTrial, iFrame, side);
        end

        function enableTrackingMode(obj, enabled, side)
            obj.gui.trackingMode.enabled = enabled;
            if enabled
                assert(ismember(side, {'contra', 'ipsi'}))
                obj.gui.trackingMode.side = side;
            else
                obj.gui.trackingMode.side = '';
            end
        end

        function onMouseDown(obj, src, evnt)
            obj.gui.isMouseDown = true;
            [~, ~, side] = obj.saveCursorPos();

            if obj.gui.isKeyDown.capslock
                obj.enableTrackingMode(true, side);
            else
                obj.enableTrackingMode(false);
            end
        end

        function onMouseUp(obj, src, evnt)
            obj.gui.isMouseDown = false;
        end

        function onMouseMove(obj, src, evnt)
            if ~obj.gui.isMouseDown && ~obj.gui.trackingMode.enabled
                return
            end

            obj.saveCursorPos();
        end

        function onMouseScroll(obj, src, evnt)
            if evnt.VerticalScrollCount > 0
                dir = 'next';
            elseif evnt.VerticalScrollCount < 0
                dir = 'prev';
            end

            if obj.gui.isKeyDown.control
                obj.showFrame(exp=dir, trial=1, frame=1);
            elseif obj.gui.isKeyDown.shift
                obj.showFrame(trial=dir, frame=1);
            else
                obj.showFrame(frame=dir);
                if obj.gui.isMouseDown || obj.gui.trackingMode.enabled
                    obj.saveCursorPos();
                end
            end
        end

        function onKeyDown(obj, src, evnt)
            switch evnt.Key
                case {'rightarrow', 'downarrow', 'period', 'd', 'space'}
                    if obj.gui.isKeyDown.shift
                        obj.showFrame(trial='next', frame=obj.gui.index.frame);
                    else
                        obj.showFrame(frame='next');
                        if obj.gui.isMouseDown || obj.gui.trackingMode.enabled
                            obj.saveCursorPos();
                        end
                    end
                case {'leftarrow', 'uparrow', 'comma', 'a', 'tab'}
                    if obj.gui.isKeyDown.shift
                        obj.showFrame(trial='prev', frame=obj.gui.index.frame);
                    else
                        obj.showFrame(frame='prev');
                        if obj.gui.isMouseDown || obj.gui.trackingMode.enabled
                            obj.saveCursorPos();
                        end
                    end
                case 'pagedown'
                    obj.showFrame(exp='next', trial=1, frame=1);
                case 'pageup'
                    obj.showFrame(exp='prev', trial=1, frame=1);
                case {'backspace', 'home', 'z'}
                    obj.showFrame(frame=1)
                    if obj.gui.isMouseDown || obj.gui.trackingMode.enabled
                        obj.saveCursorPos();
                    end
                case {'c', 'end'}
                    obj.showFrame(frame=obj.getLength('frame'))
                    if obj.gui.isMouseDown || obj.gui.trackingMode.enabled
                        obj.saveCursorPos();
                    end
                case 'alt'
                    if ~obj.gui.precisionMode.enabled
                        obj.gui.precisionMode.enabled = true;
                        obj.gui.precisionMode.x0 = obj.gui.ax.CurrentPoint(1, 1);
                        obj.gui.precisionMode.y0 = obj.gui.ax.CurrentPoint(1, 2);
                    end
                    obj.gui.isKeyDown.alt = true;
                case {'control', 'shift'}
                    obj.gui.isKeyDown.(evnt.Key) = true;
                % Capslock+click to enable tracking mode. Capslock again to
                % exit
                case 'capslock'
                    if obj.gui.trackingMode.enabled && ~obj.gui.isKeyDown.(evnt.Key)
                        obj.enableTrackingMode(false);
                    end
                    obj.gui.isKeyDown.(evnt.Key) = true;
                case 's'
                    if obj.gui.isKeyDown.control
                        obj.save();
                    end
            end
        end

        function onKeyUp(obj, src, evnt)
            switch evnt.Key
                case 'alt'
                    obj.gui.precisionMode.enabled = false;
                    obj.gui.precisionMode.x0 = NaN;
                    obj.gui.precisionMode.y0 = NaN;
                    obj.gui.isKeyDown.alt = false;
                case {'control', 'shift', 'capslock'}
                    obj.gui.isKeyDown.(evnt.Key) = false;
            end
        end

        function trials = getTrials(obj, iExp)
            if ischar(obj.clipParams.trials)
                trials = obj.exp(iExp).eu(1).Trials.(obj.clipParams.trials);
            else
%                 error('Not tested')
                if isempty(obj.clipParams.trials)
                    trials = obj.exp(iExp).eu(1).Trials.Press;
                else
                    trials = obj.exp(iExp).eu(1).Trials.Press(obj.clipParams.trials);
                end
            end
        end

        function showFrame(obj, varargin)
            p = inputParser();
            p.addParameter('exp', 'curr', @(x) isnumeric(x) || ismember(x, {'curr', 'next', 'prev'}));
            p.addParameter('trial', 'curr', @(x) isnumeric(x) || ismember(x, {'curr', 'next', 'prev'}));
            p.addParameter('frame', 'curr', @(x) isnumeric(x) || ismember(x, {'curr', 'next', 'prev'}));
            p.parse(varargin{:})
            iExp = obj.getIndex('exp', p.Results.exp);
            iTrial = obj.getIndex('trial', p.Results.trial);
            iFrame = obj.getIndex('frame', p.Results.frame);

            cla(obj.gui.ax)
            imagesc(obj.gui.ax, obj.clips{iExp}{iTrial}(:, :, :, iFrame))
            obj.gui.index.exp = iExp;
            obj.gui.index.trial = iTrial;
            obj.gui.index.frame = iFrame;
            nFrames = obj.getLength('frame');
            iFrameRef = nFrames - obj.clipParams.nFramesAfter;
            trials = obj.getTrials(iExp);
            tRef = trials(iTrial).Stop;

            % Show title text
            obj.gui.text.title.String = sprintf('%s\nExp: %i / %i (Ctrl)\nTrial: %03i / %03i (Shift)\nFrame: %02i / %02i', obj.exp(iExp).name, iExp, obj.getLength('exp'), iTrial, obj.getLength('trial'), iFrame, nFrames);

            % Show camera textD
            w = size(obj.clips{iExp}{iTrial}, 2) ./ 3;
            sidetext = ...
            { ...
                sprintf('Front\n(%02i/%02i)\n%.1fms', iFrame, nFrames, 1e3*(obj.t{iExp}.front(iTrial, iFrame) - tRef)), ...
                sprintf('Contra-%s\n(%02i/%02i)\n%.1fms', obj.convertSide(iExp, 'contra'), iFrame, nFrames, 1e3*(obj.t{iExp}.contra(iTrial, iFrame) - tRef)), ...
                sprintf('Ipsi-%s\n(%02i/%02i)\n%.1fms', obj.convertSide(iExp, 'ipsi'), iFrame, nFrames, 1e3*(obj.t{iExp}.ipsi(iTrial, iFrame) - tRef)) ...
            };
            for iCam = 1:3
                obj.gui.text.side(iCam) = text(obj.gui.ax, w*(iCam-1) + 20, 20, sidetext{iCam}, BackgroundColor='black', Color='white', FontName='Arial', FontWeight='bold', FontSize=18, VerticalAlignment='top', HorizontalAlignment='left');
            end

            % Show annotation
            hasData = ~isempty(obj.data);
            if hasData
                side = 'contra';
                for iCam = [1, 2]
                    [x, y] = obj.getData(iExp, iTrial, iFrame, side, iCam);
                    obj.drawLabel(x, y, side, exp=iExp, trial=iTrial, frame=iFrame);
                end
                side = 'ipsi';
                for iCam = [1, 3]
                    [x, y] = obj.getData(iExp, iTrial, iFrame, side, iCam);
                    obj.drawLabel(x, y, side, exp=iExp, trial=iTrial, frame=iFrame);
                end
            end

            if obj.hasEMG(iExp)
                obj.drawEMG(iExp, iTrial, iFrame);
            end
        end

        function [iExp, iTrial, iFrame, iCam] = drawLabel(obj, x, y, side, varargin)
            p = inputParser();
            p.addRequired('x', @isnumeric)
            p.addRequired('y', @isnumeric)
            p.addRequired('side', @(x) ismember(x, {'contra', 'ipsi'}))
            p.addParameter('exp', 'curr', @(x) isnumeric(x) || ismember(x, {'curr', 'prev', 'next'}))
            p.addParameter('trial', 'curr', @(x) isnumeric(x) || ismember(x, {'curr', 'prev', 'next'}))
            p.addParameter('frame', 'curr', @(x) isnumeric(x) || ismember(x, {'curr', 'prev', 'next'}))
            p.parse(x, y, side, varargin{:})
            x = p.Results.x;
            y = p.Results.y;
            side = p.Results.side;
            iExp = obj.getIndex('exp', p.Results.exp);
            iTrial = obj.getIndex('trial', p.Results.trial);
            iFrame = obj.getIndex('frame', p.Results.frame);

            iCam = 1;
            if isnan(x) || isnan(y)
                return
            end

            w = size(obj.clips{iExp}{iTrial}, 2) ./ 3;
            h = size(obj.clips{iExp}{iTrial}, 1);
            if x < 0 || x > w*3 || y < 0 || y > h
                return
            end

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

            if ~isvalid(obj.gui.label(iSide, iCam)) || ~isgraphics(obj.gui.label(iSide, iCam))
                obj.gui.label(iSide, iCam) = plot(obj.gui.ax, x, y, style, MarkerSize=25);
            else
                set(obj.gui.label(iSide, iCam), XData=x, YData=y);
            end
        end

        function drawEMG(obj, iExp, iTrial, iFrame)
            ax = obj.gui.axEMG;
            cla(ax)
            hold(ax, 'on')
            xEMG = obj.dataEMG(iExp).X(iTrial, :);
            tEMG = obj.dataEMG(iExp).T(iTrial, :);
            tFrame = obj.t{iExp}.contra(iTrial, iFrame);
            plot(ax, 1000.*(tEMG - tEMG(end)), xEMG, 'k');
            [~, sel] = min(abs(tEMG - tFrame));
            plot(ax, 1000.*(tEMG(sel) - tEMG(end)), xEMG(sel), 'r.', MarkerSize=25)
            hold(ax, 'off')
        end

        function setData(obj, x, y, iExp, iTrial, iFrame, side)
            w = size(obj.clips{iExp}{iTrial}, 2) ./ 3;
            h = size(obj.clips{iExp}{iTrial}, 1);
            if x < 0 || x > w*3 || y < 0 || y > h
                return
            end
            if x <= w
                iCam = 1;
            elseif x <= 2*w
                iCam = 2;
            else
                iCam = 3;
            end

            switch side
                case 'contra'
                    switch iCam
                        case 1
                            ix = 1;
                            iy = 2;
                        case 2
                            ix = 3;
                            iy = 4;
                        case 3
                            return
                    end
                case 'ipsi'
                    switch iCam
                        case 1
                            ix = 5;
                            iy = 6;
                        case 2
                            return
                        case 3
                            ix = 7;
                            iy = 8;
                    end
            end

            obj.data{iExp}{iTrial}(iFrame, ix) = x;
            obj.data{iExp}{iTrial}(iFrame, iy) = y;            
        end

        function [x, y] = getData(obj, iExp, iTrial, iFrame, side, iCam)
            switch side
                case 'contra'
                    switch iCam
                        case 1
                            ix = 1;
                            iy = 2;
                        case 2
                            ix = 3;
                            iy = 4;
                    end
                case 'ipsi'
                    switch iCam
                        case 1
                            ix = 5;
                            iy = 6;
                        case 3
                            ix = 7;
                            iy = 8;
                    end
            end

            x = obj.data{iExp}{iTrial}(iFrame, ix);
            y = obj.data{iExp}{iTrial}(iFrame, iy);
        end

        % x right, y forward, z up (not Maya/Unity)
        function [contra, ipsi, targetPos, t] = getTrajectories(obj, varargin)
            p = inputParser();
            p.addRequired('exp', @isnumeric)
            p.addOptional('trials', [], @isnumeric)
            p.addParameter('videoWidth', 640, @isnumeric);
            p.addParameter('zero', 'off', @(x) isnumeric(x) || ismember(x, {'off', 'start', 'end'}));
            p.parse(varargin{:})
            iExp = p.Results.exp;
            selTrials = p.Results.trials;
            videoWidth = p.Results.videoWidth;
            zeroMode = p.Results.zero;

            if isempty(selTrials)
                nTrials = obj.getLength('trial', exp=iExp);
                selTrials = 1:nTrials;
            else
                assert(isempty(obj.clipParams.trials), 'Double trial selection. Schenannigans afoot.')
            end

            rawData = cat(1, obj.data{iExp}{:});
            rawData(any(isnan(rawData), 2), :) = [];

            switch obj.exp(iExp).animalName
                case {'desmond28', 'desmond30', 'daisy23', 'daisy24'}
                    rawData(:, 3) = 2*videoWidth - rawData(:, 3);
                    rawData(:, 7) = rawData(:, 7) - 2*videoWidth;
                case {'desmond29', 'daisy25'}
                    rawData(:, 3) = rawData(:, 3) - videoWidth;
                    rawData(:, 7) = 3*videoWidth - rawData(:, 7);
                otherwise
                    error();
            end

            switch obj.exp(iExp).animalName
                case {'desmond28', 'desmond29', 'desmond30'}
                    trialType = 'Press';
                case {'daisy23', 'daisy24', 'daisy25'}
                    trialType = 'PressSpontaneous';
                otherwise
                    error();
            end

            scale = struct( ...
                contra=std(rawData(:, 4))./std(rawData(:, 2)), ...
                ipsi=std(rawData(:, 8))./std(rawData(:, 6)) ...
                );
            nTrials = length(selTrials);
            nFrames = obj.clipParams.nFramesAfter + obj.clipParams.nFramesBefore + 1;
            contraX = NaN(nTrials, nFrames);
            contraY = contraX;
            contraZ = contraX;
            ipsiX = contraX;
            ipsiY = contraX;
            ipsiZ = contraX;
            
            targetPos = VideoDirReachTrial(obj.exp(iExp), trialType).getTargetPos();
            targetPos = targetPos(selTrials);


            for iTrial = 1:nTrials
                if any(isnan(obj.data{iExp}{iTrial}), 2)
                    continue;
                end
            
                rawData = array2table(obj.data{iExp}{iTrial}, ...
                    VariableNames={'contraFrontX', 'contraFrontY', 'contraSideX', 'contraSideY', ...
                    'ipsiFrontX', 'ipsiFrontY', 'ipsiSideX', 'ipsiSideY'});

                switch obj.exp(iExp).animalName
                    case {'desmond28', 'desmond30', 'daisy23', 'daisy24'}
                        rawData.contraSideX = 2*videoWidth - rawData.contraSideX;
                        rawData.ipsiSideX = rawData.ipsiSideX - 2*videoWidth;
                    case {'desmond29', 'daisy25'}
                        rawData.contraSideX = rawData.contraSideX - videoWidth;
                        rawData.ipsiSideX = 3*videoWidth - rawData.ipsiSideX;
                    otherwise
                        error();
                end
                
                trials = obj.getTrials(iExp);
                tGlobal = (-15:0)./30 + trials(iTrial).Stop;
                tGlobal = tGlobal(:);
                t = (-15:0)./30;
                t = t(:);
                
                rawData.contraFrontX(:) = interp1(obj.t{iExp}.front(iTrial, :), rawData.contraFrontX, tGlobal, 'linear', 'extrap');
                rawData.contraFrontY(:) = interp1(obj.t{iExp}.front(iTrial, :), rawData.contraFrontY, tGlobal, 'linear', 'extrap');
                rawData.ipsiFrontX(:) = interp1(obj.t{iExp}.front(iTrial, :), rawData.ipsiFrontX, tGlobal, 'linear', 'extrap');
                rawData.ipsiFrontY(:) = interp1(obj.t{iExp}.front(iTrial, :), rawData.ipsiFrontY, tGlobal, 'linear', 'extrap');
                rawData.contraSideX(:) = interp1(obj.t{iExp}.contra(iTrial, :), rawData.contraSideX, tGlobal, 'linear', 'extrap');
                rawData.contraSideY(:) = interp1(obj.t{iExp}.contra(iTrial, :), rawData.contraSideY, tGlobal, 'linear', 'extrap');
                rawData.ipsiSideX(:) = interp1(obj.t{iExp}.ipsi(iTrial, :), rawData.ipsiSideX, tGlobal, 'linear', 'extrap');
                rawData.ipsiSideY(:) = interp1(obj.t{iExp}.ipsi(iTrial, :), rawData.ipsiSideY, tGlobal, 'linear', 'extrap');

                varNames = rawData.Properties.VariableNames;
                if isnumeric(zeroMode)
                    rawData = array2table(table2array(rawData) - table2array(rawData(zeroMode, :)), VariableNames=varNames);
                else
                    switch zeroMode
                        case 'start'
                            rawData = array2table(table2array(rawData) - table2array(rawData(1, :)), VariableNames=varNames);
                        case 'end'
                            rawData = array2table(table2array(rawData) - table2array(rawData(end, :)), VariableNames=varNames);
                    end
                end

                % x right, y forward, z up
                contraX(iTrial, :) = rawData.contraFrontX;
                contraY(iTrial, :) = rawData.contraSideX./scale.contra;
                contraZ(iTrial, :) = 0.5*(rawData.contraFrontY + rawData.contraSideY./scale.contra);
                
                ipsiX(iTrial, :) = rawData.ipsiFrontX;
                ipsiY(iTrial, :) = rawData.ipsiSideX./scale.ipsi;
                ipsiZ(iTrial, :) = 0.5*(rawData.ipsiFrontY + rawData.ipsiSideY./scale.ipsi);
            end

            contra.x = contraX;
            contra.y = contraY;
            contra.z = contraZ;
            ipsi.x = ipsiX;
            ipsi.y = ipsiY;
            ipsi.z = ipsiZ;
        end

        % Get number of exp/trial/frame available in obj.clips
        function n = getLength(obj, type, varargin)
            switch lower(type)
                case 'exp'
                    try
                        n = length(obj.clips);
                    catch
                        n = length(obj.data);
                    end
                case 'trial'
                    p = inputParser();
                    p.addOptional('exp', 'curr', @(x) isnumeric(x) ||  ismember(x, {'curr', 'prev', 'next'}))
                    p.parse(varargin{:})
                    if ischar(p.Results.exp)
                        iExp = obj.getIndex('exp', p.Results.exp);
                    else
                        iExp = p.Results.exp;
                    end
                    try
                        n = length(obj.clips{iExp});
                    catch
                        n = length(obj.data{iExp});
                    end
                case 'frame'
                    p = inputParser();
                    p.addOptional('exp', 'curr', @(x) isnumeric(x) ||  ismember(x, {'curr', 'prev', 'next'}))
                    p.addOptional('trial', 'curr', @(x) isnumeric(x) ||  ismember(x, {'curr', 'prev', 'next'}))
                    p.parse(varargin{:})
                    if ischar(p.Results.exp)
                        iExp = obj.getIndex('exp', p.Results.exp);
                    else
                        iExp = p.Results.exp;
                    end
                    if ischar(p.Results.trial)
                        iTrial = obj.getIndex('trial', p.Results.trial);
                    else
                        iTrial = p.Results.trial;
                    end
                    try
                        n = size(obj.data{iExp}{iTrial}, 1);
                        assert(n>1)
                    catch
                        n = size(obj.clips{iExp}{iTrial}, 4);
                    end
                otherwise
                    error()
            end
        end

        % Get exp/trial/frame index bounded by [1, max] available to
        % obj.clips
        function i = getIndex(obj, type, value)
            if nargin < 3
                value = 'curr';
            end
            assert(isnumeric(value) || ismember(value, {'curr', 'prev', 'next'}))

            currentValue = obj.gui.index.(type);
            switch lower(type)
                case 'exp'
                    maxValue = obj.getLength('exp');
                case 'trial'
                    maxValue = obj.getLength('trial', exp=obj.getIndex('exp', 'curr'));
                case 'frame'
                    maxValue = obj.getLength('frame', exp=obj.getIndex('exp', 'curr'), trial=obj.getIndex('trial', 'curr'));
            end
            if ischar(value)
                switch lower(value)
                    case 'curr'
                        value = currentValue;
                    case 'prev'
                        value = currentValue - 1;
                    case 'next'
                        value = currentValue + 1;
                end
            end
            i = max(min(value, maxValue), 1);
        end

        function str = convertSide(obj, iExp, side)
            switch obj.exp(iExp).animalName
                case {'desmond28', 'desmond30', 'daisy23', 'daisy24'}
                    switch side
                        case 'contra'
                            str = 'L';
                        case 'ipsi'
                            str = 'R';
                        case 'L'
                            str = 'contra';
                        case 'R'
                            str = 'ipsi';
                    end
                case {'desmond29', 'daisy25'}
                    switch side
                        case 'contra'
                            str = 'R';
                        case 'ipsi'
                            str = 'L';
                        case 'L'
                            str = 'ipsi';
                        case 'R'
                            str = 'contra';
                    end
                otherwise
                    error('Unknown animal %s, please manually specify camera side vs. contra/ipsi.', obj.exp(iExp).animalName)
            end            
        end
    end
end