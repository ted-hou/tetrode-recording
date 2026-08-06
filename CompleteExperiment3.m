classdef CompleteExperiment3 < CompleteExperiment
    properties
        vtdF = []
        tce = []
    end

    methods
        function obj = CompleteExperiment3(varargin)
            if nargin == 0
                return
            end
            p = inputParser();
            p.addRequired('eu', @(x) isa(x, 'EphysUnit'))
            p.addParameter('cameras', 'flr', @(x) all(ismember(x, 'flr')))
            p.addParameter('deeplabcutPath', '', @(x) ischar(x) || isstring(x))
            p.parse(varargin{:});
            eu = p.Results.eu;
            cameras = p.Results.cameras;
            deeplabcutPath = p.Results.deeplabcutPath;

            [uniqueExpNames, ~, expIndices] = unique({eu.ExpName});
            nExp = length(uniqueExpNames);
            obj(nExp) = CompleteExperiment3();
            for i = 1:nExp
                obj(i).name = uniqueExpNames{i};
                obj(i).eu = eu(expIndices==i);
                try
                    obj(i).tce = obj(i).eu(1).loadTwoColorExperiment();
                catch
                    obj(i).tce = [];
                end
                if ismember('f', cameras)
                    obj(i).vtdF = obj(i).readOrCreateVideoTrackingData(obj(i).name, 'f', deeplabcutPath);
                end
                if ismember('l', cameras)
                    obj(i).vtdL = obj(i).readOrCreateVideoTrackingData(obj(i).name, 'l', deeplabcutPath);
                end
                if ismember('r', cameras)
                    obj(i).vtdR = obj(i).readOrCreateVideoTrackingData(obj(i).name, 'r', deeplabcutPath);
                end
                obj(i).ac = CompleteExperiment.readArduino(obj(i).name);
            end
        end

        function vtd = readOrCreateVideoTrackingData(obj, expName, side, path)
            sidenum = obj.getCameraIndex(side);
            animalName = strsplit(expName, '_');
            animalName = animalName{1};

            if sidenum == -1
                vtd = [];
                return
            end

            if sidenum ~= 0
                if isempty(path)
                    if isempty(obj.tce)
                        deepLabCutFiles = dir(sprintf('C:\\SERVER\\%s\\%s\\%s_%g*.csv', animalName, expName, expName, sidenum));
                    else
                        deepLabCutFiles = dir(sprintf('C:\\SERVER\\%s\\%s\\%s_laser_%g*.csv', animalName, expName, expName, sidenum));
                    end
                else
                    if isempty(obj.tce)
                        deepLabCutFiles = dir(sprintf('%s\\%s_%g*.csv', path, expName, sidenum));
                    else
                        deepLabCutFiles = dir(sprintf('%s\\%s_laser_%g*.csv', path, expName, sidenum));
                    end
                end
            else
                if isempty(path)
                    deepLabCutFiles = dir(sprintf('C:\\SERVER\\%s\\%s\\%s*.csv', animalName, expName, expName));
                else
                    deepLabCutFiles = dir(sprintf('%s\\%s*.csv', path, expName));
                end
            end
            % Deep lab cut
            if ~isempty(deepLabCutFiles)
                if isempty(obj.tce)
                    vtd = obj.readVideoTrackingData(expName, side, sprintf('%s', expName), deeplabcutPath=path);
                else
                    vtd = obj.readVideoTrackingData(expName, side, sprintf('%s_laser', expName), deeplabcutPath=path);
                end
            % Pawnalyzer2 manually labeled
            else
                vidFile = dir(sprintf('C:\\SERVER\\%s\\%s\\%s*_%g.mp4', animalName, expName, expName, sidenum));
                vid = VideoReader(sprintf('%s\\%s', vidFile.folder, vidFile.name));
                nFrames = round(vid.FrameRate*vid.Duration);
                delete(vid);
                vtd = table(transpose(0:nFrames - 1), VariableNames={'FrameNumber'});
            end
        end

        function alignTimestamps(obj, varargin)
            p = inputParser();
            p.addParameter('refEventNameArduino', 'CUE_ON', @(x) ischar(x) || iscell(x))
            p.addParameter('refEventNameEphys', 'Cue', @(x) ischar(x) || iscell(x))
            p.addParameter('trialDurationTolerance', 0.1, @isnumeric)
            p.addParameter('shiftDurationTolerance', 1.5, @isnumeric)
            p.parse(varargin{:});
            refEventNameArduino = p.Results.refEventNameArduino;
            refEventNameEphys = p.Results.refEventNameEphys;
            trialDurationTolerance = p.Results.trialDurationTolerance;
            shiftDurationTolerance = p.Results.shiftDurationTolerance;
            % Use CUE_ON events because this is recorded in arduino and ephys
            if length(obj) == 1
                fprintf(1, '%s (%g units), alingning ephys & camera using %s and %s timestamps.\n', obj.name, length(obj.eu), string(refEventNameEphys).join, string(refEventNameArduino).join)

                % eventId = find(strcmp(obj.ac.EventMarkerNames, refEventNameArduino));
                % eventDateNum = obj.ac.EventMarkersUntrimmed(obj.ac.EventMarkersUntrimmed(:, 1) == eventId, 3)';
                % eventDateTime = datetime(eventDateNum, ConvertFrom='datenum', TimeZone='America/New_York');
                if ischar(refEventNameArduino)
                    eventDateTime = obj.ac.GetEventMarker(refEventNameArduino, 'datetime', Untrimmed=true);
                else
                    eventDateTime = [];
                    for iEvent = 1:length(refEventNameArduino)
                        eventDateTime = vertcat(eventDateTime, obj.ac.GetEventMarker(refEventNameArduino{iEvent}, 'datetime', Untrimmed=true));
                    end
                    eventDateTime = sort(eventDateTime, 'ascend');
                end

                if ~isempty(obj.vtdF)
                    fcamDateTime = datetime([obj.ac.Cameras(obj.getCameraIndex('f')).Camera.EventLog.Timestamp], ConvertFrom='datenum', TimeZone='America/New_York');
                    fcamFrameNum = [obj.ac.Cameras(obj.getCameraIndex('f')).Camera.EventLog.FrameNumber];
                end
                if ~isempty(obj.vtdL)
                    lcamDateTime = datetime([obj.ac.Cameras(obj.getCameraIndex('l')).Camera.EventLog.Timestamp], ConvertFrom='datenum', TimeZone='America/New_York');
                    lcamFrameNum = [obj.ac.Cameras(obj.getCameraIndex('l')).Camera.EventLog.FrameNumber];
                end
                if ~isempty(obj.vtdR)
                    sidenum = obj.getCameraIndex('r');
                    if sidenum ~= 0
                        rcamDateTime = datetime([obj.ac.Cameras(obj.getCameraIndex('r')).Camera.EventLog.Timestamp], ConvertFrom='datenum', TimeZone='America/New_York');
                        rcamFrameNum = [obj.ac.Cameras(obj.getCameraIndex('r')).Camera.EventLog.FrameNumber];
                    elseif isempty(obj.ac.Cameras)
                        rcamDateTime = datetime([obj.ac.Camera.EventLog.Timestamp], ConvertFrom='datenum', TimeZone='America/New_York');
                        rcamFrameNum = [obj.ac.Camera.EventLog.FrameNumber];
                    else
                        assert(length(obj.ac.Cameras) == 1)
                        rcamDateTime = datetime([obj.ac.Cameras(1).Camera.EventLog.Timestamp], ConvertFrom='datenum', TimeZone='America/New_York');
                        rcamFrameNum = [obj.ac.Cameras(1).Camera.EventLog.FrameNumber];
                    end
                end

                % Find event in ephystime
                if ischar(refEventNameEphys)
                    eventEphysTime = obj.eu(1).EventTimes.(refEventNameEphys);
                else
                    eventEphysTime = [];
                    for iEvent = 1:length(refEventNameEphys)
                        eventEphysTime = vertcat(eventEphysTime(:), obj.eu(1).EventTimes.(refEventNameEphys{iEvent})(:));
                    end
                    eventEphysTime = sort(eventEphysTime, 'ascend');
                end
                
                eventDateTime = eventDateTime(:);
                eventEphysTime = eventEphysTime(:);
                % eventDateTime = unique(eventDateTime(:));
                % eventEphysTime = unique(eventEphysTime(:));

                % Some assertions: 
                try
                    assert(length(eventEphysTime) == length(eventDateTime))
                catch
                    warning('Arduino has %g %s events, but ephys has %g %s events.', length(eventDateTime), string(refEventNameArduino), length(eventEphysTime), string(refEventNameEphys))
                    itiEphys = diff(eventEphysTime);
                    itiArduino = seconds(diff(eventDateTime));
                    n = length(eventDateTime) - length(eventEphysTime);
                    if n > 0
                        % case 1: remove first n arduino events
                        if all(abs(itiArduino(n+1:end) - itiEphys) < shiftDurationTolerance)
                            eventDateTime = eventDateTime(n+1:end);
                        % case 2: remove last n arduino events
                        elseif all(abs(itiArduino(1:end-n) - itiEphys) < shiftDurationTolerance)
                            eventDateTime = eventDateTime(1:end-n);
                        else
                            shifts = 0:n;
                            df = zeros(size(shifts));
                            for iShift = 1:length(shifts)
                                shift = shifts(iShift);
                                df(iShift) = mean(abs(itiArduino(1+shift:end-n+shift) - itiEphys));
                            end
                            [~, iShift] = min(df);
                            shift = shifts(iShift);
                            eventDateTime = eventDateTime(1+shift:end-n+shift);
                            if all(abs(itiArduino(1+shift:end-n+shift) - itiEphys) < shiftDurationTolerance)
                                fprintf(1, '\tShifting data by %i samples achieved a max distance of %g, which is below the min tolerance of %g.\n', shift, max(abs(itiArduino(1+shift:end-n+shift) - itiEphys)), shiftDurationTolerance)
                            else
                                warning('Shifting data by %i (n=%i) samples achieved a max distance of %g, which exceeded the min tolerance of %g.', shift, n, max(abs(itiArduino(1+shift:end-n+shift) - itiEphys)), shiftDurationTolerance)
                            end
                            % error('Arduino has %g ref events, but ephys has %g ref events.', length(eventDateTime), length(eventEphysTime));
                        end
                    else
                        n = -n;
                        if strcmpi(obj.name, 'desmond23_20220504')
                            eventEphysTime = eventEphysTime(n:end-1);
                        % case 1: remove first n ephys events
                        elseif all(abs(itiArduino - itiEphys(n+1:end)) < shiftDurationTolerance)
                            eventEphysTime = eventEphysTime(n+1:end);
                        % case 2: remove last n ephys events
                        elseif all(abs(itiArduino - itiEphys(1:end-n)) < shiftDurationTolerance)
                            eventEphysTime = eventEphysTime(1:end-n);
                        else
                            error('Arduino has %g ref events, but ephys has %g ref events.', length(eventDateTime), length(eventEphysTime));
                        end
                    end
                end
                assert(all(abs(diff(eventEphysTime(:)) - seconds(diff(eventDateTime(:)))) < trialDurationTolerance), 'Adruino trial lengths differe significantly from ephys, max different: %g.', max(abs(diff(eventEphysTime(:)) - seconds(diff(eventDateTime(:))))))

                fprintf(1, '\tInter-ref-intervals match between ephys and arduino for %g trials with a tolerance of %gs.\n', length(eventEphysTime), trialDurationTolerance);

                [uniqueEventDateTime, ia] = unique(eventDateTime);
                if length(uniqueEventDateTime) < length(eventEphysTime)
                    eventDateTime = uniqueEventDateTime;
                    eventEphysTime = eventEphysTime(ia);
                end
                assert(length(unique(eventDateTime))==length(unique(eventEphysTime)))

                % Clean up restarting framenums
                if ~isempty(obj.vtdF)
                    if nnz(fcamFrameNum == 0) > 1
                        iStart = find(fcamFrameNum == 0, 1, 'last');
                        fcamFrameNum = fcamFrameNum(iStart:end);
                        fcamDateTime = fcamDateTime(iStart:end);
                        fprintf('%s front camera had a restart. Only the last batch of framenumbers and timestamps are kept. %.2f seconds of data are useless.\n', obj.name, (iStart-1)*10/30)
                    end
                end
                if ~isempty(obj.vtdL)
                    if nnz(lcamFrameNum == 0) > 1
                        iStart = find(lcamFrameNum == 0, 1, 'last');
                        lcamFrameNum = lcamFrameNum(iStart:end);
                        lcamDateTime = lcamDateTime(iStart:end);
                        fprintf('%s left camera had a restart. Only the last batch of framenumbers and timestamps are kept. %.2f seconds of data are useless.\n', obj.name, (iStart-1)*10/30)
                    end
                end
                if ~isempty(obj.vtdR)
                    if nnz(rcamFrameNum == 0) > 1
                        iStart = find(rcamFrameNum == 0, 1, 'last');
                        rcamFrameNum = rcamFrameNum(iStart:end);
                        rcamDateTime = rcamDateTime(iStart:end);
                        fprintf('%s right camera had a restart. Only the last batch of framenumbers and timestamps are kept. %.2f seconds of data are useless.\n', obj.name, (iStart-1)*10/30)
                    end
                end

                if ~isempty(obj.vtdF)
                    assert(all(diff(fcamFrameNum) == 10))
                    fcamEphysTime = interp1(eventDateTime, eventEphysTime, fcamDateTime, 'linear', 'extrap');
                    fvtdEphysTime = interp1(fcamFrameNum, fcamEphysTime, obj.vtdF.FrameNumber, 'linear', 'extrap');
                    obj.vtdF.Timestamp = fvtdEphysTime;
                else
                    fcamEphysTime = NaN;
                end

                if ~isempty(obj.vtdL)
                    assert(all(diff(lcamFrameNum) == 10))
                    lcamEphysTime = interp1(eventDateTime, eventEphysTime, lcamDateTime, 'linear', 'extrap');
                    lvtdEphysTime = interp1(lcamFrameNum, lcamEphysTime, obj.vtdL.FrameNumber, 'linear', 'extrap');
                    obj.vtdL.Timestamp = lvtdEphysTime;
                else
                    lcamEphysTime = NaN;
                end

                if ~isempty(obj.vtdR)
                    assert(all(diff(rcamFrameNum) == 10))
                    rcamEphysTime = interp1(eventDateTime, eventEphysTime, rcamDateTime, 'linear', 'extrap');
                    rvtdEphysTime = interp1(rcamFrameNum, rcamEphysTime, obj.vtdR.FrameNumber, 'linear', 'extrap');
                    obj.vtdR.Timestamp = rvtdEphysTime;
                else
                    rcamEphysTime = NaN;
                end

                fprintf(1, '\tFrame 0 in ephys time: %.3f s, %.3f s, %.3f s\n', fcamEphysTime(1), lcamEphysTime(1), rcamEphysTime(1))
                
            else
                for i = 1:length(obj)
                    try
                        obj(i).alignTimestamps(varargin{:});
                    catch ME
                        fprintf(1, '%g: %s has error. %g EphysUnits involved.\n', i, obj(i).name, length(obj(i).eu))
                        warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
                    end
                end
            end
        end

        function i = getCameraIndex(obj, side)
            assert(length(obj) == 1)
            animalName = strsplit(obj.name, '_');
            animalName = animalName{1};
            switch lower(animalName)
                case {'daisy2', 'daisy3', 'daisy8', 'daisy9'}
                    switch lower(side)
                        case {'r', 'right'}
                            i = 0;
                        otherwise
                            i = -1;
                            warning('Unrecognized side string: ''%s'' for animal ''%s''', side, animalName)
                    end
                case {'daisy23', 'daisy24', 'daisy25', 'daisy26', 'daisy27', 'daisy28', 'daisy29', 'daisy30', 'daisy31', 'daisy32', 'daisy33', 'daisy34', 'desmond38', 'desmond39', 'desmond40', 'desmond41', 'desmond42', 'desmond43', 'desmond44', 'desmond45', 'desmond46', 'desmond47'}
                    switch lower(side)
                        case {'f', 'front'}
                            i = 2;
                        case {'l', 'left'}
                            i = 3;
                        case {'r', 'right'}
                            i = 1;
                        otherwise
                            i = -1;
                            warning('Unrecognized side string: ''%s'' for animal ''%s''', side, animalName)
                    end
                case {'desmond28', 'desmond29', 'desmond30'}
                    switch lower(side)
                        case {'f', 'front'}
                            i = 1;
                        case {'l', 'left'}
                            i = 2;
                        case {'r', 'right'}
                            i = 3;
                        otherwise
                            i = -1;
                            warning('Unrecognized side string: ''%s'' for animal ''%s''', side, animalName)
                    end
                case {'daisy14', 'daisy15', 'daisy16', 'desmond23', 'desmond24', 'desmond25', 'desmond26', 'desmond27'}
                    switch lower(side)
                        case {'l', 'left'}
                            i = 2;
                        case {'r', 'right'}
                            i = 1;
                        otherwise
                            i = -1;
                            warning('Unrecognized side string: ''%s'' for animal ''%s''', side, animalName)
                    end
                otherwise
                    error('Unrecognized animal name "%s"', animalName)
            end
        end

        function [x, y, l] = getTrajectory(obj, t, side, feature, varargin)
            p = inputParser();
            p.addRequired('t', @isnumeric)
            p.addRequired('side', @(x) ismember(x, {'l', 'r', 'f'}))
            p.addRequired('feature', @ischar)
            p.addParameter('likelihoodThreshold', 0, @isnumeric)
            p.addParameter('data', [], @istable)
            p.parse(t, side, feature, varargin{:})
            t = p.Results.t;
            side = p.Results.side;
            feature = p.Results.feature;
            likelihoodThreshold = p.Results.likelihoodThreshold;

            if isempty(p.Results.data)
                switch side
                    case 'l'
                        vtd = obj.vtdL;
                    case 'r'
                        vtd = obj.vtdR;
                    case 'f'
                        vtd = obj.vtdF;
                end
            else
                vtd = p.Results.data;
            end

            x = vtd.([feature, '_X']);
            y = vtd.([feature, '_Y']);
            l = vtd.([feature, '_Likelihood']);
            x(l < likelihoodThreshold) = NaN;
            y(l < likelihoodThreshold) = NaN;
            l(l < likelihoodThreshold) = NaN;
            x = interp1(vtd.Timestamp, x, t, 'linear');                        
            y = interp1(vtd.Timestamp, y, t, 'linear');                        
            l = interp1(vtd.Timestamp, l, t, 'linear');   
        end

        function [X, Y, L, t] = getTrajectoryByTrial(obj, side, feature, varargin)
            p = inputParser();
            p.addRequired('side', @(x) ismember(x, {'l', 'r', 'f'}))
            p.addRequired('feature', @ischar)
            p.addParameter('trialType', 'press', @(x) ismember(x, {'press', 'lick'}))
            p.addParameter('trials', [], @(x) isa(x, 'Trial'))
            p.addParameter('alignTo', 'stop', @(x) ismember(x, {'stop', 'start'}))
            p.addParameter('window', [-1, 1], @(x) isnumeric(x) && x(2) > x(1))
            p.addParameter('includeInvalid', false, @islogical)
            p.addParameter('resolution', 1/30, @isnumeric)
            p.addParameter('likelihoodThreshold', 0, @isnumeric)
            p.addParameter('data', [], @istable)
            p.addParameter('interp', 'linear', @ischar)
            p.parse(side, feature, varargin{:})
            side = p.Results.side;
            trialType = p.Results.trialType;
            trials = p.Results.trials;
            alignTo = p.Results.alignTo;
            window = p.Results.window;
            includeInvalid = p.Results.includeInvalid;
            resolution = p.Results.resolution;
            likelihoodThreshold = p.Results.likelihoodThreshold;
            interp = p.Results.interp;

            if isempty(trials)
                trials = obj.eu(1).getTrials(trialType);
            end

            if isempty(p.Results.data)
                switch side
                    case 'l'
                        vtd = obj.vtdL;
                    case 'r'
                        vtd = obj.vtdR;
                    case 'f'
                        vtd = obj.vtdF;
                end
            else
                vtd = p.Results.data;
            end
            t = window(1):resolution:window(2);
            X = NaN(length(trials), length(t));
            Y = X;
            L = X;

            switch alignTo
                case 'start'
                    alignToStart = true;
                case 'stop'             
                    alignToStart = false;
            end

            for i = 1:length(trials)
                if alignToStart
                    tt = trials(i).Start + t;                   
                    if includeInvalid
                        inTrial = true(size(t));
                    else
                        inTrial = tt <= trials(i).Stop;
                    end
                else            
                    tt = trials(i).Stop + t;
                    if includeInvalid
                        inTrial = true(size(t));
                    else
                        inTrial = tt >= trials(i).Start;
                    end
                end

                xx = vtd.([feature, '_X']);
                yy = vtd.([feature, '_Y']);
                ll = vtd.([feature, '_Likelihood']);
                xx(ll < likelihoodThreshold) = NaN;
                yy(ll < likelihoodThreshold) = NaN;
                ll(ll < likelihoodThreshold) = NaN;
                xx = interp1(vtd.Timestamp, xx, tt, interp);                        
                yy = interp1(vtd.Timestamp, yy, tt, interp);                        
                ll = interp1(vtd.Timestamp, ll, tt, interp);     
                X(i, inTrial) = xx(inTrial);
                Y(i, inTrial) = yy(inTrial);
                L(i, inTrial) = ll(inTrial);
            end
        end


        function [vtd, meta, motPos, trials] = readVideoTrackingDataShortConcatenated(obj, fdir, fname, iter, side, varargin)
            p = inputParser();
            p.addRequired('fdir', @isfolder);
            p.addRequired('fname', @ischar)
            p.addRequired('iter', @isnumeric);
            p.addRequired('side', @ischar);
            p.addParameter('assign', false, @islogical)
            p.parse(fdir, fname, iter, side, varargin{:})
            r = p.Results;


            metafile = dir(sprintf('%s\\%s_meta_iter%i.mat', r.fdir, r.fname, r.iter));
            meta = load(sprintf('%s\\%s', metafile.folder, metafile.name));

            files_vtd = sortrows(struct2table(dir(sprintf('%s\\*_iter%i.csv', r.fdir, r.iter))), 'datenum', 'descend');
            % Use the newest csv file generated by DeepLabCut if multiple matches
            % are found
            if height(files_vtd) > 1
                fname_vtd = sprintf('%s\\%s', files_vtd.folder{1}, files_vtd.name{1});
            else
                fname_vtd = sprintf('%s\\%s', files_vtd.folder, files_vtd.name);
            end
            opts = detectImportOptions(fname_vtd, 'NumHeaderLines', 3);
            opts.VariableNamesLine = 2;
            
            
            t_read = tic();
            fprintf(1, 'Reading short concatenated video tracking data from file %s...', fname_vtd);
            vtd = readtable(fname_vtd, opts);
            fprintf(1, '\nDone (%s).\n', seconds(toc(t_read)));
            
            % Set colnames
            vtd.Properties.VariableNames{1} = 'FrameNumber';
            w = length(vtd.Properties.VariableNames);
            for i = 2:w
                splitName = strsplit(vtd.Properties.VariableNames{i}, '_');
                if length(splitName) == 1
                    vtd.Properties.VariableNames{i} = [splitName{1}, '_X'];
                    spos = table2array(smoothdata(vtd(:, i:i+1), 'gaussian', 5));
                    vtd = addvars(vtd, spos(:, 1), spos(:, 2), table2array(vtd(:, i+2)), 'NewVariableNames', {[splitName{1}, '_Smooth_X'], [splitName{1}, '_Smooth_Y'], [splitName{1}, '_Smooth_Likelihood']});
                elseif splitName{2} == '1'
                    vtd.Properties.VariableNames{i} = [splitName{1}, '_Y'];
                elseif splitName{2} == '2'
                    vtd.Properties.VariableNames{i} = [splitName{1}, '_Likelihood'];
                end
            end

            vtd.Timestamp = meta.ephysTimes;
            vtdOld = vtd;
            vtd = cell(length(obj), 1);

            fprintf(1, 'Splitting into %i experiments...\n', length(obj))
            nFramesPerTrial = length(meta.ephysTimes) / length(meta.trials);
            assert(mod(nFramesPerTrial, 1) == 0)
            trialEdges = cumsum([0; meta.nTrialsInExp]);
            frameEdges = trialEdges * nFramesPerTrial;
            motPos = cell(length(obj), 1);
            trials = cell(length(obj), 1);
            for iObj = 1:length(obj)
                iExp = find(strcmpi(obj(iObj).name, meta.expNames));
                iTrialStart = trialEdges(iExp) + 1;
                iTrialEnd = trialEdges(iExp + 1);
                trials{iObj} = meta.trials(iTrialStart:iTrialEnd);
                motPos{iObj} = meta.motPos(iTrialStart:iTrialEnd);

                iFrameStart = frameEdges(iExp) + 1;
                iFrameStop = frameEdges(iExp + 1);
                if (r.assign)
                    obj(iObj).(sprintf('vtd%s', upper(side))) = vtdOld(iFrameStart:iFrameStop, :);
                end
                vtd{iObj} = vtdOld(iFrameStart:iFrameStop, :);
            end
        end        
    end
end