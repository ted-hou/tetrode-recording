classdef CompleteExperiment3 < CompleteExperiment
    properties
        vtdF = []
    end

    methods
        function obj = CompleteExperiment3(varargin)
            if nargin == 0
                return
            end
            p = inputParser();
            p.addRequired('eu', @(x) isa(x, 'EphysUnit'))
            p.parse(varargin{:});
            eu = p.Results.eu;

            [uniqueExpNames, ~, expIndices] = unique({eu.ExpName});
            nExp = length(uniqueExpNames);
            obj(nExp) = CompleteExperiment3();
            for i = 1:nExp
                obj(i).name = uniqueExpNames{i};
                obj(i).eu = eu(expIndices==i);
                obj(i).vtdF = obj.readVideoTrackingData(obj(i).name, 'f');
                obj(i).vtdL = obj.readVideoTrackingData(obj(i).name, 'l');
                obj(i).vtdR = obj.readVideoTrackingData(obj(i).name, 'r');
                obj(i).ac = CompleteExperiment.readArduino(obj(i).name);
            end
        end

        function alignTimestamps(obj)
            % Use CUE_ON events because this is recorded in arduino and ephys
            if length(obj) == 1
                fprintf(1, '%s (%g units), alingning ephys & camera using CUE timestamps.\n', obj.name, length(obj.eu))

                eventId = find(strcmp(obj.ac.EventMarkerNames, 'CUE_ON'));
                eventDateNum = obj.ac.EventMarkersUntrimmed(obj.ac.EventMarkersUntrimmed(:, 1) == eventId, 3)';
                eventDateTime = datetime(eventDateNum, ConvertFrom='datenum', TimeZone='America/New_York');
                fcamDateTime = datetime([obj.ac.Cameras(1).Camera.EventLog.Timestamp], ConvertFrom='datenum', TimeZone='America/New_York');
                fcamFrameNum = [obj.ac.Cameras(1).Camera.EventLog.FrameNumber];
                lcamDateTime = datetime([obj.ac.Cameras(2).Camera.EventLog.Timestamp], ConvertFrom='datenum', TimeZone='America/New_York');
                lcamFrameNum = [obj.ac.Cameras(2).Camera.EventLog.FrameNumber];
                rcamDateTime = datetime([obj.ac.Cameras(3).Camera.EventLog.Timestamp], ConvertFrom='datenum', TimeZone='America/New_York');
                rcamFrameNum = [obj.ac.Cameras(3).Camera.EventLog.FrameNumber];

                % Find event in ephystime
                eventEhpysTime = obj.eu(1).EventTimes.Cue;
                
                % Some assertions: 
                try
                    assert(length(eventEhpysTime) == length(eventDateTime), 'Arduino has %g cue events, but ephys has %g cue events.', length(eventDateTime), length(eventEhpysTime))
                catch
                    % The last cue event might not have been saved to
                    % arduino, if trial was interrupted before reaching
                    % resultCode, and the user did not manually save after.
                    if length(eventEhpysTime) == length(eventDateTime) + 1
                        eventEhpysTime = eventEhpysTime(1:end-1);
                        fprintf(1, '\tArduino has %g cue events, but ephys has %g+1 cue events. Now we try to ignore the last ephys cue event.\n', length(eventDateTime), length(eventEhpysTime));
                    else
                        error('Arduino has %g cue events, but ephys has %g cue events.', length(eventDateTime), length(eventEhpysTime));
                    end
                end
                assert(all(abs(diff(eventEhpysTime) - seconds(diff(eventDateTime))) < 0.1), 'Adruino trial lengths differe significantly from ephys, max different: %g.', max(abs(diff(eventEhpysTime) - seconds(diff(eventDateTime)))))
                fprintf(1, '\tInter-cue-intervals match between ephys and arduino for %g trials with a tolerance of 0.1s.\n', length(eventEhpysTime));

                % Clean up restarting framenums
                if nnz(fcamFrameNum == 0) > 1
                    iStart = find(fcamFrameNum == 0, 1, 'last');
                    fcamFrameNum = fcamFrameNum(iStart:end);
                    fcamDateTime = fcamDateTime(iStart:end);
                    warning('%s front camera had a restart. Only the last batch of framenumbers and timestamps are kept. %.2f seconds of data are useless.', obj.name, (iStart-1)*10/30)
                end
                if nnz(lcamFrameNum == 0) > 1
                    iStart = find(lcamFrameNum == 0, 1, 'last');
                    lcamFrameNum = lcamFrameNum(iStart:end);
                    lcamDateTime = lcamDateTime(iStart:end);
                    warning('%s left camera had a restart. Only the last batch of framenumbers and timestamps are kept. %.2f seconds of data are useless.', obj.name, (iStart-1)*10/30)
                end
                if nnz(rcamFrameNum == 0) > 1
                    iStart = find(rcamFrameNum == 0, 1, 'last');
                    rcamFrameNum = rcamFrameNum(iStart:end);
                    rcamDateTime = rcamDateTime(iStart:end);
                    warning('%s right camera had a restart. Only the last batch of framenumbers and timestamps are kept. %.2f seconds of data are useless.', obj.name, (iStart-1)*10/30)
                end
                assert(all(diff(fcamFrameNum) == 10))
                assert(all(diff(lcamFrameNum) == 10))
                assert(all(diff(rcamFrameNum) == 10))

                fcamEphysTime = interp1(eventDateTime, eventEhpysTime, fcamDateTime, 'linear', 'extrap');
                lcamEphysTime = interp1(eventDateTime, eventEhpysTime, lcamDateTime, 'linear', 'extrap');
                rcamEphysTime = interp1(eventDateTime, eventEhpysTime, rcamDateTime, 'linear', 'extrap');

                fvtdEphysTime = interp1(fcamFrameNum, fcamEphysTime, obj.vtdF.FrameNumber, 'linear', 'extrap');
                lvtdEphysTime = interp1(lcamFrameNum, lcamEphysTime, obj.vtdL.FrameNumber, 'linear', 'extrap');
                rvtdEphysTime = interp1(rcamFrameNum, rcamEphysTime, obj.vtdR.FrameNumber, 'linear', 'extrap');

                obj.vtdF.Timestamp = fvtdEphysTime;
                obj.vtdL.Timestamp = lvtdEphysTime;
                obj.vtdR.Timestamp = rvtdEphysTime;

                fprintf(1, '\tFrame 0 in ephys time: %.3f s, %.3f s, %.3f s\n', fcamEphysTime(1), lcamEphysTime(1), rcamEphysTime(1))
                
            else
                for i = 1:length(obj)
                    try
                        obj(i).alignTimestamps();
                    catch ME
                        fprintf(1, '%g: %s has error. %g EphysUnits involved.\n', i, obj(i).name, length(obj(i).eu))
                        warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
                    end
                end
            end
        end

        function i = getCameraIndex(obj, side)
            switch lower(side)
                case {'f', 'front'}
                    i = 1;
                case {'l', 'left'}
                    i = 2;
                case {'r', 'right'}
                    i = 3;
                otherwise
                    error('Unrecognized side string: ''%s''', side)
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
            p.parse(side, feature, varargin{:})
            side = p.Results.side;
            trialType = p.Results.trialType;
            trials = p.Results.trials;
            alignTo = p.Results.alignTo;
            window = p.Results.window;
            includeInvalid = p.Results.includeInvalid;
            resolution = p.Results.resolution;
            likelihoodThreshold = p.Results.likelihoodThreshold;

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
                xx = interp1(vtd.Timestamp, xx, tt, 'linear');                        
                yy = interp1(vtd.Timestamp, yy, tt, 'linear');                        
                ll = interp1(vtd.Timestamp, ll, tt, 'linear');     
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