classdef CompleteExperiment2 < CompleteExperiment
    % camera 2, animal's right
    % camera 3, animal's left

    methods
        function obj = CompleteExperiment2(varargin)
            if nargin == 0
                return
            end
            p = inputParser();
            p.addRequired('eu', @(x) isa(x, 'EphysUnit'))
            p.parse(varargin{:});
            eu = p.Results.eu;

            [uniqueExpNames, ~, expIndices] = unique({eu.ExpName});
            nExp = length(uniqueExpNames);
            obj(nExp) = CompleteExperiment2();
            for i = 1:nExp
                obj(i).name = uniqueExpNames{i};
                obj(i).eu = eu(expIndices==i);
                obj(i).vtdL = obj.readVideoTrackingData(obj(i).name, 'l');
                obj(i).vtdR = obj.readVideoTrackingData(obj(i).name, 'r');
                obj(i).ac = CompleteExperiment.readArduino(obj(i).name);
            end
        end        

        function i = getCameraIndex(obj, side)
            switch lower(side)
                case {'f', 'front'}
                    i = 1;
                case {'l', 'left'}
                    i = 3;
                case {'r', 'right'}
                    i = 2;
                otherwise
                    error('Unrecognized side string: ''%s''', side)
            end
        end

        function alignTimestamps(obj)
            % Use CUE_ON events because this is recorded in arduino and ephys
            if length(obj) == 1
                fprintf(1, '%s (%g units), alingning ephys & camera using CUE timestamps.\n', obj.name, length(obj.eu))

                eventId = find(strcmp(obj.ac.EventMarkerNames, 'TRIAL_START'));
                eventDateNum = obj.ac.EventMarkersUntrimmed(obj.ac.EventMarkersUntrimmed(:, 1) == eventId, 3)';
                eventDateTime = datetime(eventDateNum, ConvertFrom='datenum', TimeZone='America/New_York');
                lcamDateTime = datetime([obj.ac.Cameras(3).Camera.EventLog.Timestamp], ConvertFrom='datenum', TimeZone='America/New_York');
                lcamFrameNum = [obj.ac.Cameras(3).Camera.EventLog.FrameNumber];
                rcamDateTime = datetime([obj.ac.Cameras(2).Camera.EventLog.Timestamp], ConvertFrom='datenum', TimeZone='America/New_York');
                rcamFrameNum = [obj.ac.Cameras(2).Camera.EventLog.FrameNumber];

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
                assert(all(diff(lcamFrameNum) == 10))
                assert(all(diff(rcamFrameNum) == 10))

                lcamEphysTime = interp1(eventDateTime, eventEhpysTime, lcamDateTime, 'linear', 'extrap');
                rcamEphysTime = interp1(eventDateTime, eventEhpysTime, rcamDateTime, 'linear', 'extrap');

                lvtdEphysTime = interp1(lcamFrameNum, lcamEphysTime, obj.vtdL.FrameNumber, 'linear', 'extrap');
                rvtdEphysTime = interp1(rcamFrameNum, rcamEphysTime, obj.vtdR.FrameNumber, 'linear', 'extrap');

                obj.vtdL.Timestamp = lvtdEphysTime;
                obj.vtdR.Timestamp = rvtdEphysTime;

                fprintf(1, '\tFrame 0 in ephys time: %.3f s, %.3f s\n', lcamEphysTime(1), rcamEphysTime(1))
                
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
    end
end