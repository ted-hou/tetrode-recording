classdef VideoDirReachTrial < Trial
    properties
        Side % Implant hemisphere, 3 lever targets are contra, 1 is ipsi
        TargetPosRaw % 1, 2, 3, 4 = (0, 0), (0, 1), (1, 0), (1, 1) on the 2-bit motpos digin lines
        TargetPosRelative %in
        TargetPosName
    end

    methods
        function obj = VideoDirReachTrial(exp, trialType)
            if nargin < 2
                trialType = 'PressSpontaneous';
            end
            if nargin == 0
                return
            end

            % Copy over press trial start/stop times from experiment
            assert(length(exp) == 1)
            assert(isa(exp, 'CompleteExperiment3'))

            trials = exp.eu(1).Trials.(trialType);
            nTrials = length(trials);

            if nTrials == 0
                obj = VideoDirReachTrial.empty;
                return
            end

            switch exp.animalName
                % Right implant, press with left arm, raw 1-4 are ipsi-front, contra-in, contra-front, contra-out
                case {'desmond28', 'desmond30'}
                    side = 'R';
                    posOrder = [4, 3, 2, 1];
                    targetPosNames = {'contra-out', 'contra-front', 'contra-in', 'ipsi-front'};
                % Left implant, press with right arm
                % contra-front, contra-out
                case 'desmond29'
                    side = 'L';
                    posOrder = [1, 2, 3, 4];
                    targetPosNames = {'contra-out', 'contra-front', 'contra-in', 'ipsi-front'};
                case {'daisy23', 'daisy24'}
                    side = 'R';
                    posOrder = [2, 1];
                    targetPosNames = {'contra-out', 'contra-in'};
                case 'daisy25'
                    side = 'L';
                    posOrder = [2, 1];
                    targetPosNames = {'contra-out', 'contra-in'};
                otherwise
                    error()
            end

            obj(nTrials) = VideoDirReachTrial();
            for i = 1:nTrials
                obj(i).Start = trials(i).Start;
                obj(i).Stop = trials(i).Stop;
            end

            % Embed mot pos information
            motPos = exp.eu(1).getMotorState([obj.Stop]);
            for i = 1:nTrials
                obj(i).Side = side;
                obj(i).TargetPosRaw = motPos(i);
                obj(i).TargetPosRelative = posOrder(motPos(i));
                obj(i).TargetPosName = targetPosNames{posOrder(motPos(i))};
            end
        end

        function pos = getTargetPos(obj, varargin)
            p = inputParser();
            p.addOptional('mode', 'name', @(x) ismember(lower(x), {'name', 'relative', 'raw'})) 
            % name: 'contra-out', 'contra-front', 'contra-in', 'ipsi-front'
            % relative: 1 2 3 4 (corresponding to 'contra-out', 'contra-front', 'contra-in', 'ipsi-front')
            % raw: 1 2 3 4 (corresponding to 'right-out', 'right-front', 'right-in', 'left-front')
            p.parse(varargin{:})
            mode = p.Results.mode;

            switch lower(mode)
                case 'name'
                    pos = string({obj.TargetPosName});
                case 'relative'
                    pos = [obj.TargetPosRelative];
                case 'raw'
                    pos = [obj.TargetPosRaw];
                otherwise
                    error('Unrecognized mode: ''%s''.', mode)
            end
        end
    end
end