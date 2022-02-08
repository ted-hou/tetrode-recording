classdef Trial < handle
    properties
        Start = []
        Stop = []
    end
    
    % public methods
    methods
        function obj = Trial(start, stop, varargin)
            if nargin == 0
                return
            end
            
            p = inputParser();
            p.addRequired('start', @isnumeric); % Trial start events
            p.addRequired('stop', @isnumeric); % Trial end events
            p.addOptional('stopMode', 'first', @(x) ischar(x) && ismember(lower(x), {'first', 'last'}))
            p.addOptional('exclude', [], @isnumeric); % Trial exclusion events (trial invalid if an exclusion event is detected between start and stop)
            p.parse(start, stop, varargin{:})
            start = p.Results.start;
            stop = p.Results.stop;
            stopMode = p.Results.stopMode;
            exclude = p.Results.exclude;
            
            if isempty(start) || isempty(stop)
                obj(1) = [];
                return
            end
            
            [~, start, stop] = Trial.findEdges(start, stop, stopMode, exclude);
            nTrials = length(start);
            
            if nTrials == 0
                obj(1) = [];
                return
            end
            
            obj(nTrials) = Trial();
            for i = 1:nTrials
                obj(i).Start = start(i);
                obj(i).Stop = stop(i);
            end
        end
        
        function l = duration(obj)
            l = [obj.Stop] - [obj.Start];
        end
    end

    % public static methods
    methods (Static)
        function varargout = findEdges(start, stop, varargin)
            % TRIAL.FINDEDGES find edges of trials
            %   edges = TRIAL.FINDEDGES(start, stop) finds trial edges from a list of trial start and trial stop event timestamps.
            %   edges = TRIAL.FINDEDGES(start, stop, exclude) would not include trials if an 'exclude' event was detected between start and stop.
            %   edges = TRIAL.FINDEDGES(start, stop, exclude, 'first') uses the first stop event after start as end of trial
            %   edges = TRIAL.FINDEDGES(start, stop, exclude, 'last') uses the last stop event before the next start event as end of trial
            %   [edges, start, stop] = TRIAL.FINDEDGES(start, stop, _) also returns the updated trial start, stop timestamps, taking trial-exclusion into account
            p = inputParser();
            p.addRequired('start', @isnumeric);
            p.addRequired('stop', @isnumeric);
            p.addOptional('stopMode', 'first', @(x) ischar(x) && ismember(lower(x), {'first', 'last'}))
            p.addOptional('exclude', [], @isnumeric);
            p.parse(start, stop, varargin{:});
            start = p.Results.start;
            stop = p.Results.stop;
            stopMode = p.Results.stopMode;
            exclude = p.Results.exclude;
            
            if isempty(start) || isempty(stop)
                % warning('Could not find trial edges, empty array provided as start/stop events.')
                varargout = {[], [], []};
                return
            end
            
            % Make sure to use row vectors for event times, col vectors for
            % edges
            start = reshape(start, 1, []);
            stop = reshape(stop, 1, []);
            exclude = reshape(exclude, 1, []);
            
            j = 0;
            while nnz(diff(start) <= 0) > 0
                start = start(2:end);
                j = j + 1;
            end
            if j > 0
                warning('Skipped %i start events', j)
            end
            
            % Find trial edges using histcounts()
            edges = reshape([start(1:end-1), max(start(end), stop(end))], [], 1);
            
            [~, ~, bins] = histcounts(stop, edges);
            stop = stop(bins > 0); % Remove 0th bin, stop event cannot occur befor start
            bins = bins(bins > 0);
            
            % Update trial edges using 'first' or 'last' stop event in trial
            [iStart, iStop] = unique(bins, stopMode);
            start = start(iStart); % This exclude trials where there are no stop events
            stop = stop(iStop);
            edges = reshape([start; stop], [], 1);
            
            % Exclusion
            if ~isempty(exclude)
                [~, ~, bins] = histcounts(exclude, edges);
                inTrial = rem(bins, 2) ~= 0; % Oddbins represent an exclusion event in trial
                shouldExclude = (unique(bins(inTrial)) + 1) / 2;
                start(shouldExclude) = [];
                stop(shouldExclude) = [];
                edges = reshape([start; stop], [], 1);
            end
            
            varargout = {edges, start, stop};
        end
        
        function varargout = extendEdges(varargin)
            p = inputParser();
            if isnumeric(varargin{1}) && isnumeric(varargin{2})
                p.addRequired('start', @isnumeric)
                p.addRequired('stop', @isnumeric)
            elseif isnumeric(varargin{1})
                p.addRequired('edges', @isnumeric)
            else
                error('The first (or the first two) input(s) are expected to be numeric.')
            end
            p.addParameter('window', [-1, 1], @(x) isnumeric(x) && length(x) >= 2)
            p.addParameter('onInvalid', 'shrink', @(x) ischar(x) && ismember(lower(x), {'shrink', 'cull'}))
            p.addParameter('shrinkPriority', 'prewindow', @(x) ischar(x) && ismember(lower(x), {'prewindow', 'postwindow'}))
            p.parse(varargin{:})
            
            if isfield(p.Results, 'start')
                start = reshape(p.Results.start, 1, []);
                stop = reshape(p.Results.stop, 1, []);
                edges = reshape([start; stop], [], 1);
            else
                edges = reshape(p.Results.edges, [], 1);
                reshapedEdges = reshape(edges, 2, []);
                start = reshapedEdges(1, :);
                stop = reshapedEdges(2, :);
            end
            window = p.Results.window;
            onInvalid = lower(p.Results.onInvalid);
            shrinkPriority = lower(p.Results.shrinkPriority);
            
            % Validate original edges
            if sum(diff(edges) <= 0) > 0
                error('Input edges must be monotonic increasing.')
            end
            
            if isempty(nonzeros(window))
                varargout = {edges, [0, 0]};
                return
            end
            
            extendedStart = start + window(1);
            extendedStop = stop + window(2);
            extendedEdges = reshape([extendedStart; extendedStop], [], 1);
            
            culled = [];
            
            % Validate extended edges (by recursively shrinking the window, or by culling invalid trials)
            if sum(diff(extendedEdges) <= 0) > 0
                switch onInvalid
                    case 'shrink'
                        % Can shrink either way
                        if window(1) < -1e-3 && window(2) > 1e-3
                            switch shrinkPriority
                                case 'prewindow'
                                    window(1) = window(1) * 0.5;
                                case 'postwindow'
                                    window(2) = window(2) * 0.5;
                            end
                        % Can only shrink post window
                        elseif window(2) > 1e-3
                            window(1) = 0;
                            window(2) = window(2) * 0.5;
                        % Can only shrink pre window
                        elseif window(1) < -1e-3
                            window(1) = window(1) * 0.5;
                            window(2) = 0;
                        % Cannot shrink either
                        else
                            window = [0, 0];
                        end
                        [extendedEdges, window] = Trial.extendEdges(edges, 'window', window, 'onInvalid', 'shrink', 'shrinkPriority', shrinkPriority);                            
                    case 'cull'
                        assert(sum(diff(start) <= 0) == 0, 'start event times are not monotonic increasing (%i violatiions).', sum(diff(start) <= 0))
                        assert(sum(diff(stop) <= 0) == 0, 'stop event times are not monotonic increasing (%i violatiions).', sum(diff(stop) <= 0))
                        assert(sum(stop <= start) == 0, 'stop events must occur after corresponding start events (%i violatiions).', sum(stop <= start))
                        culled = [extendedStop(1:end-1) >= extendedStart(2:end), false];
                        extendedStop(culled) = [];
                        extendedStart(culled) = [];
                        extendedEdges = reshape([extendedStart; extendedStop], [], 1);
                end
            end
            
            varargout = {extendedEdges, window, culled};
        end
    end
    
    % private methods
    methods (Access = {})
    end
end