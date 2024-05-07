classdef ProbeMap
    properties
        site % Struct containing:
            % ml % middle of probe, i.e. craniotomy ML (negative for left, positive for right)
            % ap % middle of probe, i.e. craniotomy AP
            % dv % tip of probe (depth from dura)
            % duraOffset % dura depth relative to bregma, this will be added to ProbeMap.dv
        facing % {'front'/0, 'back'/180, 'left'/270, 'right'/90} (which way are the probe contacts/optrode optic fiber facing, from the animal's perspective)
        model % {'128DN'}
        map % final global coordinates for each channel taking probe rotation & implant coordinates/dura offset into account
        localMap % local coordinates for each channel (x, y, z) -> (horz, 0, vert)
    end

    methods
        function obj = ProbeMap(varargin)
            p = inputParser();
            p.addRequired('ml', @isnumeric);
            p.addRequired('ap', @isnumeric);
            p.addRequired('dv', @(x) isnumeric(x) && x < 0);
            p.addParameter('duraOffset', 0, @(x) isnumeric(x) && x <= 0);
            p.addParameter('facing', 'front', @(x) isnumeric(x) || ismember(x, {'front', 'back', 'left', 'right'})) % 0, 180, 270, 90
            p.addParameter('model', '128DN', @(x) ismember(x, {'128D', '128DN'}))
            p.parse(varargin{:})
            ml = p.Results.ml;
            ap = p.Results.ap;
            dv = p.Results.dv;
            duraOffset = p.Results.duraOffset;
            facing = p.Results.facing;
            model = p.Results.model;

            obj.site = struct(ml=ml, ap=ap, dv=dv, duraOffset=duraOffset);
            obj.facing = ProbeMap.parseFacing(facing);
            obj.model = model;
            obj.localMap = ProbeMap.getLocalMap(model);
            obj.map = ProbeMap.transformMap(obj.localMap, obj.facing, [obj.site.ml, obj.site.ap, obj.site.dv + obj.site.duraOffset]);
        end

        function side = getSide(obj)
            if obj.site.ml < 0
                side = 'l';
            elseif obj.site.ml > 0
                side = 'r';
            else
                error('Hemispheric side cannot be determined because ML=0.')
            end
        end
        
        function plot(obj, varargin)
            p = inputParser();
            if ~isempty(varargin) && isgraphics(varargin{1}) && isa(varargin{1}, 'Axes')
                p.addRequired('ax')
            end
            p.addOptional('sz', 25)
            p.addOptional('c', 'k')
            p.addParameter('autoMirror', true, @islogical)
            p.parse(varargin{:});
            if isfield(p.Results, 'ax')
                ax = p.Results.ax;
                fig = ax.Parent;
            else
                fig = figure();
                ax = axes(fig);
            end
            sz = p.Results.sz;
            c = p.Results.c;
            autoMirror = p.Results.autoMirror;

            x = obj.map.ml./1000;
            y = obj.map.ap./1000;
            if autoMirror
                x = abs(x);
            end
            z = obj.map.dv./1000;
            scatter3(ax, x, y, z, sz, c);
            text(ax, x, y, z, arrayfun(@(x) num2str(x), 1:height(obj.map), UniformOutput=false), HorizontalAlignment='center')
            if autoMirror && obj.getSide == 'l'
                xlabel(ax, 'ML (mirrored)')
            else
                xlabel(ax, 'ML')
            end
            ylabel(ax, 'AP')
            zlabel(ax, 'DV')
            axis(ax, 'equal')

            % select best view
            az = mod(obj.facing, 180); % prevent axis flipping
            if mod(obj.facing, 90) == 0
                el = 0;
            else
                el = 15;
            end
            view(az, el);

            xlim(ax, [min(x) - 0.1, max(x) + 0.1])
            ylim(ax, [min(y) - 0.1, max(y) + 0.1])
            zlim(ax, [min(z) - 0.1, max(z) + 0.1])
        end
    end

    methods (Static)
        function localMap = getLocalMap(model)
            switch upper(model)
                case '128D'
                    file = '128D_bottom.mat';
                case '128DN'
                    file = '128DN_bottom.mat';
                otherwise
                    error('Unknown probe model: ''%s''', model)
            end

            S = load(file);
            s = S.s;
            % z-correction: tip electrode is distance from tip of probe to bottom edge of the bottom-most electrode (we add it to s.z assuming electrode pads are infintely small, i.e. edge is center)
            % x-correction: substract median so 0 is halfway between the two center shanks
            localMap = table(s.x - median(s.x), s.y, s.z + s.tipelectrode, s.shaft, ...
                VariableNames={'x', 'y', 'z', 'shank'});
        end

        function map = transformMap(localMap, facing, varargin)
            p = inputParser();
            p.addRequired('localMap', @istable)
            p.addRequired('facing', @(x) isnumeric(x) || ismember(x, {'front', 'back', 'left', 'right'})) % 0, 180, 270, 90; think compass headings where forward is north.
            p.addOptional('offset', [0, 0, 0], @(x) length(x)==3 && isnumeric(x)) % ml, ap, dv offsets
            p.parse(localMap, facing, varargin{:})
            localMap = p.Results.localMap;
            facing = p.Results.facing;
            offset = p.Results.offset;

            facing = ProbeMap.parseFacing(facing);
            % How much clockwise rotation to apply to the probe
            
            theta = -(facing - 180)/180; % convert to counterclockwise rotation in radians, substract 180deg because default heading for localMap is 180
            R = [cospi(theta), -sinpi(theta); sinpi(theta), cospi(theta)]; % Counterclockwise rotation matrix R gives: new coords = R*[x; y]
            Rv = R*[localMap.x(:)'; localMap.y(:)'];
            map = table(Rv(1, :)' + offset(1), Rv(2, :)' + offset(2), localMap.z + offset(3), localMap.shank, VariableNames={'ml', 'ap', 'dv', 'shank'});
        end

        function facing = parseFacing(facing) % convert text to heading angle, and ensure value is between [0, 360)
            if ischar(facing)
                switch facing
                    case 'front'
                        facing = 0;
                    case 'back'
                        facing = 180;
                    case 'left'
                        facing = 270;
                    case 'right'
                        facing = 90;
                end
            end
            facing = mod(facing, 360);
            assert(isnumeric(facing) && facing < 360 && facing >= 0)
        end
    end
end