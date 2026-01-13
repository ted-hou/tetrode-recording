function [map, coords] = getNeuroPixelChannelMap(eu, varargin)
%% Read channel map
% Probe is mounted along ML, the flat bit is facing forward, the dovetail
% is facing posterior of animal (if otherwise, map needs to be mirrored).
% Looking from behind the animal is equivalent to SpikeGLX view
% Map is looking from behind animal:
%   shank 1 - 4 are left-to-right.
%   col 1 - 2 are left-to-right, on same shank
%   row 1 - n are bottom to top
%   x is horizontal coord (um, 0 is center of all 4 shanks, negative is left)
%   y is vertical coord (um, 0 is tip of shank, negative is down)
p = inputParser();
p.addRequired('eu', @(x) isa(x, 'EphysUnit'))
p.addParameter('ml', 1300, @(x) isnumeric(x) && x>=0); % Absolute value of ML (left/right is determined by animal name)
p.addParameter('ap', -3280, @isnumeric);
p.addParameter('dv', -4700, @isnumeric); % from dura (which is usually 200um below bregma for SNr)
p.parse(eu, varargin{:})
eu = p.Results.eu;
ml = p.Results.ml;
ap = p.Results.ap;
dv = p.Results.dv;


tipLength = 175; % um
shankPitch = 250;
columnPitch = 15; % vertical distance between electrodes
rowPitch = 32; % horizontal distance between electrodes

if ~exist("D:\Data\tips_out.imro", 'file')
    file = "C:\SERVER\tips_out.imro";
else
    file = "E:\Data\tips_out.imro";
end
mapStr = string(fileread(file));
tokens = regexp(mapStr, '(?<=\()[^)]*(?=\))', 'match');
assert(tokens(1) == "2013,384")
tokens = tokens(2:end);
clear map
map(length(tokens)) = struct(channel=[], shank=[], col=[], row=[], x=[], y=[]);
for i = 1:length(tokens)
    subtokens = strsplit(tokens(i), " ");
    map(i).channel = str2double(subtokens(1)) + 1;
    assert(map(i).channel == i)
    map(i).shank = str2double(subtokens(2)) + 1; %
    idOnShank = str2double(subtokens(5));
    map(i).col = 1 + (mod(idOnShank, 2) == 1);
    switch map(i).col
        case 1
            map(i).row = idOnShank/2 + 1;
        case 2
            map(i).row = (idOnShank-1)/2 + 1;
    end
    map(i).x = (map(i).shank-1)*shankPitch - 1.5*shankPitch;
    switch map(i).col
        case 1
            map(i).x = map(i).x - rowPitch/2;
        case 2
            map(i).x = map(i).x + rowPitch/2;
    end
    map(i).y = tipLength + (map(i).row-1)*columnPitch;
end

% Plot map
% assert(all(ismember(eu.getAnimalName(), {'daisy27', 'daisy28', 'desmond39'}))) % 'daisy27/28, desmond39 were recorded from left hemisphere (SNr/SC)'
% ml = -1.3*1e3;
% ap = -3.28*1e3;
% dv = -4.7*1e3; % From dura (which is about -0.2mm below bregma here)
coordsLocal = [vertcat(map.x), vertcat(map.y), zeros(size(vertcat(map.x)))];
coords = zeros(length(eu), 3);
for iEu = 1:length(eu)
    chn = eu(iEu).Channel;
    switch eu(iEu).getAnimalName()
        case {'daisy27', 'daisy28', 'desmond39'}
            coords(iEu, :) = [map(chn).x - ml, map(chn).y + dv, ap];
        case {'daisy26', 'desmond38'}
            coords(iEu, :) = [map(chn).x + ml, map(chn).y + dv, ap];
        otherwise
            error('Unrecognized animal name %s', eu(iEu).getAnimalName())
    end
end
