%% Load SNr_SCRetro
if ~exist('E:\Data\Units\TwoColor_SNr_SCRetro\SingleUnit_NonDuplicate_NonDrift_SNr', 'dir')
    eu = EphysUnit.load('C:\SERVER\Units\TwoColor_SNr_SCRetro\SingleUnit_NonDuplicate_NonDrift_SNr', waveforms=false, spikecounts=false, spikerates=false);
    metaPath = 'C:\SERVER\Units\meta_TwoColor_SNr_SCRetro_20260113.mat';
    metaSavePath = 'C:\\SERVER\\Units\\meta_TwoColor_SNr_SCRetro_%s.mat';    
else
    eu = EphysUnit.load('E:\Data\Units\TwoColor_SNr_SCRetro\SingleUnit_NonDuplicate_NonDrift_SNr', waveforms=false, spikecounts=false, spikerates=false);
    metaPath = 'E:\Data\Units\meta_TwoColor_SNr_SCRetro_20260113.mat';
    metaSavePath = 'E:\\Data\\Units\\meta_TwoColor_SNr_SCRetro_%s.mat';
end

% First time, generate meta and boostrap
% [SNr_SCRetro.eu, SNr_SCRetro.rd, SNr_SCRetro.eta, SNr_SCRetro.meta, SNr_SCRetro.p, SNr_SCRetro.c, SNr_SCRetro.boot] = read_SNr_SCRetro( ...
%     eu=eu, recalculateBootstrap=true, ...
%     metaSavePath=metaSavePath, ...
%     stimBluePowers=[100, 500, 2000]*1e-6, ...
%     stimRedPowers=[100, 500, 2000, 8000, 16000]*1e-6, ...
%     stimBluePowersMustContain=["*1e6==100", "*1e6==500", "*1e6==2000"], ...
%     stimRedPowersMustContain="*1e6>=2000", ...
%     groupRedPowersAbove=2000*1e-6 ...
% );

% Subsequent times, regenerate meta but load boostrap
% [SNr_SCRetro.eu, SNr_SCRetro.rd, SNr_SCRetro.eta, SNr_SCRetro.meta, SNr_SCRetro.p, SNr_SCRetro.c, SNr_SCRetro.boot] = read_SNr_SCRetro( ...
%     eu=eu, metaPath=metaPath, recalculateBootstrap=false, recalculateETA=true, ...
%     metaSavePath=metaSavePath, ...
%     stimBluePowers=[100, 500, 2000, 8000, 16000]*1e-6, ...
%     stimRedPowers=[100, 500, 2000, 8000, 16000]*1e-6, ...
%     stimBluePowersMustContain=["*1e6==100", "*1e6==500", "*1e6==2000"], ...
%     stimRedPowersMustContain="*1e6>=2000", ...
%     groupRedPowersAbove=2000*1e-6 ...
% );

% Just load
[SNr_SCRetro.eu, SNr_SCRetro.rd, SNr_SCRetro.eta, SNr_SCRetro.meta, SNr_SCRetro.p, SNr_SCRetro.c, SNr_SCRetro.boot] = read_SNr_SCRetro( ...
    eu=eu, metaPath=metaPath, recalculateBootstrap=false, recalculateETA=false, ...
    metaSavePath=metaSavePath, ...
    stimBluePowers=[100, 500, 2000, 8000, 16000]*1e-6, ...
    stimRedPowers=[100, 500, 2000, 8000, 16000]*1e-6, ...
    stimBluePowersMustContain=["*1e6==100", "*1e6==500", "*1e6==2000"], ...
    stimRedPowersMustContain="*1e6>=2000", ...
    groupRedPowersAbove=2000*1e-6 ...
);

clear metaPath metaSavePath

% Load SNr_SCRetro_ReverseInjection
if ~exist('E:\Data\Units\TwoColor_SNr_SCRetro\ReverseInjection\SingleUnit_NonDuplicate_NonDrift_SNr_withTrials', 'dir')
    euReverseInjection = EphysUnit.load('C:\SERVER\Units\TwoColor_SNr_SCRetro\ReverseInjection\SingleUnit_NonDuplicate_NonDrift_SNr_withTrials', waveforms=false, spikecounts=false, spikerates=false);
    metaPath = 'C:\SERVER\Units\meta_TwoColor_SNr_SCRetro_ReverseInjection_20260113.mat';
    metaSavePath = 'C:\\SERVER\\Units\\meta_TwoColor_SNr_SCRetro_ReverseInjection_%s.mat';    
else
    euReverseInjection = EphysUnit.load('E:\Data\Units\TwoColor_SNr_SCRetro\ReverseInjection\SingleUnit_NonDuplicate_NonDrift_SNr_withTrials', waveforms=false, spikecounts=false, spikerates=false);
    metaPath = 'E:\Data\Units\meta_TwoColor_SNr_SCRetro_ReverseInjection_20260113.mat';
    metaSavePath = 'E:\\Data\\Units\\meta_TwoColor_SNr_SCRetro_ReverseInjection_%s.mat';
end
% First time, generate meta and boostrap
% [SNr_SCRetro_ReverseInjection.eu, SNr_SCRetro_ReverseInjection.rd, SNr_SCRetro_ReverseInjection.eta, SNr_SCRetro_ReverseInjection.meta, SNr_SCRetro_ReverseInjection.p, SNr_SCRetro_ReverseInjection.c, SNr_SCRetro_ReverseInjection.boot] = read_SNr_SCRetro( ...
%     eu=euReverseInjection, recalculateBootstrap=true, ...
%     metaSavePath=metaSavePath, ...
%     stimBluePowers=[100, 500, 2000]*1e-6, ...
%     stimRedPowers=[100, 500, 2000, 8000, 16000]*1e-6, ...
%     stimBluePowersMustContain=["*1e6==100", "*1e6==500", "*1e6==2000"], ...
%     stimRedPowersMustContain="*1e6>=2000", ...
%     groupRedPowersAbove=2000*1e-6 ...
% );

% Subsequent times, regenerate meta but load boostrap
% [SNr_SCRetro_ReverseInjection.eu, SNr_SCRetro_ReverseInjection.rd, SNr_SCRetro_ReverseInjection.eta, SNr_SCRetro_ReverseInjection.meta, SNr_SCRetro_ReverseInjection.p, SNr_SCRetro_ReverseInjection.c, SNr_SCRetro_ReverseInjection.boot] = read_SNr_SCRetro( ...
%     eu=euReverseInjection, metaPath=metaPath, recalculateBootstrap=false, recalculateETA=true, ...
%     metaSavePath=metaSavePath, ...
%     stimBluePowers=[100, 500, 2000, 8000, 16000]*1e-6, ...
%     stimRedPowers=[100, 500, 2000, 8000, 16000]*1e-6, ...
%     stimBluePowersMustContain=["*1e6==100", "*1e6==500", "*1e6==2000"], ...
%     stimRedPowersMustContain="*1e6>=2000", ...
%     groupRedPowersAbove=2000*1e-6 ...
% );

% Just load
[SNr_SCRetro_ReverseInjection.eu, SNr_SCRetro_ReverseInjection.rd, SNr_SCRetro_ReverseInjection.eta, SNr_SCRetro_ReverseInjection.meta, SNr_SCRetro_ReverseInjection.p, SNr_SCRetro_ReverseInjection.c, SNr_SCRetro_ReverseInjection.boot] = read_SNr_SCRetro( ...
    eu=euReverseInjection, metaPath=metaPath, recalculateBootstrap=false, recalculateETA=false, ...
    metaSavePath=metaSavePath, ...
    stimBluePowers=[100, 500, 2000, 8000, 16000]*1e-6, ...
    stimRedPowers=[100, 500, 2000, 8000, 16000]*1e-6, ...
    stimBluePowersMustContain=["*1e6==100", "*1e6==500", "*1e6==2000"], ...
    stimRedPowersMustContain="*1e6>=2000", ...
    groupRedPowersAbove=2000*1e-6 ...
);

clear metaPath metaSavePath

%% Read channel map, filter out non-SNr neurons
% Probe is mounted along ML, the flat bit is facing forward, the dovetail
% is facing posterior of animal (if otherwise, map needs to be mirrored).
% Looking from behind the animal is equivalent to SpikeGLX view
% Map is looking from behind animal:
%   shank 1 - 4 are left-to-right.
%   col 1 - 2 are left-to-right, on same shank
%   row 1 - n are bottom to top
%   x is horizontal coord (um, 0 is center of all 4 shanks, negative is left)
%   y is vertical coord (um, 0 is tip of shank, negative is down)

dvDuraCorrection = -200;

[~, SNr_SCRetro.coords] = getNeuroPixelChannelMap(SNr_SCRetro.eu, ml=1300, ap=-3280, dv=-4700+dvDuraCorrection);
SNr_SCRetro.c.isSNr = SNr_SCRetro.coords(:, 2) <= -3.8*1e3;

[~, SNr_SCRetro_ReverseInjection.coords] = getNeuroPixelChannelMap(SNr_SCRetro_ReverseInjection.eu, ml=1300, ap=-3280, dv=-4700+dvDuraCorrection);
SNr_SCRetro_ReverseInjection.c.isSNr = SNr_SCRetro_ReverseInjection.coords(:, 2) <= -3.8*1e3;

%% S8. load data (slow)
% Load metadata
if exist('E:\Data\Units\meta_Lite_NonDuplicate_NonDrift.mat', 'file')
    eu = EphysUnit.load('E:\Data\Units\SNr_nonDuplicate_nonDrift_withITI', waveforms=false, spikecounts=false, spikerates=false);
    load('E:\Data\Units\meta_Lite_NonDuplicate_NonDrift.mat')
else
    eu = EphysUnit.load('C:\SERVER\Units\Lite_NonDuplicate_NonDrift', waveforms=false, spikecounts=false, spikerates=false);
    load('C:\SERVER\Units\meta_Lite_NonDuplicate_NonDrift.mat')
end

if exist('E:\Data\Units\PressVsLick_ArtifactsRemoved_Full\FixedEventsAndTrials', 'dir')
    euArtiFree = EphysUnit.load('E:\Data\Units\PressVsLick_ArtifactsRemoved_Full\FixedEventsAndTrials');
    load('E:\Data\Units\meta_PressVsLick_ArtifactsRemoved_Full_20260107.mat');
    etaArtiFree = metaArtiFree.eta;
else
    euArtiFree = EphysUnit.load('C:\SERVER\Units\PressVsLick_ArtifactsRemoved_Full\FixedEventsAndTrials');
    load('C:\SERVER\Units\meta_PressVsLick_ArtifactsRemoved_Full_20260107.mat');
    etaArtiFree = metaArtiFree.eta;
end

%% Calculate some metadata for lickVsReach_artifactRemoved
% Subselect euPos
[Lia, Locb] = ismember(euArtiFree.getName(), eu.getName());
assert(all(Lia))
euPosArtiFree = euPos(Locb, :);

% Calculate circlick amplitude (z-scored)
sdArtiFree = vertcat(etaArtiFree.pressNorm.stats.sd) ./ diff(etaArtiFree.pressNorm.t(1:2));

zAmpNormArtiFree = metaArtiFree.circlick.amp(:) ./ sdArtiFree(:);

%% Merge euPos for Intan and Neuropixel
euPosMerge = vertcat(euPos, SNr_SCRetro.coords);
additionalDVOffset = -100; % To account for additional dura/cortex that was scraped off. (i.e. probe was probably lower than recorded)
euPosMerge(:, 2) = euPosMerge(:, 2) + additionalDVOffset;

cMerge = struct();
for fn = ["isPressUp", "isPressDown", "isPressResponsive", "isLickUp", "isLickDown", "isLickResponsive"]
    cMerge.(fn) = vertcat(reshape(c.(fn), [], 1), reshape(SNr_SCRetro.c.(fn), [], 1));
end
cMerge.name = vertcat(string(reshape(eu.getName(), [], 1)), string(reshape(SNr_SCRetro.eu.getName(), [], 1)));
cMerge.hasPress = vertcat(c.hasPress(:), true(length(SNr_SCRetro.eu), 1));
cMerge.hasLick = vertcat(c.hasLick(:), true(length(SNr_SCRetro.eu), 1));
cMerge.hasPos = vertcat(c.hasPos(:), true(length(SNr_SCRetro.eu), 1));

metaMerge = struct();
for fn = ["press", "lick"]
    metaMerge.(fn) = [reshape(meta.(fn), [], 1); reshape(SNr_SCRetro.meta.(fn), [], 1)];
end

clear fn



%% S9. SC Reverse injection and Somatotopy (combined Neuropixel and Intan)
close all

layout.w = 5/3*4;
layout.h = 4;
fig = figure(Units='inches', Position=[0, 0, layout.w, layout.h], DefaultAxesFontSize=p.fontSize);
layout.tl = tiledlayout(fig, 2, 3, TileSpacing='compact', Padding='loose');

TITLE = ["all", "reach-modulated", "lick-modulated"; ...
    "reach-dec, lick-inc", "reach-inc, lick-dec", "lick-entrained"];

selCommon = cMerge.hasPress & cMerge.hasLick & cMerge.hasPos & euPosMerge(:, 2) < -3800;
SEL = { ...
    selCommon, selCommon & cMerge.isPressResponsive, selCommon & cMerge.isLickResponsive; ...
    selCommon & cMerge.isPressDown & cMerge.isLickUp, selCommon & cMerge.isPressUp & cMerge.isLickDown, metaArtiFree.cc.isLick(:) & metaArtiFree.cc.isIntan(:), ...
    };
SELBACKGROUND = { ...
    selCommon, selCommon, selCommon; ...
    selCommon, selCommon, metaArtiFree.cc.isIntan(:), ...
    };
STATS = { ...
    repmat(0.5, size(SEL{2})), metaMerge.press, metaMerge.lick; ...
    repmat(0.5, size(SEL{4})), repmat(0.5, size(SEL{5})), zAmpNormArtiFree, ...
    };
SRANGE = { ...
    [0, 7], [0, 7], [0, 7]; ...
    [0, 7], [0, 7], [0, 3], ...
    };
COLOR = { ...
    [0.15, 0.15, 0.15], [], []; ...
    [0.15, 0.15, 0.15], [0.15, 0.15, 0.15], [0.15, 0.15, 0.15], ...
    };
POS = { ...
    euPosMerge, euPosMerge, euPosMerge; ...
    euPosMerge, euPosMerge, euPosArtiFree, ...
    };

ALPHA = repmat(0.25, 2, 3);

AX = gobjects(2, 3);
LETTERS = 'abcdef';

for i = 1:2
    for j = 1:3
        ax = nexttile(layout.tl, 3*(i-1) + j);
        sel = SEL{i, j};
        selBackground = SELBACKGROUND{i, j};
        coords = POS{i, j}(sel, :);
        stats = STATS{i, j}(sel);

        hold(ax, 'on')
        AcuteRecording.plotMap(ax, coords, stats, SRANGE{i, j}, 0, UseSignedML=false, BubbleSize=[1, 5], MarkerAlpha=ALPHA(i, j), ...
            MarkerEdgeAlpha=0.8, Color=COLOR{i, j});
        % plot(snrPoly, FaceAlpha=0, EdgeColor='black', LineStyle='--', LineWidth=1)
    
        if i == 1 && j == 2
            hDummy = gobjects(2, 1);
            hDummy(1) = scatter(ax, 0, 0, 1, [1, 0, 0], 'filled', DisplayName='inc');
            hDummy(2) = scatter(ax, 0, 0, 1, [0, 0, 1], 'filled', DisplayName='dec');            
        end

        title(ax, sprintf("%s\n(n=%i/%i)", TITLE(i, j), nnz(sel), nnz(selBackground)))
        axis(ax, 'image')
        xlim(ax, [0.7, 2])
        ylim(ax, [-5.1, -3.8])
        xticks(ax, [0.9, 1.3, 1.7])
        yticks(ax, [-4.9, -3.9])
    
        xlabel(ax, 'ML')
        if j == 1
            ylabel(ax, 'DV')
        else
            ylabel(ax, '')
        end
    
        fontsize(ax, p.fontSize, 'points');
        fontname(ax, 'Arial')
        fprintf('%i/%i %s modulated units\n', nnz(sel), nnz(selBackground), TITLE{i, j})

        hLetter = text(ax, 0, 0, LETTERS(3*(i-1)+j), FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
        ax.Units = 'inches';
        hLetter.HorizontalAlignment = 'right';
        hLetter.VerticalAlignment = 'top';
        hLetter.Position = [-0.4, ax.Position(4) + 0.3, 0];
    end
end

lgd = legend(hDummy, Orientation='horizontal');
fontsize(lgd, p.fontSize, 'points');
fontname(lgd, 'Arial')
lgd.Layout.Tile = 'north';
lgd.Position(1) = lgd.Position(1) + 4;

euMerge = [eu, SNr_SCRetro.eu];
nSessions = length(unique({euMerge(selCommon).ExpName}));
nAnimals = length(unique(euMerge(selCommon).getAnimalName()));
clear euMerge
fprintf('nAnimals=%i, nSessions=%i\n', nAnimals, nSessions)


nSessions = length(unique({euArtiFree(metaArtiFree.cc.isIntan).ExpName}));
nAnimals = length(unique(euArtiFree(metaArtiFree.cc.isIntan).getAnimalName()));
fprintf('Lick-entrained: nAnimals=%i, nSessions=%i\n', nAnimals, nSessions)

copygraphics(fig, BackgroundColor='none', ContentType='vector')


% %% Draw an SNr ROI based on paxinos
% img = imread("C:\GIT\tetrode-recording\Publication\SNr_Paxinos.jpg");
% img = imresize(img, [3000, 2000]);
% image(img);
% axis('image');
% 
% h = drawpolygon;
% %%
% posIJ = h.Position;
% posML = (h.Position(:, 1) - 0) * (2000-0) / (2000-0) + 0;
% posDV = (h.Position(:, 2) - 0) * (-6000-(-3000)) / (3000-0) + (-3000);
% posAP = -3280;
% pos = [posML, posDV];
% 
% save('E:\Data\SNr_Paxinos_3280.mat', 'img', 'posIJ', 'posML', 'posDV', 'posAP', 'pos')
% 
% % (x, y) coordinates. 
% % x: (0, 2000) -> ml(0, 2000) um
% % y: (0, 3000) -> dv(-3000, -6000) um