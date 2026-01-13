read_SC_opto_trajectories_DLC;

%% Load SNr_SCRetro
if ~exist('E:\Data\Units\TwoColor_SNr_SCRetro\SingleUnit_NonDuplicate_NonDrift_SNr', 'dir')
    eu = EphysUnit.load('C:\SERVER\Units\TwoColor_SNr_SCRetro\SingleUnit_NonDuplicate_NonDrift_SNr', waveforms=false, spikecounts=false, spikerates=false);
    metaPath = 'C:\SERVER\Units\meta_TwoColor_SNr_SCRetro_20260112.mat';
    metaSavePath = 'C:\\SERVER\\Units\\meta_TwoColor_SNr_SCRetro_%s.mat';    
else
    eu = EphysUnit.load('E:\Data\Units\TwoColor_SNr_SCRetro\SingleUnit_NonDuplicate_NonDrift_SNr', waveforms=false, spikecounts=false, spikerates=false);
    metaPath = 'E:\Data\Units\meta_TwoColor_SNr_SCRetro_20260112.mat';
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
[SNr_SCRetro.eu, SNr_SCRetro.rd, SNr_SCRetro.eta, SNr_SCRetro.meta, SNr_SCRetro.p, SNr_SCRetro.c, SNr_SCRetro.boot] = read_SNr_SCRetro( ...
    eu=eu, metaPath=metaPath, recalculateBootstrap=false, recalculateETA=true, ...
    metaSavePath=metaSavePath, ...
    stimBluePowers=[100, 500, 2000, 8000, 16000]*1e-6, ...
    stimRedPowers=[100, 500, 2000, 8000, 16000]*1e-6, ...
    stimBluePowersMustContain=["*1e6==100", "*1e6==500", "*1e6==2000"], ...
    stimRedPowersMustContain="*1e6>=2000", ...
    groupRedPowersAbove=2000*1e-6 ...
);

clear metaPath metaSavePath

%% Load SNr_SCRetro_ReverseInjection
if ~exist('E:\Data\Units\TwoColor_SNr_SCRetro\ReverseInjection\SingleUnit_NonDuplicate_NonDrift_SNr_withTrials', 'dir')
    euReverseInjection = EphysUnit.load('C:\SERVER\Units\TwoColor_SNr_SCRetro\ReverseInjection\SingleUnit_NonDuplicate_NonDrift_SNr_withTrials', waveforms=false, spikecounts=false, spikerates=false);
    metaPath = 'C:\SERVER\Units\meta_TwoColor_SNr_SCRetro_ReverseInjection_20260112.mat';
    metaSavePath = 'C:\\SERVER\\Units\\meta_TwoColor_SNr_SCRetro_ReverseInjection_%s.mat';    
else
    euReverseInjection = EphysUnit.load('E:\Data\Units\TwoColor_SNr_SCRetro\ReverseInjection\SingleUnit_NonDuplicate_NonDrift_SNr_withTrials', waveforms=false, spikecounts=false, spikerates=false);
    metaPath = 'E:\Data\Units\meta_TwoColor_SNr_SCRetro_ReverseInjection_20260112.mat';
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
[SNr_SCRetro_ReverseInjection.eu, SNr_SCRetro_ReverseInjection.rd, SNr_SCRetro_ReverseInjection.eta, SNr_SCRetro_ReverseInjection.meta, SNr_SCRetro_ReverseInjection.p, SNr_SCRetro_ReverseInjection.c, SNr_SCRetro_ReverseInjection.boot] = read_SNr_SCRetro( ...
    eu=euReverseInjection, metaPath=metaPath, recalculateBootstrap=false, recalculateETA=true, ...
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

[~, SNr_SCRetro.coords] = getNeuroPixelChannelMap(SNr_SCRetro.eu, ml=1300, ap=-3280, dv=-4700);
SNr_SCRetro.c.isSNr = SNr_SCRetro.coords(:, 2) <= -3.8*1e3;

[~, SNr_SCRetro_ReverseInjection.coords] = getNeuroPixelChannelMap(SNr_SCRetro_ReverseInjection.eu, ml=1300, ap=-3280, dv=-4700);
SNr_SCRetro_ReverseInjection.c.isSNr = SNr_SCRetro_ReverseInjection.coords(:, 2) <= -3.8*1e3;


%% Fig 8b SC Stim causes movements (medial SC stim vs. lateral SC stim)
close all

useYYAxis = false;
windowPreStim = [-0.25, 0];
windowPostStim = [0, 0.45];
WAVELENGTHS = [470, 635];
COLORS = ["blue", "red"];
BODYPARTS = ["HandCameraSide", "Jaw"];
BODYPARTDISPNAMES = ["Forepaw", "Jaw"];
YYAXIS = ["left", "right"];
BODYPARTCOLORS = arrayfun(@(i) getColor(i, 3, 0.7), [1, 3], UniformOutput=false);
% YLIMS = {[-2, 10], [-1, 5]};
YLIMS = {[-1, 5], [-1, 5]};
% YTICKS = {[0, 6], [0, 3]};
YTICKS = {[0, 3], [0, 3]};
p.fontSize = 9;

iExp = length(trajectoriesSC);
fig = figure(Units="inches", Position=[1, 1, 3.5, 1.75]);
tl = tiledlayout(fig, 1, 2);
h = gobjects(2, 1);
AX = gobjects(length(BODYPARTS), 1);
for iColor = 1:length(COLORS)
    color = COLORS(iColor);
    switch color
        case "blue"
            mwPower = pSC.stimBluePowers*1e3;
        case "red"
            mwPower = pSC.stimRedPowers*1e3;
    end

    ax = nexttile(tl);
    AX(iColor) = ax;
    hold(ax, 'on')
    colororder(ax, getColor([1, 3], 3, 0.7))
    for iBodypart = 1:length(BODYPARTS)
        if useYYAxis
            yyaxis(ax, YYAXIS(iBodypart));
        end
        bodypart = BODYPARTS(iBodypart);
        t = trajectoriesSC(iExp).(bodypart).(color).t;
        X = trajectoriesSC(iExp).(bodypart).(color).X;
        Y = -trajectoriesSC(iExp).(bodypart).(color).Y;
        nTrials = size(X, 1);
        velX = diff(X, 1, 2)./diff(t);
        velX = [NaN([size(velX, 1), 1]), velX];
        velY = diff(Y, 1, 2)./diff(t);
        velY = [NaN([size(velY, 1), 1]), velY];
        spd = sqrt(velX.^2 + velY.^2);
        spd = (spd - mean(spd(:, t<0), 'all', 'omitnan')) ./ std(spd(:, t<0), 0, 'all', 'omitnan');
        mu = mean(spd, 1, 'omitnan');
        sd = std(spd, 0, 1, 'omitnan');
        mu(isnan(mu)) = 0;
        sd(isnan(sd)) = 0;
        assert(isscalar(mwPower))
        col = BODYPARTCOLORS{iBodypart};
        h(iBodypart) = plot(ax, t*1e3, mu, Color=col, DisplayName=BODYPARTDISPNAMES(iBodypart), LineWidth=1.5);
        fprintf("%s %s\n", BODYPARTDISPNAMES(iBodypart), num2str(col));
        patch(ax, [t, flip(t)]*1e3, [mu-sd, flip(mu+sd)], col, LineStyle='-', FaceAlpha=0.075, FaceColor=col, EdgeAlpha=0.075, EdgeColor=col)

        ylim(ax, YLIMS{iBodypart})
        yticks(ax, YTICKS{iBodypart})
        if useYYAxis
            ylabel(ax, sprintf("%s", BODYPARTDISPNAMES(iBodypart)))
        end
    end
    
    xline(ax, [0, 20], ':')
    xlim(ax, [windowPreStim(1), windowPostStim(2)] * 1e3)
    % fprintf("%i nm, %g mW, %i trials, %i sessions\n", WAVELENGTHS(iColor), mwPower, nTrials, length(exp));
    title(ax, sprintf("%i nm, %g mW", WAVELENGTHS(iColor), mwPower))
    fontsize(ax, p.fontSize, 'points')
    set(ax.YAxis, TickLength=[0.04, 0.025])
end
xlabel(tl, 'Time from laser on (ms)', FontSize=p.fontSize);
ylabel(tl, 'Speed (a.u.)', FontSize=p.fontSize)
if ~useYYAxis
    lgd = legend(h, Orientation='horizontal');
    lgd.Layout.Tile = 'north';
end

copygraphics(fig, BackgroundColor='none', ContentType='vector')
clear windowPreStim windowPostStim WAVELENGTHS COLORS BODYPARTS BODYPARTDISPNAMES YYAXIS BODYPARTCOLORS YLIMS YTICKS
clear iExp tl h AX iColor color mwPower ax iBodypart bodypart X Y t nTrials velX velY spd mu sd col h lgd useYYAxis
clear fig

%% Fig 8d. Example rasters of optotagging SNr neurons (left: blue, right: red)
% close all
% path = 'E:\Figures\Rasters_SNr_SCRetro';
% fig = figure(Units='inches', Position=[1 1 7 7]);
% ax = axes(fig);
% for iEu = 1:length(SNr_SCRetro.rd.stim)
%     if ~SNr_SCRetro.c.isSNr(iEu)
%         continue
%     end
%     cla(ax);
%     EphysUnit.plotRaster(ax, SNr_SCRetro.rd.stim(iEu), xlim=[-0.1, 0.4]);
%     print(fig, sprintf("%s\\%s.png", path, SNr_SCRetro.rd.stim(iEu).name), '-dpng', '-r0')
% end 
% close all
% path = 'E:\Figures\Rasters_SNr_SCRetro_ReverseInjection';
% fig = figure(Units='inches', Position=[1 1 7 7]);
% ax = axes(fig);
% for iEu = 1:length(SNr_SCRetro_ReverseInjection.rd.stim)
%     if ~SNr_SCRetro_ReverseInjection.c.isSNr(iEu)
%         continue
%     end
%     cla(ax);
%     EphysUnit.plotRaster(ax, SNr_SCRetro_ReverseInjection.rd.stim(iEu), xlim=[-0.1, 0.4]);
%     print(fig, sprintf("%s\\%s.png", path, SNr_SCRetro_ReverseInjection.rd.stim(iEu).name), '-dpng', '-r0')
% end
% clear path fig ax iEu 

egUnitNames = { ...
    % 'daisy26_20250425_Channel83_Unit1', ... Only blue
    % 'daisy26_20250425_Channel143_Unit1', ... Only blue, latency scales
    % 'desmond38_20250403_Channel362_Unit1', ... Only blue,  25uW -> 8mW
    "desmond39_20250423_Channel380_Unit1", ... Only blue (100uW->2mW), not red (100uW->16mW)
    "daisy26_20250425_Channel185_Unit1", ... blue (100uW->2mW), also red (2mW->16mW)
    };

close all
p.sz = 1.5;
fig = figure(Units='inches', Position=[1, 1, 5, 2]);
tl = tiledlayout(fig, length(egUnitNames), 2, TileSpacing='tight');
ax = gobjects(2, 2);
for i = 1:2
    for j = 1:2
        ax(i, j) = nexttile(tl);
    end
end
for iUnit = 1:length(egUnitNames)
    iEu = find(string({SNr_SCRetro.rd.stim.name}) == egUnitNames{iUnit});
    EphysUnit.plotRaster(ax(1, iUnit), SNr_SCRetro.rd.stim(iEu), xlim=[-0.02, 0.05], filterByTwoColorConditions=struct(wavelength=[470, 473], power=["*1e6>25"], duration=[20e-3]), mergeStimTrains=true, sz=p.sz);
    EphysUnit.plotRaster(ax(2, iUnit), SNr_SCRetro.rd.stim(iEu), xlim=[-0.02, 0.05], filterByTwoColorConditions=struct(wavelength=[590, 593, 635], power=["*1e6>25"], duration=[20e-3]), mergeStimTrains=true, sz=p.sz);
    if iUnit == 2
        ax(1, iUnit).Legend.FontSize = 7;
        ax(2, iUnit).Legend.FontSize = 7;
    else
        delete(ax(1, iUnit).Legend);
        delete(ax(2, iUnit).Legend);
    end
    title(ax(1, iUnit), sprintf('SNr unit %i', iUnit))
    title(ax(2, iUnit), '')
    % title(tl, SNr_SCRetro.rd.stim(iEu).name, Interpreter='none')
end
xlim(ax, [-0.02, 0.05])
xticks(ax, [0, 20, 50]*1e-3)
xticklabels(ax, ["0", "20", "50"])
xlabel(ax, '')
ylabel(ax, '')
xticklabels(ax(1, :), ["", "", ""])
fontsize(ax, p.fontSize, 'points')
xlabel(tl, 'Time from opto onset (ms)', FontSize=p.fontSize)
ylabel(tl, 'Trial', FontSize=p.fontSize)

ax(1, 2).Legend.Location = 'eastoutside';
ax(2, 2).Legend.Location = 'eastoutside';
copygraphics(fig, BackgroundColor='none', ContentType='vector')
% clear i j fig tl ax egUnitNames iUnit iEu

%% Fig 8e. Line-point plot Laser pwr vs. response (\deltaSR)
%% Fig 8f. Line-point plot Laser pwr vs. response latency
%% Fig 8g. Heatmap, All optotagged (thus SC-projecting) SNr neurons (lick vs. reach)
%% Fig 8h. Heatmap, SNr->LickSC (lick vs. reach)
%% Fig 8i. Heatmap, SNr->ReachSC (lick vs. reach)
