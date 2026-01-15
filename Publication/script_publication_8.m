read_SC_opto_trajectories_DLC;

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

[~, SNr_SCRetro.coords] = getNeuroPixelChannelMap(SNr_SCRetro.eu, ml=1300, ap=-3280, dv=-4700);
SNr_SCRetro.c.isSNr = SNr_SCRetro.coords(:, 2) <= -3.6*1e3;

[~, SNr_SCRetro_ReverseInjection.coords] = getNeuroPixelChannelMap(SNr_SCRetro_ReverseInjection.eu, ml=1300, ap=-3280, dv=-4700);
SNr_SCRetro_ReverseInjection.c.isSNr = SNr_SCRetro_ReverseInjection.coords(:, 2) <= -3.6*1e3;


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
% BODYPARTCOLORS = arrayfun(@(i) getColor(i, 7, 0.7), [1, 3], UniformOutput=false);
BODYPARTCOLORS = {hsl2rgb([170/360, 0.5, 0.35]); hsl2rgb([308/360, 0.5, 0.5])};
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
xlabel(tl, 'Time from opto onset (ms)', FontSize=p.fontSize);
ylabel(tl, 'Speed (a.u.)', FontSize=p.fontSize)
if ~useYYAxis
    lgd = legend(h, Orientation='horizontal');
    lgd.Layout.Tile = 'north';
end
ax = AX(1);
hLetter = text(ax, 0, 0, 'b', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.25, ax.Position(4) + 0.5, 0];

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
fig = figure(Units='inches', Position=[1, 1, 4.7, 2]);
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
ax(1, 2).Legend.Position = [0.704491725768321,0.604253468207187,0.234027777777778,0.223958333333333];

hLetter = text(ax(1, 1), 0, 0, 'd', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax(1).Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.4, ax(1).Position(4) + 0.25, 0];

copygraphics(fig, BackgroundColor='none', ContentType='vector')
clear i j fig tl ax egUnitNames iUnit iEu hLetter

%% Fig 8e/f. Bar chart, modulation strength/latency scales with light power
p.stimThreshold = 2;

% Calculate onset timing
etaRed = SNr_SCRetro.eta.stimRed(["100", "500", "2000+"]);
etaBlue = SNr_SCRetro.eta.stimBlue(["100", "500", "2000"]);
clear tOnset
tOnset.red = NaN(length(eu), 3);
tOnset.blue = NaN(length(eu), 3);
for iPower = 1:length(etaRed)
    t = etaRed(iPower).t;
    X = etaRed(iPower).X(:, t>=0);
    t = t(t>=0);
    for iEu = 1:size(X, 1)
        iStart = strfind([false, X(iEu, :)>2], [0, 1, 1]);
        if isempty(iStart)
            tOnset.red(iEu, iPower) = NaN;
        else
            iStart = iStart(1) - 1;
            if iStart > 0
                tOnset.red(iEu, iPower) = t(iStart);
            else
                tOnset.red(iEu, iPower) = 0;
            end
        end
    end
end

for iPower = 1:length(etaBlue)
    t = etaBlue(iPower).t;
    X = etaBlue(iPower).X(:, t>=0);
    t = t(t>=0);
    for iEu = 1:size(X, 1)
        iStart = strfind([false, X(iEu, :)>p.stimThreshold], [0, 1, 1]);
        if isempty(iStart)
            tOnset.blue(iEu, iPower) = NaN;
        else
            iStart = iStart(1) - 1;
            if iStart > 0
                tOnset.blue(iEu, iPower) = t(iStart);
            else
                tOnset.blue(iEu, iPower) = 0;
            end
        end
    end
end

% Calculate stim response
metaRed = cat(2, SNr_SCRetro.meta.stimRed.values{:});
metaBlue = cat(2, SNr_SCRetro.meta.stimBlue.values{:});
metaRed = array2table(metaRed, VariableNames=SNr_SCRetro.meta.stimRed.keys);
metaBlue = array2table(metaBlue, VariableNames=SNr_SCRetro.meta.stimBlue.keys);

cStim.isRed = metaRed{:, "2000+"}>p.stimThreshold & SNr_SCRetro.c.isSNr;
cStim.isBlue100 = metaBlue{:, "100"}>p.stimThreshold & SNr_SCRetro.c.isSNr;
cStim.isBlue = metaBlue{:, "500"}>p.stimThreshold & SNr_SCRetro.c.isSNr;
cStim.isBlue500 = metaBlue{:, "500"}>p.stimThreshold & SNr_SCRetro.c.isSNr;
cStim.isBlue2000 = metaBlue{:, "2000"}>p.stimThreshold & SNr_SCRetro.c.isSNr;
cStim.isBlueNotRed = cStim.isBlue500 & ~cStim.isRed & SNr_SCRetro.c.isSNr;
cStim.isRedNotBlue = cStim.isRed & ~cStim.isBlue2000 & SNr_SCRetro.c.isSNr;
fprintf("%i (of %i) neurons with DV < -3.6\n", nnz(SNr_SCRetro.c.isSNr), length(SNr_SCRetro.c.isSNr));
fprintf("isRed = %i, isBlue = %i, isBlueNotRed = %i, isRedNotBlue = %i\n", nnz(cStim.isRed), nnz(cStim.isBlue500), nnz(cStim.isBlueNotRed), nnz(cStim.isRedNotBlue));


% 8e. Plot stim response
close all
fig = figure(Units='inches', Position=[1, 1, 6.5, 1.5]);
tlp = tiledlayout(fig, 1, 2, TileSpacing='loose', Padding='compact');
tl = gobjects(1, 2);
tl(1) = tiledlayout(tlp, 1, 2, TileSpacing='compact', Padding='compact'); tl(1).Layout.Tile = 1;
tl(2) = tiledlayout(tlp, 1, 2, TileSpacing='compact', Padding='compact'); tl(2).Layout.Tile = 2;
ax = gobjects(1, 2);
for i = 1:2
    ax(i) = nexttile(tl(1));
end
hold(ax, 'on')

h = gobjects(2, 2);
h(1, :) = bar(ax(1), 1:3, [mean(metaBlue{cStim.isBlueNotRed, 1:3}, 1, 'omitnan'); mean(metaRed{cStim.isBlueNotRed, 1:3}, 1, 'omitnan')], 1, FaceAlpha=0.33, Clipping='off');
errorbar(ax(1), (1:3)-0.15, mean(metaBlue{cStim.isBlueNotRed, 1:3}, 1, 'omitnan'), std(metaBlue{cStim.isBlueNotRed, 1:3}, 0, 1, 'omitnan')*0.1, LineStyle='none', Color='black', CapSize=3, Clipping='off');
errorbar(ax(1), (1:3)+0.15, mean(metaRed{cStim.isBlueNotRed, 1:3}, 1, 'omitnan'), std(metaRed{cStim.isBlueNotRed, 1:3}, 0, 1, 'omitnan')*0.1, LineStyle='none', Color='black', CapSize=3, Clipping='off');
title(ax(1), sprintf('CoChR^+\n(n=%i)', nnz(cStim.isBlueNotRed)))

h(2, :) = bar(ax(2), 1:3, [mean(metaBlue{cStim.isRed, 1:3}, 1, 'omitnan'); mean(metaRed{cStim.isRed, 1:3}, 1, 'omitnan')], 1, FaceAlpha=0.33, Clipping='off');
errorbar(ax(2), (1:3)-0.15, mean(metaBlue{cStim.isRed, 1:3}, 1, 'omitnan'), std(metaBlue{cStim.isRed, 1:3}, 0, 1, 'omitnan')*0.1, LineStyle='none', Color='black', CapSize=3, Clipping='off')
errorbar(ax(2), (1:3)+0.15, mean(metaRed{cStim.isRed, 1:3}, 1, 'omitnan'), std(metaRed{cStim.isRed, 1:3}, 0, 1, 'omitnan')*0.1, LineStyle='none', Color='black', CapSize=3, Clipping='off')
title(ax(2), sprintf('ChrimsonR^+\n(n=%i)', nnz(cStim.isRed)))

set(h(:, 1), FaceColor='blue', DisplayName='470nm')
set(h(:, 2), FaceColor='red', DisplayName='635nm')

xticks(ax, 1:3)
xticklabels(ax(1), ["0.1", "0.5", "2+"])
xticklabels(ax(2), ["0.1", "0.5", "2+"])
xtickangle(ax, 0)
xlim(ax, [0.5, 3.5])
yl = vertcat(ax.YLim);
% ylim(ax, [0, max(yl(:, 2))])
ylim(ax, [0, 25])
yticks(ax(2), [])

xlabel(tl(1), 'Light power (mW)')
ylabel(tl(1), ["Response", "(a.u.)"])
fontsize(tl(1), p.fontSize, 'points')

hLetter = text(ax(1), 0, 0, 'e', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax(1).Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.45, ax(1).Position(4) + 0.45, 0];

% 8f. Plot stim latency
ax = gobjects(1, 2);
for i = 1:2
    ax(i) = nexttile(tl(2));
end
hold(ax, 'on')

h = gobjects(2, 2);
h(1, :) = bar(ax(1), 1:3, 1e3*[mean(tOnset.blue(cStim.isBlueNotRed, 1:3), 1, 'omitnan'); mean(tOnset.red(cStim.isBlueNotRed, 1:3), 1, 'omitnan')], 1, FaceAlpha=0.33);
errorbar(ax(1), (1:3)-0.15, 1e3*mean(tOnset.blue(cStim.isBlueNotRed, 1:3), 1, 'omitnan'), std(1e3*tOnset.blue(cStim.isBlueNotRed, 1:3), 0, 1, 'omitnan')*0.1, LineStyle='none', Color='black', CapSize=3, Clipping='off')
errorbar(ax(1), (1:3)+0.15, 1e3*mean(tOnset.red(cStim.isBlueNotRed, 1:3), 1, 'omitnan'), std(1e3*tOnset.red(cStim.isBlueNotRed, 1:3), 0, 1, 'omitnan')*0.1, LineStyle='none', Color='black', CapSize=3, Clipping='on')
title(ax(1), sprintf('CoChR^+\n(n=%i)', nnz(cStim.isBlueNotRed)))

h(2, :) = bar(ax(2), 1:3, 1e3*[mean(tOnset.blue(cStim.isRed, 1:3), 1, 'omitnan'); mean(tOnset.red(cStim.isRed, 1:3), 1, 'omitnan')], 1, FaceAlpha=0.33);
errorbar(ax(2), (1:3)-0.15, 1e3*mean(tOnset.blue(cStim.isRed, 1:3), 1, 'omitnan'), std(1e3*tOnset.blue(cStim.isRed, 1:3), 0, 1, 'omitnan')*0.1, LineStyle='none', Color='black', CapSize=3, Clipping='off')
errorbar(ax(2), (1:3)+0.15, 1e3*mean(tOnset.red(cStim.isRed, 1:3), 1, 'omitnan'), std(1e3*tOnset.red(cStim.isRed, 1:3), 0, 1, 'omitnan')*0.1, LineStyle='none', Color='black', CapSize=3, Clipping='on')
title(ax(2), sprintf('ChrimsonR^+\n(n=%i)', nnz(cStim.isRed)))

set(h(:, 1), FaceColor='blue', DisplayName='470nm')
set(h(:, 2), FaceColor='red', DisplayName='635nm')

xticks(ax, 1:3)
xticklabels(ax(1), ["0.1", "0.5", "2+"])
xticklabels(ax(2), ["0.1", "0.5", "2+"])
xtickangle(ax, 0)
xlim(ax, [0.5, 3.5])
yl = vertcat(ax.YLim);
% ylim(ax, [0, max(yl(:, 2))])
ylim(ax, [0, 100])
yticks(ax(2), [])

xlabel(tl(2), 'Light power (mW)')
ylabel(tl(2), ["Latency", "(ms)"])
fontsize(tl(2), p.fontSize, 'points')


hLetter = text(ax(1), 0, 0, 'f', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax(1).Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.6, ax(1).Position(4) + 0.45, 0];

l = legend(h(2, :), Orientation='vertical', FontSize=7);
l.Layout.Tile = 'east';
copygraphics(fig, BackgroundColor='none', ContentType='vector')

clear metaRed metaBlue fig tl tlp ax i yl h l

%% Fig 8e/f/g. Heatmap, {SNr->SC, SNr->LickSC, SNr->ReachSC} SNr neurons (lick vs. reach)
cStim.isSCProjecting = (cStim.isRed | cStim.isBlue);
cStim.isReachSCProjecting = cStim.isBlueNotRed;
cStim.isLickSCProjecting = cStim.isRed;

xPress = SNr_SCRetro.meta.press;
xLick = SNr_SCRetro.meta.lick;

groupVar = NaN(length(eu), 1);
groupVar(xPress<0 & xLick>0) = 0;
groupVar(xPress>0 & xLick<0) = 1;
groupVar(xPress<0 & xLick<0) = 2;
groupVar(xPress>0 & xLick>0) = 3;

groupVarPress = NaN(length(eu), 1);
groupVarPress(xPress<0 & xLick>0) = 0;
groupVarPress(xPress<0 & xLick<0) = 1;
groupVarPress(xPress>0 & xLick>0) = 2;
groupVarPress(xPress>0 & xLick<0) = 3;

groupVarLick = NaN(length(eu), 1);
groupVarLick(xLick<0 & xPress>0) = 0;
groupVarLick(xLick<0 & xPress<0) = 1;
groupVarLick(xLick>0 & xPress>0) = 2;
groupVarLick(xLick>0 & xPress<0) = 3;

SEL = {cStim.isSCProjecting, cStim.isReachSCProjecting, cStim.isLickSCProjecting};
TITLE = {"SC-projecting", "SC^{arm}-projecting", "SC^{orofacial}-projecting"};
GROUPVAR = {groupVar, groupVarPress, groupVarLick};
LETTER = 'ghi';

close all
fig = figure(Units='inches', Position=[1, 1, 6.5, 3.25]);
tlp = tiledlayout(fig, 1, 3, TileSpacing='loose', Padding='compact');
tl = gobjects(1, 3);
tl(1) = tiledlayout(tlp, 1, 2, TileSpacing='tight', Padding='tight'); tl(1).Layout.Tile = 1;
tl(2) = tiledlayout(tlp, 1, 2, TileSpacing='tight', Padding='tight'); tl(2).Layout.Tile = 2;
tl(3) = tiledlayout(tlp, 1, 2, TileSpacing='tight', Padding='tight'); tl(3).Layout.Tile = 3;

for iTl = 1:3
    sel = SEL{iTl};
    ax = gobjects(1, 2);
    for i = 1:2
        ax(i) = nexttile(tl(iTl));
    end
    [~, order] = EphysUnit.plotETA(ax(1), SNr_SCRetro.eta.press, sel, hideColorbar=true, xlim=[-2, 0], ...
        sortGroup=GROUPVAR{iTl}(sel), sortWindow=[-2.5, 0], signWindow=[-0.3, 0], sortThreshold=0.25, negativeSortThreshold=0.25);
    EphysUnit.plotETA(ax(2), SNr_SCRetro.eta.lick, sel, hideColorbar=true, xlim=[-2, 0], ...
        sortGroup=GROUPVAR{iTl}(sel), order=order);
    applyCustomColormap(ax(1), [-1.5, 1.5], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);
    applyCustomColormap(ax(2), [-1.5, 1.5], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);

    N = histcounts(GROUPVAR{iTl}(sel), [-0.5, 0.5, 1.5, 2.5, 3.5]);
    yline(ax(1), cumsum(N(1:end-1)) + 0.5, 'k--');
    yline(ax(2), cumsum(N(1:end-1)) + 0.5, 'k--');
    xline(ax(1), 0, 'k--')
    xline(ax(2), 0, 'k--')

    yt = unique([cumsum(N(1:end-1)) + 0.5, nnz(sel)]);
    ytl = string(unique([cumsum(N(1:end-1)), nnz(sel)]));
    yticks(ax(1), yt)
    yticklabels(ax(1), ytl)
    yticks(ax(2), [])
    ylim(ax, [0.5, nnz(sel)+0.5])
    set(ax(1).YAxis, TickLength=[0, 0])

    % Manual yticks for the last panel
    if iTl == 3
        yt = unique([cumsum(N(1:end-1)) + 0.5, nnz(sel)]);
        ytl = string(unique([cumsum(N(1:end-1)), nnz(sel)]));
        yt(end-1) = yt(end-1) - 1.5;
        yt(end) = yt(end) + 0.5;
        yticks(ax(1), yt)
        yticklabels(ax(1), ytl)
    end

    title(ax(1), 'reach')
    title(ax(2), 'lick')
    xlabel(ax, '')
    ylabel(ax, '')
    title(tl(iTl), sprintf("%s (n=%i)", TITLE{iTl}, nnz(sel)), FontWeight='bold', FontSize=p.fontSize)

    fontsize(ax, p.fontSize, 'points')

    if iTl == 1
        hLetter = text(ax(1), 0, 0, LETTER(iTl), FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
        ax(1).Units = 'inches';
        hLetter.HorizontalAlignment = 'right';
        hLetter.VerticalAlignment = 'top';
        hLetter.Position = [-0.35, ax(1).Position(4) + 0.25, 0];
    end
end
xlabel(tl, 'Time to bar/spout contact (s)', FontSize=p.fontSize)
ylabel(tlp, 'Unit', FontSize=p.fontSize)

h = colorbar(ax(2)); 
h.Layout.Tile = 'east';
h.Label.String = 'Normalized spike rate (a.u.)';

copygraphics(fig, BackgroundColor='none', ContentType='vector')