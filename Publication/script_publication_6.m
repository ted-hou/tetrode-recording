%%% Figure 6. Lick vs. Reach ephys + Osci lick and somatotopy

%% Load all units
load_ephysunits;
boot_response_dir;


%% plot params
p.fontSize = 8;
p.width = 8;
p.height = 10;
p.heightFirstRow = 3;
p.heightSecondRow = 3;
p.heightThirdRow = 1.5;
p.rowHeight = [3, 3, 2.25, 1.75];

p.etaSortWindow = [-3, 0];
p.etaSignWindow = [-0.3, 0];
p.etaLatencyThresholdPos = 0.5;
p.etaLatencyThresholdNeg = 0.25;

%% 6a. Double rasters (reach vs lick) 4 example untis
unitNames = { ...
    'desmond24_20220510_Channel44_Unit1'; ...
    'Daisy2_20180420_Channel14_Unit1'; ...
    'daisy13_20220106_Electrode97_Unit1'; ...
    'daisy8_20210709_Channel7_Unit1'; ...
    };
files = cellfun(@(name) sprintf('C:\\SERVER\\Units\\Lite_NonDuplicate\\%s.mat', name), unitNames, UniformOutput=false);
euEg_lickVsPress = EphysUnit.load(files);

close all
for iEu = 1:length(unitNames)
    thisRdPress = euEg_lickVsPress(iEu).getRasterData('press', window=[0, 2], sort=true);
    thisRdLick = euEg_lickVsPress(iEu).getRasterData('lick', window=[0, 2], sort=true);
    ax = EphysUnit.plotMultiRaster([thisRdPress, thisRdLick], label=["Reach", "Lick"], ...
        xlim=[-3, 2], iti=false, sz=0.75, maxTrials=80, ...
        figUnits='inches', figPosition=[0, 0, p.width/length(unitNames), p.rowHeight(1)], figMargins=1.2*[0.16, 0.09]);
    l = findobj(Type='Legend', Parent=ax(1).Parent);
    if iEu == 1
        delete(l(1))
        l(2).Orientation = 'horizontal';
        l(2).Position = [0.052677548682703,0.947314293123628,0.906249988824129,0.0442708323244];
    else
        delete(l)
        ylabel(ax, '')
    end
    xlabel(ax(1), '')
    xlabel(ax(2), 'Time to contact (s)')
    l = findobj(Type='Legend', Parent=ax(1).Parent);
    fontsize([ax; l], p.fontSize, 'points')
    fontname([ax; l], 'Arial')
end
clear thisRd ax iEu

%% 6b. PETH Reach vs. Lick vs. Osci Lick

% Reach, Lick, Osci Lick, sorted amongst themselves
close all
clear ax
fig = figure(Units='inches', Position=[0, 0, p.width, p.rowHeight(2)], DefaultAxesFontSize=p.fontSize);
ax(1) = axes(fig, Position=[0.135507244803911,0.11,0.207898552297538,0.815], FontSize=p.fontSize);
ax(2) = axes(fig, Position=[0.416304346253187,0.11,0.207898552297538,0.815], FontSize=p.fontSize);
ax(3) = axes(fig, Position=[0.697101447702462,0.11,0.207898552297538,0.815], FontSize=p.fontSize);
[~, ~] = EphysUnit.plotETA(ax(1), eta.press, c.hasPress & c.hasLick, ...
    clim=[-2, 2], xlim=[-4, 0], sortWindow=p.etaSortWindow, signWindow=p.etaSignWindow, ...
    sortThreshold=p.etaLatencyThresholdPos, negativeSortThreshold=p.etaLatencyThresholdNeg, hidecolorbar=true);
[~, ~] = EphysUnit.plotETA(ax(2), eta.lick, c.hasPress & c.hasLick, ...
    clim=[-2, 2], xlim=[-4, 0], sortWindow=p.etaSortWindow, signWindow=p.etaSignWindow, ...
    sortThreshold=p.etaLatencyThresholdPos, negativeSortThreshold=p.etaLatencyThresholdNeg, hidecolorbar=true);
[~, ~] = EphysUnit.plotETA(ax(3), eta.anyLickNorm, c.isLick, ...
    clim=[-2, 2], xlim=[-0.25, 0.25], sortWindow=[-0.13, -0.01], signWindow=[-0.13, -0.01], hidecolorbar=true);
title(ax(1), 'Pre-reach')
title(ax(2), 'Pre-lick')
title(ax(3), 'Osci-lick')
ylabel(ax(2:3), '')
xlabel(ax(1:3), 'Time to spout-contact (s)')
h = colorbar(ax(3)); 
h.Position = [0.913242151752656,0.109479305740988,0.013611111111111,0.815754339118825];
h.Label.String = 'z-scored spike rate (a.u.)';
fontsize(ax, p.fontSize, 'points')
fontname(ax, 'Arial')

%%% Compare lick vs osci lick
% close all
% clear ax
% fig = figure(Units='inches', Position=[0, 0, p.width, p.heightSecondRow], DefaultAxesFontSize=p.fontSize);
% ax(1) = axes(fig, Position=[0.135507244803911,0.11,0.207898552297538,0.815], FontSize=p.fontSize);
% ax(2) = axes(fig, Position=[0.416304346253187,0.11,0.207898552297538,0.815], FontSize=p.fontSize);
% ax(3) = axes(fig, Position=[0.697101447702462,0.11,0.207898552297538,0.815], FontSize=p.fontSize);
% [~, order] = EphysUnit.plotETA(ax(1), eta.lick, c.hasLick & c.isLick, ...
%     clim=[-2, 2], xlim=[-4, 0], sortWindow=p.etaSortWindow, signWindow=p.etaSignWindow, ...
%     sortThreshold=p.etaLatencyThresholdPos, negativeSortThreshold=p.etaLatencyThresholdNeg, hidecolorbar=true);
% [~, ~] = EphysUnit.plotETA(ax(2), eta.anyLickNorm, c.hasLick & c.isLick, order=order, ...
%     clim=[-2, 2], xlim=[-0.25, 0.25], hidecolorbar=true);
% [~, ~] = EphysUnit.plotETA(ax(3), eta.anyLickNorm, c.hasLick & c.isLick, ...
%     clim=[-2, 2], xlim=[-0.25, 0.25], sortWindow=[-0.13, -0.01], signWindow=[-0.13, -0.01], ...
%     sortThreshold=p.etaLatencyThresholdPos, negativeSortThreshold=p.etaLatencyThresholdNeg, hidecolorbar=true);
% title(ax(1), 'Pre-lick')
% title(ax(2), 'Osci-lick')
% title(ax(3), 'Osci-lick (phase-sorted)')
% ylabel(ax(2:3), '')
% xlabel(ax(1:3), 'Time to spout-contact (s)')
% h = colorbar(ax(3)); 
% h.Position = [0.913242151752656,0.109479305740988,0.013611111111111,0.815754339118825];
% h.Label.String = 'z-scored spike rate (a.u.)';
% fontsize(ax, p.fontSize, 'points')
% fontname(ax, 'Arial')

%% 6c. Reach vs. Lick (scatter), Reach vs. Lick (osci subset, scatter), Pie-chart
close all
sz = 2;
SEL = { ...
    c.hasLick & c.hasPress, ...
    c.hasLick & c.hasPress & c.isLick, ...
    c.hasLick & c.hasPress, ...
    };

XDATA = { ...
    meta.lickRaw*10 - msr, ...
    meta.lickRaw*10 - msr, ...
    };

YDATA = { ...
    meta.pressRaw*10 - msr, ...
    meta.pressRaw*10 - msr, ...
    };

XDATA = { ...
    meta.lick, ...
    meta.lick, ...
    };

YDATA = { ...
    meta.press, ...
    meta.press, ...
    };

XNAME = { ...
    'Pre-lick (\Deltasp/s)', ...
    'Pre-lick (\Deltasp/s)'; ...
    };
YNAME = { ...
    'Pre-reach (\Deltasp/s)', ...
    'Pre-reach (\Deltasp/s)'; ...
    };
TITLE = { ...
    'All units', ...
    'Lick-osci units', ...
    sprintf('N = %g', nnz(SEL{3}))
    };

fig = figure(Units='inches', Position=[0, 0, p.width, p.rowHeight(3)], DefaultAxesFontSize=p.fontSize);
ax = arrayfun(@(i) subplot(1, 3, i), 1:3);
ax(1).Position = [0.075,0.151851851851852,0.213405797101449,0.769753086419753];
ax(2).Position = [0.345,0.151851851851852,0.213405797101449,0.769753086419753];
ax(3).Position = [0.58, 0.151851851851852,0.213405797101449,0.769753086419753];
% ax(1).Position(1) = 0.075;
% ax(2).Position(1) = 0.345;
% ax(3).Position(1) = 0.58;
h = gobjects(2, 2);
hold(ax, 'on')

% 1. Reach vs. Lick scatter
for i = 1:2
    sel = SEL{i};
    x = XDATA{i}(sel);
    y = YDATA{i}(sel);
    h(i, 1) = scatter(ax(i), x, y, sz, 'black', 'filled', DisplayName=sprintf('%g units', nnz(sel)));
    
    if i == 1
        xl = ax(i).XLim; yl = ax(i).YLim;
        set(ax(i), XLimMode='manual', YLimMode='manual');
    else
        set(ax(i), XLimMode='manual', YLimMode='manual');
        xlim(ax(i), xl);
        ylim(ax(i), yl);
    end
    plot(ax(i), xl, [0, 0], 'k:');
    plot(ax(i), [0, 0], yl, 'k:')
    
    mdl = fitlm(x, y);
    h(i, 2) = plot(ax(i), xl', mdl.predict(xl'), 'k--', LineWidth=1, DisplayName=sprintf('R^2 = %.2f', mdl.Rsquared.Ordinary));
    xlabel(ax(i), XNAME{i})
    ylabel(ax(i), YNAME{i})
    legend(ax(i), h(i, :), Location='northwest')
    title(ax(i), TITLE{i})
end

% Pie chart to count units
i = 3;
selCommon = SEL{i};
selCats = { ...
        ~c.isPressResponsive & c.isLickResponsive, ...         % No response
        c.isPressDown & c.isLickDown, ...         % Same direction
        c.isPressUp & c.isLickUp, ...
        c.isPressDown & ~c.isLickResponsive, ...
        c.isPressUp & ~c.isLickResponsive, ...        % Simple selective
        c.isLickDown & ~c.isPressResponsive, ...
        c.isLickUp & ~c.isPressResponsive, ...
        c.isPressDown & c.isLickUp, ...         % Super selective
        c.isPressUp & c.isLickDown, ...
    };

catNames = { ...
        'Reach~ Lick~', ...
        'Reach\downarrow Lick\downarrow', ...
        'Reach\uparrow Lick\uparrow', ...
        'Reach\downarrow Lick~', ...
        'Reach\uparrow Lick~', ...
        'Reach~ Lick\downarrow', ...
        'Reach~ Lick\uparrow', ...
        'Reach\downarrow Lick\uparrow', ...
        'Reach\uparrow Lick\downarrow', ...
    };

catHue = [linspace(270, 215, 3), linspace(140, 90, 2), linspace(55, 35, 2), linspace(15, 0, 2)] ./ 360;
catLightness =  [0.4, 0.4, 0.4, 0.3, 0.4, 0.4, 0.4, 0.3, 0.4];

% explode = [0 0 0 1 1 1 1 1 1];
explode = [0 0 0 0 0 0 0 0 0];

pie(ax(i), cellfun(@(sel) nnz(selCommon & sel), selCats), explode, '%.0f%%')
axis(ax(i), 'equal', 'off')
legend(ax(i), catNames, Position=[0.825883408717677,0.091446975129934,0.167100694444444,0.847222222222221], ...
    FontSize=p.fontSize, FontName='Courier');
hLabels = findobj(Parent=ax(i), Type='Text');
set(hLabels, FontSize=p.fontSize, FontName='Arial');
% hLabels(1).Position = [-0.033332041944592,1.186850167563696,0];
% hLabels(2).Position = [0.415221112669196,1.171033576756159,0];
hLabels(1).Position = [-0.028964898552372,1.067127375229603,0];
hLabels(2).Position = [0.321524523997384,1.040929743283167,0];

text(ax(i), -0.196325626228295, -1.323064002842862, TITLE{i}, FontSize=p.fontSize*1.1, FontName='Arial')

% Recolor
hPatches = flip(findobj(Parent=ax(i), Type='Patch'));
for iPatch = 1:length(hPatches)
    hPatches(iPatch).FaceColor = hsl2rgb([catHue(iPatch), 1, catLightness(iPatch)]);
    % text(ax(i), mean(hPatches(iPatch).XData), mean(hPatches(iPatch).YData), num2str(iPatch))
end

set(ax, FontSize=p.fontSize, FontName='Arial');
hold(ax, 'off')

% clear ax sel x y i mdl h sz SEL XDATA YDATA XNAME XNAME TITLE

%% 6d. Salt and pepper map (lick vs. reach vs. osci lick)
close all
MOVETYPE = {'press', 'lick', 'lickosci'};
SEL = {c.isPressResponsive; c.isLickResponsive; c.isLick};
STATS = {meta.press, meta.lick, meta.anyLickNorm};
TITLE = {'Pre-reach', 'Pre-lick', 'Lick-osci'};
COLOR = {[], [], hsl2rgb([50/360, 1, 0.4])};
ALPHA = {0.125, 0.125, 0.5};
POS = { ...
        [0.09,0.187301712531548,0.213405797101449,0.721428477506101], ...
        [0.324,0.187301712531548,0.213405797101449,0.721428477506101], ...
        [0.56,0.187301712531548,0.213405797101449,0.721428477506101], ...
    };

fig = figure(Units='inches', Position=[0, 0, p.width*3/4, p.rowHeight(4)], DefaultAxesFontSize=p.fontSize);
for iMove = 1:length(MOVETYPE)
    % ax = subplot(1, length(MOVETYPE), iMove);
    ax = axes(fig, Position=POS{iMove});
    sel = SEL{iMove};
    coords = euPos(sel, :);
    stats = STATS{iMove}(sel);
    AcuteRecording.plotMap(ax, coords, stats, [0 5], 0, UseSignedML=false, BubbleSize=[1, 5], MarkerAlpha=ALPHA{iMove}, Color=COLOR{iMove});
    title(ax, TITLE{iMove})
    axis(ax, 'image')
    xlim(ax, [0.9, 1.7])
    ylim(ax, [-4.8, -3.7])
    xticks(ax, [1, 1.6])
    yticks(ax, [-4.7, -3.8])
    if iMove > 1
        ylabel(ax, "");
    end
    xlabel(ax, 'ML');
    fontsize(ax, p.fontSize, 'points');
    fontname(ax, 'Arial')
    fontsize(fig, p.fontSize, 'points');
    fontname(fig, 'Arial')
end

% Make a dummy axes for the legends (red, blue, yellow circles)
ax = axes(fig, Position=[0.5, 0.8, 0, 0]); hold(ax, 'on')
ax.Visible = 'off';
scatter(ax, 0, 0, 1, [1, 0, 0], 'filled', DisplayName='excited')
scatter(ax, 0, 0, 1, [0, 0, 1], 'filled', DisplayName='suppressed')
scatter(ax, 0, 0, 1, COLOR{3}, 'filled', DisplayName='oscillatory')
h = legend(ax, Position=[0.783564814814813,0.474206349206349,0.169560185185185,0.232142857142857], Orientation='vertical');
fontsize(h, p.fontSize, 'points');
fontname(h, 'Arial')