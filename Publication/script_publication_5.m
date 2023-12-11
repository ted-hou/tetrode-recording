%%% Figure 6. Lick vs. Reach ephys + Osci lick and somatotopy

%% Load all units
load_ephysunits;
boot_response_dir;
find_osci_lick_circ;

%% plot params
p.fontSize = 9;
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
% close all
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
% [~, I] = sort(phase(freq==8, c.isLick));
% [~, ~] = EphysUnit.plotETA(ax(3), eta.firstLickNorm, c.isLick, order=I, ...
%     clim=[-2, 2], xlim=[0.01, 0.5], hidecolorbar=true);
[~, I] = sort(angle(meanZ(c.hasPress & c.hasLick & magH')));
[~, ~] = EphysUnit.plotETA(ax(3), eta.lickBoutNorm, c.hasPress & c.hasLick & magH', order=I, ...
    clim=[-2, 2], xlim=[0, 2*pi*maxBoutCycles], hidecolorbar=true);
xticks (ax(3), (0:2:8).*pi);
xticklabels(ax(3), arrayfun(@(x) sprintf('%i\\pi', x), 0:2:8, UniformOutput=false));
title(ax(1), 'Pre-reach')
title(ax(2), 'Pre-lick')
title(ax(3), 'Osci-lick')
ylabel(ax(2:3), '')
xlabel(ax(1), 'Time to bar contact (s)')
xlabel(ax(2), 'Time to spout contact (s)')
xlabel(ax(3), 'Lick phase')
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
sz = 5;
SEL = { ...
    c.hasLick & c.hasPress, ...
    c.hasLick & c.hasPress & magH', ...
    };

% XDATA = { ...
%     meta.lickRaw*10 - msr, ...
%     meta.lickRaw*10 - msr, ...
%     };
% 
% YDATA = { ...
%     meta.pressRaw*10 - msr, ...
%     meta.pressRaw*10 - msr, ...
%     };

XDATA = { ...
    meta.lick, ...
    meta.lick, ...
    };

YDATA = { ...
    meta.press, ...
    meta.press, ...
    };

XNAME = { ...
    'Pre-lick response (a.u.)', ...
    'Pre-lick response (a.u.)'; ...
    };
YNAME = { ...
    'Pre-reach response (a.u.)', ...
    'Pre-reach response (a.u.)'; ...
    };
TITLE = { ...
    'All units', ...
    'Lick-osci units'
    };
LEGENDPOS = {...
    [0.244956331030447,0.173453205894066,0.117187498486601,0.122685182149764], ...
    [0.513185497503088,0.182712465153325,0.11067708201396,0.122685182149764], ...
    };

% TITLE = { ...
%     '', ...
%     ''
%     };

fig = figure(Units='inches', Position=[0, 0, p.width, p.rowHeight(3)], DefaultAxesFontSize=p.fontSize);
ax = arrayfun(@(i) subplot(1, 3, i), 1:3);
% ax(1).Position = [0.075,0.151851851851852,0.213405797101449,0.769753086419753];
% ax(2).Position = [0.345,0.151851851851852,0.213405797101449,0.769753086419753];
% ax(3).Position = [0.58, 0.151851851851852,0.213405797101449,0.769753086419753];
% ax(1).Position(1) = 0.075;
% ax(2).Position(1) = 0.345;
% ax(3).Position(1) = 0.58;
h = gobjects(2, 3);
hold(ax, 'on')

% 1. Reach vs. Lick scatter
for i = 1:2
    sel = SEL{i};
    x = XDATA{i}(sel);
    y = YDATA{i}(sel);
    nnz(sel)
    subselResp = sel & c.isPressResponsive & c.isLickResponsive;
    subselNone = sel & (~c.isPressResponsive | ~c.isLickResponsive);
    h(i, 1) = scatter(ax(i), XDATA{i}(subselNone), YDATA{i}(subselNone), sz, 'black', 'filled', DisplayName=sprintf('%g units', nnz(subselNone)));   
    h(i, 2) = scatter(ax(i), XDATA{i}(subselResp), YDATA{i}(subselResp), sz, 'red', 'filled', DisplayName=sprintf('%g units', nnz(subselResp)));

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
%     h(i, 3) = plot(ax(i), xl', mdl.predict(xl'), 'k--', LineWidth=1, DisplayName=sprintf('R^2 = %.2f', mdl.Rsquared.Ordinary));
    if i == 1
        xlabel(ax(i), XNAME{i}, Units='normalized', Position=[1.18902486151029,-0.115261042292962,0])
        ylabel(ax(i), YNAME{i})
    else
%         yticks(ax(i), [])
    end
    legend(ax(i), h(i, 1:2), Units='normalized', Position=LEGENDPOS{i});

    title(ax(i), TITLE{i})
end

% 2. Peri-lick lick prob histogram
%
histogram(ax(3), BinEdges=lickHistEdges(2:end), BinCounts=lickHistCounts(2:end), Normalization='probability')
yl = ax(3).YLim;
histogram(ax(3), BinEdges=lickHistEdges(1:end), BinCounts=lickHistCounts(1:end), Normalization='probability', EdgeColor='black', FaceColor='black', FaceAlpha=1)
xlim(ax(3), [0, 0.5]);
xticks(ax(3), 0:1/4:0.5)
yticks(ax(3), [])
xlabel(ax(3), 'Time from first lick (s)')
ylabel(ax(3), 'Lick probability')
ylim(ax(3), yl);

% clear I expIndices osciEuIndices expEuIndices iExp iEu osciLickTimes iLick

set(ax, FontSize=p.fontSize, FontName='Arial');
hold(ax, 'off')


%% Pie chart to count units
% i = 3;
% selCommon = SEL{i};
% selCats = { ...
%         ~c.isPressResponsive & c.isLickResponsive, ...         % No response
%         c.isPressDown & c.isLickDown, ...         % Same direction
%         c.isPressUp & c.isLickUp, ...
%         c.isPressDown & ~c.isLickResponsive, ...
%         c.isPressUp & ~c.isLickResponsive, ...        % Simple selective
%         c.isLickDown & ~c.isPressResponsive, ...
%         c.isLickUp & ~c.isPressResponsive, ...
%         c.isPressDown & c.isLickUp, ...         % Super selective
%         c.isPressUp & c.isLickDown, ...
%     };
% 
% catNames = { ...
%         'Reach~ Lick~', ...
%         'Reach\downarrow Lick\downarrow', ...
%         'Reach\uparrow Lick\uparrow', ...
%         'Reach\downarrow Lick~', ...
%         'Reach\uparrow Lick~', ...
%         'Reach~ Lick\downarrow', ...
%         'Reach~ Lick\uparrow', ...
%         'Reach\downarrow Lick\uparrow', ...
%         'Reach\uparrow Lick\downarrow', ...
%     };
% 
% catHue = [linspace(270, 215, 3), linspace(140, 90, 2), linspace(55, 35, 2), linspace(15, 0, 2)] ./ 360;
% catLightness =  [0.4, 0.4, 0.4, 0.3, 0.4, 0.4, 0.4, 0.3, 0.4];
% 
% % explode = [0 0 0 1 1 1 1 1 1];
% explode = [0 0 0 0 0 0 0 0 0];
% 
% pie(ax(i), cellfun(@(sel) nnz(selCommon & sel), selCats), explode, '%.0f%%')
% axis(ax(i), 'equal', 'off')
% legend(ax(i), catNames, Position=[0.825883408717677,0.091446975129934,0.167100694444444,0.847222222222221], ...
%     FontSize=p.fontSize, FontName='Courier');
% hLabels = findobj(Parent=ax(i), Type='Text');
% set(hLabels, FontSize=p.fontSize, FontName='Arial');
% % hLabels(1).Position = [-0.033332041944592,1.186850167563696,0];
% % hLabels(2).Position = [0.415221112669196,1.171033576756159,0];
% hLabels(1).Position = [-0.028964898552372,1.067127375229603,0];
% hLabels(2).Position = [0.321524523997384,1.040929743283167,0];
% 
% text(ax(i), -0.196325626228295, -1.323064002842862, TITLE{i}, FontSize=p.fontSize*1.1, FontName='Arial')
% 
% % Recolor
% hPatches = flip(findobj(Parent=ax(i), Type='Patch'));
% for iPatch = 1:length(hPatches)
%     hPatches(iPatch).FaceColor = hsl2rgb([catHue(iPatch), 1, catLightness(iPatch)]);
%     % text(ax(i), mean(hPatches(iPatch).XData), mean(hPatches(iPatch).YData), num2str(iPatch))
% end

% clear ax sel x y i mdl h sz SEL XDATA YDATA XNAME XNAME TITLE
