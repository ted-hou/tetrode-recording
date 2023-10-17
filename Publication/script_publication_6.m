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
p.heightThirdRow = 3;

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
        figUnits='inches', figPosition=[0, 0, p.width/length(unitNames), p.heightFirstRow], figMargins=1.2*[0.16, 0.09]);
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
fig = figure(Units='inches', Position=[0, 0, p.width, p.heightSecondRow], DefaultAxesFontSize=p.fontSize);
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

fig = figure(Units='inches', Position=[0, 0, p.width, p.height - p.heightFirstRow - p.heightSecondRow], DefaultAxesFontSize=p.fontSize);
ax = arrayfun(@(i) subplot(1, 3, i), 1:3);
hold(ax, 'on')

% 


hold(ax, 'off')

%% 6d. Salt and pepper map (lick vs. reach vs. osci lick)