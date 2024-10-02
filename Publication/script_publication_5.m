%%% Figure 6. Lick vs. Reach ephys + Osci lick and somatotopy

%% Load all units
load_ephysunits;
% boot_response_dir; % included in load('C:\SERVER\Units\meta_Lite_NonDuplicate_NonDrift.mat')
% load('C:\SERVER\bootMoveResponse_20240830.mat') 
find_osci_lick_circ;
[bootPressVsLick.h, bootPressVsLick.p, bootPressVsLick.ci, bootPressVsLick.obs] = bootstrapPressVsLick(eu(c.hasLick & c.hasPress), ...
    pressWindow=p.metaWindowPress, lickWindow=p.metaWindowLick);
c.isPressVsLickSelective = false(1, length(eu));
c.isPressVsLickSelective(c.hasLick & c.hasPress) = bootPressVsLick.h;

%% plot params
nEgUnits = 3;
p.fontSize = 9;
p.etaSortWindow = [-3, 0];
p.etaSignWindow = [-0.3, 0];
p.etaLatencyThresholdPos = 0.25;
p.etaLatencyThresholdNeg = 0.25;

layout.w = 7;
layout.h = 8;
layout.top.h = 3;
layout.middle.h = 3;
layout.middle.left.w = 3;
layout.middle.right.w = 2;
layout.middle.right.top.h = 8;
layout.middle.right.bottom.h = 3;
layout.bottom.h = 2;

close all
fig = figure(Units='inches', Position=[0, 0, layout.w, layout.h]);
layout.tl = tiledlayout(fig, layout.top.h + layout.middle.h + layout.bottom.h, 1, TileSpacing='loose');

layout.top.tl = tiledlayout(layout.tl, 2, nEgUnits, TileSpacing='compact', Padding='compact');
l = layout.top.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [layout.top.h, 1];

layout.middle.tl = tiledlayout(layout.tl, 1, layout.middle.left.w + layout.middle.right.w, TileSpacing='loose', Padding='loose');
l = layout.middle.tl; l.Layout.Tile = 1 + layout.top.h; l.Layout.TileSpan = [layout.middle.h, 1];

layout.middle.left.tl = tiledlayout(layout.middle.tl, 1, 2, TileSpacing='compact', Padding='compact');
l = layout.middle.left.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [1, layout.middle.left.w];

layout.middle.right.tl = tiledlayout(layout.middle.tl, layout.middle.right.top.h + layout.middle.right.bottom.h, 1, TileSpacing='compact', Padding='compact');
l = layout.middle.right.tl; l.Layout.Tile = 1 + layout.middle.left.w; l.Layout.TileSpan = [1, layout.middle.right.w];

layout.bottom.tl = tiledlayout(layout.tl, 1, 3, TileSpacing='compact', Padding='compact');
l = layout.bottom.tl; l.Layout.Tile = 1 + layout.top.h + layout.middle.h; l.Layout.TileSpan = [layout.bottom.h, 1];


% 5a. Double rasters (reach vs lick) 4 example units
unitNames = { ...
    'desmond24_20220510_Channel44_Unit1'; ...
    'Daisy2_20180420_Channel14_Unit1'; ...
    'daisy13_20220106_Electrode97_Unit1'; ...
%     'daisy8_20210709_Channel7_Unit1'; ...
    };
files = cellfun(@(name) sprintf('C:\\SERVER\\Units\\Lite_NonDuplicate\\%s.mat', name), unitNames, UniformOutput=false);
euEg = EphysUnit.load(files);
assert(nEgUnits == length(unitNames))

AX = gobjects(2, nEgUnits);
for i = 1:2
    for iEu = 1:nEgUnits
        AX(i, iEu) = nexttile(layout.top.tl);
    end
end
TRIALTYPE = {'press'; 'lick'};
TITLE = ["Reach trials"; "Lick trials"];
XLABEL = ["Time to bar-contact(s)"; "Time to spout-contact (s)"];
MINTRIALDURATION = [...
        2, 2, 2; ...
        2, 2, 2; ...
    ];
MAXTRIALDURATION = [...
        Inf, Inf, Inf; ...
        Inf, Inf, Inf; ...
    ];
EVERYNTH = [4, 4, 4];
YLIM = {[0, 30], [20, 140], [20, 80]};
for iEu = 1:nEgUnits
    for i = 1:2
        ax = AX(i, iEu);
        trials = euEg(iEu).getTrials(TRIALTYPE{i});
        trials = trials(trials.duration() >= MINTRIALDURATION(i, iEu) & trials.duration() <= MAXTRIALDURATION(i, iEu));
        thisRD = euEg(iEu).getRasterData(TRIALTYPE{i}, window=[0, 2], sort=true, trials=trials);
        thisETA = euEg(iEu).getETA('count', TRIALTYPE{i}, p.etaWindow, normalize='none', trials=trials, includeInvalid=false);
        yyaxis(ax, 'right')
        EphysUnit.plotRaster(ax, thisRD, xlim=[-4,2], iti=false, sz=1, maxTrials=40, maxTrialsMethod='uniformsample', ...
            sz=1, everyNth=EVERYNTH(iEu), timingCriterion=4);
        ylabel(ax, '')
        yticks(ax, [])
        ax.YAxis(2).Direction = 'reverse';
        yyaxis(ax, 'left')
        plot(ax, thisETA.t, thisETA.X./0.1, LineWidth=1.5, Color='black')
        hold(ax, 'on')
        set(ax.YAxis, FontSize=p.fontSize, Color=[0.15, 0.15, 0.15]);
        ylabel(ax, 'Spike rate (sp/s)')
        delete(ax.Legend)
        title(ax, TITLE(i))
        xlabel(ax, XLABEL{i})
        plot(ax, [0, 0], [0, 100], 'k--')
        ylim(ax, YLIM{iEu})
        fontsize(ax, p.fontSize, 'points')
    end
end
% xlabel(AX, '')
ylabel(AX, '')
% xlabel(layout.top.tl, 'Time to bar/spout contact (s)', FontSize=p.fontSize)
ylabel(layout.top.tl, 'Spike rate (sp/s)', FontSize=p.fontSize)
clear thisRD ax iEu AX

% 5b/c. PETH Reach vs. Lick vs. Osci Lick
% Reach, Lick, Osci Lick, sorted amongst themselves
% close all
ax = gobjects(1, 3);
ax(1) = nexttile(layout.middle.left.tl);
ax(2) = nexttile(layout.middle.left.tl);
ax(3) = nexttile(layout.middle.right.tl, [layout.middle.right.top.h, 1]);
EphysUnit.plotETA(ax(1), eta.press, c.hasPress & c.hasLick, ...
    clim=[-2, 2], xlim=[-4, 0], sortWindow=p.etaSortWindow, signWindow=p.etaSignWindow, ...
    sortThreshold=p.etaLatencyThresholdPos, negativeSortThreshold=p.etaLatencyThresholdNeg, hidecolorbar=true);
EphysUnit.plotETA(ax(2), eta.lick, c.hasPress & c.hasLick, ...
    clim=[-2, 2], xlim=[-4, 0], sortWindow=p.etaSortWindow, signWindow=p.etaSignWindow, ...
    sortThreshold=p.etaLatencyThresholdPos, negativeSortThreshold=p.etaLatencyThresholdNeg, hidecolorbar=true);

% Visualize onset
% hold(ax(1:2), 'on')
% plot(ax(1), onset.press(order.press), 1:length(order.press))
% plot(ax(2), onset.lick(order.lick), 1:length(order.lick))

sel = c.hasPress & c.hasLick & c.isLick; 
meanLickInterval = mean(cat(2, durations{:}));
meanBinWidth = meanLickInterval / length(eta.circLick.t);
eta.circLick.Z = eta.circLick.X.*exp(eta.circLick.t*1i)./meanBinWidth;
eta.circLick.Z(:, 1) = mean(eta.circLick.X(:, [2, 30]), 2).*exp(eta.circLick.t(1)*1i)./meanBinWidth;
meanZ = mean(eta.circLick.Z, 2);
eta.lickBoutNorm = eta.lickBout;
eta.lickBoutNorm.X = normalize(eta.lickBout.X, 2, 'zscore', 'robust');

[~, I] = sort(angle(meanZ(sel)));
[~, ~] = EphysUnit.plotETA(ax(3), eta.lickBoutNorm, sel, order=I, ...
    clim=[-2, 2], xlim=[0, 2*pi*maxBoutCycles], hidecolorbar=true);
xticks (ax(3), (0:2:8).*pi);
xticklabels(ax(3), [{'0'}, arrayfun(@(x) sprintf('%i\\pi', x), 2:2:8, UniformOutput=false)]);
title(ax(1), 'Pre-reach')
title(ax(2), 'Pre-lick')
title(ax(3), 'Lick-entrained')
ylabel(layout.middle.left.tl, 'Unit', FontSize=p.fontSize)
ylabel(ax(3), 'Unit')
ylabel(ax(1:2), '')
xlabel(ax(1), 'Time to bar contact (s)')
xlabel(ax(2), 'Time to spout contact (s)')
xlabel(ax(3), 'Lick phase')
fontsize(ax, p.fontSize, 'points')
fontname(ax, 'Arial')
axc = ax(2);

%5d. Peri-lick lick prob histogram
ax = nexttile(layout.middle.right.tl, [layout.middle.right.bottom.h, 1]);
histogram(ax, BinEdges=lickHistEdges(2:end), BinCounts=lickHistCounts(2:end), Normalization='probability')
yl = ax.YLim;
histogram(ax, BinEdges=lickHistEdges(1:end), BinCounts=lickHistCounts(1:end), Normalization='probability', EdgeColor='none', FaceColor='black', FaceAlpha=1)
xlim(ax, [0, 1/9*4]); % Assuming 9 Hz, 5 licks (4 cycles)
xticks(ax, 0:0.2:0.5)
yticks(ax, [])
xlabel(ax, 'Time from first lick (s)')
ylabel(ax, 'lick')
set(ax, FontSize=p.fontSize, FontName='Arial');
hold(ax, 'off')

% 5e. Reach vs. Lick (scatter), Reach vs. Lick (osci subset, scatter), Pie-chart
sz = 3;
SEL = { ...
    c.hasLick & c.hasPress, ...
    c.hasLick & c.hasPress & c.isLick, ...
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
    'Pre-lick response (a.u.)', ...
    'Pre-lick response (a.u.)'; ...
    };
YNAME = { ...
    'Pre-reach response (a.u.)', ...
    'Pre-reach response (a.u.)'; ...
    };
TITLE = { ...
    'All units', ...
    'Lick-entrained units'
    };
LEGENDPOS = {...
    [0.244956331030447,0.173453205894066,0.117187498486601,0.122685182149764], ...
    [0.513185497503088,0.182712465153325,0.11067708201396,0.122685182149764], ...
    };

h = gobjects(2, 3);

AX = gobjects(1, 2);
% 1. Reach vs. Lick scatter
for i = 1:2
    ax = nexttile(layout.bottom.tl);
    hold(ax, 'on')
    sel = SEL{i};
    x = XDATA{i}(sel);
    y = YDATA{i}(sel);
    subselResp = sel & c.isPressVsLickSelective;
    subselNone = sel & ~c.isPressVsLickSelective;
    h(i, 1) = scatter(ax, XDATA{i}(subselNone), YDATA{i}(subselNone), sz, 'black', DisplayName=sprintf('%g units', nnz(subselNone)));   
    h(i, 2) = scatter(ax, XDATA{i}(subselResp), YDATA{i}(subselResp), sz, 'red', DisplayName=sprintf('%g units', nnz(subselResp)));

    plot(ax, [-10, 10], [0, 0], 'k:');
    plot(ax, [0, 0], [-10, 10], 'k:');
    plot(ax, [-10, 10], [-10, 10], 'k:')

    if i == 1
        xlabel(ax, XNAME{i}, Units='normalized', Position=[1.18902486151029,-0.115261042292962,0])
        ylabel(ax, YNAME{i})
    else
    end
    axis(ax, 'equal')
    xlim(ax, [-2, 6])
    ylim(ax, [-2, 6])

    title(ax, TITLE{i})
end

% Reach vs. Lick onset times
ax = nexttile(layout.bottom.tl);
hold(ax, 'on')
scatter(ax, onset.press(c.hasPress & c.hasLick), onset.lick(c.hasPress & c.hasLick), sz, 'black')
scatter(ax, onset.press(c.hasPress & c.hasLick & c.isPressVsLickSelective), onset.lick(c.hasPress & c.hasLick & c.isPressVsLickSelective), sz, 'red')
plot(ax, [-2, 0], [-2, 0], 'k:')
hold(ax, 'off')
axis(ax, 'equal')
xlim(ax, [-2, 0])
ylim(ax, [-2, 0])
xlabel('Reach resp onset (s)')
ylabel('Lick resp onset (s)')

h = colorbar(axc); 
h.Layout.Tile = 'east';
h.Label.String = 'normalized spike rate (a.u.)';
copygraphics(fig, ContentType='vector')

%% Supplement S5
close all
% All sign-changing units
p.fontSize = 9;
signChangingUnitIndices = [find(c.isPressUp & c.isLickDown), find(c.isPressDown & c.isLickUp)];
nRows = 7;
nCols = 5;
assert(length(signChangingUnitIndices) <= nRows*nCols)
fig = figure(Units='inches', Position=[1 1 4 6]);
tl = tiledlayout(fig, nRows, nCols, TileSpacing='tight', Padding='tight');
iAx = 0;
for i = signChangingUnitIndices
    ax = nexttile(tl);
    iAx = iAx + 1;
    hold(ax, 'on')
    plot(ax, eta.pressRaw.t, eta.pressRaw.X(i, :)./0.1, 'r', LineWidth=1.5)
    plot(ax, eta.lickRaw.t, eta.lickRaw.X(i, :)./0.1, 'b', LineWidth=1.5)
    plot(ax, [0, 0], [0, 100], 'k:')
    hold(ax, 'off')
    xlim(ax, [-4, 0])
    ymedian = median([mean(eta.pressRaw.X(i, :), 1, 'omitnan')./0.1, mean(eta.lickRaw.X(i, :), 1, 'omitnan')./0.1]);
    ylim(ax, ymedian + [-40, 40])
    yticks(ax, round(ymedian/10)*10 + [-30, 0, 30])
    fontsize(ax, p.fontSize, 'points')
end
xlabel(tl, 'Time to bar/spout contact (s)', FontSize=p.fontSize)
ylabel(tl, 'Spike rate (sp/s)', FontSize=p.fontSize)

% Mean ETA grouped by press/reach up/down
SEL_HASTRIAL = c.hasPress & c.hasLick;
SEL_I = {c.isPressUp; ~c.isPressResponsive; c.isPressDown};
SEL_J = {c.isLickDown; ~c.isLickResponsive; c.isLickUp};

fig = figure(Units='inches', Position=[1, 1, 4, 3]);
tl = tiledlayout(fig, length(SEL_I), length(SEL_J), TileSpacing='compact', Padding='compact');
h = gobjects(2, 1);
for i = 1:length(SEL_I)
    for j = 1:length(SEL_J)
        ax = nexttile(tl);
        sel = SEL_HASTRIAL & SEL_I{i} & SEL_J{j};
        assert(numel(sel) == numel(eu))
        hold(ax, 'on')
        h(1) = plot(ax, eta.pressRaw.t, mean(eta.pressRaw.X(sel, :), 1, 'omitnan')./0.1, 'r', LineWidth=1.5, DisplayName='reach');
        h(2) = plot(ax, eta.lickRaw.t, mean(eta.lickRaw.X(sel, :), 1, 'omitnan')./0.1, 'b', LineWidth=1.5, DisplayName='lick');
        plot(ax, [0, 0], [0, 100], 'k:')
        hold(ax, 'off')
        title(ax, sprintf('n = %i', nnz(sel)))
        xlim(ax, [-4, 0.5])
        ylim(ax, [0, 80])
        fontsize(ax, p.fontSize, 'points')
        if i < length(SEL_I)
%             xticks(ax, [])
        else
%             xticks(ax, [-2, 0])
        end
        ax.XGrid = 'on';
        ax.XMinorTick = 'on';
        ax.XMinorGrid = 'on';
        ax.XAxis.MinorTickValues = -0.5:0.1:0.5;
        if j > 1
            yticks(ax, [])
        else
            yticks(ax, [10, 40, 70])
        end
    end
end
lgd = legend(h, Orientation='horizontal');
lgd.Layout.Tile = 'north';
xlabel(tl, 'Time to bar/spout contact (s)', FontSize=p.fontSize)
ylabel(tl, 'Spike rate (sp/s)', FontSize=p.fontSize)