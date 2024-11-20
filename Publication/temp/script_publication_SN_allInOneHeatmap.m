eta.lickBoutEnd = euMixed.getETA('count', 'lickboutend', window=[0, 2], resolution=[2*pi/5, 0.025], normalize='none', ...
    maxInterval=0.2, minInterval=0.05, minBoutCycles=6);
eta.correctLickBout = euMixed.getETA('count', 'lick+lickbout', window=[-4, 0], resolution=[0.025, 2*pi/5], normalize='none', minTrialDuration=4, ...
    maxInterval=0.2, minInterval=0.05, minBoutCycles=2, maxBoutCycles=6);
eta.correctPressBout = euMixed.getETA('count', 'press+lickbout', window=[-4, 0], resolution=[0.025, 2*pi/32, 2*pi/5], normalize='none', minTrialDuration=4, ...
    maxInterval=0.2, minInterval=0.05, minBoutCycles=2, maxBoutCycles=6);
% eta.incorrectLickBout = euMixed.getETA('count', 'lick+lickbout', window=[-4, 0], resolution=[0.025, 2*pi/5], normalize='none', minTrialDuration=2, maxTrialDuration=4, ...
%     maxInterval=0.2, minInterval=0.05, minBoutCycles=2, maxBoutCycles=6);
% eta.pressCueRaw = eu.getETA('count', 'press', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize='none', includeInvalid=true, resolution=0.025);
% eta.lickCueRaw = eu.getETA('count', 'lick', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize='none', includeInvalid=true, resolution=0.025);


% etaFine.lickCorrectRaw = eu.getETA('count', 'lick', p.etaWindow, minTrialDuration=4, normalize='none', resolution=0.025);
% etaFine.lickIncorrectRaw = eu.getETA('count', 'lick', p.etaWindow, minTrialDuration=2, maxTrialDuration=4, normalize='none', resolution=0.025);
% etaFine.p.minTrialDuration = p.minTrialDuration;
% etaFine.p.norm = p.etaNorm;
% etaFine.p.resolution = 0.025;



% Normalized first lick/press + lick bout
eta.correctPressBoutNorm = euMixed.getETA('count', 'press+lickbout', window=[-4, 0], resolution=[0.025, 2*pi/32, 2*pi/5], normalize=[-4, -2], minTrialDuration=4, ...
    maxInterval=0.2, minInterval=0.05, minBoutCycles=2, maxBoutCycles=6);
eta.correctLickBoutNorm = euMixed.getETA('count', 'lick+lickbout', window=[-4, 0], resolution=[0.025, 2*pi/5], normalize=[-4, -2], minTrialDuration=4, ...
    maxInterval=0.2, minInterval=0.05, minBoutCycles=2, maxBoutCycles=6);

stats.correctPressBout_025 = arrayfun(@(stats) struct(mean=stats.mean .* 0.025, sd=stats.sd .* 0.025), eta.correctPressBoutNorm.stats);
stats.correctLickBout_025 = arrayfun(@(stats) struct(mean=stats.mean .* 0.025, sd=stats.sd .* 0.025), eta.correctLickBoutNorm.stats);

stats.correctPressBout = eta.correctPressBoutNorm.stats;
stats.correctLickBout = eta.correctLickBoutNorm.stats;

eta.pressCueNorm = eu.getETA('count', 'press', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize=stats.correctPressBout_025, includeInvalid=true, resolution=0.025);
eta.lickCueNorm = eu.getETA('count', 'lick', p.cueEtaWindow, alignTo='start', minTrialDuration=p.minTrialDuration, normalize=stats.correctLickBout_025, includeInvalid=true, resolution=0.025);

eta.lickBoutEndNorm = euMixed.getETA('count', 'lickboutend', window=[0, 2], resolution=[2*pi/5, 0.025], normalize=eta.correctLickBoutNorm.stats, ...
    maxInterval=0.2, minInterval=0.05, minBoutCycles=6);


%% 2x2 PETH
close all
% YL = {[-0.5, 2], [-2, 0.5]};
yl = [0, 120];
selCommon = c.hasPress & c.hasLick;
TASK = ["press", "lick"];
ETACUE = {eta.pressCueRaw, eta.lickCueRaw};
ETAPERIMOVE = {eta.correctPressBout, eta.correctLickBout};
ETABOUTEND = {eta.lickBoutEnd, eta.lickBoutEnd};
SEL = {selCommon & c.isPressUp, selCommon & c.isPressUp; ...
    selCommon & c.isPressDown, selCommon & c.isPressDown};
SUBSEL = {c.isLickUp, c.isLickDown};
TITLE = ["Self-timed Reach", "Self-timed Lick"];
LABEL1 = ["Reach Inc", "Reach Inc"; "Reach Dec", "Reach Dec"];
LABEL2 = ["lick inc", "lick dec"];
COLORS = 'rb';

nRows = size(SEL, 1);
nCols = size(SEL, 2);
clear layout
layout.wp = [0, sum([1.8, 2+0.8+6/8, 1+6/8]*100), sum([1.8, 2+6/8, 1+6/8]*100)];

% Plot ETA aligned to cue (i.e. flinchiness) vs. aligned to lever touch
% close all
fig = figure(Units='inches', Position=[1 1 7 3.5], AutoResizeChildren='off');
tlp = tiledlayout(fig, nRows, sum(layout.wp), TileSpacing='compact', Padding='compact');

AX = gobjects(nRows, nCols, 3);
hLine = gobjects(2, 2, 2);
hPatch = gobjects(2, 2, 2);
hPatchDummy = gobjects(2, 2, 2);
TL = gobjects(2, 2);
for iRow = 1:nRows
    for iCol = 1:nCols
        switch TASK(iCol)
            case "press"
                layout.w = [1.8, 2+0.8+6/8, 1+6/8]*100;
            case "lick"
                layout.w = [1.8, 2+6/8, 1+6/8]*100;
        end
        layout.gap = [2, 2]*10;

        tl = tiledlayout(tlp, 1, sum(layout.w) + sum(layout.gap), TileSpacing='none', Padding='tight');
        tl.Layout.Tile = (iRow-1)*sum(layout.wp) + 1 + sum(layout.wp(1:iCol));
        tl.Layout.TileSpan = [1, layout.wp(iCol + 1)];
        TL(iRow, iCol) = tl;

        bgAx = axes(tl, XTick=[], YTick=[], Box='off');
        bgAx.Layout.Tile = 1;
        bgAx.Layout.TileSpan = [1, sum(layout.w)];
        if iRow == 1
            title(bgAx, TITLE(iCol));
        end

        % 1. Peri-cue
        ax = axes(tl); AX(iRow, iCol, 1) = ax;
        ax.Layout.Tile = 1;
        ax.Layout.TileSpan = [1, layout.w(1)];
        sel = SEL{iRow, iCol};

        hold(ax, 'on')
        t = ETACUE{iCol}.t;
        tt = [t, flip(t)];

        for iSubSel = 1:2
            subSel = sel & SUBSEL{iSubSel};
            X = ETACUE{iCol}.X(subSel, :)./0.025;
            err = std(X, 0, 1, 'omitnan') ./ sqrt(size(X, 1));
            xx = [mean(X, 1, 'omitnan') + err, flip(mean(X, 1, 'omitnan') - err)];

            plot(ax, t, mean(X, 1, 'omitnan'), COLORS(iSubSel), LineWidth=1.5);
            patch(ax, tt(~isnan(xx)), xx(~isnan(xx)), COLORS(iSubSel), FaceAlpha=0.25, EdgeColor='none') 
        end

        switch TASK(iCol)
            case "press"
                name = 'bar deploy';
            case "lick"
                name = 'spout deploy';
        end
        hPatch(iRow, iCol, 1) = patch(ax, XData=[-0.8, 0, 0, -0.8], YData=[yl(1), yl(1), yl(2), yl(2)], FaceColor=[0.15, 0.15, 0.15], FaceAlpha=0.1, EdgeColor='none', DisplayName=name);
        hPatch(iRow, iCol, 2) = patch(ax, XData=[0, 0.1, 0.1, 0], YData=[yl(1), yl(1), yl(2), yl(2)], FaceColor=[0.8, 0.2, 0.2], FaceAlpha=0.3, EdgeColor='none', DisplayName='tone');

        ax.Box = 'off';
        xlim(ax, [-0.95, 1])
        xticks(ax, [-0.8, 0, 1])
        xline(ax, 0, '-')
        xlabel(ax, 'start cue')
        ax.XAxis.TickLabelRotation = 0;
        
        % 2. Peri-move
        ax = axes(tl); AX(iRow, iCol, 2) = ax;
        ax.Layout.Tile = sum(layout.w(1)) + 1 + layout.gap(1);
        ax.Layout.TileSpan = [1, layout.w(2)];
        hold(ax, 'on')
        t = ETAPERIMOVE{iCol}.t;
        switch TASK(iCol)
            case 'press'
                t(t>0 & t<=2*pi) = t(t>0 & t<=2*pi) / (2*pi) * 0.8;
                t(t>2*pi) = (t(t>2*pi) - 2*pi) ./ (2*pi) / 8 + 0.8;
            case 'lick'
                t(t>0) = t(t>0) ./ (2*pi) / 8;
        end
        tt = [t, flip(t)];

        for iSubSel = 1:2
            subSel = sel & SUBSEL{iSubSel};
            X = ETAPERIMOVE{iCol}.X(subSel, :);
            err = std(X, 0, 1, 'omitnan') ./ sqrt(size(X, 1));
            xx = [mean(X, 1, 'omitnan') + err, flip(mean(X, 1, 'omitnan') - err)];

            hLine(iRow, iCol, iSubSel) = plot(ax, t, mean(X, 1, 'omitnan'), COLORS(iSubSel), LineWidth=1.5, DisplayName=sprintf('%s (n=%i)', LABEL2(iSubSel), nnz(subSel)));
            patch(ax, tt(~isnan(xx)), xx(~isnan(xx)), COLORS(iSubSel), FaceAlpha=0.25, EdgeColor='none')
        end
        hPatchDummy(iRow, iCol, 1) = patch(ax, XData=[NaN], YData=[NaN], FaceColor=[0.15, 0.15, 0.15], FaceAlpha=0.1, EdgeColor='none', DisplayName=name);
        hPatchDummy(iRow, iCol, 2) = patch(ax, XData=[NaN], YData=[NaN], FaceColor=[0.8, 0.2, 0.2], FaceAlpha=0.3, EdgeColor='none', DisplayName='tone');
        

        switch TASK(iCol)
            case 'press'
                xline(ax, [0.8,  0.8+(1:5)/8], LineStyle=':')
                xlim(ax, [-2, 0.8+5/8])
                xticks(ax, [-2, -1, 0, 0.8,  0.8+(1:5)/8])
                xticklabels(ax, {'-2', '-1', '        0      \newlinebar contact', '2\pi', '', '', '', '', '12\pi'})
            case 'lick'        
                xline(ax, (1:6)/8, LineStyle=':')       
                xlim(ax, [-2, 6/8])
                xticks(ax, [-2:1:0, (1:6)/8])
                xticklabels(ax, {'-2', '-1', '          0      \newlinespout contact', '', '', '', '', '', '12\pi'}) 
        end
        ax.XAxis.TickLabelRotation = 0;
        xline(ax, 0, '-')
        ax.Box = 'off';
        ax.YAxis.Visible = 'off';
    
        % 3. End of bout
        ax = axes(tl); AX(iRow, iCol, 3) = ax;
        ax.Layout.Tile = sum(layout.w(1:2)) + 1 + sum(layout.gap(1:2));
        ax.Layout.TileSpan = [1, layout.w(3)];
        hold(ax, 'on')
        
        t = eta.lickBoutEnd.t;
        t(t<0) = t(t<0) ./ (2*pi) / 8;
        tt = [t, flip(t)];

        for iSubSel = 1:2
            subSel = sel & SUBSEL{iSubSel};
            X = ETABOUTEND{iCol}.X(subSel, :);
            err = std(X, 0, 1, 'omitnan') ./ sqrt(size(X, 1));
            xx = [mean(X, 1, 'omitnan') + err, flip(mean(X, 1, 'omitnan') - err)];

            plot(ax, t, mean(X, 1, 'omitnan'), COLORS(iSubSel), LineWidth=1.5);
            patch(ax, tt(~isnan(xx)), xx(~isnan(xx)), COLORS(iSubSel), FaceAlpha=0.25, EdgeColor='none')
        end

        yline(ax, 0, 'k')
        xlim(ax, [-6/8, 1])
        xticks(ax, [(-6:-1)/8, 0:1])
        xline(ax, (-6:-1)/8, LineStyle=':')       
        xticklabels(ax, {'', '', '', '-6\pi', '', '', '0', '1'})
        ax.XAxis.TickLabelRotation = 0;
        xlabel(ax, 'last lick')
        ax.Box = 'off';
        ax.YAxis.Visible = 'off';
        xline(ax, 0, '-')

        ylim(AX(iRow, iCol, :), yl)
        linkaxes(AX(iRow, iCol, :), 'y')
    end
end

xlabel(tlp, 'Time (s) or lick phase', FontSize=p.fontSize)
ylabel(tlp, 'Spike rate (sp/s)', FontSize=p.fontSize)
fontsize(AX(:), p.fontSize, 'points')


hLineDummy = gobjects(2, 1);
for iRow = 1:2
    for iCol = 1:2
        axDummy = axes(fig, Position=TL(iRow, iCol).Position);
        hold(axDummy, 'on')
        subSel = SEL{iRow, 1} & SUBSEL{1};
        hLineDummy(1) = plot(axDummy, NaN, NaN, 'red', LineWidth=1.5, DisplayName=sprintf('%s (n=%i)', LABEL2(1), nnz(subSel)));
        subSel = SEL{iRow, 1} & SUBSEL{2};
        hLineDummy(2) = plot(axDummy, NaN, NaN, 'blue', LineWidth=1.5, DisplayName=sprintf('%s (n=%i)', LABEL2(2), nnz(subSel)));
        axDummy.Visible = 'off';
        legend(axDummy, hLineDummy, Location='northwest')
    end
end

copygraphics(fig, ContentType='vector', BackgroundColor='none')


%% Heatmap
% close all

CLIM = [-2.5, 2.5];
ETACUE = {eta.pressCueNorm, eta.lickCueNorm};
ETAPERIMOVE = {eta.correctPressBoutNorm, eta.correctLickBoutNorm};
ETABOUTEND = {eta.lickBoutEndNorm, eta.lickBoutEndNorm};
CONTACTNAME = ["bar contact", "spout contact"];

fig = figure(Units='inches', Position=[1, 1, 14, 5]);
sel = c.hasPress(:) & c.hasLick(:) & c.isPressResponsive(:) & c.isLickResponsive(:);

clear layout
layout.wp = [0, sum([1.8, 2+0.8+6/8, 1+6/8]*100), sum([1.8, 2+6/8, 1+6/8]*100)];
layout.w{1} = [1.8, 2+0.8+6/8, 1+6/8]*100;
layout.w{2} = [1.8, 2+0.8+6/8, 1+6/8]*100;
layout.gap = [2, 2]*10;

tlp = tiledlayout(fig, 1, sum(layout.wp), TileSpacing='compact', Padding='compact');

AX = gobjects(iCol, 3);
clear order
for iCol = 1:2
    tl = tiledlayout(tlp, 1, sum(layout.w{iCol}) + sum(layout.gap), TileSpacing='none', Padding='tight');
    tl.Layout.Tile = 1 + sum(layout.wp(1:iCol));
    tl.Layout.TileSpan = [1, layout.wp(iCol + 1)];
    
    bgAx = axes(tl, XTick=[], YTick=[], Box='off');
    bgAx.Layout.Tile = 1;
    bgAx.Layout.TileSpan = [1, sum(layout.w{iCol})];
    if iRow == 1
        title(bgAx, TITLE(iCol));
    end
    
    % 2. Peri-move
    ax = axes(tl); AX(iCol, 2) = ax;
    ax.Layout.Tile = sum(layout.w{iCol}(1)) + 1 + layout.gap(1);
    ax.Layout.TileSpan = [1, layout.w{iCol}(2)];
    hold(ax, 'on')

    thisETA = ETAPERIMOVE{iCol};
    t = thisETA.t;
    switch iCol
        case 1
            t(t>0 & t<=2*pi) = t(t>0 & t<=2*pi) / (2*pi) * 0.8;
            t(t>2*pi) = (t(t>2*pi) - 2*pi) ./ (2*pi) / 8 + 0.8;
        case 2
            t(t>0) = t(t>0) ./ (2*pi) / 8;
    end
    thisETA.t = t;
    thisETA.X = thisETA.X(sel, :);
    thisETA.N = thisETA.N(sel);
    thisETA.D = thisETA.D(sel);
    thisETA.stats = thisETA.stats(sel);
    groupVar = NaN(length(sel), 1);
    groupVar(c.isPressDown & c.isLickDown) = 1;
    groupVar(c.isPressUp & c.isLickDown) = 0;
    groupVar(c.isPressDown & c.isLickUp) = 3;
    groupVar(c.isPressUp & c.isLickUp) = 2;
    groupVar = groupVar(sel);
    if ~exist('order', 'var')
        [~, order] = EphysUnit.plotETA(ax, thisETA, clim=CLIM, signWindow=[-0.3, 0], sortWindow=[-3, -0.1], sortThreshold=0.25, sortGroup=groupVar);
    else
        EphysUnit.plotETA(ax, thisETA, clim=[-2, 2], order=order);
    end

    ax.Box = 'off';
    ax.YAxis.Visible = 'off';
    colorbar(ax, 'off')


    switch iCol
        case 1
            xline(ax, [0.8,  0.8+(1:5)/8], LineStyle=':')
            xlim(ax, [-2, 0.8+5/8])
            xticks(ax, [-2, -1, 0, 0.8,  0.8+(1:5)/8])
            if iRow == 1
                xticklabels(ax, {'-2', '-1', '0', '2\pi', '', '', '', '', '12\pi'})
            else
                xticklabels(ax, {'-2', '-1', '        0      \newlinebar contact', '2\pi', '', '', '', '', '12\pi'})
            end
        case 2        
            xline(ax, (1:6)/8, LineStyle='--')       
            xlim(ax, [-2, 6/8])
            xticks(ax, [-2:1:0, (1:6)/8])
            if iRow == 1
                xticklabels(ax, {'-2', '-1', '0', '', '', '', '', '', '12\pi'}) 
            else
                xticklabels(ax, {'-2', '-1', '          0      \newlinespout contact', '', '', '', '', '', '12\pi'}) 
            end
    end
    
    xlabel(ax, 'Time (s) or lick phase')
    xline(ax, 0, '-')
    ax.XAxis.TickLabelRotation = 0;
    ax.YAxis.Direction = 'reverse';
    colorbar(ax, 'off')

    % 1. Peri-cue
    ax = axes(tl); AX(iCol, 1) = ax;
    ax.Layout.Tile = 1;
    ax.Layout.TileSpan = [1, layout.w{iCol}(1)];
    
    thisETA = ETACUE{iCol};
    thisETA.X = thisETA.X(sel, :);
    thisETA.N = thisETA.N(sel);
    thisETA.D = thisETA.D(sel);
    EphysUnit.plotETA(ax, thisETA, clim=CLIM, order=order)

    ax.Box = 'off';
    xlim(ax, [-0.95, 1])
    xticks(ax, [-0.8, 0, 1])
    xline(ax, 0, '-')
    xlabel(ax, 'start cue')
    ax.XAxis.TickLabelRotation = 0;
    ax.YAxis.Direction = 'reverse';
    colorbar(ax, 'off')
    yticks(ax, 0:50:300)

    % 3. End of bout
    ax = axes(tl); AX(iCol, 3) = ax;
    ax.Layout.Tile = sum(layout.w{iCol}(1:2)) + 1 + sum(layout.gap(1:2));
    ax.Layout.TileSpan = [1, layout.w{iCol}(3)];
    hold(ax, 'on')
    
    thisETA = ETABOUTEND{iCol};
    t = thisETA.t;
    t(t<0) = t(t<0) ./ (2*pi) / 8;
    tt = [t, flip(t)];
    thisETA.t = t;
    thisETA.X = thisETA.X(sel, :);
    thisETA.N = thisETA.N(sel);
    thisETA.D = thisETA.D(sel);
    EphysUnit.plotETA(ax, thisETA, clim=CLIM, order=order)

    yline(ax, 0, 'k')
    xlim(ax, [-6/8, 1])
    xticks(ax, [(-6:-1)/8, 0:1])
    xline(ax, (-6:-1)/8, LineStyle=':')       
    xticklabels(ax, {'', '', '', '-6\pi', '', '', '0', '1'})
    ax.XAxis.TickLabelRotation = 0;
    ax.YAxis.Direction = 'reverse';
    xlabel(ax, 'last lick')
    ax.Box = 'off';
    ax.YAxis.Visible = 'off';
    xline(ax, 0, '-')

    ylim(AX(iCol, :), [0, nnz(sel)])
    linkaxes(AX(iCol, :), 'y')
    if iCol == 1
        colorbar(ax, 'off')
    else
        ax.Colorbar.Layout.Tile = 'east';
    end
end

title(AX, '')
title(AX(1, 2), 'Self-timed Reach')
title(AX(2, 2), 'Self-timed Lick (same order)')

%% Heatmap (w/o time warping/lick bouts)
close all
clear order

for iMoveType = 1:3
    CLIM = [-2, 2];
    ETACUE = {eta.pressCueNorm, eta.lickCueNorm};
    XLCUE = {[-0.95, 0.5], [-0.95, 0.5]};
    XLPERIMOVE = {[-2, 2], [-2, 2]};
    switch iMoveType
        case 1
            ETAPERIMOVE = {etaFine.press, etaFine.lick};
            TITLE = ["Self-timed reach", "Self-timed lick"];
        case 2
            ETAPERIMOVE = {etaFine.correctPress, etaFine.correctLick};
            TITLE = ["Self-timed reach (rewarded)", "Self-timed lick (rewarded)"];
        case 3
            ETAPERIMOVE = {etaFine.incorrectPress, etaFine.incorrectLick};
            TITLE = ["Self-timed reach (unrewarded)", "Self-timed lick (unrewarded)"];
    end
    CONTACTNAME = ["bar contact", "spout contact"];
    
    fig = figure(Units='normalized', Position=[0, 0.025 + (iMoveType-1)*0.3, 0.5, 0.3]);
    sel = c.hasPress(:) & c.hasLick(:);
    
    clear layout
    layout.wp = [0, sum([sum(abs(XLCUE{1})), sum(abs(XLPERIMOVE{1}))]*100), sum([sum(abs(XLCUE{2})), sum(abs(XLPERIMOVE{2}))]*100)];
    layout.w{1} = [sum(abs(XLCUE{1})), sum(abs(XLPERIMOVE{1}))]*100;
    layout.w{2} = [sum(abs(XLCUE{2})), sum(abs(XLPERIMOVE{2}))]*100;
    layout.gap = 2*10;
    
    tlp = tiledlayout(fig, 1, sum(layout.wp), TileSpacing='compact', Padding='compact');
    
    AX = gobjects(iCol, 2);
    for iCol = 1:2
        tl = tiledlayout(tlp, 1, sum(layout.w{iCol}) + sum(layout.gap), TileSpacing='none', Padding='tight');
        tl.Layout.Tile = 1 + sum(layout.wp(1:iCol));
        tl.Layout.TileSpan = [1, layout.wp(iCol + 1)];
        
        bgAx = axes(tl, XTick=[], YTick=[], Box='off');
        bgAx.Layout.Tile = 1;
        bgAx.Layout.TileSpan = [1, sum(layout.w{iCol})];
        if iRow == 1
            title(bgAx, TITLE(iCol));
        end
        
        % 2. Peri-move
        ax = axes(tl); AX(iCol, 2) = ax;
        ax.Layout.Tile = sum(layout.w{iCol}(1)) + 1 + layout.gap(1);
        ax.Layout.TileSpan = [1, layout.w{iCol}(2)];
        hold(ax, 'on')
    
        thisETA = ETAPERIMOVE{iCol};
        t = thisETA.t;
        thisETA.t = t;
        thisETA.X = thisETA.X(sel, :);
        thisETA.N = thisETA.N(sel);
        thisETA.D = thisETA.D(sel);
        thisETA.stats = thisETA.stats(sel);
        groupVar = NaN(length(sel), 1);
        groupVar(c.isPressUnresponsiveButDown & c.isLickUp) = 0;
        groupVar(c.isPressUnresponsiveButDown & c.isLickUnresponsiveButUp) = 0;
        groupVar(c.isPressDown & c.isLickUp) = 0;
        groupVar(c.isPressDown & c.isLickUnresponsiveButUp) = 0;
        groupVar(c.isPressUp & c.isLickUnresponsiveButDown) = 1;
        groupVar(c.isPressUp & c.isLickDown) = 1;
        groupVar(c.isPressUnresponsiveButUp & c.isLickUnresponsiveButDown) = 1;
        groupVar(c.isPressUnresponsiveButUp & c.isLickDown) = 1;

        groupVar(c.isPressUnresponsiveButDown & c.isLickUnresponsiveButDown) = 2;
        groupVar(c.isPressUnresponsiveButDown & c.isLickDown) = 2;
        groupVar(c.isPressDown & c.isLickUnresponsiveButDown) = 2;
        groupVar(c.isPressDown & c.isLickDown) = 2;
        groupVar(c.isPressUp & c.isLickUp) = 3;
        groupVar(c.isPressUp & c.isLickUnresponsiveButUp) = 3;
        groupVar(c.isPressUnresponsiveButUp & c.isLickUp) = 3;
        groupVar(c.isPressUnresponsiveButUp & c.isLickUnresponsiveButUp) = 3;
        groupVar = groupVar(sel);
        N = histcounts(groupVar, -0.5:2:3.5);
        yline(ax, cumsum(N(1:end-1)) + 1, '--');
        if ~exist('order', 'var')
            [~, order] = EphysUnit.plotETA(ax, thisETA, clim=CLIM, signWindow=[-0.3, 0], sortWindow=[-3, -0.1], sortThreshold=0.25, sortGroup=groupVar);
        else
            EphysUnit.plotETA(ax, thisETA, clim=CLIM, order=order);
        end
    
        ax.Box = 'off';
        ax.YAxis.Visible = 'off';
        colorbar(ax, 'off')
        xlim(ax, XLPERIMOVE{iCol})
    
        xlabel(ax, 'Time (s) or lick phase')
        xline(ax, 0, '-')
        ax.XAxis.TickLabelRotation = 0;
        ax.YAxis.Direction = 'reverse';
        colorbar(ax, 'off')
    
        % 1. Peri-cue
        ax = axes(tl); AX(iCol, 1) = ax;
        ax.Layout.Tile = 1;
        ax.Layout.TileSpan = [1, layout.w{iCol}(1)];
        
        thisETA = ETACUE{iCol};
        thisETA.X = thisETA.X(sel, :);
        thisETA.N = thisETA.N(sel);
        thisETA.D = thisETA.D(sel);
        EphysUnit.plotETA(ax, thisETA, clim=CLIM, order=order)
    
        ax.Box = 'off';
        xlim(ax, XLCUE{iCol})
        xticks(ax, [-0.8, 0, 1])
        xline(ax, 0, '-')
        xlabel(ax, 'start cue')
        ax.XAxis.TickLabelRotation = 0;
        ax.YAxis.Direction = 'reverse';
        yticks(ax, [1, 50:50:427, 427])
    
    
        ylim(AX(iCol, :), [1, nnz(sel)+1])
        linkaxes(AX(iCol, :), 'y')
        if iCol == 1
            colorbar(ax, 'off')
        else
            ax.Colorbar.Layout.Tile = 'east';
        end
    end
    
    title(AX, '')
    title(AX(1, 2), TITLE(1))
    title(AX(2, 2), TITLE(2))
end

% Heatmap (w/o time warping/lick bouts)
% close all
clear order

for iMoveType = 1:3
    CLIM = [-2, 2];
    ETACUE = {eta.pressCueNorm, eta.lickCueNorm};
    XLCUE = {[-0.95, 0.5], [-0.95, 0.5]};
    XLPERIMOVE = {[-2, 2], [-2, 2]};
    switch iMoveType
        case 1
            ETAPERIMOVE = {etaFine.press, etaFine.lick};
            TITLE = ["Self-timed reach", "Self-timed lick"];
        case 2
            ETAPERIMOVE = {etaFine.correctPress, etaFine.correctLick};
            TITLE = ["Self-timed reach (rewarded)", "Self-timed lick (rewarded)"];
        case 3
            ETAPERIMOVE = {etaFine.incorrectPress, etaFine.incorrectLick};
            TITLE = ["Self-timed reach (unrewarded)", "Self-timed lick (unrewarded)"];
    end
    CONTACTNAME = ["bar contact", "spout contact"];
    
    fig = figure(Units='normalized', Position=[0.5, 0.025 + (iMoveType-1)*0.3, 0.5, 0.3]);
    sel = c.hasPress(:) & c.hasLick(:);
    
    clear layout
    layout.wp = [0, sum([sum(abs(XLCUE{1})), sum(abs(XLPERIMOVE{1}))]*100), sum([sum(abs(XLCUE{2})), sum(abs(XLPERIMOVE{2}))]*100)];
    layout.w{1} = [sum(abs(XLCUE{1})), sum(abs(XLPERIMOVE{1}))]*100;
    layout.w{2} = [sum(abs(XLCUE{2})), sum(abs(XLPERIMOVE{2}))]*100;
    layout.gap = 2*10;
    
    tlp = tiledlayout(fig, 1, sum(layout.wp), TileSpacing='compact', Padding='compact');
    
    AX = gobjects(iCol, 2);
    for iCol = 2:-1:1
        tl = tiledlayout(tlp, 1, sum(layout.w{iCol}) + sum(layout.gap), TileSpacing='none', Padding='tight');
        tl.Layout.Tile = 1 + sum(layout.wp(1:iCol));
        tl.Layout.TileSpan = [1, layout.wp(iCol + 1)];
        
        bgAx = axes(tl, XTick=[], YTick=[], Box='off');
        bgAx.Layout.Tile = 1;
        bgAx.Layout.TileSpan = [1, sum(layout.w{iCol})];
        if iRow == 1
            title(bgAx, TITLE(iCol));
        end
        
        % 2. Peri-move
        ax = axes(tl); AX(iCol, 2) = ax;
        ax.Layout.Tile = sum(layout.w{iCol}(1)) + 1 + layout.gap(1);
        ax.Layout.TileSpan = [1, layout.w{iCol}(2)];
        hold(ax, 'on')
    
        thisETA = ETAPERIMOVE{iCol};
        t = thisETA.t;
        thisETA.t = t;
        thisETA.X = thisETA.X(sel, :);
        thisETA.N = thisETA.N(sel);
        thisETA.D = thisETA.D(sel);
        thisETA.stats = thisETA.stats(sel);
        groupVar = NaN(length(sel), 1);
        groupVar(c.isLickUnresponsiveButDown & c.isPressUp) = 0;
        groupVar(c.isLickUnresponsiveButDown & c.isPressUnresponsiveButUp) = 0;
        groupVar(c.isLickDown & c.isPressUp) = 0;
        groupVar(c.isLickDown & c.isPressUnresponsiveButUp) = 0;
        groupVar(c.isLickUp & c.isPressUnresponsiveButDown) = 1;
        groupVar(c.isLickUp & c.isPressDown) = 1;
        groupVar(c.isLickUnresponsiveButUp & c.isPressUnresponsiveButDown) = 1;
        groupVar(c.isLickUnresponsiveButUp & c.isPressDown) = 1;

        groupVar(c.isLickUnresponsiveButDown & c.isPressUnresponsiveButDown) = 2;
        groupVar(c.isLickUnresponsiveButDown & c.isPressDown) = 2;
        groupVar(c.isLickDown & c.isPressUnresponsiveButDown) = 2;
        groupVar(c.isLickDown & c.isPressDown) = 2;
        groupVar(c.isLickUp & c.isPressUp) = 3;
        groupVar(c.isLickUp & c.isPressUnresponsiveButUp) = 3;
        groupVar(c.isLickUnresponsiveButUp & c.isPressUp) = 3;
        groupVar(c.isLickUnresponsiveButUp & c.isPressUnresponsiveButUp) = 3;
        groupVar = groupVar(sel);
        N = histcounts(groupVar, -0.5:2:3.5);
        yline(ax, cumsum(N(1:end-1)) + 1, '--');
        if ~exist('order', 'var')
            [~, order] = EphysUnit.plotETA(ax, thisETA, clim=CLIM, signWindow=[-0.3, 0], sortWindow=[-3, -0.1], sortThreshold=0.25, sortGroup=groupVar);
        else
            EphysUnit.plotETA(ax, thisETA, clim=CLIM, order=order);
        end
    
        ax.Box = 'off';
        ax.YAxis.Visible = 'off';
        colorbar(ax, 'off')
        xlim(ax, XLPERIMOVE{iCol})
    
        xlabel(ax, 'Time (s) or lick phase')
        xline(ax, 0, '-')
        ax.XAxis.TickLabelRotation = 0;
        ax.YAxis.Direction = 'reverse';
        colorbar(ax, 'off')
    
        % 1. Peri-cue
        ax = axes(tl); AX(iCol, 1) = ax;
        ax.Layout.Tile = 1;
        ax.Layout.TileSpan = [1, layout.w{iCol}(1)];
        
        thisETA = ETACUE{iCol};
        thisETA.X = thisETA.X(sel, :);
        thisETA.N = thisETA.N(sel);
        thisETA.D = thisETA.D(sel);
        EphysUnit.plotETA(ax, thisETA, clim=CLIM, order=order)
    
        ax.Box = 'off';
        xlim(ax, XLCUE{iCol})
        xticks(ax, [-0.8, 0, 1])
        xline(ax, 0, '-')
        xlabel(ax, 'start cue')
        ax.XAxis.TickLabelRotation = 0;
        ax.YAxis.Direction = 'reverse';
        yticks(ax, [1, 50:50:427, 427])
    
    
        ylim(AX(iCol, :), [1, nnz(sel)+1])
        linkaxes(AX(iCol, :), 'y')
        if iCol == 1
            colorbar(ax, 'off')
        else
            ax.Colorbar.Layout.Tile = 'east';
        end
    end
    
    title(AX, '')
    title(AX(1, 2), TITLE(1))
    title(AX(2, 2), TITLE(2))
end