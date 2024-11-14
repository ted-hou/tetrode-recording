clear etaSpontaneous
euSpontaneousCombined = ([euSpontaneous, euReachDir2Tgt]);
etaSpontaneousCombined.press = euSpontaneousCombined.getETA('count', 'press_spontaneous', [-4, 2], minTrialDuration=4, normalize=[-4, -2], resolution=0.025);
etaSpontaneousCombined.pressCorrect = euSpontaneousCombined.getETA('count', 'press_spontaneous_correct', [-4, 2], minTrialDuration=4, normalize=[-4, -2], resolution=0.025);
etaSpontaneousCombined.pressIncorrect = euSpontaneousCombined.getETA('count', 'press_spontaneous_incorrect', [-4, 2], minTrialDuration=4, normalize=[-4, -2], resolution=0.025);

metaSpontaneousCombined.press = mean(etaSpontaneousCombined.press.X(:, etaSpontaneousCombined.press.t <= 0 & etaSpontaneousCombined.press.t >= -0.3), 2, 'omitnan');

%% Plot correct vs. incorrect
close all
p.fontSize = 9;
ETA = {etaSpontaneousCombined.press, etaSpontaneousCombined.pressCorrect, etaSpontaneousCombined.pressIncorrect};
TITLE = ["All", "Rewarded", "Incorrect"];
COLORS = {[0.8, 0.2, 0.2], [0.2, 0.8, 0.2], [0.4, 0.4, 0.4]};

fig = figure(Units='inches', Position=[1 1 8 4]);
tl = tiledlayout(fig, 2, 2, TileSpacing='compact', Padding='compact', TileIndexing='columnmajor');

ax = nexttile(tl);
hold(ax, 'on')
h = gobjects(1, 3);
for iLine = 2:3
    h(iLine) = plot(ax, ETA{iLine}.t, mean(ETA{iLine}.X(metaSpontaneousCombined.press < 0, :), 1, 'omitnan'), LineWidth=1.5, DisplayName=TITLE(iLine), Color=COLORS{iLine});
    plot(ax, ETA{iLine}.t, mean(ETA{iLine}.X(metaSpontaneousCombined.press > 0, :), 1, 'omitnan'), LineWidth=1.5, DisplayName=TITLE(iLine), Color=COLORS{iLine})
end
yline(ax, 0, 'k--')
xline(ax, [0, 0.8], 'k--')
legend(h(2:3), Location='northwest')
% xlabel(ax, 'Time to bar contact (s)')
ylabel(ax, 'Norm spike rate (a.u.)')
fontsize(ax, p.fontSize, 'points')

for iRow = 1:3
    ax = nexttile(tl);
    if iRow == 1
        [~, order] = EphysUnit.plotETA(ax, ETA{iRow}, clim=[-1.5, 1.5], signWindow=[-0.3, 0], sortWindow=[-3, 0], sortThreshold=0.25);
    else
        EphysUnit.plotETA(ax, ETA{iRow}, clim=[-1.5, 1.5], order=order);
    end
    xline(ax, [0, 0.8], '--')
    xlim(ax, [-4, 2])
    title(ax, TITLE(iRow))
    xlabel(ax, '')
    ylabel(ax, 'Unit')
    yticks(ax, [1, 100:100:length(ETA{iRow}.N), length(ETA{iRow}.N)])
    fontsize(ax, p.fontSize, 'points')
end
xlabel(tl, 'Time to bar contact (s)', FontSize=p.fontSize)

copygraphics(fig, ContentType='vector')

clear ETA TITLE fig tl iRow order iLine ax

%% Plot correct vs. incorrect
close all
p.fontSize = 9;
ETA = {etaFine.press, etaFine.correctPress, etaFine.incorrectPress};
TITLE = ["All", "Rewarded", "Incorrect"];
COLORS = {[0.8, 0.2, 0.2], [0.2, 0.8, 0.2], [0.4, 0.4, 0.4]};

fig = figure(Units='inches', Position=[1 1 8 4]);
tl = tiledlayout(fig, 2, 2, TileSpacing='compact', Padding='compact', TileIndexing='columnmajor');

ax = nexttile(tl);
hold(ax, 'on')
h = gobjects(1, 3);
for iLine = 2:3
    h(iLine) = plot(ax, ETA{iLine}.t, mean(ETA{iLine}.X(meta.press < 0 & c.hasPress, :), 1, 'omitnan'), LineWidth=1.5, DisplayName=TITLE(iLine), Color=COLORS{iLine});
    plot(ax, ETA{iLine}.t, mean(ETA{iLine}.X(meta.press > 0 & c.hasPress, :), 1, 'omitnan'), LineWidth=1.5, DisplayName=TITLE(iLine), Color=COLORS{iLine})
end
yline(ax, 0, 'k--')
xline(ax, [0, 0.8], 'k--')
legend(h(2:3), Location='northwest')
% xlabel(ax, 'Time to bar contact (s)')
ylabel(ax, 'Norm spike rate (a.u.)')
fontsize(ax, p.fontSize, 'points')

for iRow = 1:3
    ax = nexttile(tl);
    if iRow == 1
        [~, order] = EphysUnit.plotETA(ax, ETA{iRow}, c.hasPress, clim=[-1.5, 1.5], signWindow=[-0.3, 0], sortWindow=[-3, 0], sortThreshold=0.25);
    else
        EphysUnit.plotETA(ax, ETA{iRow}, c.hasPress, clim=[-1.5, 1.5], order=order);
    end
    xline(ax, [0, 0.8], '--')
    xlim(ax, [-4, 2])
    title(ax, TITLE(iRow))
    xlabel(ax, '')
    ylabel(ax, 'Unit')
    yticks(ax, [1, 200:200:nnz(c.hasPress), nnz(c.hasPress)])
%     ylim(ax, [1, length(ETA{iRow}.N) + 1])
    fontsize(ax, p.fontSize, 'points')
end
xlabel(tl, 'Time to bar contact (s)', FontSize=p.fontSize)

copygraphics(fig, ContentType='vector')

clear ETA TITLE fig tl iRow order iLine ax

%% Plot correct vs. incorrect (lick x press)
close all
p.fontSize = 9;
ETA = {etaFine.press, etaFine.correctPress, etaFine.incorrectPress};
TITLE = ["All", "Rewarded", "Incorrect"];
COLORS = { ...
    [], hsl2rgb([0.7, 1, 0.5]), [hsl2rgb([0.7, 0.5, 0.5]), 0.5]; ...
    [], hsl2rgb([0.0, 1, 0.5]), [hsl2rgb([0.0, 0.5, 0.5]), 0.5]; ...
    [], hsl2rgb([0.7, 1, 0.5]), [hsl2rgb([0.7, 0.5, 0.5]), 0.5]; ...
    [], hsl2rgb([0.0, 1, 0.5]), [hsl2rgb([0.0, 0.5, 0.5]), 0.5]; ...
    };

fig = figure(Units='inches', Position=[1 1 6 4]);
tl = tiledlayout(fig, 1, 1, TileSpacing='compact', Padding='compact', TileIndexing='columnmajor');

ax = nexttile(tl);
hold(ax, 'on')
h = gobjects(4, 3);
sel = c.hasPress & c.hasLick;
for iLine = 2:3
    h(1, iLine) = plot(ax, ETA{iLine}.t, mean(ETA{iLine}.X(c.isPressDown & c.isLickDown & sel, :), 1, 'omitnan'), LineWidth=1.5, DisplayName=sprintf('%s (lick dec)', TITLE(iLine)), Color=COLORS{1, iLine});
    h(2, iLine) = plot(ax, ETA{iLine}.t, mean(ETA{iLine}.X(c.isPressDown & c.isLickUp & sel, :), 1, 'omitnan'), LineWidth=1.5, DisplayName=sprintf('%s (lick inc)', TITLE(iLine)), Color=COLORS{2, iLine});
    h(3, iLine) = plot(ax, ETA{iLine}.t, mean(ETA{iLine}.X(c.isPressUp & c.isLickDown & sel, :), 1, 'omitnan'), LineWidth=1.5, DisplayName=TITLE(iLine), Color=COLORS{3, iLine});
    h(4, iLine) = plot(ax, ETA{iLine}.t, mean(ETA{iLine}.X(c.isPressUp & c.isLickUp & sel, :), 1, 'omitnan'), LineWidth=1.5, DisplayName=TITLE(iLine), Color=COLORS{4, iLine});
end
yline(ax, 0, 'k--')
xline(ax, [0, 0.8], 'k--')
h = h(1:2, 2:3);
h = h(:);
legend(h, Location='northwest')
% xlabel(ax, 'Time to bar contact (s)')
ylabel(ax, 'Norm spike rate (a.u.)')
fontsize(ax, p.fontSize, 'points')

xlabel(tl, 'Time to bar contact (s)', FontSize=p.fontSize)

copygraphics(fig, ContentType='vector')

clear ETA TITLE fig tl iRow order iLine ax sel
