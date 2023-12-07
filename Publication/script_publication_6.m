p.fontSize = 9;
p.width = 5;
p.height = 2;

%% Load all units
load_ephysunits;
boot_response_dir;

%% 6d. Salt and pepper map (lick vs. reach vs. osci lick)
close all
MOVETYPE = {'press', 'lick', 'lickosci'};
SEL = {c.isPressResponsive; c.isLickResponsive; c.isLick};
STATS = {meta.press, meta.lick, pi - abs(phase(freq==8, :))};%meta.anyLickNorm};
TITLE = {'Pre-reach', 'Pre-lick', 'Lick-osci'};
COLOR = {[], [], hsl2rgb([50/360, 1, 0.4])};
ALPHA = {0.125, 0.125, 0.5};
POS = { ...
        [0.05,0.18,0.8/3,0.6], ...
        [0.10+0.8/3,0.18,0.8/3,0.6], ...
        [0.15+0.8/3*2,0.18,0.8/3,0.6], ...
    };

fig = figure(Units='inches', Position=[0, 0, p.width, p.height], DefaultAxesFontSize=p.fontSize);
for iMove = 1:length(MOVETYPE)
%     ax = subplot(1, length(MOVETYPE), iMove);
    ax = axes(fig, Position=POS{iMove});
    sel = SEL{iMove};
    coords = euPos(sel, :);
    stats = STATS{iMove}(sel);
    AcuteRecording.plotMap(ax, coords, stats, [0, pi], 0, UseSignedML=false, BubbleSize=[1, 5], MarkerAlpha=ALPHA{iMove}, Color=COLOR{iMove});
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
%     fontsize(fig, p.fontSize, 'points');
%     fontname(fig, 'Arial')
    fontsize(ax, p.fontSize, 'points');
    fontname(ax, 'Arial')
end

% Make a dummy axes for the legends (red, blue, yellow circles)
ax = axes(fig, Position=[0.5, 0.8, 0, 0]); hold(ax, 'on')
ax.Visible = 'off';
scatter(ax, 0, 0, 1, [1, 0, 0], 'filled', DisplayName='excited')
scatter(ax, 0, 0, 1, [0, 0, 1], 'filled', DisplayName='suppressed')
scatter(ax, 0, 0, 1, COLOR{3}, 'filled', DisplayName='oscillatory')
h = legend(ax, Position=[0.181886578916949,0.868923612275264,0.643749990314244,0.098958331005027], Orientation='horizontal');
fontsize(h, p.fontSize, 'points');
fontname(h, 'Arial')

copygraphics(fig, ContentType='vector')