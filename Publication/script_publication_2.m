%%% Figure 2c and Figure S2
p.fontSize = 8;
%% Only run once (make bs.mat and save to C:\SERVER_PRIVATE)
% make_behavior_sessions

%% Load bs.mat, do some processing
load_behavior_sessions

%% Fig 2c. Movement time histogram (self-timed reach)

close all
clear fig ax
edges = 0:0.5:10;

% Plot aggregate histograms as line plots
fig = figure(Units='inches', Position=[0, 0, 7.5, 2.5], DefaultAxesFontSize=p.fontSize, DefaultAxesFontName='Arial');
ax = axes(fig);
centers = 0.5*(edges(2:end) + edges(1:end-1));
hold(ax, 'on')
ndayshown = 18;
for id = 1:ndayshown
    N = histcounts(ptcat{id}, edges, Normalization='probability');
    % plot(ax, centers, N, Color=hsl2rgb([0.7*(id-1)/(ndayshown-1), 1, 0.5]), LineWidth=1, DisplayName=sprintf('Day %g (%g trials)', daysPress(id), nTrialsPress(id)))
    plot(ax, centers, N, Color=hsl2rgb([0.7*(id-1)/(ndayshown-1), 1, 0.4]), LineWidth=1, DisplayName=sprintf('Day %g', daysPress(id)))
    xlabel(ax, 'Contact time (s)')
    ylabel(ax, 'Probability')
end
hold(ax, 'off')
legend(ax, Location='northeast', NumColumns=2)
title(ax, sprintf('Reach task performance (%g animals)', nAnimalsPress), FontSize=p.fontSize)
fontsize(fig, p.fontSize, 'points')
fontname(fig, 'Arial')
print(fig, 'Fig 2c aggregate reach time histogram across training sessions.fig');
clear ax fig id

