load_behavior_sessions

%% Fig S1. Movement time histograms for each animal (self-timed reach)
close all
animalNames = cellfun(@(bs) bs(1).animalName, bs, UniformOutput=false);


edges = 1:1:10;
centers = 0.5*(edges(2:end) + edges(1:end-1));

ncols = 5;
nrows = ceil(nAnimalsPress/ncols);
fig = figure(Units='inches', Position=[0, 0, 6.5, 6], DefaultAxesFontSize=p.fontSize, DefaultAxesFontName='Arial', Name='Reach task training progress');
tl = tiledlayout(fig, nrows, ncols, TileSpacing='compact');

hasPress = nSessionsPress > 0;
hasLick = nSessionsLick > 0;

[~, bestDayPress] = sort(-cellfun(@(t) nnz(t >= 3 & t <= 7) / nnz(t >= 1), pt), 2);
[~, bestDayLick] = sort(-cellfun(@(t) nnz(t >= 3 & t <= 7) / nnz(t >= 1), lt), 2);


for ia = find(hasPress(:)')
    ax = nexttile(tl);
    hold(ax, 'on')
    ptsel = pt(ia, bestDayPress(ia, 1:3));
    ptsel = cat(1, ptsel{:});
    ptsel = ptsel(ptsel >= edges(1));
    N = histcounts(ptsel, edges, Normalization='probability');
    plot(ax, centers, N, 'r', LineWidth=2)
    ylim(ax, [0, max(N) + 0.01])
    
%     ndays = nSessionsPress(ia);
%     for id = 1:ndays
%         N = histcounts(pt{ia, id}, edges, Normalization='probability');
%         plot(ax, centers, N, Color=hsl2rgb([0.7*(id-1)/(ndays-1), 0.1, 0.5]), LineWidth=0.1, DisplayName=sprintf('Day %g', daysPress(id)))
%     end
    hold(ax, 'off')
    title(ax, ai(ia).displayName);
    fontsize(ax, p.fontSize, 'points');
end

% for ia = find(hasLick(:)')
%     ax = subplot(nrows, ncols, ia);
%     hold(ax, 'on')
%     ltsel = lt(ia, bestDayLick(ia, 1:5));
%     ltsel = cat(1, ltsel{:});
%     N = histcounts(ltsel, edges, Normalization='probability');
%     plot(ax, centers, N, 'b', LineWidth=2)
% end
% 
% annotation(fig, 'textbox', [0.064397435897436,0.030638888888889,0.9,0.05], String='Time to contact (s)', ...
%     HorizontalAlignment='center', LineStyle='none', FontSize=11);
% annotation(fig, 'textbox', [0.093576923076923,0.264583333333332,0.45,0.05], String='Probability', ...
%     HorizontalAlignment='center', LineStyle='none', FontSize=11, Rotation=90);

xlabel(tl, 'Bar-contact time relative to cue (s)', fontSize=p.fontSize)
ylabel(tl, 'Probability', fontSize=p.fontSize)

copygraphics(fig, ContentType='vector')

clear ax fig id

