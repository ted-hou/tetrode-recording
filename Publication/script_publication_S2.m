load_behavior_sessions

%% S2. Lick vs. Reach by animal
edges = 1:1:10;
centers = 0.5*(edges(2:end) + edges(1:end-1));

hasPress = nSessionsPress > 0;
hasLick = nSessionsLick > 0;

ncols = 4;
nrows = ceil(nnz(hasPress & hasLick)/ncols);
fig = figure(Units='inches', Position=[0, 0, 5, 3], DefaultAxesFontSize=p.fontSize, DefaultAxesFontName='Arial', Name='Reach task training progress');

hasPress = nSessionsPress > 0;
hasLick = nSessionsLick > 0;

[~, bestDay] = sort(-cellfun(@(t) nnz(t >= 3 & t <= 7) / nnz(t >= 1), pt) - cellfun(@(t) nnz(t >= 3 & t <= 7) / nnz(t >= 1), lt), 2);

i = 0;
for ia = find(hasPress(:)' & hasLick(:)')
    i = i + 1;
    ax = subplot(nrows, ncols, i);
    hold(ax, 'on')
    ptsel = pt(ia, bestDay(ia, 1:3));
    ptsel = cat(1, ptsel{:});
    ptsel(ptsel < edges(1)) = [];
    N1 = histcounts(ptsel, edges, Normalization='probability');
    plot(ax, centers, N1, 'r', LineWidth=2, DisplayName='Reach')

    ltsel = lt(ia, bestDay(ia, 1:3));
    ltsel = cat(1, ltsel{:});
    ltsel(ltsel < edges(1)) = [];
    N2 = histcounts(ltsel, edges, Normalization='probability');
    plot(ax, centers, N2, 'b', LineWidth=2, DisplayName='Lick')

    ylim(ax, [0, max([N1, N2]) + 0.01])
    
%     ndays = nSessionsPress(ia);
%     for id = 1:ndays
%         N = histcounts(pt{ia, id}, edges, Normalization='probability');
%         plot(ax, centers, N, Color=hsl2rgb([0.7*(id-1)/(ndays-1), 0.1, 0.5]), LineWidth=0.1, DisplayName=sprintf('Day %g', daysPress(id)))
%     end
    hold(ax, 'off')
    set(ax, FontSize=p.fontSize, FontName='Arial');
    if contains(lower(animalNames{ia}), 'daisy')
        title(ax, sprintf('%d (F)', ia))
    else
        title(ax, sprintf('%d (M)', ia))
    end

end
legend(ax, Position=[0.353156636288793,0.946397568363208,0.308617231363285,0.046874998913457], Orientation='horizontal')

% for ia = find(hasLick(:)')
%     ax = subplot(nrows, ncols, ia);
%     hold(ax, 'on')
%     ltsel = lt(ia, bestDayLick(ia, 1:5));
%     ltsel = cat(1, ltsel{:});
%     N = histcounts(ltsel, edges, Normalization='probability');
%     plot(ax, centers, N, 'b', LineWidth=2)
% end

annotation(fig, 'textbox', [0.062393427881404,0.028034722222222,0.899999999999997,0.05], String='Time to contact (s)', ...
    HorizontalAlignment='center', LineStyle='none', FontSize=11);
annotation(fig, 'textbox', [0.06151279482041,0.202083333333332,0.449999999999998,0.05], String='Probability', ...
    HorizontalAlignment='center', LineStyle='none', FontSize=11, Rotation=90);


print(fig, 'Fig S2 lick time histogram per animal best 3 sessions.fig');


clear ax fig id