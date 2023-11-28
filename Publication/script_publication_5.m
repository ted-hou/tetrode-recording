p.fontSize=9;
% p.width = 7;
% p.height = 5.5;
% p.leftPanelWidth = 3.5;
% p.rightPanelWidth = p.width - p.leftPanelWidth;

clear layout
layout.w = 7;
layout.h = 5;
layout.left.w = 7; % 3.5
layout.right.w = 7; % 3.5
layout.left.top.h = 4;
layout.left.bottom.h = 5;
layout.right.top.h = 2;
layout.right.bottom.h = 4;

p.lineWidth = 1.5;


%% Load bs.mat, do some processing
load_behavior_sessions
%% Read DLC data (SLOW!)
read_DLC_data;


%% Figure 5 
close all
fig = figure(Units='inches', Position=[0, 0, layout.w, layout.h]);
layout.tl = tiledlayout(fig, 1, layout.left.w + layout.right.w, TileSpacing='compact', Padding='compact');
layout.left.tl = tiledlayout(layout.tl, layout.left.top.h + layout.left.bottom.h, 1, TileSpacing='loose');
l = layout.left.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [1, layout.left.w];

layout.right.tl = tiledlayout(layout.tl, layout.right.top.h + layout.right.bottom.h, 1, TileSpacing='loose');
l = layout.right.tl; l.Layout.Tile = layout.left.w + 1; l.Layout.TileSpan = [1, layout.right.w];

layout.right.top.tl = tiledlayout(layout.right.tl, 1, 2, TileSpacing='compact');
l = layout.right.top.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [layout.right.top.h, 1];

layout.right.bottom.tl = tiledlayout(layout.right.tl, 2, 2, TileSpacing='compact');
l = layout.right.bottom.tl; l.Layout.Tile = layout.right.top.h + 1; l.Layout.TileSpan = [layout.right.bottom.h, 1];


% 5b. Lick time histogram across training session
edges = 0:0.5:10;

% Plot aggregate histograms as line plots
ax = nexttile(layout.left.tl, layout.left.top.h + 1, [layout.left.bottom.h, 1]);
centers = 0.5*(edges(2:end) + edges(1:end-1));
hold(ax, 'on')
ndayshown = 18;
for id = 1:ndayshown
    N = histcounts(ltcat{id}, edges, Normalization='probability');
    % plot(ax, centers, N, Color=hsl2rgb([0.7*(id-1)/(ndayshown-1), 1, 0.5]), LineWidth=1, DisplayName=sprintf('Day %g (%g trials)', daysPress(id), nTrialsPress(id)))
    plot(ax, centers, N, Color=hsl2rgb([0.7*(id-1)/(ndayshown-1), 1, 0.4]), LineWidth=1.5, DisplayName=sprintf('Day %g', daysPress(id)))
    xlabel(ax, 'Spout contact time (s)')
    ylabel(ax, 'Probability')
end
hold(ax, 'off')


% Legends (we're gonna be in this city)
assert(ndayshown==18)
ticks = [1, 18];
cmap = arrayfun(@(id) hsl2rgb([0.7*(id-1)/(ndayshown-1), 1, 0.4]), 1:ndayshown, 'UniformOutput', false);
cmap = cat(1, cmap{:});
colormap(ax, cmap)
h = colorbar(ax);
h.Ticks = (1:ndayshown) / ndayshown;
ticklabels = arrayfun(@(x) sprintf('%i', x), 1:ndayshown, UniformOutput=false);
for i = 1:ndayshown
    if ~ismember(i, ticks)
        ticklabels{i} = '';
    end
end
h.TickLabels = ticklabels;
h.Label.String = 'session';
h.Label.Position = [1.103333312471708,0.505319625773328,0];

% title(ax, sprintf('Lick task performance (%g animals)', nAnimalsLick), FontSize=p.fontSize)
xticks(ax, [0, 4, 10])
yticks(ax, ax.YLim(2))
ax.YLabel.Position = [-0.5,0.125,-1];

fontsize(ax, p.fontSize, 'points')
fontname(ax, 'Arial')
clear ax id ticks cmap h ticklabels i

% 5c. Timing comparison (average lick time vs. average press time)
euIndices = find(c.hasPress & c.hasLick);
[~, ia, ~] = unique({eu(euIndices).ExpName});
euIndices = euIndices(ia);

medianPressTime = arrayfun(@(iEu) median([eu(iEu).Trials.Press.Stop] - [eu(iEu).Trials.Press.Start]), euIndices);
medianLickTime = arrayfun(@(iEu) median([eu(iEu).Trials.Lick.Stop] - [eu(iEu).Trials.Lick.Start]), euIndices);

aggrPressTime = arrayfun(@(iEu) [eu(iEu).Trials.Press.Stop] - [eu(iEu).Trials.Press.Start], euIndices, UniformOutput=false);
aggrPressTime = cat(2, aggrPressTime{:});

aggrLickTime = arrayfun(@(iEu) [eu(iEu).Trials.Lick.Stop] - [eu(iEu).Trials.Lick.Start], euIndices, UniformOutput=false);
aggrLickTime = cat(2, aggrLickTime{:});
edges = 1:0.5:10;
centers = (edges(1:end-1) + edges(2:end))*0.5;
histPress = histcounts(aggrPressTime, edges, Normalization='probability');
histLick = histcounts(aggrLickTime, edges, Normalization='probability');

ax = nexttile(layout.right.top.tl);
hold(ax, 'on')
plot(ax, centers, histPress, 'r', LineWidth=1.5, DisplayName=sprintf('Reach', length(aggrPressTime)))
plot(ax, centers, histLick, 'b', LineWidth=1.5, DisplayName=sprintf('Lick', length(aggrLickTime)))
hold(ax, 'off')
legend(ax, Location='northoutside', Orientation='horizontal')
xlabel(ax, 'Contact time (s)')
ylabel(ax, 'Probability')
xticks(ax, [edges(1), 4, 10])
xlim(ax, [0, 10])
% yticks(ax, ax.YLim(2))
yticks(ax, [0, 0.1])
% ax.YLabel.Position = [-0.5,0.125,-1];
fontsize(ax, p.fontSize, 'points');
fontname(ax, 'Arial');

ax = nexttile(layout.right.top.tl);
hold(ax, 'on')
h = scatter(ax, medianPressTime, medianLickTime, 5, 'k', 'filled', DisplayName=sprintf('%i sessions', length(euIndices)));
plot(ax, [0, 10], [0, 10], 'k:')
xlim(ax, [0, 10])
ylim(ax, [0, 10])
% legend(ax, h, Location='northoutside', Orientation='horizontal')
xlabel(ax, 'Median reach time (s)')
ylabel(ax, 'Median lick time (s)')
hold(ax, 'off')
xticks(ax, [0, 4, 10])
yticks(ax, [0, 4, 10])
fontsize(ax, p.fontSize, 'points');
fontname(ax, 'Arial');

% 5e,5f. Lick vs. Reach DLC traces
resultNames = {'incorrect', 'correct'};
trialTypes = {'press', 'lick'};

%     fig = figure(Units='inches', Position=[0+(iExp-1)*4.5, -1, 2.5, 3], DefaultAxesFontSize=p.fontSize);
for iTrialType = 1:length(trialTypes)
    for iExp = 1:length(fstats)
        trialTypeName = trialTypes{iTrialType};
        ax = nexttile(layout.right.bottom.tl);
        hold(ax, 'on')

        switch trialTypeName
            case {'press', 'lick'}
                mu = table2array(fstats{iExp}.(trialTypeName).mean(:, fnames));
                sd = table2array(fstats{iExp}.(trialTypeName).sd(:, fnames));
                n = fstats{iExp}.(trialTypeName).nTrials;
        end

        ftNames = {'handContra', 'tongue'};
        sdNames = {'handContra', 'tongue'};
        axisSide = {'left', 'right'};

        h = gobjects(1, length(ftNames));
        colororder(ax, getColor(1:length(ftNames), length(ftNames), 0.7));
        for iFt = 1:length(ftNames)
            yyaxis(ax, axisSide{iFt});
            iVar = find(strcmpi(ftNames{iFt}, fnames));
            col = getColor(iFt, length(ftNames), 0.7);
            h(iFt) = plot(ax, t, mu(:, iVar), Color=col, LineWidth=1.5, DisplayName=fnamesDisp{iVar});
            if ismember(ftNames{iFt}, sdNames)
                sel = ~isnan(mu(:, iVar)+sd(:, iVar));
                patch(ax, [t(sel)'; flip(t(sel)')], [mu(sel, iVar)-sd(sel, iVar); flip(mu(sel, iVar)+sd(sel, iVar))], 'r', ...
                    LineStyle='none', FaceAlpha=0.075, FaceColor=col)
            end
        end

        switch trialTypeName
            case 'press'
                trialTypeDispName = 'Reach';
            case 'lick'
                trialTypeDispName = 'Lick';
        end
        plot(ax, [0, 0], [-100, 100], 'k--')
%         title(ax, sprintf('%s trials (%s, N=%d)', trialTypeDispName, resultNames{iExp}, n));
%         title(ax, sprintf('%s trials (%s)', trialTypeDispName, resultNames{iExp}));
        title(ax, sprintf('%s (%s)', trialTypeDispName, resultNames{iExp}));
        hold(ax, 'off')
        xlim(ax, [-2, 2])

        yyaxis(ax, 'left')
        if iExp == 1 && strcmp(trialTypeName, 'press')
            ylabel(ax, 'Contra hand velocity (a.u.)', Units='normalized', Position=[-0.20,-0,0])
        end
        if iExp == 2
            yticks(ax, [])
        end
        ylim(ax, [-5, 5])

        yyaxis(ax, 'right')
        if iExp == 2 && strcmp(trialTypeName, 'press')
            ylabel(ax, 'Lick probability', Units='normalized', VerticalAlignment='bottom', Position=[1.23,-0,0])
        end
        if iExp == 1
            yticks(ax, [])
        else
            yticks(ax, [0, 1])
        end
        ylim(ax, [-1.5, 1.5])
        xlabel(layout.right.bottom.tl, 'Time to bar/spout contact (s)', FontSize=p.fontSize, FontName='Arial')
        if strcmp(trialTypeName, 'press')
            xticks(ax, [])
        end

        fontsize(ax, p.fontSize, 'points');
        fontname(ax, 'Arial');
    end
end


copygraphics(fig, ContentType='vector')

clear iExp ax fig trialTypeName h iVar axAll n
