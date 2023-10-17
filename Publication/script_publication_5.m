p.fontSize=8;

%% Load bs.mat, do some processing
load_behavior_sessions

%% 5b. Lick time histogram across training session
clear fig ax
edges = 0:0.5:10;
fig = figure(Units='inches', Position=[0, 0, 5, 2.5], DefaultAxesFontSize=p.fontSize, DefaultAxesFontName='Arial', Name='Reach task training progress');
ax = axes(fig);
centers = 0.5*(edges(2:end) + edges(1:end-1));
hold(ax, 'on')
ndayshown = 18;
for id = 1:ndayshown
    N = histcounts(ltcat{id}, edges, Normalization='probability');
    % plot(ax, centers, N, Color=hsl2rgb([0.7*(id-1)/(ndayshown-1), 1, 0.5]), LineWidth=1, DisplayName=sprintf('Day %g (%g trials)', daysPress(id), nTrialsLick(id)))
    plot(ax, centers, N, Color=hsl2rgb([0.7*(id-1)/(ndayshown-1), 1, 0.5]), LineWidth=1, DisplayName=sprintf('Day %g', daysPress(id)))
    xlabel(ax, 'Contact time (s)')
    ylabel(ax, 'Probability')
end
hold(ax, 'off')
legend(ax, Location='northeast', NumColumns=2)
title(ax, sprintf('Lick task performance (%g animals)', nAnimalsLick))
set(ax, FontSize=p.fontSize);

print(fig, 'Fig 5c aggregate lick time histogram across training sessions.fig');

clear ax fig id

%% 5c. Timing comparison (average lick time vs. average press time)
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

fig = figure(Units='inches', Position=[0, 0, 5, 1.75]);
ax = gobjects(2, 1);

ax(1) = subplot(1, 2, 1);
hold(ax(1), 'on')
plot(ax(1), centers, histPress, 'r', LineWidth=1.5, DisplayName=sprintf('Reach', length(aggrPressTime)))
plot(ax(1), centers, histLick, 'b', LineWidth=1.5, DisplayName=sprintf('Lick', length(aggrLickTime)))
hold(ax(1), 'off')
legend(ax(1))
xlabel(ax(1), 'Time from cue (s)')
ylabel(ax(1), 'Probability')
xticks(ax(1), [1, 4, 7, 10])
xlim(ax(1), [0, 10])

ax(2) = subplot(1, 2, 2);
hold(ax(2), 'on')
h = scatter(ax(2), medianPressTime, medianLickTime, 5, 'k', 'filled', DisplayName=sprintf('%i sessions', length(euIndices)));
plot(ax(2), [0, 10], [0, 10], 'k:')
xlim(ax(2), [0, 10])
ylim(ax(2), [0, 10])
legend(ax(2), h)
xlabel(ax(2), 'Median reach time (s)')
ylabel(ax(2), 'Median lick time (s)')
hold(ax(2), 'off')
xticks(ax(2), [1, 4, 7, 10])
yticks(ax(2), [1, 4, 7, 10])

fontsize(ax, p.fontSize, 'points');
fontname(ax, 'Arial');
%% Read DLC data (SLOW!)
read_DLC_data;

%% 5c,5d. Lick vs. Reach DLC traces
close all
resultNames = {'incorrect', 'correct'};

for iExp = 1:length(fstats)
    fig = figure(Units='inches', Position=[0+(iExp-1)*4.5, -1, 2.5, 3], DefaultAxesFontSize=p.fontSize);
    axAll = gobjects(1, 2);
    iTrialType = 0;
    for trialTypeName = {'press', 'lick'}
        trialTypeName = trialTypeName{1};
        iTrialType = iTrialType + 1;
        ax = subplot(2, 1, iTrialType);
        axAll(iTrialType) = ax;
        hold(ax, 'on')

        switch trialTypeName
            case {'press', 'lick'}
                mu = table2array(fstats{iExp}.(trialTypeName).mean(:, fnames));
                sd = table2array(fstats{iExp}.(trialTypeName).sd(:, fnames));
                n = fstats{iExp}.(trialTypeName).nTrials;
        end

        sdNames = {'handContra', 'handIpsi', 'tongue'};
        ftNames = {'handContra', 'handIpsi', 'tongue'};


        h = gobjects(1, length(ftNames));
        ii = 0;
        for iVar = 1:length(fnames)
            col = hsl2rgb([0.8*(iVar-1)/(length(fnames)-1), 1, 0.5]);
            if ismember(fnames{iVar}, ftNames)
                ii = ii + 1;
                h(ii) = plot(ax, t, mu(:, iVar), Color=col, LineWidth=1.5, DisplayName=fnamesDisp{iVar});
            end
            if ismember(fnames{iVar}, sdNames)
                sel = ~isnan(mu(:, iVar)+sd(:, iVar));
                patch(ax, [t(sel)'; flip(t(sel)')], [mu(sel, iVar)-sd(sel, iVar); flip(mu(sel, iVar)+sd(sel, iVar))], 'r', ...
                    LineStyle='none', FaceAlpha=0.075, FaceColor=col)
            end
        end
        switch trialTypeName
            case 'press'
                xlabel(ax, 'time to bar contact (s)')
                trialTypeDispName = 'Reach';
            case 'lick'
                xlabel(ax, 'time to spout contact (s)')
                trialTypeDispName = 'Lick';
        end
%         if iExp == 1
            ylabel(ax, 'z-scored speed (a.u.)')
%         end
        plot(ax, [0, 0], [-100, 100], 'k--')
%         if iExp == 1 && iTrialType == 1
        if iTrialType == 1
            legend(ax, h, Location='northwest')
        end
        title(ax, sprintf('%s trials (%s, N=%d)', trialTypeDispName, resultNames{iExp}, n));
        % title(ax, sprintf('%s trials (%s)', trialTypeDispName, resultNames{iExp}));
        hold(ax, 'off')
        fontsize(ax, p.fontSize, 'points');
        fontname(ax, 'Arial');
    end

    set(axAll(1:2), YLim=[0, 12])
    set(axAll, XLim=[-3, 1])

    fontsize(fig, p.fontSize, 'points');
    fontname(fig, 'Arial');
end
clear iTrialType iExp ax fig trialTypeName h iVar axAll n
