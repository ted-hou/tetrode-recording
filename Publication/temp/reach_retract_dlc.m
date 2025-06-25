%%
load_ephysunits;
read_DLC_data;
%%
expIndices = zeros(size(euAcute));
for iExp = 1:length(expAcute)
    expIndices(ismember(euAcute, expAcute(iExp).eu)) = iExp;
end

%%

clear rd
rd.press = eu.getRasterData('press', [-4, 4], minTrialDuration=2, maxTrialDuration=Inf);
rd.lick = eu.getRasterData('lick', [-4, 4], minTrialDuration=2, maxTrialDuration=Inf);

eta.correctPress = eu.getETA('count', 'press', [-4, 4], minTrialDuration=4, normalize=p.etaNorm, resolution=0.1);
eta.incorrectPress = eu.getETA('count', 'press', [-4, 4], minTrialDuration=p.minTrialDuration, maxTrialDuration=4, normalize=p.etaNorm, resolution=0.1);

eta.correctLick = eu.getETA('count', 'lick', [-4, 4], minTrialDuration=4, normalize=p.etaNorm, resolution=0.1);
eta.incorrectLick= eu.getETA('count', 'lick', [-4, 4], minTrialDuration=p.minTrialDuration, maxTrialDuration=4, normalize=p.etaNorm, resolution=0.1);

% %% Plot reach and lick raster for reach-dec (among 427) cells.
% sel = c.hasPress & c.hasLick & c.isPressDown;
% tl = tiledlayout(figure, 2, 1);
% ax = gobjects(2, 1);
% ax(1) = nexttile(tl);
% ax(2) = nexttile(tl);
% for i = find(sel)
%     cla(ax(1));
%     cla(ax(2));
%     EphysUnit.plotRaster(ax(1), rd.press(i), xlim=[-4, 4]);
%     xline(ax(1), 0, 'r-', LineWidth=1.5)
%     yline(ax(1), find(rd.press(i).duration>4, 1), 'g-', LineWidth=1.5)
% 
%     EphysUnit.plotRaster(ax(2), rd.lick(i), xlim=[-4, 4]);
%     xline(ax(2), 0, 'r-', LineWidth=1.5)
%     yline(ax(2), find(rd.lick(i).duration>4, 1), 'g-', LineWidth=1.5)
% 
%     title(ax(1), 'Reach')
%     title(ax(2), 'Lick')
%     title(tl, rd.press(i).name, Interpreter='none')
% 
%     waitforbuttonpress();
% end
% 
% %% Plot reach raster/peth for reach-dec (among 784) cells.
% sel = c.hasPress & c.isPressDown;
% tl = tiledlayout(figure(Units='inches', OuterPosition=[1, 1, 8, 12]), 2, 1);
% ax = gobjects(2, 1);
% ax(1) = nexttile(tl);
% ax(2) = nexttile(tl);
% for i = find(sel)
%     cla(ax(1));
%     cla(ax(2));
%     EphysUnit.plotRaster(ax(1), rd.press(i), xlim=[-4, 4]);
%     set(ax(1).Legend, AutoUpdate=false)
%     xline(ax(1), 0, 'r-', LineWidth=1.5)
%     yline(ax(1), find(rd.press(i).duration>4, 1), 'g-', LineWidth=1.5)
% 
%     hold(ax(2), 'on')
%     h = gobjects(2, 1);
%     h(1) = plot(ax(2), eta.incorrectPress.t, eta.incorrectPress.X(i, :), Color=hsl2rgb([0, 0.6, 0.6]), DisplayName='Incorrect Reach', LineWidth=1.5);
%     h(2) = plot(ax(2), eta.correctPress.t, eta.correctPress.X(i, :), Color=hsl2rgb([0.3, 0.4, 0.4]), DisplayName='Correct Reach', LineWidth=1.5);
%     legend(h, AutoUpdate=false)
%     xline(ax(2), 0, 'k--')
%     yline(ax(2), 0, 'k--')
% 
%     xlim(ax(2), [-4, 4])
% 
%     title(ax(1), 'Reach Raster')
%     title(ax(2), 'Reach PETH')
%     title(tl, rd.press(i).name, Interpreter='none')
% 
%     waitforbuttonpress();
% end

%% Plot reach raster, peth, DLC trajectories for reach-dec DLC cells (among 163)
% close all
if ~exist('C:\SERVER\Figures\reach_retract_dlc\reach_decrease', 'dir')
    mkdir('C:\SERVER\Figures\reach_retract_dlc\reach_decrease')
end
if ~exist('C:\SERVER\Figures\reach_retract_dlc\reach_increase', 'dir')
    mkdir('C:\SERVER\Figures\reach_retract_dlc\reach_increase')
end
if ~exist('C:\SERVER\Figures\reach_retract_dlc\reach_flat', 'dir')
    mkdir('C:\SERVER\Figures\reach_retract_dlc\reach_flat')
end

DATA = {fIncorrect, fCorrect};
RESULTNAMES = {'incorrect', 'correct'};
FTNAMES = {'handContra_xPos', 'tongue'};
FTDISPNAMES = {'Contra forepaw', 'Lick'};
COLORS = arrayfun(@(i) getColor(i, 3, 0.7), [1, 3], UniformOutput=false);
YYAXIS = {'left', 'right'};

% fig = figure(Units='inches', OuterPosition=[1, 1, 24, 8]);
% tl = tiledlayout(fig, 2, 7, TileIndexing='columnmajor');
% for iExpAcute = 1:length(expAcute)
%     for iResult = 1:2
%         ax = nexttile(tl);
%         hold(ax, 'on')
%         h = gobjects(length(FTNAMES), 1);
%         colororder(ax, [getColor(1, 3, 0.7); getColor(3, 3, 0.7)])
%         t = fAll.press.t;
%         for iFt = 1:length(FTNAMES)
%             yyaxis(ax, YYAXIS{iFt})
%             X = DATA{iResult}{iExpAcute}.press(:, find(strcmpi(fnames, FTNAMES{iFt})), :); % t, feature, trial
%             if contains(FTNAMES{iFt}, {'yPos', 'yVel'})
%                 X = -X;
%             end
%             n = size(X, 3);
%             mu = mean(X, 3, 'omitnan');
%             sd = std(X, 0, 3, 'omitnan');
%             col = COLORS{iFt};
%             h(iFt) = plot(ax, t, mu, LineStyle='-', Color=col, LineWidth=1.5, DisplayName=FTDISPNAMES{iFt});
%             selErr = ~isnan(mu + sd);
%             patch(ax, [t(selErr), flip(t(selErr))], [mu(selErr)-sd(selErr); flip(mu(selErr)+sd(selErr))], 'r', ...
%                 LineStyle='-', FaceAlpha=0.075, FaceColor=col, EdgeAlpha=0.075, EdgeColor=col);
%             xline(ax, 0, 'k--')
%             yline(ax, 0, 'k--')
%         end
%         xlim(ax, [-2, 3])
% 
%         yyaxis(ax, 'left')
%         ylabel(ax, 'Contra forepaw AP-pos (a.u.)')
%         ylim(ax, [-6, 6])
%         yticks(ax, [-3, 0, 3])
%         yyaxis(ax, 'right')
%         ylabel(ax, 'Lick probability (a.u.)')
%         ylim(ax, [-1.5, 1.5])
%         yticks(ax, [0, 1])
% 
%         title(ax, sprintf('Reach (%s, n=%i)', RESULTNAMES{iResult}, n))
%     end
%     xlabel(tl, 'Time to bar-contact (s)')
% end


[sel, ~] = ismember(eu.getName, euAcute.getName);


fig = figure(Units='normalized', OuterPosition=[0.1, 0.1, 0.8, 0.8]);
tl = tiledlayout(fig, 4, 2, TileIndexing='columnmajor');
ax = gobjects(4, 1);
for iAx = 1:8
    ax(iAx) = nexttile(tl);
end
for i = find(sel)
    iEuAcute = find(euAcute == eu(i));
    iExpAcute = expIndices(iEuAcute);


    % 1: Raster
    cla(ax(1), 'reset')
    EphysUnit.plotRaster(ax(1), rd.press(i), xlim=[-4, 4]);
    set(ax(1).Legend, AutoUpdate=false)
    xline(ax(1), 0, 'k-', LineWidth=1.5)
    yline(ax(1), find(rd.press(i).duration>4, 1), 'k-', LineWidth=1.5)
    delete(ax(1).Legend)

    % 2: PETH
    cla(ax(2), 'reset')
    hold(ax(2), 'on')
    h = gobjects(2, 1);
    h(1) = plot(ax(2), eta.incorrectPress.t, eta.incorrectPress.X(i, :), 'k-.', DisplayName='Incorrect', LineWidth=1.5);
    h(2) = plot(ax(2), eta.correctPress.t, eta.correctPress.X(i, :), 'k-', DisplayName='Correct', LineWidth=1.5);
    legend(h, AutoUpdate=false)
    xline(ax(2), 0, 'k--')
    yline(ax(2), 0, 'k--')
    xlim(ax(2), [-4, 4])
    ylim(ax(2), [-3, 3])

    % 3/4: DLC incorrect/correct
    for iResult = 1:2
        iAx = iResult + 2;
        cla(ax(iAx), 'reset')
        hold(ax(iAx), 'on')
        h = gobjects(length(FTNAMES), 1);
        colororder(ax(iAx), [getColor(1, 3, 0.7); getColor(3, 3, 0.7)])
        t = fAll.press.t;
        for iFt = 1:length(FTNAMES)
            yyaxis(ax(iAx), YYAXIS{iFt})
            X = DATA{iResult}{iExpAcute}.press(:, find(strcmpi(fnames, FTNAMES{iFt})), :); % t, feature, trial
            n = size(X, 3);
            mu = mean(X, 3, 'omitnan');
            sd = std(X, 0, 3, 'omitnan');
            col = COLORS{iFt};
            h(iFt) = plot(ax(iAx), t, mu, LineStyle='-', Color=col, LineWidth=1.5, DisplayName=FTDISPNAMES{iFt});
            selErr = ~isnan(mu + sd);
            patch(ax(iAx), [t(selErr), flip(t(selErr))], [mu(selErr)-sd(selErr); flip(mu(selErr)+sd(selErr))], 'r', ...
                LineStyle='-', FaceAlpha=0.075, FaceColor=col, EdgeAlpha=0.075, EdgeColor=col);
            xline(ax(iAx), 0, 'k--')
            yline(ax(iAx), 0, 'k--')
        end
        xlim(ax(iAx), [-4, 4])

        yyaxis(ax(iAx), 'left')
        ylabel(ax(iAx), 'Contra forepaw AP-pos (a.u.)')
        ylim(ax(iAx), [-6, 6])
        yticks(ax(iAx), [-3, 0, 3])
        yyaxis(ax(iAx), 'right')
        ylabel(ax(iAx), 'Lick probability (a.u.)')
        ylim(ax(iAx), [-1.5, 1.5])
        yticks(ax(iAx), [0, 1])

        title(ax(iAx), sprintf('DeepLabCut, Reach (%s, n=%i)', RESULTNAMES{iResult}, n))
        xlabel(ax(iAx), 'Time to bar-contact (s)')

    end

    title(ax(1), 'Reach Raster')
    title(ax(2), 'Reach PETH')


    % 5: Raster
    cla(ax(5), 'reset')
    EphysUnit.plotRaster(ax(5), rd.lick(i), xlim=[-4, 4]);
    set(ax(5).Legend, AutoUpdate=false)
    xline(ax(5), 0, 'k-', LineWidth=1.5)
    yline(ax(5), find(rd.lick(i).duration>4, 1), 'k-', LineWidth=1.5)
    delete(ax(5).Legend)

    % 6 PETH
    cla(ax(6), 'reset')
    hold(ax(6), 'on')
    h = gobjects(2, 1);
    h(1) = plot(ax(6), eta.incorrectLick.t, eta.incorrectLick.X(i, :), 'k-.', DisplayName='Incorrect', LineWidth=1.5);
    h(2) = plot(ax(6), eta.correctLick.t, eta.correctLick.X(i, :), 'k-', DisplayName='Correct', LineWidth=1.5);
    legend(h, AutoUpdate=false)
    xline(ax(6), 0, 'k--')
    yline(ax(6), 0, 'k--')
    xlim(ax(6), [-4, 4])
    ylim(ax(6), [-3, 3])

    % 7/8: DLC incorrect/correct
    for iResult = 1:2
        iAx = iResult + 6;
        cla(ax(iAx), 'reset')
        hold(ax(iAx), 'on')
        h = gobjects(length(FTNAMES), 1);
        colororder(ax(iAx), [getColor(1, 3, 0.7); getColor(3, 3, 0.7)])
        t = fAll.lick.t;
        for iFt = 1:length(FTNAMES)
            yyaxis(ax(iAx), YYAXIS{iFt})
            X = DATA{iResult}{iExpAcute}.lick(:, find(strcmpi(fnames, FTNAMES{iFt})), :); % t, feature, trial
            n = size(X, 3);
            mu = mean(X, 3, 'omitnan');
            sd = std(X, 0, 3, 'omitnan');
            col = COLORS{iFt};
            h(iFt) = plot(ax(iAx), t, mu, LineStyle='-', Color=col, LineWidth=1.5, DisplayName=FTDISPNAMES{iFt});
            selErr = ~isnan(mu + sd);
            patch(ax(iAx), [t(selErr), flip(t(selErr))], [mu(selErr)-sd(selErr); flip(mu(selErr)+sd(selErr))], 'r', ...
                LineStyle='-', FaceAlpha=0.075, FaceColor=col, EdgeAlpha=0.075, EdgeColor=col);
            xline(ax(iAx), 0, 'k--')
            yline(ax(iAx), 0, 'k--')
        end
        xlim(ax(iAx), [-4, 4])

        yyaxis(ax(iAx), 'left')
        ylabel(ax(iAx), 'Contra forepaw AP-pos (a.u.)')
        ylim(ax(iAx), [-6, 6])
        yticks(ax(iAx), [-3, 0, 3])
        yyaxis(ax(iAx), 'right')
        ylabel(ax(iAx), 'Lick probability (a.u.)')
        ylim(ax(iAx), [-1.5, 1.5])
        yticks(ax(iAx), [0, 1])

        title(ax(iAx), sprintf('DeepLabCut, Lick (%s, n=%i)', RESULTNAMES{iResult}, n))
        xlabel(ax(iAx), 'Time to bar-contact (s)')

    end

    title(ax(1), 'Reach Raster')
    title(ax(2), 'Reach PETH')
    title(ax(5), 'Lick Raster')
    title(ax(6), 'Lick PETH')



    title(tl, rd.press(i).name, Interpreter='none')

    xlim(ax, [-2, 3])

    % waitforbuttonpress();
    if (c.isPressDown(i))
        print(fig, sprintf('C:\\SERVER\\Figures\\reach_retract_dlc\\reach_decrease\\%s.png', eu(i).getName()), '-dpng')
    elseif (c.isPressUp(i))
        print(fig, sprintf('C:\\SERVER\\Figures\\reach_retract_dlc\\reach_increase\\%s.png', eu(i).getName()), '-dpng')
    else
        print(fig, sprintf('C:\\SERVER\\Figures\\reach_retract_dlc\\reach_flat\\%s.png', eu(i).getName()), '-dpng')
    end
end


%%