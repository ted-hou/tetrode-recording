%%
load_ephysunits;
%% Align ArduinoConnection events to ephys time using a common event (try CueOn)
clear ac
[ac.ac, ac.euIndicesFirstInSession, ac.expIndices, ac.uniqueExpNames] = eu.loadArduinoConnection();
%%
tEU = alignTimestamps(eu, ["REWARD_ON", "LEVER_RELEASED", "LEVER_RETRACT_START", "LEVER_RETRACT_END", "TUBE_RETRACT_START", "TUBE_RETRACT_END"], ac={ac.ac, ac.euIndicesFirstInSession, ac.expIndices, ac.uniqueExpNames});

for iEu = 1:length(eu)
    try
        leverReleaseTimes = eu(iEu).EventTimes.LEVER_RELEASED;
        leverRetractTimes = eu(iEu).EventTimes.LEVER_RETRACT_START;
        eu(iEu).Trials.Press = Trial(eu(iEu).EventTimes.Cue, eu(iEu).EventTimes.Press, 'first');
        eu(iEu).Trials.PressIncorrect = eu(iEu).Trials.Press(eu(iEu).Trials.Press.duration() < 4 & eu(iEu).Trials.Press.duration() >= 2);
        eu(iEu).Trials.PressCorrect = eu(iEu).Trials.Press(eu(iEu).Trials.Press.duration() >= 4);
        [~, leverRetractTimesIncorrect] = eu(iEu).Trials.PressIncorrect.inTrial(leverRetractTimes, [0, 2], windowMode='stop');
        assert(nnz(leverRetractTimesIncorrect) > 0)
        [~, leverRetractTimesCorrect] = eu(iEu).Trials.PressCorrect.inTrial(leverRetractTimes, [-0.1, 8], windowMode='stop');
        eu(iEu).Trials.RetractReleaseIncorrect = Trial([leverRetractTimesIncorrect, Inf], leverReleaseTimes, 'first');
        eu(iEu).Trials.RetractReleaseCorrect = Trial([leverRetractTimesCorrect, Inf], leverReleaseTimes, 'first');
        eu(iEu).Trials.PressReleaseIncorrect = Trial([eu(iEu).Trials.PressIncorrect.Start, Inf], [eu(iEu).Trials.RetractReleaseIncorrect.Stop], 'first');
        eu(iEu).Trials.PressReleaseCorrect = Trial([eu(iEu).Trials.PressCorrect.Stop, Inf], [eu(iEu).Trials.RetractReleaseCorrect.Stop], 'first');
        % 
        % fprintf(['press=%i\npressIncorrect=%i, pressCorrect=%i;\nleverRetractTimesIncorrect=%i, leverRetractTimesCorrect=%i;\n' ...
        %     'eu(iEu).Trials.RetractReleaseIncorrect=%i, eu(iEu).Trials.RetractReleaseCorrect=%i;\n' ...
        %     'eu(iEu).Trials.PressReleaseIncorrect=%i, eu(iEu).Trials.PressReleaseCorrect=%i\n'], length(eu(iEu).Trials.Press), length(eu(iEu).Trials.PressIncorrect), length(eu(iEu).Trials.PressCorrect), ...
        %     length(leverRetractTimesIncorrect), length(leverRetractTimesCorrect), ...
        %     length(eu(iEu).Trials.RetractReleaseIncorrect), length(eu(iEu).Trials.RetractReleaseCorrect), ...
        %     length(eu(iEu).Trials.PressReleaseIncorrect), length(eu(iEu).Trials.PressReleaseCorrect));
        % disp(1)
    catch
        warning('Counld not process iEu = %i, "%s"', iEu, eu(iEu).getName())
    end
end
%%
read_DLC_data;
%% Assign exp to each eu
expIndices = zeros(size(euAcute));
for iExp = 1:length(expAcute)
    expIndices(ismember(euAcute, expAcute(iExp).eu)) = iExp;
end


%%

clear rd
rd.press = eu.getRasterData('press', [-4, 4], minTrialDuration=2, maxTrialDuration=Inf);
rd.releaseCorrect = eu.getRasterData('press_release_correct', [-4, 4], minTrialDuration=0, maxTrialDuration=Inf);
rd.releaseIncorrect = eu.getRasterData('press_release_incorrect', [-4, 4], minTrialDuration=0, maxTrialDuration=Inf);
rd.lick = eu.getRasterData('lick', [-4, 4], minTrialDuration=2, maxTrialDuration=Inf);

eta.correctPress = eu.getETA('count', 'press', [-4, 4], minTrialDuration=4, normalize='none', resolution=0.1);
eta.incorrectPress = eu.getETA('count', 'press', [-4, 4], minTrialDuration=p.minTrialDuration, maxTrialDuration=4, normalize='none', resolution=0.1);

eta.correctRelease = eu.getETA('count', 'press_release_correct', [-4, 4], alignTo='stop', normalize='none', resolution=0.1);
eta.incorrectRelease = eu.getETA('count', 'press_release_incorrect', [-4, 4], alignTo='stop', normalize='none', resolution=0.1);

eta.correctLick = eu.getETA('count', 'lick', [-4, 4], minTrialDuration=4, normalize='none', resolution=0.1);
eta.incorrectLick= eu.getETA('count', 'lick', [-4, 4], minTrialDuration=p.minTrialDuration, maxTrialDuration=4, normalize='none', resolution=0.1);

for field = ["correctPress", "incorrectPress", "correctLick", "incorrectLick", "correctRelease", "incorrectRelease"]
    eta.(field).X = eta.(field).X ./ 0.1;
end

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
if ~exist('E:\DATA\Figures\reach_retract_dlc\reach_decrease', 'dir')
    mkdir('E:\DATA\Figures\reach_retract_dlc\reach_decrease')
end
if ~exist('E:\DATA\Figures\reach_retract_dlc\reach_increase', 'dir')
    mkdir('E:\DATA\Figures\reach_retract_dlc\reach_increase')
end
if ~exist('E:\DATA\Figures\reach_retract_dlc\reach_flat', 'dir')
    mkdir('E:\DATA\Figures\reach_retract_dlc\reach_flat')
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
tl = tiledlayout(fig, 4, 3, TileIndexing='columnmajor');
ax = gobjects(12, 1);
for iAx = 1:12
    ax(iAx) = nexttile(tl);
end
for i = find(sel)
    try
        iEuAcute = find(euAcute == eu(i));
        iExpAcute = expIndices(iEuAcute);
    
    
        % 1: Raster
        cla(ax(1), 'reset')
        EphysUnit.plotRaster(ax(1), rd.press(i), xlim=[-4, 4]);
        set(ax(1).Legend, AutoUpdate=false)
        xline(ax(1), 0, 'k-', LineWidth=1.5)
        yline(ax(1), find(rd.press(i).duration>4, 1), 'k-', LineWidth=3)
        delete(ax(1).Legend)
        legend(ax(1), ["", "trial start"])
    
        % 2: PETH
        cla(ax(2), 'reset')
        hold(ax(2), 'on')
        h = gobjects(2, 1);
        h(1) = plot(ax(2), eta.incorrectPress.t, eta.incorrectPress.X(i, :), Color=hsl2rgb([0.4, 0.4, 0.4]), LineStyle='-', DisplayName='Incorrect', LineWidth=2);
        h(2) = plot(ax(2), eta.correctPress.t, eta.correctPress.X(i, :), 'k-', DisplayName='Correct', LineWidth=2);
        legend(h, AutoUpdate=false)
        xline(ax(2), 0, 'k--')
        yline(ax(2), 0, 'k--')
        xlim(ax(2), [-4, 4])
        % ylim(ax(2), [0, 100])
    
        % 3/4: DLC incorrect/correct
        for iResult = 1:2
            iAx = iResult + 2;
            cla(ax(iAx), 'reset')
            hold(ax(iAx), 'on')
            h = gobjects(length(FTNAMES), 1);
            colororder(ax(iAx), [getColor(1, 3, 0.7); getColor(3, 3, 0.7)])
            tEU = fAll.press.t;
            for iFt = 1:length(FTNAMES)
                yyaxis(ax(iAx), YYAXIS{iFt})
                X = DATA{iResult}{iExpAcute}.press(:, find(strcmpi(fnames, FTNAMES{iFt})), :); % t, feature, trial
                n = size(X, 3);
                mu = mean(X, 3, 'omitnan');
                sd = std(X, 0, 3, 'omitnan');
                col = COLORS{iFt};
                h(iFt) = plot(ax(iAx), tEU, mu, LineStyle='-', Color=col, LineWidth=1.5, DisplayName=FTDISPNAMES{iFt});
                selErr = ~isnan(mu + sd);
                patch(ax(iAx), [tEU(selErr), flip(tEU(selErr))], [mu(selErr)-sd(selErr); flip(mu(selErr)+sd(selErr))], 'r', ...
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
        EphysUnit.plotRaster(ax(5), rd.releaseCorrect(i), xlim=[-4, 4]);
        set(ax(5).Legend, AutoUpdate=false)
        xline(ax(5), 0, 'k-', LineWidth=1.5)
        % yline(ax(5), find(rd.lick(i).duration>4, 1), 'k-', LineWidth=1.5)
        delete(ax(5).Legend)
        legend(ax(5), ["", "bar-contact"])
    
        % 6 PETH
        cla(ax(6), 'reset')
        hold(ax(6), 'on')
        h = gobjects(2, 1);
        h(1) = plot(ax(6), eta.incorrectRelease.t, eta.incorrectRelease.X(i, :), Color=hsl2rgb([0.4, 0.4, 0.4]), LineStyle='-', DisplayName='Incorrect', LineWidth=2);
        h(2) = plot(ax(6), eta.correctRelease.t, eta.correctRelease.X(i, :), 'k-', DisplayName='Correct', LineWidth=2);
        legend(h, AutoUpdate=false)
        xline(ax(6), 0, 'k--')
        yline(ax(6), 0, 'k--')
        xlim(ax(6), [-4, 4])
        % ylim(ax(6), [0, 100])
    
        % 7/8: DLC incorrect/correct
        for iResult = 1:2
            iAx = iResult + 6;
            cla(ax(iAx), 'reset')
            hold(ax(iAx), 'on')
            h = gobjects(length(FTNAMES), 1);
            colororder(ax(iAx), [getColor(1, 3, 0.7); getColor(3, 3, 0.7)])
            tEU = fAll.lick.t;
            for iFt = 1:length(FTNAMES)
                yyaxis(ax(iAx), YYAXIS{iFt})
                X = DATA{iResult}{iExpAcute}.press_release(:, find(strcmpi(fnames, FTNAMES{iFt})), :); % t, feature, trial
                n = size(X, 3);
                mu = mean(X, 3, 'omitnan');
                sd = std(X, 0, 3, 'omitnan');
                col = COLORS{iFt};
                h(iFt) = plot(ax(iAx), tEU, mu, LineStyle='-', Color=col, LineWidth=1.5, DisplayName=FTDISPNAMES{iFt});
                selErr = ~isnan(mu + sd);
                patch(ax(iAx), [tEU(selErr), flip(tEU(selErr))], [mu(selErr)-sd(selErr); flip(mu(selErr)+sd(selErr))], 'r', ...
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
    
            title(ax(iAx), sprintf('DeepLabCut, Retract (%s, n=%i)', RESULTNAMES{iResult}, n))
            xlabel(ax(iAx), 'Time to bar-release (s)')
    
        end
    
        % 9: Raster
        cla(ax(9), 'reset')
        EphysUnit.plotRaster(ax(9), rd.lick(i), xlim=[-4, 4]);
        set(ax(9).Legend, AutoUpdate=false)
        xline(ax(9), 0, 'k-', LineWidth=1.5)
        yline(ax(9), find(rd.lick(i).duration>4, 1), 'k-', LineWidth=3)
        delete(ax(9).Legend)
        legend(ax(9), ["", "trial start"])
    
        % 10 PETH
        cla(ax(10), 'reset')
        hold(ax(10), 'on')
        h = gobjects(2, 1);
        h(1) = plot(ax(10), eta.incorrectLick.t, eta.incorrectLick.X(i, :), Color=hsl2rgb([0.4, 0.4, 0.4]), LineStyle='-', DisplayName='Incorrect', LineWidth=2);
        h(2) = plot(ax(10), eta.correctLick.t, eta.correctLick.X(i, :), 'k-', DisplayName='Correct', LineWidth=2);
        legend(h, AutoUpdate=false)
        xline(ax(10), 0, 'k--')
        yline(ax(10), 0, 'k--')
        xlim(ax(10), [-4, 4])
        % ylim(ax(6), [0, 100])
    
        % 11/12: DLC incorrect/correct
        for iResult = 1:2
            iAx = iResult + 10;
            cla(ax(iAx), 'reset')
            hold(ax(iAx), 'on')
            h = gobjects(length(FTNAMES), 1);
            colororder(ax(iAx), [getColor(1, 3, 0.7); getColor(3, 3, 0.7)])
            tEU = fAll.lick.t;
            for iFt = 1:length(FTNAMES)
                yyaxis(ax(iAx), YYAXIS{iFt})
                X = DATA{iResult}{iExpAcute}.lick(:, find(strcmpi(fnames, FTNAMES{iFt})), :); % t, feature, trial
                n = size(X, 3);
                mu = mean(X, 3, 'omitnan');
                sd = std(X, 0, 3, 'omitnan');
                col = COLORS{iFt};
                h(iFt) = plot(ax(iAx), tEU, mu, LineStyle='-', Color=col, LineWidth=1.5, DisplayName=FTDISPNAMES{iFt});
                selErr = ~isnan(mu + sd);
                patch(ax(iAx), [tEU(selErr), flip(tEU(selErr))], [mu(selErr)-sd(selErr); flip(mu(selErr)+sd(selErr))], 'r', ...
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
            xlabel(ax(iAx), 'Time to spout-contact (s)')
    
        end
    
        title(ax(1), 'Reach Raster')
        title(ax(2), 'Reach PETH')
        title(ax(5), 'Retract Raster (correct trials)')
        title(ax(6), 'Retract PETH')
        title(ax(11), 'Lick Raster')
        title(ax(12), 'Lick PETH')
        xlabel(ax, '')
        xlabel(ax(4), 'Time to bar-contact (s)')
        xlabel(ax(8), 'Time to bar-release (s)')
        xlabel(ax(12), 'Time to spout-contact (s)')
    
    
    
        title(tl, rd.press(i).name, Interpreter='none')
    
        xlim(ax, [-3, 3])
    
        ylim(ax([2, 6]), [min([ax(2).YLim, ax(6).YLim]), max([ax(2).YLim, ax(6).YLim])]);
    
        % waitforbuttonpress();
        if (c.isPressDown(i))
            print(fig, sprintf('E:\\DATA\\Figures\\reach_retract_dlc\\reach_decrease\\%s.png', eu(i).getName()), '-dpng')
        elseif (c.isPressUp(i))
            print(fig, sprintf('E:\\DATA\\Figures\\reach_retract_dlc\\reach_increase\\%s.png', eu(i).getName()), '-dpng')
        else
            print(fig, sprintf('E:\\DATA\\Figures\\reach_retract_dlc\\reach_flat\\%s.png', eu(i).getName()), '-dpng')
        end
    catch
        error()
    end
end

