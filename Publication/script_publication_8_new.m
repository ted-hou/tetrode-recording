read_SC_opto_trajectories_DLC;
[SNr_SCRetro.eu, SNr_SCRetro.rd, SNr_SCRetro.eta, SNr_SCRetro.meta, SNr_SCRetro.p, SNr_SCRetro.c, SNr_SCRetro.boot] = read_SNr_SCRetro();

%% Fig 8a SC Stim causes movements (medial SC stim vs. lateral SC stim)
close all

useYYAxis = false;
windowPreStim = [-0.25, 0];
windowPostStim = [0, 0.45];
WAVELENGTHS = [470, 635];
COLORS = ["blue", "red"];
BODYPARTS = ["HandCameraSide", "Jaw"];
BODYPARTDISPNAMES = ["Forepaw", "Jaw"];
YYAXIS = ["left", "right"];
BODYPARTCOLORS = arrayfun(@(i) getColor(i, 3, 0.7), [1, 3], UniformOutput=false);
% YLIMS = {[-2, 10], [-1, 5]};
YLIMS = {[-1, 5], [-1, 5]};
% YTICKS = {[0, 6], [0, 3]};
YTICKS = {[0, 3], [0, 3]};
p.fontSize = 9;

iExp = length(trajectoriesSC);
fig = figure(Units="inches", Position=[1, 1, 3.5, 1.75]);
tl = tiledlayout(fig, 1, 2);
h = gobjects(2, 1);
AX = gobjects(length(BODYPARTS), 1);
for iColor = 1:length(COLORS)
    color = COLORS(iColor);
    switch color
        case "blue"
            mwPower = pSC.stimBluePowers*1e3;
        case "red"
            mwPower = pSC.stimRedPowers*1e3;
    end

    ax = nexttile(tl);
    AX(iColor) = ax;
    hold(ax, 'on')
    colororder(ax, getColor([1, 3], 3, 0.7))
    for iBodypart = 1:length(BODYPARTS)
        if useYYAxis
            yyaxis(ax, YYAXIS(iBodypart));
        end
        bodypart = BODYPARTS(iBodypart);
        t = trajectoriesSC(iExp).(bodypart).(color).t;
        X = trajectoriesSC(iExp).(bodypart).(color).X;
        Y = -trajectoriesSC(iExp).(bodypart).(color).Y;
        nTrials = size(X, 1);
        velX = diff(X, 1, 2)./diff(t);
        velX = [NaN([size(velX, 1), 1]), velX];
        velY = diff(Y, 1, 2)./diff(t);
        velY = [NaN([size(velY, 1), 1]), velY];
        spd = sqrt(velX.^2 + velY.^2);
        spd = (spd - mean(spd(:, t<0), 'all', 'omitnan')) ./ std(spd(:, t<0), 0, 'all', 'omitnan');
        mu = mean(spd, 1, 'omitnan');
        sd = std(spd, 0, 1, 'omitnan');
        mu(isnan(mu)) = 0;
        sd(isnan(sd)) = 0;
        assert(isscalar(mwPower))
        col = BODYPARTCOLORS{iBodypart};
        h(iBodypart) = plot(ax, t*1e3, mu, Color=col, DisplayName=BODYPARTDISPNAMES(iBodypart), LineWidth=1.5);
        fprintf("%s %s\n", BODYPARTDISPNAMES(iBodypart), num2str(col));
        patch(ax, [t, flip(t)]*1e3, [mu-sd, flip(mu+sd)], col, LineStyle='-', FaceAlpha=0.075, FaceColor=col, EdgeAlpha=0.075, EdgeColor=col)

        ylim(ax, YLIMS{iBodypart})
        yticks(ax, YTICKS{iBodypart})
        if useYYAxis
            ylabel(ax, sprintf("%s", BODYPARTDISPNAMES(iBodypart)))
        end
    end
    
    xline(ax, [0, 20], ':')
    xlim(ax, [windowPreStim(1), windowPostStim(2)] * 1e3)
    % fprintf("%i nm, %g mW, %i trials, %i sessions\n", WAVELENGTHS(iColor), mwPower, nTrials, length(exp));
    title(ax, sprintf("%i nm, %g mW", WAVELENGTHS(iColor), mwPower))
    fontsize(ax, p.fontSize, 'points')
    set(ax.YAxis, TickLength=[0.04, 0.025])
end
xlabel(tl, 'Time from laser on (ms)', FontSize=p.fontSize);
ylabel(tl, 'Speed (a.u.)', FontSize=p.fontSize)
if ~useYYAxis
    lgd = legend(h, Orientation='horizontal');
    lgd.Layout.Tile = 'north';
end

clear windowPreStim windowPostStim WAVELENGTHS COLORS BODYPARTS BODYPARTDISPNAMES YYAXIS BODYPARTCOLORS YLIMS YTICKS
clear iExp tl h AX iColor color mwPower ax iBodypart bodypart X Y t nTrials velX velY spd mu sd col h lgd


copygraphics(fig, BackgroundColor='none', ContentType='vector')

%% Fig 8b. Rasters of optotagging SNr neurons
