
euSC = EphysUnit.load('C:\SERVER\Units\TwoColor_SC\SingleUnit_NonDuplicate_NonDrift_SC');
%
expSC = CompleteExperiment3(euSC, cameras='r');
expSC.alignTimestamps(refEventNameArduino={'LEVER_PRESSED', 'OPTO1_ON', 'OPTO2_ON'}, refEventNameEphys={'PressOn', 'StimOn'}, trialDurationTolerance=0.6);
% exp.alignTimestamps(refEventNameArduino={'OPTO1_ON', 'OPTO2_ON'}, refEventNameEphys={'StimOn'}, trialDurationTolerance=0.15);

%%
pSC.stimBluePowers = [2000]*1e-6; 
pSC.stimRedPowers = [16000]*1e-6;
pSC.stimBlueDurations = [20]*1e-3;
pSC.stimRedDurations = [20]*1e-3;

clear groupsSC trajectoriesSC
groupsSC(length(expSC)) = struct(blue=[], red=[]);
trajectoriesSC(length(expSC)) = struct(HandCameraSide=[], HandOppositeSide=[], Jaw=[]);
for iExp = 1:length(expSC)
    groupsSC(iExp).blue = expSC(iExp).eu(1).groupTwoColorStimTrials({'wavelength', 'power', 'duration'}, selectBy=struct(power=pSC.stimBluePowers, duration=pSC.stimBlueDurations, location=[], wavelength=[470, 473]));
    groupsSC(iExp).red = expSC(iExp).eu(1).groupTwoColorStimTrials({'wavelength', 'power', 'duration'}, selectBy=struct(power=pSC.stimRedPowers, duration=pSC.stimRedDurations, location=[], wavelength=635));
    
    for bodypart = ["HandCameraSide", "HandOppositeSide", "Jaw"]
        for color = ["blue", "red"]
            [trajectoriesSC(iExp).(bodypart).(color).X, trajectoriesSC(iExp).(bodypart).(color).Y, trajectoriesSC(iExp).(bodypart).(color).L, trajectoriesSC(iExp).(bodypart).(color).t] = expSC(iExp).getTrajectoryByTrial('r', char(bodypart), trials=[groupsSC(iExp).(color).trials], window=[-0.5, 0.5], likelihoodThreshold=0.9, includeInvalid=true, alignTo='start', interp='previous');
        end
    end
end

clear iExp bodypart color

% Make session average
for bodypart = ["HandCameraSide", "HandOppositeSide", "Jaw"]
    for color = ["blue", "red"]
        for field = ["X", "Y", "L"]
            fielddata = arrayfun(@(traj) traj.(bodypart).(color).(field), trajectoriesSC(1:length(expSC)), UniformOutput=false);
            trajectoriesSC(length(expSC)+1).(bodypart).(color).(field) = cat(1, fielddata{:});
        end
        trajectoriesSC(length(expSC)+1).(bodypart).(color).t = trajectoriesSC(1).(bodypart).(color).t;
    end
end
clear bodypart color field fielddata

%% Get video clips around time of stim (for manual validation of video/ephys alignment)
clear clips
for color = ["blue", "red"]
    switch color
        case "blue"
            trials = [groupsBlue.trials];
        case "red"
            trials = [groupsRed.trials];
    end
    [clips.(color), t] = expSC.getVideoClip([trials.Start], 'r', numFramesBefore=15, numFramesAfter=15, bodyParts={'HandCameraSide', 'Jaw'}, file='C:\SERVER\desmond39\desmond39_20250522\desmond39_20250522_laser_1.mp4');
end

%% Session average figure, just spd (Alternative comparison)
close all

windowPreStim = [-0.25, 0];
windowPostStim = [0, 0.45];
WAVELENGTHS = [470, 635];
COLORS = ["blue", "red"];
BODYPARTS = ["HandCameraSide", "Jaw"];
BODYPARTDISPNAMES = ["Forepaw", "Jaw"];
YYAXIS = ["left", "right"];
BODYPARTCOLORS = arrayfun(@(i) getColor(i, 3, 0.7), [1, 3], UniformOutput=false);
YLIMS = {[-2, 10], [-1, 5]};
YTICKS = {[0, 6], [0, 3]};
p.fontSize = 9;

iExp = length(trajectoriesSC);
fig = figure(Units="inches", Position=[1, 1, 5, 2]);
tl = tiledlayout(fig, 1, 2);
iTile = 0;
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
        yyaxis(ax, YYAXIS(iBodypart));
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
        ylabel(ax, sprintf("%s speed (a.u.)", BODYPARTDISPNAMES(iBodypart)))
    end
    
    xline(ax, [0, 20], ':')
    xlim(ax, [windowPreStim(1), windowPostStim(2)] * 1e3)
    % fprintf("%i nm, %g mW, %i trials, %i sessions\n", WAVELENGTHS(iColor), mwPower, nTrials, length(exp));
    title(ax, sprintf("%i nm, %g mW", WAVELENGTHS(iColor), mwPower))
    fontsize(ax, p.fontSize, 'points')
end
xlabel(tl, 'Time from laser on (ms)', FontSize=p.fontSize);
% ylabel(tl, 'Speed (a.u.)', FontSize=p.fontSize)
% lgd = legend(h, Orientation='horizontal');
% lgd.Layout.Tile = 'north';

copygraphics(fig, BackgroundColor='none', ContentType='vector')

%% Session average figure, just spd
close all

BODYPARTS = ["HandCameraSide", "Jaw"];
BODYPARTDISPNAMES = ["Forepaw", "Jaw"];
YLIMS = {[-2, 9], [-2, 9]};
p.fontSize = 9;

iExp = length(trajectoriesSC);
fig = figure(Units="inches", Position=[1, 1, 5, 2]);
tl = tiledlayout(fig, 1, 2);
iTile = 0;
h = gobjects(2, 1);
AX = gobjects(length(BODYPARTS), 1);
for iBodypart = 1:length(BODYPARTS)
    bodypart = BODYPARTS(iBodypart);
    ax = nexttile(tl);
    AX(iBodypart) = ax;
    hold(ax, 'on')
    windowPreStim = [-0.25, 0];
    windowPostStim = [0, 0.5];
    WAVELENGTHS = [470, 635];
    COLORS = ["blue", "red"];
    for iColor = 1:length(COLORS)
        color = COLORS(iColor);
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
        switch color
            case "blue"
                col = hsl2rgb([232.48/360, 1, 0.5]);
                mwPower = pSC.stimBluePowers*1e3;
            case "red"
                col = hsl2rgb([355.71/360, 1, 0.5]);
                mwPower = pSC.stimRedPowers*1e3;
        end
        assert(isscalar(mwPower))
        h(iColor) = plot(ax, t*1e3, mu, Color=col, DisplayName=sprintf("%i nm, %g mW", WAVELENGTHS(iColor), mwPower), LineWidth=1.5);
        fprintf("%i nm, %g mW, %i trials, %i sessions\n", WAVELENGTHS(iColor), mwPower, nTrials, length(expSC));
        patch(ax, [t, flip(t)]*1e3, [mu-sd, flip(mu+sd)], col, LineStyle='-', FaceAlpha=0.075, FaceColor=col, EdgeAlpha=0.075, EdgeColor=col)
    end

    xline(ax, 0, '--')
    xline(ax, 20, '--')
    xlim(ax, [windowPreStim(1), windowPostStim(2)] * 1e3)
    ylim(ax, YLIMS{iBodypart})
    yticks(ax, [0, 5])
    title(ax, BODYPARTDISPNAMES(iBodypart))
    fontsize(ax, p.fontSize, 'points')
end
xlabel(tl, 'Time from laser on (ms)', FontSize=p.fontSize);
ylabel(tl, 'Speed (a.u.)', FontSize=p.fontSize)
lgd = legend(h, Orientation='horizontal');
lgd.Layout.Tile = 'north';

copygraphics(fig, BackgroundColor='none', ContentType='vector')

%% Figure with all individual sessions plus session average, showing xvel, yvel, spd
close all
for iExp = 1:length(trajectoriesSC)
    fig = figure(Units="inches", Position=[1, 1, 6, 6]);
    tlp = tiledlayout(fig, 2, 1);
    iTile = 0;
    h = gobjects(2, 1);
    for bodypart = ["HandCameraSide", "Jaw"]
        ax = gobjects(1, 3);
        tl = tiledlayout(tlp, 1, 3);
        iTile = iTile + 1;
        tl.Layout.Tile = iTile;
        ax(1) = nexttile(tl);
        ax(2) = nexttile(tl);
        ax(3) = nexttile(tl);
        hold(ax, 'on')
        windowPreStim = [-0.25, 0];
        windowPostStim = [0, 0.5];
        clear nTrials
        WAVELENGTHS = [470, 635];
        for color = ["blue", "red"]
            t = trajectoriesSC(iExp).(bodypart).(color).t;
            X = trajectoriesSC(iExp).(bodypart).(color).X;
            Y = -trajectoriesSC(iExp).(bodypart).(color).Y;
            % X = (X - mean(X, 'all', 'omitnan')) ./ std(X, 1, 'all', 'omitnan');
            % Y = (Y - mean(Y, 'all', 'omitnan')) ./ std(Y, 1, 'all', 'omitnan');
            nTrials.(color) = size(X, 1);
            velX = diff(X, 1, 2)./diff(t);
            velX = [NaN([size(velX, 1), 1]), velX];
            velY = diff(Y, 1, 2)./diff(t);
            velY = [NaN([size(velY, 1), 1]), velY];
            spd = sqrt(velX.^2 + velY.^2);
            % plot(ax, mean(X(:, t>=windowPostStim(1) & t<=windowPostStim(2)), 1, 'omitnan'), mean(Y(:, t>=windowPostStim(1) & t<=windowPostStim(2)), 1, 'omitnan'), color)
            % plot(ax, mean(X(:, t>=windowPreStim(1) & t<=windowPreStim(2)), 1, 'omitnan'), mean(Y(:, t>=windowPreStim(1) & t<=windowPreStim(2)), 1, 'omitnan'), color, LineStyle=':')
            % plot(ax, X(:, t>=windowPostStim(1) & t<=windowPostStim(2)), Y(:, t>=windowPostStim(1) & t<=windowPostStim(2)), color)
            % plot(ax, X(:, t>=windowPreStim(1) & t<=windowPreStim(2)), Y(:, t>=windowPreStim(1) & t<=windowPreStim(2)), color, LineStyle=':')
            iColor = find(["blue", "red"] == color);
            h(iColor) = plot(ax(1), t, mean(velX, 1, 'omitnan'), color, DisplayName=sprintf("%i nm (%i trials)", WAVELENGTHS(iColor), nTrials.(color)));
            plot(ax(2), t, mean(velY, 1, 'omitnan'), color)
            plot(ax(3), t, mean(spd, 1, 'omitnan'), color)
            % axis(ax, 'equal')
            % axis(ax, 'image')
        end

        for i = 1:3
            xline(ax(i), 0, '--')
            yline(ax(i), 0, '--')
            ylim(ax(i), [-60, 60])
            xlim(ax(i), [windowPreStim(1), windowPostStim(2)])
        end
        ylim(ax(3), [-60, 60])

        % title(ax(1), 'velX')
        % title(ax(2), 'velY')
        % title(ax(3), 'spd')
        ylabel(ax(1), 'velX')
        ylabel(ax(2), 'velY')
        ylabel(ax(3), 'spd')
        xlabel(ax, 't')

        title(tl, char(bodypart), FontSize=9)
    end
    lgd = legend(h);
    lgd.Layout.Tile = 'south';

    if iExp <= length(expSC)
        title(tlp, expSC(iExp).name, Interpreter='none')
    else
        title(tlp, sprintf('%i sessions', length(expSC)))
    end
end
