% %%
% load_ephysunits
p.kernel = struct(type='exponential', params=struct(lambda1=5, lambda2=10, resolution=0.01, width=1));
%%
doPlot = questdlg('Generate and save plots?');

switch lower(doPlot)
    case 'yes'
        doPlot = true;
    case 'no'
        doPlot = false;
    otherwise
        error('Operation cancelled by user.')
end

TRIALTYPES = {'lick', 'lick', 'press', 'press'};
ISUP = [true, false, true, false];
SELUNITS = {c.hasPress & c.hasLick & c.isLickUp, c.hasPress & c.hasLick & c.isLickDown, c.hasPress & c.isPressUp, c.hasPress & c.isPressDown};
TRAINWINDOW = {[-1.5, 0], [-1.5, 0], [-1.5, 0], [-1.5, 0]};
PATHNAME = ["LickUp", "LickDown", "PressUp", "PressDown"];
YLIM = {[-2.5, 5], [-5, 2.5], [-2.5, 5], [-5, 2.5]};

PARAM = ["tau", "a", "c"];
PARAMNAME = ["\tau", "a", "c"];
PARAMYLIM = {[0, 1], [0.5, 1.5], [-0.75, 0.75]};
minTrialLength = 1.5;
minNumTrials = 20;
cueBlankWindow = 0.5;

% Plot
close all
clear substats stats
substats = struct(slope=NaN, p=NaN);
stats(length(eu)) = struct(press=struct(tau=substats, a=substats, c=substats, mdl=[]), lick=struct(tau=substats, a=substats, c=substats, mdl=[]));
for iGroup = 1:4
    trialType = TRIALTYPES{iGroup};
    euIndices = find(SELUNITS{iGroup});
    isUp = ISUP(iGroup);
    trainWindow = TRAINWINDOW{iGroup};

    for i = 1:length(euIndices)
        iEu = euIndices(i); 
        trials = eu(iEu).getTrials(trialType);
        trials = trials(trials.duration() >= minTrialLength);
        [~, I] = sort(trials.duration(), 'ascend');
        trials = trials(I);
        fprintf('%i trials longer than %gs.\n', length(trials), minTrialLength)
        if length(trials) < minNumTrials
            continue
        end

        % Get x, t, z
        [xx, t] = eu(iEu).getTrialAlignedData('rate', window=[-4, 2], trialType=trialType, resolution=0.01, trials=trials, allowedTrialDuration=[0, Inf], includeInvalid=false, startBlankWindow=[-Inf, cueBlankWindow], kernel=p.kernel);
        zz = trials.duration();

        % Normalize x
        inNormWindow = t >= -4 & t <= -2;
        inRespWindow = t >= -2 & t <= 0;
        mu = mean(xx(:, inNormWindow), 'all', 'omitnan');
        sigma = std(xx(:, inNormWindow), 0, 'all', 'omitnan');
        xx = (xx - mu) ./ sigma;

        % Do a binned average of x and z
        k = 5;
        nTrials = size(xx, 1);
        nBins = ceil(nTrials/k);
        x = NaN(nBins, size(xx, 2));
        z = NaN(nBins, 1);
        for iBin = 1:nBins
            iTrial = 1 + (iBin-1)*k;
            selTrials = iTrial:min(iTrial+k, nTrials);
            x(iBin, :) = mean(xx(selTrials, :), 1, 'omitnan');
            z(iBin) = mean(zz(selTrials));
        end
        nTrials = nBins;

        clear iTrial xSel
        % Fit exponential to movement aligned spike rate
        selT = t >= trainWindow(1) & t <= trainWindow(2);
        tTrain = t(selT);
        xTrain = x(:, selT);
        
        clear mdl
        mdl = struct(a=NaN(size(xTrain, 1), 1), b=NaN(size(xTrain, 1), 1), c=NaN(size(xTrain, 1), 1), xHat=NaN(size(xTrain)), tau=[], z=[]);
        if isUp
            g = fittype(@(a, b, c, t) a*exp(b*t)+c, dependent='x', independent='t');
        else
            g = fittype(@(a, b, c, t) -a*exp(b*t)+c, dependent='x', independent='t');
        end
            
        for iTrial = 1:size(xTrain, 1)
            try
                isValid = ~isnan(xTrain(iTrial, :));
                thisMdl = fit(tTrain(isValid)', xTrain(iTrial, isValid)', g, StartPoint=[1, 10, 0], Lower=[0.1, 0.1, -5], Upper=[10, 100, 5]);
                mdl.a(iTrial) = thisMdl.a;
                mdl.b(iTrial) = thisMdl.b;
                mdl.c(iTrial) = thisMdl.c;
                mdl.xHat(iTrial, :) = thisMdl(tTrain');
            end
        end
        mdl.tau = 1./mdl.b;
        mdl.z = z;

        if doPlot
            if exist('fig', 'var') && isvalid(fig)
                close(fig)
            end
            fig = figure(Units='inches', Position=[1 1 10 7]);
            tl = tiledlayout(fig, 2, 6);

            ax = nexttile(tl, [1, 3]); hold(ax, 'on')
            for iTrial = 1:nTrials
                plot(ax, t, x(iTrial, :), Color=[getColor(iTrial, nTrials, 0.8), 0.35]);
            end
            xlabel('time to move (inv, s)')
            xlim(ax, [-2, 0.5])
            ylabel('norm spike rate (a.u.)')
            xline(ax, trainWindow(1))
            xline(ax, trainWindow(2))
            ylim(ax, YLIM{iGroup})

            ax = nexttile(tl, [1, 3]); hold(ax, 'on')
            for iTrial = 1:nTrials
                plot(ax, tTrain, mdl.xHat(iTrial, :), Color=[getColor(iTrial, nTrials, 0.8), 0.35]);
            end
            xlabel('time to move (inv, s)')
            xlim(ax, [-2, 0.5])
            ylabel('norm spike rate (a.u.)')
            xline(ax, trainWindow(1))
            xline(ax, trainWindow(2))
            ylim(ax, YLIM{iGroup})
            title(ax, 'y = a*exp(-x/\tau) + c')
            print(fig, sprintf('C:\\SERVER\\Figures\\SingleTrialExponential\\%s\\%s.png', PATHNAME(iGroup), eu(iEu).getName()), '-dpng')
        end

        % Fit LM to params of exponential (tau, a, c) vs. movement time (z)
        clear thisMdl
        stats(iEu).(trialType).mdl = mdl;
        for iParam = 1:3
            thisMdl = fitlm(mdl.z, mdl.(PARAM(iParam)), VarNames=["z", PARAM(iParam)]);
            stats(iEu).(trialType).(PARAM(iParam)).slope = table2array(thisMdl.Coefficients('z', 'Estimate'));
            stats(iEu).(trialType).(PARAM(iParam)).p = table2array(thisMdl.Coefficients('z', 'pValue'));

            if doPlot
                ax = nexttile(tl, [1, 2]); hold(ax, 'on')
                plot(ax, [min(mdl.z); max(mdl.z)], thisMdl.predict([min(mdl.z); max(mdl.z)]), 'k--', LineWidth=1.5)
                scatter(ax, mdl.z, mdl.(PARAM(iParam)), 10, getColor(1:length(mdl.z), length(mdl.z), 0.8))
                xlabel(ax, 'z = move time (s)')
                ylabel(ax, PARAMNAME(iParam), Interpreter='tex')
                title(ax, sprintf('%s: slope = %.3f, p < %.3f', PARAMNAME(iParam), table2array(thisMdl.Coefficients('z', 'Estimate')), table2array(thisMdl.Coefficients('z', 'pValue'))), Interpreter='tex')
                xlim(ax, [0, 10])
                if min(mdl.(PARAM(iParam))) >= PARAMYLIM{iParam}(1) && max(mdl.(PARAM(iParam))) <= PARAMYLIM{iParam}(2)
                    ylim(ax, PARAMYLIM{iParam})
                end
            end
        end
    end
end

save('C:\SERVER\Figures\SingleTrialExponential\SingleTrialExponential.mat', 'stats', 'p', 'c')

%% Plot params
TRIALTYPES = {'press', 'press', 'lick', 'lick'};
UNITTYPES = {'press increase', 'press decrease', 'lick increase', 'lick decrease'};
SELUNITS = {c.hasPress & c.isPressUp; c.hasPress & c.isPressDown; c.hasPress & c.hasLick & c.isLickUp; c.hasPress & c.hasLick & c.isLickDown};
fig = figure(Units='inches', Position=[1 1 9 6]);
tlp = tiledlayout(4, 1);

for iTrialType = 1:4
    trialType = TRIALTYPES{iTrialType};
    tl = tiledlayout(tlp, 1, 3);
    tl.Layout.Tile = iTrialType;
    for paramName = ["tau", "a", "c"]
        ax = nexttile(tl);
        hold(ax, 'on')
        histogram(arrayfun(@(p) p.(paramName).slope, [stats(SELUNITS{iTrialType}).(trialType)]), -2:0.1:2, FaceColor='blue');
        hold(ax, 'off')
        title(ax, paramName)
        xlabel(ax, 'slope')
        ylabel(ax, 'count')
    end
    title(tl, UNITTYPES{iTrialType})
end

%% Plot params (new)
% S2b (reach)
% 1 row (inc neurons): PETH sig; PETH insig;
% 2 row (dec neurons): PETH sig; PETH insig;
% 3 row: distribution of slope: inc vs. dec (color by p value)
% 4 row (inc neurons): PETH sig, {ETH insig;
% 5 row (dec neurons): PETH sig, {ETH insig;
% S6b
clear layout
layout.w = 6;
layout.h = 6;
layout.top.h = 3;
layout.middle.h = 1;
layout.bottom.h = 3;

close all
fig = figure(Units='inches', Position=[1 1 layout.w layout.h]);

layout.tl = tiledlayout(fig, layout.top.h + layout.middle.h + layout.bottom.h, 1);
layout.top.tl = tiledlayout(layout.tl, 2, 2);
l = layout.top.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [layout.top.h, 1];
layout.middle.tl = tiledlayout(layout.tl, 1, 2);
l = layout.middle.tl; l.Layout.Tile = 1 + layout.top.h; l.Layout.TileSpan = [layout.middle.h, 1];
layout.bottom.tl = tiledlayout(layout.tl, 2, 2);
l = layout.bottom.tl; l.Layout.Tile = 1 + layout.top.h + layout.middle.h; l.Layout.TileSpan = [layout.bottom.h, 1];

ax = nexttile(layout.middle.tl);


%%
% X = spike rate (nTrials, nTimestampsPerTrial), t = timestamp (1, nTimestampsPerTrial), z = movement time (nTrials, 1)
function boot = bootParamSlope(X, t, z, varargin)
    p = inputParser();
    p.addRequired('X', @isnumeric)
    p.addRequired('t', @isnumeric)
    p.addRequired('z', @isnumeric)
    p.addParameter('nboot', 10000, @(n) n>=10000);
    p.addParameter('isUp', true, @islogical);
    p.addParameter('k', 5, @isnumeric);
    p.parse(X, t, z, varargin{:})
    r = p.Results;
    k = r.k;

    nTrials = size(X, 1);
    nTimestamps = size(X, 2);
    assert(length(t) == nTimestamps);
    assert(length(z) == nTrials);
    t = reshape(t, 1, []);
    z = reshape(z, [], 1);

    boot.args = r;
    boot.tau = struct(obs=[], ci95=[], ci99=[], p=[]);
    boot.a = struct(obs=[], ci95=[], ci99=[], p=[]);
    boot.c = struct(obs=[], ci95=[], ci99=[], p=[]);

    % Do a binned average of x and z
    k = 5;
    nTrials = size(xx, 1);
    nBins = ceil(nTrials/k);
    x = NaN(nBins, size(xx, 2));
    z = NaN(nBins, 1);
    for iBin = 1:nBins
        iTrial = 1 + (iBin-1)*k;
        selTrials = iTrial:min(iTrial+k, nTrials);
        x(iBin, :) = mean(xx(selTrials, :), 1, 'omitnan');
        z(iBin) = mean(zz(selTrials));
    end
    nTrials = nBins;

    if r.isUp
        g = fittype(@(a, b, c, t) a*exp(b*t)+c, dependent='x', independent='t');
    else
        g = fittype(@(a, b, c, t) -a*exp(b*t)+c, dependent='x', independent='t');
    end

    for iboot = 1:r.nboot
        % shuffle z

        a = NaN(nTrials, 1);
        b = a;
        c = a;
        xHat = NaN(nTrials, nTimestamps);
        for iTrial = nTrials
            isValid = ~isnan(X(iTrial, :));
            thisMdl = fit(t(isValid)', X(iTrial, isValid)', g, StartPoint=[1, 10, 0], Lower=[0.1, 0.1, -5], Upper=[10, 100, 5]);
            a(iTrial) = thisMdl.a;
            b(iTrial) = thisMdl.b;
            c(iTrial) = thisMdl.c;
            xHat(iTrial, :) = thisMdl(t');
        end
        tau = 1./b;
    end
end