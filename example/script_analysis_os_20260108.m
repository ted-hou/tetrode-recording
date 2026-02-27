%% THIS SCRIPTS IS BASICALLY FOR RASTERS 

%% LOAD DATA

% Load sorted EphysUnits (takes a few minutes)
eu = EphysUnit.load('C:\SERVER\Units\TwoColor_Striatonigral\SingleUnit_NonDuplicate_NonDrift_SNr', waveforms=false, spikecounts=false, spikerates=false);

% ETA Stim (binned spike counts)

% saved earlier as rates with etaRes 0.025
load('C:/GIT/mouse-behavior-arduino/MATLAB/AutoOpto/decoder/Olivia/analysis/eta_groups.mat');

%% raw data peek - rasters

% a) basic intro

% b) low power vs high power

% c) blue light vs red light


%%

expName = "desmond42_20251121";  % <-- whatever ExpName string is
cond = 79;    % 34, 74      ; 39, 79           % 1..80
bandColor = [1 0 0];

[spikesByUnit, unitIdx] = plotPopRasterFirstTrialByExp(eu, expName, cond);

win = [-0.5 0.5];

% ----- population PSTH parameters -----
binW = 0.005;                % 5 ms bins (change as you like)
tEdges = win(1):binW:win(2);
tCenters = tEdges(1:end-1) + binW/2;

smoothSigma = 0.015;         % 15 ms Gaussian smoothing (seconds)
smoothBins = max(1, round(smoothSigma/binW));
x = (-4*smoothBins):(4*smoothBins);
g = exp(-(x.^2)/(2*smoothBins^2));
g = g / sum(g);

% ----- build population rate (mean across units) -----
nUnits = numel(spikesByUnit);
counts = zeros(1, numel(tCenters));

for ii = 1:nUnits
    sp = spikesByUnit{ii};
    if isempty(sp), continue; end
    counts = counts + histcounts(sp, tEdges);
end

% mean firing rate per unit (spikes/s/unit)
popRate = (counts / binW) / max(1, nUnits);

% smooth
popRateSm = conv(popRate, g, 'same');

% ----- figure layout: top trace, bottom raster -----
figure;

% ========== TOP: population trace ==========
ax1 = subplot(2,1,1); hold(ax1, 'on');

% optional highlight band in top trace too
bandX = [0.00 0.01];          % example

bandAlpha = 0.2;            % 0..0.01

yl1 = [0 60];    % fixed population-rate scale

patch(ax1, [bandX(1) bandX(2) bandX(2) bandX(1)], ...
           [yl1(1) yl1(1) yl1(2) yl1(2)], ...
           bandColor, 'EdgeColor','none', 'FaceAlpha', bandAlpha);

plot(ax1, tCenters, popRateSm, 'k-', 'LineWidth', 1.25);

title(sprintf('%s — condition %d (first trial, %d units)', ...
        expStr, cond, numel(unitIdx)));

xline(ax1, 0, '--k');
xlim(ax1, win);
ylim(ax1, yl1);
ylabel(ax1, 'Pop. rate (sp/s/unit)');
set(ax1, 'XTickLabel', []);   % hide x labels on top panel
box(ax1, 'off');

% ========== BOTTOM: raster ==========
ax2 = subplot(2,1,2); hold(ax2, 'on');

% draw band spanning full raster height (robust)
yl2 = [0 nUnits+1];
patch(ax2, [bandX(1) bandX(2) bandX(2) bandX(1)], ...
           [yl2(1) yl2(1) yl2(2) yl2(2)], ...
           bandColor, 'EdgeColor','none', 'FaceAlpha', bandAlpha);

for ii = 1:nUnits
    sp = spikesByUnit{ii};
    if isempty(sp), continue; end
    plot(ax2, sp, ii*ones(size(sp)), 'k.', 'MarkerSize', 4);
end

xline(ax2, 0, '--k');
xlim(ax2, win);
ylim(ax2, yl2);
xlabel(ax2, 'Time (s)');
ylabel(ax2, 'Unit');

expStr = strrep(string(expName), '_', ' ');

pbaspect(ax1, [1 0.3 1]);
pbaspect(ax2, [1 0.3 1]);     % raster panel only


% [spikesByUnit, unitIdx] = plotPopRasterFirstTrialByExp(eu, expName, cond1);
% cond2 = 38;                       % 1..80
% [spikesByUnit, unitIdx] = plotPopRasterFirstTrialByExp(eu, expName, cond2);

%% Trash

% --- choose unit and settings ---
iEuPlot = 1;                 % pick your unit
condIdx = 1:80;              % conditions
tWin = [-0.5 0.5];           % display window
baseWin = [-0.5 -0.1];       % baseline window

% --- time vector ---
t = eta{iEuPlot, condIdx(1)}.t(:)';

maskT = (t >= tWin(1)) & (t <= tWin(2));
maskBase = (t >= baseWin(1)) & (t <= baseWin(2));

tPlot = t(maskT);

% --- allocate ---
Z = nan(numel(condIdx), numel(tPlot));

for c = 1:numel(condIdx)
    e = eta{iEuPlot, condIdx(c)};
    if isempty(e) || isempty(e.X)
        continue
    end

    x = e.X(:)';

    % baseline statistics for this condition
    mu = mean(x(maskBase), 'omitnan');
    sd = std(x(maskBase), 0, 'omitnan');

    % avoid divide-by-zero
    if sd == 0 || isnan(sd)
        continue
    end

    % z-score and restrict to plotting window
    Z(c,:) = (x(maskT) - mu) ./ sd;
end

% --- plot ---
figure;
imagesc(tPlot, 1:numel(condIdx), Z);
set(gca, 'YDir','normal');
xlabel('Time (s)');
ylabel('Condition');
title(sprintf('Unit %d: z-scored response (baseline = [%.1f %.1f] s)', ...
    iEuPlot, baseWin(1), baseWin(2)));
colorbar;
xline(0,'--k');

% horizontal separators every 8 conditions
hold on;
for r = 8:8:numel(condIdx)
    yline(r + 0.5, 'k-', 'LineWidth', 0.75);
end
hold off;

% --- diverging colormap: blue < 0 < red ---
n = 256;
cmap = [ ...
    linspace(0, 1, n/2)', linspace(0, 1, n/2)', ones(n/2,1); ...
    ones(n/2,1), linspace(1, 0, n/2)', linspace(1, 0, n/2)' ];
colormap(cmap);

% symmetric color scaling (optional but recommended)
maxAbs = max(abs(Z(:)), [], 'omitnan');
caxis([-maxAbs maxAbs]);

%%

function [spikesByUnit, unitIdx] = plotPopRasterFirstTrialByExp(eu, expName, cond)

%   note nPulses=10, pulseWidth=0.010, ipi=0.990

    win = [-0.5 0.5];

    % Filter units by ExpName
    allNames = string(arrayfun(@(x) x.ExpName, eu, 'UniformOutput', false));
    idxExp = find(allNames == string(expName));
    assert(~isempty(idxExp), 'No units found with ExpName == "%s"', string(expName));

    spikesByUnit = {};
    unitIdx = [];

    for k = 1:numel(idxExp)
        iEu = idxExp(k);

        g = eu(iEu).groupTwoColorStimTrials({'wavelength','power','duration','location'}, ...
            selectBy=struct(power=[], duration=10e-3, location=[], wavelength=[]));

        if numel(g) < cond || isempty(g(cond).trials)
            continue
        end

        firstTrial = g(cond).trials(1);  % <-- Trial object or ID, whatever it is

        % Use the same event name you used elsewhere; if yours is 'stimtwocolor', swap it in
        rd = eu(iEu).getRasterData( ...
            'stim', win, ...
            trials=firstTrial, ...
            alignTo='start', ...
            shutterDelay=0, sort=false, ...
            photoelectricBlankDuration=0.5e-3);

        % Extract spikes (robust to common rd formats)
        sp = [];
        if isstruct(rd) && isfield(rd,'spikeTimes')
            if iscell(rd.spikeTimes) && ~isempty(rd.spikeTimes), sp = rd.spikeTimes{1};
            elseif isnumeric(rd.spikeTimes), sp = rd.spikeTimes(:)'; end
        elseif isstruct(rd) && isfield(rd,'t')
            if iscell(rd.t) && ~isempty(rd.t), sp = rd.t{1};
            elseif isnumeric(rd.t), sp = rd.t(:)'; end
        end

        unitIdx(end+1,1) = iEu; %#ok<AGROW>
        spikesByUnit{end+1,1} = sp; %#ok<AGROW>
    end

    fprintf('Selected %d units (ExpName=%s, cond=%d)\n', numel(unitIdx), string(expName), cond);
    assert(~isempty(unitIdx), 'No units had condition %d with at least 1 trial.', cond);


    win = [-0.5 0.5];

    bandX = [0.00 0.01];            % <--- start/end time (sec) of highlighted box
    bandColor = [0 0 1];            % <--- RGB (red). Examples: [0 0 1]=blue, [0 1 0]=green
    bandAlpha = 0.2; 
    nUnits = numel(spikesByUnit);

    % Plot
    figure; hold on;

    xlim(win);
    ylim([0 nUnits+1]);           % nUnits = numel(spikesByUnit)
    yl = ylim;                    % [ymin ymax]
    x0 = bandX(1);
    x1 = bandX(2);
    hBand = patch([x0 x1 x1 x0], [yl(1) yl(1) yl(2) yl(2)], bandColor, ...
        'EdgeColor','none', 'FaceAlpha', bandAlpha);
    uistack(hBand, 'bottom');

    for ii = 1:numel(spikesByUnit)
        sp = spikesByUnit{ii};
        if isempty(sp), continue; end
        plot(sp, ii*ones(size(sp)), 'k.', 'MarkerSize', 4);
    end
    xlim(win);
    ylim([0 numel(spikesByUnit)+1]);
    xlabel('Time (s)');
    ylabel('Unit');
    expStr = strrep(string(expName), '_', ' ');
    title(sprintf('%s — condition %d (first trial, %d units)', ...
        expStr, cond, numel(unitIdx)));
    pbaspect([1 0.5 1])
    xline(0,'--k');
    hold off;
end
