% eu = EphysUnit.load('C:\SERVER\Units\ReachVsLick_1225');
% load('C:\SERVER\Units\meta_ReachVsLick_1225_20260728.mat') % 'boot', 'c', 'eta'

% Try running script_optrodeSNr_CoChR_VGATCre_20260729 again;
load('C:\SERVER\Units\meta_SNr_CoChR_VGATCre_ValidVideos.mat')
eu = EphysUnit.load('C:\SERVER\Units\SNr_CoChR_VGATCre\SingleUnit_NonDuplicate_NonDrift_SNr_ValidVideos');


%%
% clearvars -except boot c eta eu exp expIndices
%% Calculate ETA/META
clear artifacts
artifacts(1) = struct(event='StimOn', length=0.5, lengthUnit='ms', direction='right');
artifacts(2) = struct(event='StimOff', length=0.5, lengthUnit='ms', direction='right');
eta.stim = eu.getETA('count', 'stim', [-4, 4], resolution=0.05, alignTo='start', normalize=[-4, -2], artifacts=artifacts);
eta.press = eu.getETA('count', 'press', [-4, 4], resolution=0.1, alignTo='stop', normalize=[-4, -2], minTrialDuration=1, artifacts=artifacts);
eta.lick = eu.getETA('count', 'lick', [-4, 4], resolution=0.1, alignTo='stop', normalize=[-4, -2], minTrialDuration=1, artifacts=artifacts);

meta.stim = mean(eta.stim.X(:, isin(eta.stim.t, [0.05, 0.2])), 2, 'omitnan');
meta.press = mean(eta.press.X(:, isin(eta.press.t, [-0.3, 0])), 2, 'omitnan');
meta.lick = mean(eta.lick.X(:, isin(eta.lick.t, [-0.3, 0])), 2, 'omitnan');

%% Plot ETAs (per session)
close all
for iExp = 1:length(uniqueExpNames)% + 1
        if iExp > length(uniqueExpNames)
        selUnits = true(1, length(eu));
    else
        selUnits = reshape(euToExpIndices == iExp, size(c.isPressUp));
    end

    fig = figure(Units='inches', Position=[1, 1, 15, 10]);
    tl = tiledlayout(fig, 2, 2); 
    
    clear h
    iLine = 1;
    ax = nexttile(tl); hold(ax, 'on')
    plot(ax, eta.press.t, eta.press.X(selUnits & c.isPressUp, :), Color=[1, 0, 0, 0.6])
    if any(selUnits & c.isPressDown)
        plot(ax, eta.press.t, eta.press.X(selUnits & c.isPressDown, :), Color=[0, 0, 1, 0.6])
    end
    plot(ax, eta.press.t, eta.press.X(selUnits & ~c.isPressResponsive, :), Color=[0, 0, 0, 0.6])
    % h(iLine) = plot(ax, eta.press.t, mean(eta.press.X(selUnits & c.isPressUp, :), 1, 'omitnan'), Color=[0.8, 0.2, 0.2], LineWidth=1.5, DisplayName=sprintf('reach-inc (n=%i)', nnz(selUnits & c.isPressUp)));
    % iLine = iLine + 1;
    % if any(selUnits & c.isPressDown)
    %     h(iLine) = plot(ax, eta.press.t, mean(eta.press.X(selUnits & c.isPressDown, :), 1, 'omitnan'), Color=[0.2, 0.2, 0.8], LineWidth=1.5, DisplayName=sprintf('reach-dec (n=%i)', nnz(selUnits & c.isPressDown)));
    %     iLine = iLine + 1;
    % end
    % h(iLine) = plot(ax, eta.press.t, mean(eta.press.X(selUnits & ~c.isPressResponsive, :), 1, 'omitnan'), Color=[0.2, 0.2, 0.2], LineWidth=1.5, DisplayName=sprintf('reach-flat (n=%i)', nnz(selUnits & ~c.isPressResponsive)));
    % iLine = iLine + 1;
    xline(ax, 0, 'k:')
    % lgd = legend(h(isvalid(h)), Location='northwest');
    % lgd.ItemTokenSize = [9, 9];
    xlim(ax, [-3, 2])
    ylim(ax, [-2, 4])
    title(ax, 'reach')
    xlabel(ax, 'time to bar contact (s)')
    
    clear h
    iLine = 1;
    ax = nexttile(tl); hold(ax, 'on')
    plot(ax, eta.lick.t, eta.lick.X(selUnits & c.isLickUp, :), Color=[1, 0, 0, 0.6])
    if any(selUnits & c.isLickDown)
        plot(ax, eta.lick.t, eta.lick.X(selUnits & c.isLickDown, :), Color=[0, 0, 1, 0.6])
    end
    plot(ax, eta.lick.t, eta.lick.X(selUnits & ~c.isLickResponsive, :), Color=[0, 0, 0, 0.6])
    % h(iLine) = plot(ax, eta.lick.t, mean(eta.lick.X(selUnits & c.isLickUp, :), 1, 'omitnan'), Color=[0.8, 0.2, 0.2], LineWidth=1.5, DisplayName=sprintf('lick-inc (n=%i)', nnz(selUnits & c.isLickUp)));
    % iLine = iLine + 1;
    % if any(selUnits & c.isLickDown)
    %     h(iLine) = plot(ax, eta.lick.t, mean(eta.lick.X(selUnits & c.isLickDown, :), 1, 'omitnan'), Color=[0.2, 0.2, 0.8], LineWidth=1.5, DisplayName=sprintf('lick-dec (n=%i)', nnz(selUnits & c.isLickDown)));
    %     iLine = iLine + 1;
    % end
    % h(iLine) = plot(ax, eta.lick.t, mean(eta.lick.X(selUnits & ~c.isLickResponsive, :), 1, 'omitnan'), Color=[0.2, 0.2, 0.2], LineWidth=1.5, DisplayName=sprintf('lick-flat (n=%i)', nnz(selUnits & ~c.isLickResponsive)));
    % iLine = iLine + 1;
    xline(ax, 0, 'k:')
    % lgd = legend(h(isvalid(h)), Location='northwest');
    % lgd.ItemTokenSize = [9, 9];
    xlim(ax, [-3, 2])
    ylim(ax, [-2, 4])
    title(ax, 'lick')
    xlabel(ax, 'time to spout contact (s)')
    
    clear h
    iLine = 1;
    ax = nexttile(tl); hold(ax, 'on')
    plot(ax, eta.stim.t, eta.stim.X(selUnits & c.isPressUp, :), Color=[1, 0, 0, 0.6])
    if any(selUnits & c.isPressDown)
        plot(ax, eta.stim.t, eta.stim.X(selUnits & c.isPressDown, :), Color=[0, 0, 1, 0.6])
    end
    plot(ax, eta.stim.t, eta.stim.X(selUnits & ~c.isPressResponsive, :), Color=[0, 0, 0, 0.6])
    % h(iLine) = plot(ax, eta.stim.t, mean(eta.stim.X(selUnits & c.isPressUp, :), 1, 'omitnan'), Color=[0.8, 0.2, 0.2], LineWidth=1.5, DisplayName=sprintf('reach-inc (n=%i)', nnz(selUnits & c.isPressUp)));
    % iLine = iLine + 1;
    % if any(selUnits & c.isPressDown)
    %     h(iLine) = plot(ax, eta.stim.t, mean(eta.stim.X(selUnits & c.isPressDown, :), 1, 'omitnan'), Color=[0.2, 0.2, 0.8], LineWidth=1.5, DisplayName=sprintf('reach-dec (n=%i)', nnz(selUnits & c.isPressDown)));
    %     iLine = iLine + 1;
    % end
    % h(iLine) = plot(ax, eta.stim.t, mean(eta.stim.X(selUnits & ~c.isPressResponsive, :), 1, 'omitnan'), Color=[0.2, 0.2, 0.2], LineWidth=1.5, DisplayName=sprintf('reach-flat (n=%i)', nnz(selUnits & ~c.isPressResponsive)));
    % iLine = iLine + 1;
    % h(iLine) = plot(ax, eta.stim.t, mean(eta.stim.X(selUnits, :), 1, 'omitnan'), Color=[0.2, 0.2, 0.2], LineStyle='-.', LineWidth=2, DisplayName=sprintf('all (n=%i)', nnz(selUnits)));
    iLine = iLine + 1;
    xline(ax, [0, 1, 3], 'k:')
    % lgd = legend(h(isvalid(h)), Location='northeast');
    % lgd.ItemTokenSize = [18, 9];
    xlim(ax, [-1, 4])
    ylim(ax, [-2, 4])
    title(ax, 'opto')
    xlabel(ax, 'time to opto onset (s)')
    
    clear h
    iLine = 1;
    ax = nexttile(tl); hold(ax, 'on')
    plot(ax, eta.stim.t, eta.stim.X(selUnits & c.isLickUp, :), Color=[1, 0, 0, 0.6])
    if any(selUnits & c.isLickDown)
        plot(ax, eta.stim.t, eta.stim.X(selUnits & c.isLickDown, :), Color=[0, 0, 1, 0.6])
    end
    plot(ax, eta.stim.t, eta.stim.X(selUnits & ~c.isLickResponsive, :), Color=[0, 0, 0, 0.6])
    % h(iLine) = plot(ax, eta.stim.t, mean(eta.stim.X(selUnits & c.isLickUp, :), 1, 'omitnan'), Color=[0.8, 0.2, 0.2], LineWidth=1.5, DisplayName=sprintf('lick-inc (n=%i)', nnz(selUnits & c.isLickUp)));
    % iLine = iLine + 1;
    % if any(selUnits & c.isLickDown)
    %     h(iLine) = plot(ax, eta.stim.t, mean(eta.stim.X(selUnits & c.isLickDown, :), 1, 'omitnan'), Color=[0.2, 0.2, 0.8], LineWidth=1.5, DisplayName=sprintf('lick-dec (n=%i)', nnz(selUnits & c.isLickDown)));
    %     iLine = iLine + 1;
    % end
    % h(iLine) = plot(ax, eta.stim.t, mean(eta.stim.X(selUnits & ~c.isLickResponsive, :), 1, 'omitnan'), Color=[0.2, 0.2, 0.2], LineWidth=1.5, DisplayName=sprintf('lick-flat (n=%i)', nnz(selUnits & ~c.isLickResponsive)));
    iLine = iLine + 1;
    % h(iLine) = plot(ax, eta.stim.t, mean(eta.stim.X(selUnits, :), 1, 'omitnan'), Color=[0.2, 0.2, 0.2], LineStyle='-.', LineWidth=2, DisplayName=sprintf('all (n=%i)', nnz(selUnits)));
    iLine = iLine + 1;
    xline(ax, [0, 1, 3], 'k:')
    % lgd = legend(h(isvalid(h)), Location='northeast');
    % lgd.ItemTokenSize = [18, 9];
    xlim(ax, [-1, 4])
    ylim(ax, [-2, 4])
    title(ax, 'opto')
    xlabel(ax, 'time to opto onset (s)')
    
    ylabel(tl, 'z-scored spike rate (a.u.)')

    % ttl = sprintf("VGAT-Cre x SNr(AAV-flex-CoChR) %i %s", iExp, eu(expToEuIndices(iExp)).ExpName);
    ttl = sprintf("WT x SNr(AAV-syn-CoChR) %i %s", iExp, eu(expToEuIndices(iExp)).ExpName);
    title(tl, ttl, Interpreter='none');
    print(fig, fullfile("C:\Users\AssadLab\Pictures\SNrOpto", sprintf("%s.png", ttl)), '-dpng')
end


%% Plot scattered METAs
tlp = tiledlayout(figure, 1, 2); 
tl = gobjects(1, 2);
tl(1) = tiledlayout(tlp, 1, 1);
tl(2) = tiledlayout(tlp, 1, 1); tl(2).Layout.Tile = 2;

AX = gobjects(1, 4);

sz = 15;
ax = nexttile(tl(1)); hold(ax, 'on'); AX(1) = ax;
h(1) = scatter(ax, meta.press(c.isPressUp), meta.stim(c.isPressUp), sz, [.8,.2,.2], 'filled', DisplayName=sprintf('reach-inc (n=%i)', nnz(c.isPressUp)));
h(2) = scatter(ax, meta.press(c.isPressDown), meta.stim(c.isPressDown), sz, [.2,.2,.8], 'filled', DisplayName=sprintf('reach-dec (n=%i)', nnz(c.isPressDown)));
h(3) = scatter(ax, meta.press, meta.stim, sz, [.2,.2,.2], DisplayName=sprintf('all (n=%i)', length(eu)));
xline(ax, 0, 'k:')
yline(ax, 0, 'k:')
xlabel(ax, 'reach')
ylabel(ax, 'opto')
axis(ax, 'equal')
legend(h, Location='northeast')

ax = nexttile(tl(1), 'east'); hold(ax, 'on'); AX(2) = ax;
edges = -2.5:0.5:6.5;
histogram(ax, meta.stim(c.isPressUp), edges, FaceColor=[.8,.2,.2], EdgeColor=[.8,.2,.2], EdgeAlpha=0.5, Orientation='horizontal')
histogram(ax, meta.stim(c.isPressDown), edges, FaceColor=[.2,.2,.8], EdgeColor=[.2,.2,.8], EdgeAlpha=0.5, Orientation='horizontal')
histogram(ax, meta.stim, edges, DisplayStyle='stairs', EdgeColor=[.2,.2,.2], EdgeAlpha=0.5, Orientation='horizontal')

ax = nexttile(tl(2)); hold(ax, 'on'); AX(3) = ax;
h(1) = scatter(ax, meta.lick(c.isLickUp), meta.stim(c.isLickUp), sz, [.8,.2,.2], 'filled', DisplayName=sprintf('lick-inc (n=%i)', nnz(c.isLickUp)));
h(2) = scatter(ax, meta.lick(c.isLickDown), meta.stim(c.isLickDown), sz, [.2,.2,.8], 'filled', DisplayName=sprintf('lick-dec (n=%i)', nnz(c.isLickDown)));
h(3) = scatter(ax, meta.lick, meta.stim, sz, [.2,.2,.2], DisplayName=sprintf('all (n=%i)', length(eu)));
xline(ax, 0, 'k:')
yline(ax, 0, 'k:')
xlabel(ax, 'lick')
ylabel(ax, 'opto')
axis(ax, 'equal')
legend(h, Location='northeast')
ax = nexttile(tl(2), 'east'); AX(4) = ax;

ax = nexttile(tl(2), 'east'); hold(ax, 'on'); AX(2) = ax;
edges = -2.5:0.5:6.5;
histogram(ax, meta.stim(c.isLickUp), edges, FaceColor=[.8,.2,.2], EdgeColor=[.8,.2,.2], EdgeAlpha=0.5, Orientation='horizontal')
histogram(ax, meta.stim(c.isLickDown), edges, FaceColor=[.2,.2,.8], EdgeColor=[.2,.2,.8], EdgeAlpha=0.5, Orientation='horizontal')
histogram(ax, meta.stim, edges, DisplayStyle='stairs', EdgeColor=[.2,.2,.2], EdgeAlpha=0.5, Orientation='horizontal')

xlim(AX([1, 3]), [-3, 7])
ylim(AX, [-3, 7])

%% Count units
expNames = string({eu.ExpName}');
[uniqueExpNames, expToEuIndices, euToExpIndices] = unique(expNames);
clc
for iExp = 1:length(uniqueExpNames)
    clear n ei
    isInExp = euToExpIndices == iExp;
    n.press.inc = nnz(isInExp & c.isPressUp(:));
    n.press.dec = nnz(isInExp & c.isPressDown(:));
    n.lick.inc = nnz(isInExp & c.isLickUp(:));
    n.lick.dec = nnz(isInExp & c.isLickDown(:));
    ei.press.inc = find(isInExp & c.isPressUp(:));
    ei.press.dec = find(isInExp & c.isPressDown(:));
    ei.lick.inc = find(isInExp & c.isLickUp(:));
    ei.lick.dec = find(isInExp & c.isLickDown(:));
    fprintf('iExp=%i %s \t %i units \t reach: %i inc \t %i dec \t lick: %i inc \t %i dec\n', ...
        iExp, eu(expToEuIndices(iExp)).getName(), nnz(isInExp), ...
        nnz(isInExp & c.isPressUp(:)), nnz(isInExp & c.isPressDown(:)), ...
        nnz(isInExp & c.isLickUp(:)), nnz(isInExp & c.isLickDown(:)))
end


%% Cull spikes that are too adjacent
% for i = 1:length(eu)
%     isi = [Inf, diff(eu(i).SpikeTimes)];
%     eu(i).SpikeTimes = eu(i).SpikeTimes(isi>1e-4);
% end

%% Increase vs. decrease population PETHs: are they always anti-correlated, even in ITI/pre-movement quiescent period?
close all
p.plotIndividualUnits = false;
p.plotPopulationMeans = true;
p.kinematicDataSource = "vel"; % pos, vel
p.spikeDataSource = "rate"; % rate, count
p.spikeRes = 0.001;
p.spikeKernelType = 'gaussian';
p.correctSpikeRateDrift = true; % Subtract a smoothed baseline spike rate
if isfield(p, 'blank')
    p = rmfield(p, 'blank');
end
p.blank(1).event = "StimOn";
p.blank(1).window = [-0.5, 0.5]*1e-3;
p.blank(1).event = "StimOff";
p.blank(1).window = [-0.5, 0.5]*1e-3;
% p.blank(1).event = "FirstPress";
% p.blank(1).window = [-2, 2];
% p.blank(1).event = "FirstLick";
% p.blank(1).window = [-2, 2];
if isfield(p, 'artifacts')
    p = rmfield(p, 'artifacts');
end
p.artifacts(1) = struct(event='StimOn', length=0.5, lengthUnit='ms', direction='right');
p.artifacts(2) = struct(event='StimOff', length=0.5, lengthUnit='ms', direction='right');

switch p.spikeKernelType
    case 'gaussian'
        p.spikeKernelSigma = 0.025; % 0.1;
        p.spikeKernelWidth = 0.1; % 0.5;
        % p.spikeKernelSigma = 0.2; % 0.1;
        % p.spikeKernelWidth = 1; % 0.5;
        [~, ~, p.spikeKernel] = eu(1).getSpikeRates('gaussian', p.spikeKernelSigma, p.spikeRes, kernelWidth=p.spikeKernelWidth, artifacts=p.artifacts);
    case 'exponential'
        p.spikeKernelLambda1 = 10;
        p.spikeKernelLambda2 = 100;
        p.spikeKernelWidth = 0.5;
        [~, ~, p.spikeKernel] = eu(1).getSpikeRates('exponential', p.spikeKernelLambda1, p.spikeKernelLambda2, p.spikeRes, kernelWidth=p.spikeKernelWidth, artifacts=p.artifacts);
end


% for trialType = ["press", "lick"]
trialType = "press";
iExp = 5;
clear n ei
isInExp = euToExpIndices == iExp;
n.press.inc = nnz(isInExp & c.isPressUp(:));
n.press.dec = nnz(isInExp & c.isPressDown(:));
n.lick.inc = nnz(isInExp & c.isLickUp(:));
n.lick.dec = nnz(isInExp & c.isLickDown(:));
ei.press.inc = find(isInExp & c.isPressUp(:));
ei.press.dec = find(isInExp & c.isPressDown(:));
ei.lick.inc = find(isInExp & c.isLickUp(:));
ei.lick.dec = find(isInExp & c.isLickDown(:));
if n.(trialType).inc == 0 || n.(trialType).dec == 0
    error("Cannot plot exp %i for trialtype %s", iExp, trialType)
end
X = struct(dec=[], inc=[]);
x = X;
XBoot = struct(dec=[], inc=[]);
nBoot = 10;
for dir = ["dec", "inc"]
    res = p.spikeRes;
    maxT = max(arrayfun(@(eu) eu.SpikeTimes(end), eu(ei.(trialType).inc)));
    edges = 0:res:maxT;
    X.(dir) = NaN(length(ei.(trialType).(dir)), length(edges) - 1, 'single');
    XBoot.(dir) = NaN(length(ei.(trialType).(dir)), length(edges) - 1, nBoot, 'single');
    ll = 0;
    for i = 1:length(ei.(trialType).(dir))
        fprintf(repmat('\b', [1, ll]))
        ll = fprintf("%s unit=%i of %i\n", dir, i, length(ei.(trialType).(dir)));
        iEu = ei.(trialType).(dir)(i);

        % Calculate observed spike rates
        switch p.spikeDataSource
            case "count"
                [xx, t] = eu(iEu).getSpikeCounts(edges);
                xx = single(xx);
                xBaseline = eu(iEu).getTrialAlignedData('count', [-4, -2], char(trialType), alignTo='stop', allowedTrialDuration=[1, Inf], ...
                    resolution=res, artifacts=p.artifacts);
            case "rate"
                switch p.spikeKernelType
                    case 'gaussian'
                        [xx, t, kernel] = eu(iEu).getSpikeRates('gaussian', p.spikeKernelSigma, edges, kernelWidth=p.spikeKernelWidth, artifacts=p.artifacts);
                    case 'exponential'
                        [xx, t, kernel] = eu(iEu).getSpikeRates('exponential', p.spikeKernelLambda1, p.spikeKernelLambda2, edges, kernelWidth=p.spikeKernelWidth, artifacts=p.artifacts);
                end
                xBaseline = eu(iEu).getTrialAlignedData('rate', [-4, -2], char(trialType), alignTo='stop', allowedTrialDuration=[1, Inf], ...
                    resolution=res, kernel=kernel);
            otherwise
                error("Unknown p.spikeDataSource=%s", p.spikeDataSource)
        end
        xx = (xx - mean(xBaseline, 'all', 'omitnan')) ./ std(xBaseline, 0, 'all', 'omitnan');
        if p.correctSpikeRateDrift
            xx = xx - smoothdata(xx, 2, 'movmedian', 1000/res);
        end

        if isfield(p, 'blank')
            for iEvent = 1:length(p.blank)
                tEvent = eu(iEu).EventTimes.(p.blank(iEvent).event);
                windows = tEvent(:) + p.blank(iEvent).window;
                for ii = 1:length(tEvent)
                    [a, b] = isin(t, windows(ii, :), true, true);
                    xx(a:b) = NaN;
                end
            end
            clear iEvent tEvent windows ii a b
        end

        X.(dir)(i, :) = xx;

        % Calculate bootstrapped spike rates (shuffle them ISIs)
        % st0 = eu(iEu).SpikeTimes(1);
        % isi = diff(eu(iEu).SpikeTimes);
        % rng(42)
        % for iBoot = 1:nBoot
        %     % Shuffle ISIs to generate new spike train
        %     st = cumsum([st0, isi(randperm(length(isi), length(isi)))]);
        % 
        %     % Calculate shuffled spike rates
        %     switch p.spikeDataSource
        %         case "count"
        %             [x, t] = eu(iEu).getSpikeCounts(edges, spikeTimes=st);
        %             x = single(x);
        %             xBaseline = eu(iEu).getTrialAlignedData('count', [-4, -2], char(trialType), alignTo='stop', allowedTrialDuration=[1, Inf], ...
        %                 resolution=res, spikeTimes=st);
        %         case "rate"
        %             switch p.spikeKernelType
        %                 case 'gaussian'
        %                     [x, t, kernel] = eu(iEu).getSpikeRates('gaussian', p.spikeKernelSigma, edges, kernelWidth=p.spikeKernelWidth, spikeTimes=st);
        %                 case 'exponential'
        %                     [x, t, kernel] = eu(iEu).getSpikeRates('exponential', p.spikeKernelLambda1, p.spikeKernelLambda2, edges, kernelWidth=p.spikeKernelWidth, spikeTimes=st);
        %             end
        %             xBaseline = eu(iEu).getTrialAlignedData('rate', [-4, -2], char(trialType), alignTo='stop', allowedTrialDuration=[1, Inf], ...
        %                 resolution=res, kernel=kernel, spikeTimes=st);
        %         otherwise
        %             error("Unknown p.spikeDataSource=%s", p.spikeDataSource)
        %     end
        %     x = (x - mean(xBaseline, 'all', 'omitnan')) ./ std(xBaseline, 0, 'all', 'omitnan');
        %     XBoot.(dir)(i, :, iBoot) = x;
        % end
    end
    x.(dir) = mean(X.(dir), 1, 'omitnan')';
    % XBoot.(dir) = squeeze(mean(XBoot.(dir), 1, 'omitnan'));
end
X.dec = X.dec';
X.inc = X.inc';

% [rho, pval] = corr(x.dec, x.inc, Rows='complete');
% sel = isfinite(x.dec) & isfinite(x.inc);
% [r, lags] = xcorr(x.dec(sel), x.inc(sel), 100/res, 'none');

% rhoBoot = NaN(1, nBoot);
% rBoot = NaN(size(r, 1), nBoot);
% for iBoot = 1:nBoot
%     rhoBoot(iBoot) = corr(XBoot.dec(:, iBoot), XBoot.inc(:, iBoot), Rows='complete');
%     rBoot(:, iBoot) = xcorr(fillmissing(XBoot.dec(:, iBoot), 'linear'), fillmissing(XBoot.inc(:, iBoot), 'linear'), 100/res, 'none');
% end
% alpha = 0.05;
% rhoBootCI = quantile(rhoBoot, [alpha/2, 1-alpha/2]);
% rBootCI = quantile(rBoot, [alpha/2, 1-alpha/2], 2);

% Do a PCA on population activity
XMerge = double([X.dec, X.inc]);
XMerge = XMerge - mean(XMerge, 1, 'omitnan');
[coeff, score, ~, ~, explained, mu] = pca(XMerge);

% fig = figure(Units='inches', Position=[1,1,6,3], DefaultAxesFontSize=9);
% tl = tiledlayout(fig, 1, 2);
% 
% ax = nexttile(tl); hold(ax, 'on')
% scatter(ax, x.dec, x.inc, 2, 'k.')
% mdl = fitlm(x.dec, x.inc);
% lims = [min(min(x.dec), min(x.inc)), max(max(x.dec), max(x.inc))];
% plot(ax, lims', mdl.predict(lims'), 'r', LineWidth=2)
% xline(ax, 0, 'k:')
% yline(ax, 0, 'k:')
% % plot(ax, mdl);
% xlabel(ax, 'decrease')
% ylabel(ax, 'increase')
% title(ax, 'z-scored spike counts')
% axis(ax, 'equal')
% ax = nexttile(tl);
% plot(ax, lags.*res, r, 'k')
% % patch(ax, [lags, flip(lags)].*res, [rBootCI(:, 1)', flip(rBootCI(:, 2)')], [0.15, 0.15, 0.15], FaceAlpha=0.1, EdgeAlpha=0.1)
% xlabel(ax, 'lag applied to increase population spike counts (s)')
% ylabel(ax, 'xcorr')
% title(tl, sprintf("exp%i - %s (%i dec vs. %i inc)\nrho=%g ([%g, %g]), pval=%g", iExp, trialType, n.(trialType).dec, n.(trialType).inc, rho, rhoBootCI(1), rhoBootCI(2), pval))
% fprintf("Exp %i, %i dec vs. %i inc: rho=%g, pval=%g;\n", iExp, n.(trialType).dec, n.(trialType).inc, rho, pval);

threshold = 3;
exclusionWindow = [0, 1]; % NON-INCLUSIVE
features = ["Jaw", "HandR", "HandL"];
for ifn = 1:length(features)
    fn = features(ifn);
    switch p.kinematicDataSource
        case "pos"
            vel.(fn).t = kinematics(iExp).(fn).t;
            vel.(fn).x = kinematics(iExp).(fn).X(:);
        case "vel"
            vel.(fn).x = diff([NaN; kinematics(iExp).(fn).X(:)]) ./ diff([NaN; kinematics(iExp).(fn).t(:)]);
            vel.(fn).st = vel.(fn).t(strfind(vel.(fn).x' >= threshold, [0, 1]) + 1); % spike time, duh
        otherwise
            error("Unknown argument p.kinematicDataSource=%s", p.kinematicDataSource)
    end
    [~, vel.(fn).st] = findpeaks(vel.(fn).x, vel.(fn).t, MinPeakHeight=1, MinPeakProminence=1);
    if isfield(p, 'blank')
        for iEvent = 1:length(p.blank)
            tEvent = eu(iEu).EventTimes.(p.blank(iEvent).event);
            windows = tEvent(:) + p.blank(iEvent).window;
            for ii = 1:length(tEvent)
                [a, b] = isin(vel.(fn).t, windows(ii, :), true, true);
                vel.(fn).x(a:b) = NaN;
            end
        end
        clear iEvent tEvent windows ii a b
    end
end
stMerge = [vel.Jaw.st(:)', vel.HandR.st(:)', vel.HandL.st(:)'];
stMerge = sort(unique(stMerge), 'ascend');
stMergeUnfiltered = stMerge;
i = 1;
while i < length(stMerge)
    [a, b] = isin(stMerge, [stMerge(i), stMerge(i) + 0.5], false, true);
    if ~isempty(a)
        stMerge(a:b) = [];
    end
    i = i + 1;
end
if isfield(p, 'blank')
    for iEvent = 1:length(p.blank)
        tEvent = eu(iEu).EventTimes.(p.blank(iEvent).event);
        windows = tEvent(:) + p.blank(iEvent).window;
        for ii = 1:length(tEvent)
            [a, b] = isin(stMerge, windows(ii, :), true, true);
            stMerge(a:b) = [];
        end
    end
end

% features = ["Jaw", "HandR", "HandL", "SpikeRate", "PCAScore", "PC1Angle"];
% featureDispName = ["jaw", "r.hand", "l.hand", "spike rate", "pca scores", "pc1 angle"];
features = ["Jaw", "HandR", "HandL", "SpikeRateWithOpto", "SpikeRateDiffWithOpto"];
featureDispName = ["jaw", "r.hand", "l.hand", "spike rate", "divergence"];
nPCs = 1;

l.h = [1, 1, 1, 2, 1];
l.ch = cumsum([1, l.h]);
fig = figure(Units='inches', Position=[1, 1, 10, 6], DefaultAxesFontSize=9);
tl = tiledlayout(fig, sum(l.h), 1, TileSpacing='none', Padding='tight');

% h = gobjects(3+2, 1);
clear h
iLine = 0;
ax = gobjects(length(features), 1);
for i = 1:length(features)
    ax(i) = nexttile(tl, l.ch(i), [l.h(i), 1]);
    hold(ax(i), 'on')
    fn = features(i);
    switch fn
        case "Jaw"
            iLine = iLine + 1;
            h(iLine) = plot(ax(i), vel.(fn).t, vel.(fn).x, 'k-', DisplayName=featureDispName(i), Clipping='on');
            tLick = eu(expToEuIndices(iExp)).EventTimes.FirstLick;
            for iLick = 1:length(tLick)
                patch(ax(i), tLick(iLick) + [0, 0.1, 0.1, 0], [-5, -5, 5, 5], [0.2, 0.8, 0.2], FaceAlpha=0.1, EdgeColor=[0.2, 0.8, 0.2], EdgeAlpha=0.8)
            end
        case {"HandR", "HandL"}
            iLine = iLine + 1;
            h(iLine) = plot(ax(i), vel.(fn).t, vel.(fn).x, 'k-', DisplayName=featureDispName(i), Clipping='on');                    
            tPress = eu(expToEuIndices(iExp)).EventTimes.FirstPress;
            for iPress = 1:length(tPress)
                patch(ax(i), tPress(iPress) + [0, 0.1, 0.1, 0], [-5, -5, 5, 5], [0.8, 0.2, 0.2], FaceAlpha=0.1, EdgeColor=[0.8, 0.2, 0.2], EdgeAlpha=0.8)
            end                    
        case "SpikeRateWithOpto"
            if p.plotPopulationMeans
                iLine = iLine + 1;
                h(iLine) = plot(ax(i), t, x.dec, 'b', DisplayName=sprintf('dec (n=%i)', n.(trialType).dec), Clipping='on');
                iLine = iLine + 1;
                h(iLine) = plot(ax(i), t, x.inc, 'r', DisplayName=sprintf('inc (n=%i)', n.(trialType).inc), Clipping='on');
            end
            if p.plotIndividualUnits
                plot(ax(i), t, X.dec, Color=[.2, .2, .8, .2])
                plot(ax(i), t, X.inc, Color=[.8, .2, .2, .2])
            end
            yline(ax(i), 0, 'k:')
            tOn = eu(expToEuIndices(iExp)).EventTimes.LaserModBlueOn;
            tOff = eu(expToEuIndices(iExp)).EventTimes.LaserModBlueOff;
            for iStim = 1:length(tOn)
                patch(ax(i), [tOn(iStim), tOff(iStim), tOff(iStim), tOn(iStim)], [-5, -5, 5, 5], [0.2, 0.2, 0.8], FaceAlpha=0.1, EdgeColor=[0.2, 0.2, 0.8], EdgeAlpha=0.5)
            end
        case "SpikeRateDiffWithOpto"
            if p.plotPopulationMeans
                iLine = iLine + 1;
                h(iLine) = plot(ax(i), t, x.inc - x.dec, 'k', DisplayName=sprintf('inc-dec (n=%i,%i)', n.(trialType).inc, n.(trialType).dec), Clipping='on');
            end
            yline(ax(i), 0, 'k:')
            tOn = eu(expToEuIndices(iExp)).EventTimes.LaserModBlueOn;
            tOff = eu(expToEuIndices(iExp)).EventTimes.LaserModBlueOff;
            for iStim = 1:length(tOn)
                patch(ax(i), [tOn(iStim), tOff(iStim), tOff(iStim), tOn(iStim)], [-5, -5, 5, 5], [0.2, 0.2, 0.8], FaceAlpha=0.1, EdgeColor=[0.2, 0.2, 0.8], EdgeAlpha=0.5)
            end
        case "PCAScore"
            for iPC = 1:nPCs
                iLine = iLine + 1;
                h(iLine) = plot(ax(i), t, score(:, iPC), Color=[getColor(iPC, 3, 0.7, s=0.5, l=0.5), 0.3], DisplayName=sprintf('PC%i (%.0f%%)', iPC, explained(iPC)), Clipping='on');
            end
            yline(ax(i), 0, 'k--')
            ylim(ax(i), [-10, 10])
        case "PC1Angle"
            for iPC = 1:nPCs
                iLine = iLine + 1;
                theta = acos(((XMerge)*coeff(:, iPC)) ./ (vecnorm(XMerge, 2, 2)));
                h(iLine) = plot(ax(i), t, theta, LineStyle='-', Color=[getColor(iPC, 3, 0.7, s=0.5, l=0.5), 0.3], DisplayName=sprintf('angle relative to PC%i', iPC), Clipping='on');
                ylim(ax(i), [0, pi])
                yticks(ax(i), [0, pi/2, pi])
                yticklabels(ax(i), ["0", "0.5\pi", "\pi"])
            end
    end
    ylabel(ax(i), featureDispName(i))
    ax(i).InteractionOptions.LimitsDimensions = "x";
end
ylim(ax(1:3), [-5, 5]);
yticks(ax(1:3), [-3, 0, 3])
ylim(ax(4), [-5, 5]);
yticks(ax(4), [-3, 0, 3])
ylim(ax(5), [-5, 5]);
yticks(ax(5), [-3, 0, 3])
% xlim(ax, eu(expToEuIndices(iExp)).EventTimes.LaserModBlueOn(1) + [-30, 30])
linkaxes(ax, 'x');
% legend(h);
xlabel(tl, 'time (s)')
box(ax, 'off')

% xlim(ax, eu(expToEuIndices(iExp)).EventTimes.LaserModBlueOn(1) + [-10, 10])

% Thanks AI for drawing pretty buttons! Would've been nicer if they were functional
uicontrol(fig, Style='pushbutton', String='|<', Units='normalized', Position=[0.80, 0.02, 0.04, 0.03], ...
    Callback=@(src, ~) showStim(src, eu(expToEuIndices(iExp)), 1));
uicontrol(fig, Style='pushbutton', String='<', Units='normalized', Position=[0.85, 0.02, 0.04, 0.03], ...
    Callback=@(src, ~) showStim(src, eu(expToEuIndices(iExp)), fig.UserData.StimIndex - 1));
uicontrol(fig, Style='pushbutton', String='>', Units='normalized', Position=[0.90, 0.02, 0.04, 0.03], ...
    Callback=@(src, ~) showStim(src, eu(expToEuIndices(iExp)), fig.UserData.StimIndex + 1));
uicontrol(fig, Style='pushbutton', String='>|', Units='normalized', Position=[0.95, 0.02, 0.04, 0.03], ...
    Callback=@(src, ~) showStim(src, eu(expToEuIndices(iExp)), length(eu(expToEuIndices(iExp)).EventTimes.LaserModBlueOn)));

showStim(tl, eu(expToEuIndices(iExp)), 1)

function showStim(src, eu, index)
    fig = src.Parent;
    ax = fig.Children(5).Children;
    index = max(1, min(index, length(eu.EventTimes.LaserModBlueOn)));
    fig.UserData.StimIndex = index;
    xlim(ax, eu.EventTimes.LaserModBlueOn(index) + [-10, 10])
    drawnow
end
