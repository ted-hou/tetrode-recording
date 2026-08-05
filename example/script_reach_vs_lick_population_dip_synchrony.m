eu = EphysUnit.load('C:\SERVER\Units\ReachVsLick_1225');
load('C:\SERVER\Units\meta_ReachVsLick_1225_20260728.mat') % 'boot', 'c', 'eta'
%%
% clearvars -except boot c eta eu exp expIndices
%% Make pairwise(same-session) rasters of increase/decrease cells
expNames = string({eu.ExpName}');
[uniqueExpNames, ~, euToExpIndices] = unique(expNames);

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
    if (n.press.inc > 0 && n.press.dec > 0) || (n.lick.inc > 0 && n.lick.dec > 0)
        fprintf('iExp=%i, %i units; %i reach-inc, %i reach-dec; %i lick-inc, %i lick-dec\n', ...
            iExp, nnz(isInExp), ...
            nnz(isInExp & c.isPressUp(:)), nnz(isInExp & c.isPressDown(:)), ...
            nnz(isInExp & c.isLickUp(:)), nnz(isInExp & c.isLickDown(:)))
    else
        continue
    end
end

%% Cull spikes that are too adjacent
for i = 1:length(eu)
    isi = [Inf, diff(eu(i).SpikeTimes)];
    eu(i).SpikeTimes = eu(i).SpikeTimes(isi>1e-4);
end

%%
for iExp = 9%1:length(uniqueExpNames) %9
    close all
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
    if (n.press.inc > 0 && n.press.dec > 0) && (n.lick.inc > 0 && n.lick.dec > 0)
        fprintf('iExp=%i, %i units; %i reach-inc, %i reach-dec; %i lick-inc, %i lick-dec\n', ...
            iExp, nnz(isInExp), ...
            nnz(isInExp & c.isPressUp(:)), nnz(isInExp & c.isPressDown(:)), ...
            nnz(isInExp & c.isLickUp(:)), nnz(isInExp & c.isLickDown(:)))
    else
        continue
    end
    for trialType = ["press", "lick"]
        fig = figure(Units='inches', Position=[1,1,14,3], DefaultAxesFontSize=9);
        layout.h = [50, 100];
        tl = tiledlayout(fig, sum(layout.h), 1, TileSpacing='none', Padding='compact');
        ax = gobjects(3, 1);
        switch trialType
            case "press"
                implementName = "bar";
                trialTypeName = "reach";
            case "lick"
                implementName = "spout";
                trialTypeName = "lick";
        end
        xlabel(tl, sprintf('time to %s contact', implementName), FontSize=9)
        for iAx = 1:3
            ax(iAx) = nexttile(tl, 1+(iAx-1)*25, [25, 1]);
        end
        rd = struct(dec=[], inc=[]);
        nTrials0 = NaN;
        for dir = ["dec", "inc"]
            rd.(dir) = eu(ei.(trialType).(dir)).getRasterData(char(trialType), window=[-0, 1], sort=true, alignTo='stop', minTrialDuration=2, maxTrialDuration=Inf);
            for iUnit = 1:length(rd.(dir))
                nTrials = max(unique(rd.(dir)(iUnit).I)); assert(max(unique(rd.(dir)(iUnit).I)) == length(unique(rd.(dir)(iUnit).I)))
                if isnan(nTrials0)
                    nTrials0 = nTrials;
                else
                    assert(nTrials == nTrials0)
                end
            end
        end

        % Plot one figure per trial
        for iTrial = 1:nTrials
            title(tl, sprintf("%s - trial %i", trialTypeName, iTrial), FontSize=9)
            cla(ax, 'reset')
            clear rdAlt
            iDir = 0;
            colors = [0.2, 0.2, 0.8; 0.8, 0.2, 0.2];
            colororder(ax(1), [0.15, 0.15, 0.15; 0.15, 0.15, 0.15])
            for dir = ["dec", "inc"]
                iDir = iDir + 1;
                iAx = iDir + 1;
                rdAlt.(dir).t = arrayfun(@(rd) rd.t(rd.I==iTrial), rd.(dir), UniformOutput=false);
                for iUnit = 1:length(rdAlt.(dir).t)
                    rdAlt.(dir).I{iUnit} = repmat(iUnit, size(rdAlt.(dir).t{iUnit}));
                end
                rdAlt.(dir).t = cat(2, rdAlt.(dir).t{:});
                rdAlt.(dir).I = cat(2, rdAlt.(dir).I{:});

                edges = -10:0.001:5;
                st = rdAlt.(dir).t;
                stJittered = st + (rand(size(rdAlt.(dir).t))-0.5)*1e-3*100;
                sr = histcounts(st, edges);
                sr = sr./length(rd.(dir));
                srJittered = histcounts(stJittered, edges);
                srJittered = srJittered./length(rd.(dir));

                resCoarse = 0.05;
                edgesCoarse = -10:resCoarse:5;
                centersCoarse = 0.5*(edgesCoarse(1:end-1) + edgesCoarse(2:end));
                srCoarse = histcounts(st, edgesCoarse);
                srCoarse = srCoarse./length(rd.(dir))./resCoarse;

                % Smooth (optional)
                sr = smoothdata(sr, 'gaussian', 10);
                srJittered = smoothdata(srJittered, 'gaussian', 10);

                hold(ax(iAx), 'on'), hold(ax(1), 'on')
                scatter(ax(iAx), st.*1e3, rdAlt.(dir).I, 3, colors(iDir, :), 'filled', DisplayName=dir);
                % scatter(ax(iAx), stJittered.*1e3, rdAlt.(dir).I+0.1, 3, [0.2, 0.2, 0.2], 'filled', DisplayName=dir);
                yyaxis(ax(1), 'left')
                if dir == "inc"
                    histogram(ax(1), BinEdges=edges.*1e3, BinCounts=sr, EdgeColor=colors(iDir, :), EdgeAlpha=0.5, DisplayStyle='stairs')
                    histogram(ax(1), BinEdges=edges.*1e3, BinCounts=srJittered, EdgeColor=[0.2, 0.2, 0.2], EdgeAlpha=0.3, DisplayStyle='stairs')
                end
                yyaxis(ax(1), 'right')
                plot(ax(1), centersCoarse.*1e3, srCoarse, Color=[colors(iDir, :), 0.5], LineStyle=':')
                hold(ax(iAx), 'off'), hold(ax(1), 'off')
                xline(ax(iAx), 0, 'k--')
                ylim(ax(iAx), [0.5, length(rd.(dir))+0.5])
                yticks(ax(iAx), [])
                ax(iAx).InteractionOptions.LimitsDimensions = "x";
            end
            xline(ax(1), 0, 'k--')
            xticks(ax(1:2), [])
            ax(1).XAxis.TickLength = [0, 0];
            ax(2).XAxis.TickLength = [0, 0];
            ax(3).XAxis.TickLength = [0, 0];
            ax(1).YAxis(1).TickLength = [0, 0];
            ax(1).YAxis(2).TickLength = [0, 0];
            ax(2).YAxis.TickLength = [0, 0];
            ax(3).YAxis.TickLength = [0, 0];
            ax(1).XAxis.Visible = 'off';
            ax(2).XAxis.Visible = 'off';
            xlim(ax, [-1500, 500])

            % Resize axes based on number of dec vs. inc units
            hDec = round(100*length(rd.dec)/(length(rd.dec) + length(rd.inc)));
            hInc = 100 - hDec;
            ax(1).Layout.Tile = 1; ax(1).Layout.TileSpan = [layout.h(1), 1];
            ax(2).Layout.Tile = layout.h(1) + 1; ax(2).Layout.TileSpan = [hDec, 1];
            ax(3).Layout.Tile = layout.h(1) + 1 + hDec; ax(3).Layout.TileSpan = [hInc, 1];
            linkaxes(ax, "x")

            yyaxis(ax(1), 'left'), ylabel(ax(1), 'sp/cell'), yticks(ax(1), [0, 0.2]), ylim(ax(1), [-0.1, 0.3])
            yyaxis(ax(1), 'right'), ylabel(ax(1), 'sp/s'), yticks(ax(1), 'auto')
            ylabel(ax(2), 'dec')
            ylabel(ax(3), 'inc')

            if ~exist(fullfile("E:\Figures\synchrony", sprintf("exp%i", iExp)), 'dir')
                mkdir(fullfile("E:\Figures\synchrony", sprintf("exp%i", iExp)))
            end
            print(fig, fullfile("E:\Figures\synchrony", sprintf("exp%i", iExp), sprintf("exp%i_%s_trial_%i.jpg", iExp, trialTypeName, iTrial)), '-dpng', '-r0')
        end
    end
end


%% Increase vs. decrease population PETHs: are they always anti-correlated, even in ITI/pre-movement quiescent period?
close all
p.blank(1).event = "StimOn";
p.blank(1).window = [-1, 1];
p.blank(1).event = "FirstPress";
p.blank(1).window = [-2, 2];
p.blank(1).event = "FirstLick";
p.blank(1).window = [-2, 2];

p.spikeDataSource = "rate"; % rate, count
p.spikeRes = 0.001;
p.spikeKernelType = 'gaussian';
switch p.spikeKernelType
    case 'gaussian'
        p.spikeKernelSigma = 0.1;
        p.spikeKernelWidth = 0.5;
        [~, ~, p.spikeKernel] = eu(1).getSpikeRates('gaussian', p.spikeKernelSigma, p.spikeRes, kernelWidth=p.spikeKernelWidth);
    case 'exponential'
        p.spikeKernelLambda1 = 10;
        p.spikeKernelLambda2 = 100;
        p.spikeKernelWidth = 0.5;
        [~, ~, p.spikeKernel] = eu(1).getSpikeRates('exponential', p.spikeKernelLambda1, p.spikeKernelLambda2, p.spikeRes, kernelWidth=p.spikeKernelWidth);
end


% for trialType = ["press", "lick"]
for trialType = "press"
    for iExp = 9%1:length(exp)
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
            continue
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
                            resolution=res);
                    case "rate"
                        switch p.spikeKernelType
                            case 'gaussian'
                                [xx, t, kernel] = eu(iEu).getSpikeRates('gaussian', p.spikeKernelSigma, edges, kernelWidth=p.spikeKernelWidth);
                            case 'exponential'
                                [xx, t, kernel] = eu(iEu).getSpikeRates('exponential', p.spikeKernelLambda1, p.spikeKernelLambda2, edges, kernelWidth=p.spikeKernelWidth);
                        end
                        xBaseline = eu(iEu).getTrialAlignedData('rate', [-4, -2], char(trialType), alignTo='stop', allowedTrialDuration=[1, Inf], ...
                            resolution=res, kernel=kernel);
                    otherwise
                        error("Unknown p.spikeDataSource=%s", p.spikeDataSource)
                end
                xx = (xx - mean(xBaseline, 'all', 'omitnan')) ./ std(xBaseline, 0, 'all', 'omitnan');

                for iEvent = 1:length(p.blank)
                    tEvent = eu(iEu).EventTimes.(p.blank(iEvent).event);
                    windows = tEvent(:) + p.blank(iEvent).window;
                    for ii = 1:length(tEvent)
                        [a, b] = isin(t, windows(ii, :), true, true);
                        xx(a:b) = NaN;
                    end
                end
                clear iEvent tEvent windows ii a b

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

        [rho, pval] = corr(x.dec, x.inc, Rows='complete');
        sel = isfinite(x.dec) & isfinite(x.inc);
        [r, lags] = xcorr(x.dec(sel), x.inc(sel), 100/res, 'none');

        % rhoBoot = NaN(1, nBoot);
        % rBoot = NaN(size(r, 1), nBoot);
        % for iBoot = 1:nBoot
        %     rhoBoot(iBoot) = corr(XBoot.dec(:, iBoot), XBoot.inc(:, iBoot), Rows='complete');
        %     rBoot(:, iBoot) = xcorr(fillmissing(XBoot.dec(:, iBoot), 'linear'), fillmissing(XBoot.inc(:, iBoot), 'linear'), 100/res, 'none');
        % end
        % alpha = 0.05;
        % rhoBootCI = quantile(rhoBoot, [alpha/2, 1-alpha/2]);
        % rBootCI = quantile(rBoot, [alpha/2, 1-alpha/2], 2);
        
        %% Do a PCA on population activity
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
            vel.(fn).t = kinematics(iExp).(fn).t;
            vel.(fn).x = diff([NaN; kinematics(iExp).(fn).X(:)]) ./ diff([NaN; kinematics(iExp).(fn).t(:)]);
            % vel.(fn).st = vel.(fn).t(strfind(vel.(fn).x' >= threshold, [0, 1]) + 1); % spike time, duh
            [~, vel.(fn).st] = findpeaks(vel.(fn).x, vel.(fn).t, MinPeakHeight=1, MinPeakProminence=1);

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
        for iEvent = 1:length(p.blank)
            tEvent = eu(iEu).EventTimes.(p.blank(iEvent).event);
            windows = tEvent(:) + p.blank(iEvent).window;
            for ii = 1:length(tEvent)
                [a, b] = isin(stMerge, windows(ii, :), true, true);
                stMerge(a:b) = [];
            end
        end

        features = ["Jaw", "HandR", "HandL", "SpikeRate", "PCAScore", "PC1Angle"];
        featureDispName = ["jaw", "r.hand", "l.hand", "spike rate", "pca scores", "pc1 angle"];
        nPCs = 1;

        l.h = [1, 1, 1, 3, 3, 3];
        l.ch = cumsum([1, l.h]);
        fig = figure(Units='normalized', Position=[0.05, 0.1, 0.9, 0.7], DefaultAxesFontSize=9);
        tl = tiledlayout(fig, sum(l.h), 1, TileSpacing='none', Padding='tight');

        h = gobjects(3+2+2*nPCs, 1);
        iLine = 0;
        ax = gobjects(length(features), 1);
        for i = 1:length(features)
            ax(i) = nexttile(tl, l.ch(i), [l.h(i), 1]);
            hold(ax(i), 'on')
            fn = features(i);
            switch fn
                case {"Jaw", "HandR", "HandL"}
                    iLine = iLine + 1;
                    h(iLine) = plot(ax(i), vel.(fn).t, vel.(fn).x, 'k-', DisplayName=featureDispName(i), Clipping='off');
                    plot(ax(i), vel.(fn).st, 0, 'bo')
                    % yline(ax(i), threshold, 'r--')
                case "SpikeRate"
                    iLine = iLine + 1;
                    h(iLine) = plot(ax(i), t, x.dec, 'b', DisplayName=sprintf('dec (n=%i)', n.(trialType).dec), Clipping='off');
                    iLine = iLine + 1;
                    h(iLine) = plot(ax(i), t, x.inc, 'r', DisplayName=sprintf('inc (n=%i)', n.(trialType).inc), Clipping='off');
                case "PCAScore"
                    for iPC = 1:nPCs
                        iLine = iLine + 1;
                        h(iLine) = plot(ax(i), t, score(:, iPC), Color=[getColor(iPC, 3, 0.7, s=0.5, l=0.5), 0.3], DisplayName=sprintf('PC%i (%.0f%%)', iPC, explained(iPC)), Clipping='off');
                    end
                    yline(ax(i), 0, 'k--')
                    plot(ax(i), stMerge, 0, 'bo')
                    ylim(ax(i), [-10, 10])
                case "PC1Angle"
                    for iPC = 1:nPCs
                        iLine = iLine + 1;
                        theta = acos(((XMerge)*coeff(:, iPC)) ./ (vecnorm(XMerge, 2, 2)));
                        h(iLine) = plot(ax(i), t, theta, LineStyle='-', Color=[getColor(iPC, 3, 0.7, s=0.5, l=0.5), 0.3], DisplayName=sprintf('angle relative to PC%i', iPC), Clipping='off');
                        ylim(ax(i), [0, pi])
                        yticks(ax(i), [0, pi/2, pi])
                        yticklabels(ax(i), ["0", "0.5\pi", "\pi"])
                    end
            end
            ylabel(ax(i), featureDispName(i))
            ax(i).InteractionOptions.LimitsDimensions = "x";
            % ax(i).XAxis.Visible = "off";
        end

        % xticklabels(ax(1:length(features)-1), [])
        ylim(ax(1:4), [-15, 15]);
        yticks(ax(1:4), [-10, 0, 10])
        ylim(ax(5), [-5, 5]);
        yticks(ax(5), [-3, 0, 3])
        xlim(ax, [36, 50])
        grid(ax, 'on')
        linkaxes(ax, 'x');
        legend(h);
        xlabel(tl, 'time (s)')
        box(ax, 'off')
    end
end


%%
% for windowWidth = 0.1:0.1:1
for windowWidth = 0.1

    periMovementScores = NaN(1, length(stMerge));
    baselineScores = NaN(1, length(stMerge));
    periMovementThetas = NaN(1, length(stMerge));
    baselineThetas = NaN(1, length(stMerge));
    iter = 0;
    rng(42)
    for iEvent = 1:length(stMerge)
        % ll = 0;
        while true
            iter = iter + 1;
            tRand = rand(1)*t(end);
            isPerimove = false;
            for t0 = stMergeUnfiltered(:)'
                if isin(tRand, t0 + [-1, 1])
                    isPerimove = true;
                    break
                end
            end
            if ~isPerimove
                break
            end
            % fprintf(repmat('\b', [1, ll]))
            % ll = fprintf('%i\n', iter);            
        end
        [a, b] = isin(t, stMerge(iEvent) + [-windowWidth/2, windowWidth/2], true, true);
        [ar, br] = isin(t, tRand + [-windowWidth/2, windowWidth/2], true, true);
        if ~isempty(a)
            periMovementThetas(iEvent) = mean(theta(a:b), 'omitnan');
            periMovementScores(iEvent) = mean(score(a:b), 'omitnan');
        else
            warning('boo!')
        end
        if ~isempty(ar)
            baselineThetas(iEvent) = mean(theta(ar:br), 'omitnan');
            baselineScores(iEvent) = mean(score(ar:br), 'omitnan');
        else
            warning('boo ar!')
        end
    end
    fig = figure;
    tl = tiledlayout(fig, 2, 1);
    ax = nexttile(tl); hold(ax, 'on')
    histogram(periMovementScores, -5:0.5:15, DisplayStyle='stairs', EdgeColor='red', DisplayName='peri-movement')
    histogram(baselineScores, -5:0.5:15, DisplayStyle='stairs', EdgeColor='black', DisplayName='baseline')
    xlabel('score on PC1')
    
    ax = nexttile(tl); hold(ax, 'on')
    histogram(periMovementThetas, linspace(0, pi, 50), DisplayStyle='stairs', EdgeColor='red', DisplayName='peri-movement')
    histogram(baselineThetas, linspace(0, pi, 50), DisplayStyle='stairs', EdgeColor='black', DisplayName='baseline')
    xlabel('\theta')
    xticks([0, pi/2, pi])
    xticklabels(["0", "0.5\pi", "pi"])
    legend

    title(tl, sprintf("windowWidth=%g", windowWidth))
%%
    fig = figure();
    pax = polaraxes(); hold(pax, 'on')
    polarhistogram(periMovementThetas, linspace(0, pi, 50), DisplayStyle='stairs', EdgeColor='red', DisplayName='peri-movement')
    polarhistogram(baselineThetas, linspace(0, pi, 50), DisplayStyle='stairs', EdgeColor='black', DisplayName='baseline')

    fig = figure();
    ax = axes(fig); hold(ax, 'on');
    for iEvent = 1:length(stMerge)
        sel = isin(t, stMerge(iEvent) + [-0.05, 0.05]);
        % scatter3(mean(score(sel, 1)), mean(score(sel, 2)), mean(score(sel, 3)), 'k.', MarkerEdgeAlpha=0.5)
        plot3(score(sel, 1), score(sel, 2), score(sel, 3), Color=[0.2, 0.2, 0.2, 0.1])
        scatter3(score(sel, 1), score(sel, 2), score(sel, 3), 'ko', MarkerEdgeAlpha=0.1)
    end
    plot3([-10, 10], [0, 0], [0, 0], Color='r', LineWidth=5)
    plot3([0, 0], [-10, 10], [0, 0], Color='g', LineWidth=5)
    plot3([0, 0], [0, 0], [-10, 10], Color='b', LineWidth=5)
    xlabel('PC1')
    ylabel('PC2')
    axis(ax, 'equal')
%%
    fig = figure(Units='inches', Position=[1, 1, 6, 6], DefaultAxesFontSize=9);
    tl = tiledlayout(fig, 4, 1);
    ax = nexttile(tl); hold(ax, 'on')
    metaPerimove = NaN(length(stMerge), size(score, 2));
    for iEvent = 1:length(stMerge)
        sel = isin(t, stMerge(iEvent) + [-0.05, 0.05]);
        metaPerimove(iEvent, :) = mean(score(sel, :), 1, 'omitnan');
        plot(ax, 1:size(score, 2), metaPerimove(iEvent, :), Color=[0.8, 0.2, 0.2, 0.05], LineStyle='-')
    end
    plot(ax, 1:size(score, 2), mean(metaPerimove, 1, 'omitnan'), Color=[0.8, 0.2, 0.2, 1], Marker='o', LineStyle='-', LineWidth=1.5)
    ylim(ax, [-6, 6])
    title(ax, sprintf('peri-move (n=%i)', length(stMerge)))

    ax = nexttile(tl); hold(ax, 'on')
    metaBaseline = NaN(length(stMerge), size(score, 2));
    rng(42);
    for iEvent = 1:length(stMerge)
        while true
            tRand = rand(1)*t(end);
            isPerimove = false;
            for t0 = stMergeUnfiltered(:)'
                if isin(tRand, t0 + [-1, 1])
                    isPerimove = true;
                    break
                end
            end
            if ~isPerimove
                break
            end   
        end
        sel = isin(t, tRand + [-0.05, 0.05]);
        metaBaseline(iEvent, :) = mean(score(sel, :), 1, 'omitnan');
        plot(ax, 1:size(score, 2), metaBaseline(iEvent, :), Color=[0.2, 0.2, 0.2, 0.05], LineStyle='-')
    end
    plot(ax, 1:size(score, 2), mean(metaBaseline, 1, 'omitnan'), Color=[0.2, 0.2, 0.2, 1], Marker='o', LineStyle='-', LineWidth=1.5)
    ylim(ax, [-6, 6])
    title(ax, sprintf('baseline (n=%i); no-movement within 2s window', length(stMerge)))

    ax = nexttile(tl); hold(ax, 'on')
    edges = -10:0.5:10;
    dispEdges = [edges(1:end-1); edges(2:end)];
    dispEdges = dispEdges(:)';
    h = gobjects(4, 1);
    for iPC = 1:size(score, 2)
        N = histcounts(metaPerimove(:, iPC), edges);
        N = N ./ max(N) .* 0.3;
        dispN = [N; N];
        dispN = dispN(:)';
        h(1) = patch(ax, iPC + dispN, dispEdges, [0.8, 0.2, 0.2], FaceAlpha=0.67, EdgeAlpha=0, DisplayName='peri-move');

        N = histcounts(metaBaseline(:, iPC), edges);
        N = N ./ max(N) .* 0.3;
        dispN = [N; N];
        dispN = dispN(:)';
        h(2) = patch(ax, iPC - dispN, dispEdges, [0.2, 0.2, 0.2], FaceAlpha=0.67, EdgeAlpha=0, DisplayName='baseline');
    end
    h(3) = scatter(ax, (1:size(score, 2))+0.15, mean(metaPerimove, 1, 'omitnan'), 15, [0.8, 0.2, 0.2], 'filled', Marker='o', MarkerEdgeAlpha=1, DisplayName='peri-move (mean)');
    h(4) = scatter(ax, (1:size(score, 2))-0.15, mean(metaBaseline, 1, 'omitnan'), 15, [0.2, 0.2, 0.2], 'filled', Marker='o', MarkerEdgeAlpha=1, DisplayName='baseline (mean)');
    ylim(ax, [-6, 6])
    xticks(ax, 1:size(score, 2))
    xlabel(ax, 'PC')
    ylabel(ax, 'PC scores')
    lgd = legend(ax, h, Location='northeast', AutoUpdate=false);
    lgd.ItemTokenSize = [9, 9];
    yline(ax, 0, 'k:')
    title(ax, 'peri-move vs. baseline')

    ax = nexttile(tl); hold(ax, 'on')
    plot(ax, 1:size(score, 2), mean(metaPerimove, 1, 'omitnan'), Color=[0.8, 0.2, 0.2, 1], Marker='o', MarkerSize=3, LineStyle='-', LineWidth=1.5, DisplayName='peri-move')
    plot(ax, 1:size(score, 2), mean(metaBaseline, 1, 'omitnan'), Color=[0.2, 0.2, 0.2, 1], Marker='o', MarkerSize=3, LineStyle='-', LineWidth=1.5, DisplayName='baseline')
    xticks(ax, 1:size(score, 2))
    ylim(ax, [-2, 2])
    xlabel(ax, 'PC')
    ylabel(ax, 'PC scores')
    legend(ax, Location='northeast', AutoUpdate=false)
    yline(ax, 0, 'k:')
    title(ax, 'peri-move vs. baseline (average across trials)')

end


%% Find a good place for example (has hand, no jaw)
tGood = [];
for t0 = vel.HandR.st(:)'
    if any(isin(vel.Jaw.st(:)', t0+[-0.5, 0.5]))
        continue
    end
    tGood = [tGood, t0];
end

for t0 = tGood(randperm(length(tGood), length(tGood)))
    xlim(gca, t0+[-15, 15])
    disp(t0) % Put a breakpoint here and keep on scrolling
end

%%
close all
% features = ["Jaw", "HandR", "HandL", "SpikeRate"];
% featureDispName = ["jaw", "r.hand", "l.hand", "spike rate"];
% lineStyles = ["-", "-.", ":", "-"];
% l.h = [1, 1, 1, 3];
features = ["Jaw", "HandR", "SpikeRate"];
featureDispName = ["jaw\nvelocity", "r.hand\nvelocity", "spike rate (a.u.)"];
lineStyles = ["-", "-", "-"];
l.h = [3, 3, 10];
l.ch = cumsum([1, l.h]);
fig = figure(Units='inches', Position=[1, 1, 4, 2.5], DefaultAxesFontSize=8);
tl = tiledlayout(fig, sum(l.h), 1, TileSpacing='none', Padding='tight');
t0 = 4244;
window = [-20, 15];

h = gobjects(length(features)+1, 1);
iLine = 0;
ax = gobjects(length(features), 1);
for i = 1:length(features)
    ax(i) = nexttile(tl, l.ch(i), [l.h(i), 1]);
    hold(ax(i), 'on')
    fn = features(i);
    switch fn
        case {"Jaw", "HandR", "HandL"}
            iLine = iLine + 1;
            h(iLine) = plot(ax(i), vel.(fn).t, vel.(fn).x, Color='k', LineWidth=1, LineStyle=lineStyles(i), DisplayName=featureDispName(i));
        case "SpikeRate"
            iLine = iLine + 1;
            h(iLine) = plot(ax(i), t, x.dec, 'b', LineWidth=1, DisplayName=sprintf('dec (n=%i)', n.(trialType).dec));
            iLine = iLine + 1;
            h(iLine) = plot(ax(i), t, x.inc, 'r', LineWidth=1, DisplayName=sprintf('inc (n=%i)', n.(trialType).inc));
            % yline(ax(i), 0, 'k--')
    end
    ylabel(ax(i), strsplit(featureDispName(i), "\\n"))
    ax(i).InteractionOptions.LimitsDimensions = "x";
end

plot(ax(length(features)), t0+window(2)+[-2, 0]-1, [-3.5, -3.5], 'k-', LineWidth=2)
text(ax(length(features)), t0+window(2)-1-1, -3.5, "1 s", VerticalAlignment="bottom", HorizontalAlignment='center', Clipping='off', FontSize=8, FontName='Arial')
for i = 1:length(features)-1
    ax(i).XAxis.Visible = "off";
end
ylim(ax(1:length(features)-1), [-10, 10]);
yticks(ax(1:length(features)-1), [])
ylim(ax(length(features)), [-4, 4]);
yticks(ax(length(features)), [-2, 0, 2])
xlim(ax, t0+window)
xticks(ax(length(features)), [])
% grid(ax(4), 'on')
linkaxes(ax, 'x');
lgd = legend(h([4, 3]), Location='northeast');
lgd.ItemTokenSize = [12, 8];
xlabel(tl, 'time')
box(ax, 'off')
fontsize(fig, 8, 'points')
fontname(fig, 'Arial')