
clear l
l.h = ones(1, length(features));
l.ch = cumsum([1, l.h]);
l.w = 10*cellfun(@diff, xl);
l.cw = cumsum([1, l.w]);
exportPath = "C:\SERVER\Figures\SNr_VGAT-Cre-CoChR";
if p.useWeakIncDec
    exportPath = fullfile(exportPath, "weakIncDec");
else
    exportPath = fullfile(exportPath, "bootIncDec");
end
if exist(exportPath, 'dir')
    rmdir(exportPath, 's')
end
mkdir(exportPath)

if ~exist('trialTypes', 'var')
    trialTypes = ["press", "lick"];
end
for iExp = selExp % nSessions+1 will plot session average
    x0 = 0;
    for trialType = trialTypes
        fig = figure(Units='normalized', Position=[x0, 0, 0.35, length(features)*0.08]);
        x0 = x0 + 0.5;
        tl = tiledlayout(fig, sum(l.h), sum(l.w), TileSpacing='compact');
        clear lgd
        lgd = struct(hist=gobjects(1, 2), line=gobjects(1, 2)); % One legend for each duration (col)
        AX = gobjects(length(features), 2);
    
        clear ax
        for ifn = 1:length(features)
            fn = char(features(ifn));
            % Plot STA
            for iDuration = 1:2
                if isnan(breakthroughThreshold)
                    breakthroughWindow = [0, p.pulseDurations(iDuration)];
                else
                    breakthroughWindow = [0, breakthroughThreshold];
                end
                ax = nexttile(tl, sum(l.w)*(ifn-1) + l.cw(iDuration), [l.h(ifn), l.w(iDuration)]);
                AX(ifn, iDuration) = ax;
                if ifn == 1
                    title(ax, sprintf('%gs opto', p.pulseDurations(iDuration)))
                end
                clear h
                iLine = 1;
                hold(ax, 'on')
    
                % Ctrl traces
                iPower = 0;
                selCtrl = stimData(iExp).trialTypeCtrl == trialType;
                clear fnx
                switch fn
                    case {'Jaw', 'HandL', 'HandR', 'JawInhibited', 'HandLInhibited', 'HandRInhibited', 'JawBreakthrough', 'HandLBreakthrough', 'HandRBreakthrough'}
                        if contains(fn, {'Inhibited', 'Breakthrough'})
                            fnx = strsplit(fn, {'Inhibited', 'Breakthrough'});
                            fnx = fnx{1};
                        else
                            fnx = fn;
                        end
                        t = stimData(iExp).(fnx).t;
                        mu = mean(stimData(iExp).(fnx).XCtrl(selCtrl, :), 1, 'omitnan');
                        err = 0.1*std(stimData(iExp).(fnx).XCtrl(selCtrl, :), 0, 1, 'omitnan');
                        ci = quantile(stimData(iExp).(fnx).XCtrl(selCtrl, :), [0.25, 0.75], 1);
                        h(iLine) = plot(ax, t, mu, Color=colorsByPower(iPower+1, :), LineWidth=1.5, DisplayName="ctrl");%sprintf("ctrl (n=%i)", nnz(selCtrl)));
                        iLine = iLine + 1;
                        patch(ax, [t, flip(t)], [ci(1, :), flip(ci(2, :))], colorsByPower(iPower+1, 1:3), FaceAlpha=0.05, EdgeAlpha=p.varianceEdgeAlhpa, LineStyle=':', LineWidth=0.5);
                    case {'psmh', 'psmhCumulative', 'psmhCumulativeInhibited', 'psmhCumulativeBreakthrough'}
                        edges = stimData(iExp).psmh.ctrl.(trialType).edges;
                        N = mean(stimData(iExp).psmh.ctrl.(trialType).N(selCtrl, :), 1, 'omitnan');
                        N(~isfinite(N)) = 0;
                        if ismember(fn, {'psmhCumulative', 'psmhCumulativeInhibited', 'psmhCumulativeBreakthrough'})
                            N = cumsum(N);
                        end
                        h(iLine) = histogram(ax, BinEdges=edges, BinCounts=N, DisplayStyle='bar', EdgeColor='none', FaceColor=colorsByPower(iPower+1, 1:3), FaceAlpha=0.5, DisplayName="ctrl");%sprintf("ctrl (n=%i)", nnz(selCtrl)));
                        iLine = iLine + 1;
                    case {'psmhFirst', 'psmhSecond', 'psmhFirstCumulative', 'psmhSecondCumulative', 'psmhFirstCumulativeInhibited', 'psmhSecondCumulativeInhibited', 'psmhFirstCumulativeBreakthrough', 'psmhSecondCumulativeBreakthrough'}
                        if contains(fn, {'Inhibited', 'Breakthrough'})
                            NMoves = stimData(iExp).psmh.ctrl.(trialType).N; % trials x t
                            centers = stimData(iExp).psmh.ctrl.(trialType).edges; centers = 0.5*(centers(1:end-1) + centers(2:end));
                        end
                        if contains(fn, 'Inhibited')
                            subsel = sum(NMoves(:, isin(centers, [-Inf, breakthroughWindow(2)])), 2) == 0;
                        elseif contains(fn, 'Breakthrough')
                            subsel = sum(NMoves(:, isin(centers, breakthroughWindow)), 2) > 0;
                        else
                            subsel = true(size(selCtrl));
                        end
                        edges = stimData(iExp).psmh.ctrl.(trialType).edges;
                        N = stimData(iExp).psmh.ctrl.(trialType).N(selCtrl & subsel, :);
                        N(~isfinite(N)) = 0;
                        N1 = zeros(size(N));
                        N2 = zeros(size(N));
                        assert(all(N<=1, 'all'))
                        for iTrial = 1:size(N, 1)
                            ii = find(N(iTrial, :)>0, 2, 'first');
                            if length(ii) >= 1
                                N1(iTrial, ii(1)) = 1;
                            end
                            if length(ii) >= 2
                                N2(iTrial, ii(2)) = 1;
                            end
                        end
                        if contains(fn, 'First')
                            N = mean(N1, 1, 'omitnan');
                        elseif contains(fn, 'Second')
                            N = mean(N2, 1, 'omitnan');
                        end
                        if contains(fn, 'Cumulative')
                            N = cumsum(N);
                        end
                        h(iLine) = histogram(ax, BinEdges=edges, BinCounts=N, DisplayStyle='bar', EdgeColor='none', FaceColor=colorsByPower(iPower+1, 1:3), FaceAlpha=0.5, DisplayName="ctrl (1st)");%sprintf("ctrl (n=%i)", nnz(selCtrl)));
                        iLine = iLine + 1;
                    case {'X', 'XInc', 'XDec', 'XFlt', 'XInhibited', 'XBreakthrough'}
                        t = stimData(iExp).psth.ctrl.t;
                        switch fn
                            case {'X', 'XInhibited', 'XBreakthrough'}
                                fnx = 'X';
                            case {'XInc', 'XDec', 'XFlt'}
                                fnx = char(string(fn) + trialType);
                                fnx(5) = upper(fnx(5)); % XIncPress
                        end
                        mumu = mean(stimData(iExp).psth.ctrl.(fnx)(selCtrl, :, :), [1, 3], 'omitnan');
                        if showIndividualTracesByPower(1)
                            switch p.showVarianceFor
                                case "trials"
                                    mu = mean(stimData(iExp).psth.ctrl.(fnx)(selCtrl, :, :), 3, 'omitnan'); % avg across units: trials x time x units, then average across 3rd dim -> trials x time
                                case "units"
                                    mu = mean(permute(stimData(iExp).psth.ctrl.(fnx)(selCtrl, :, :), [3, 2, 1]), 3, 'omitnan'); % avg across trials: units x time x trials, then average across 3rd dim -> units x time
                            end
                            if ~isempty(mu)
                                switch p.showVarianceAs
                                    case "ci"
                                        ci = quantile(mu, [0.25, 0.75], 1);
                                        patch(ax, [t, flip(t)], [ci(1, :), flip(ci(2, :))], colorsByPower(iPower+1, 1:3), FaceAlpha=0.05, EdgeColor=colorsByPower(iPower+1, 1:3), EdgeAlpha=p.varianceEdgeAlhpa, LineStyle=':', LineWidth=0.5);
                                    case "traces"
                                        plot(ax, t, mu', Color=[colorsByPower(iPower+1, 1:3), individualTracesAlpha], LineWidth=0.5);
                                    case "sd"
                                        sd = std(mu, 0, 1, 'omitnan');
                                        patch(ax, [t, flip(t)], [mumu-sd, flip(mumu+sd)], colorsByPower(iPower+1, 1:3), FaceAlpha=0.05, EdgeColor=colorsByPower(iPower+1, 1:3), EdgeAlpha=p.varianceEdgeAlhpa, LineStyle=':', LineWidth=0.5);
                                end
                            end
                        end

                        h(iLine) = plot(ax, t, mumu, Color=colorsByPower(iPower+1, :), LineWidth=1.5, DisplayName="ctrl");%sprintf("ctrl (n=%i)", nnz(selCtrl)));
                        iLine = iLine + 1;
                end
                % Stim traces
                [uniqueHash, ia] = unique(stimData(iExp).hash);
                durations = [0]; % For drawing xline at opto onset/offset
                for iHash = 1:length(uniqueHash)
                    sel = stimData(iExp).hash == uniqueHash(iHash) & stimData(iExp).trialType == trialType;
                    iPower = stimData(iExp).iPower(ia(iHash));
                    if stimData(iExp).iDuration(ia(iHash)) ~= iDuration
                        continue
                    end
                    durations = [durations, p.pulseDurations(iDuration)];
                    label = sprintf("%gmw", p.laserPowers(iPower)*1e3);
                    switch fn
                        case {'Jaw', 'HandL', 'HandR', 'JawInhibited', 'HandLInhibited', 'HandRInhibited', 'JawBreakthrough', 'HandLBreakthrough', 'HandRBreakthrough'}
                            if contains(fn, {'Inhibited', 'Breakthrough'})
                                fnx = strsplit(fn, {'Inhibited', 'Breakthrough'});
                                fnx = fnx{1};
                                NMoves = stimData(iExp).psmh.stim.(trialType).N; % trials x t
                                centers = stimData(iExp).psmh.stim.(trialType).edges; centers = 0.5*(centers(1:end-1) + centers(2:end));
                                xline(ax, breakthroughWindow(2), 'k:')
                            end
                            if contains(fn, 'Inhibited')
                                subsel = sum(NMoves(:, isin(centers, [-Inf, breakthroughWindow(2)])), 2) == 0;
                            elseif contains(fn, 'Breakthrough')
                                subsel = sum(NMoves(:, isin(centers, breakthroughWindow)), 2) > 0;
                            else
                                subsel = true(size(sel));
                            end
                            t = stimData(iExp).(fnx).t;
                            mu = mean(stimData(iExp).(fnx).X(sel & subsel, :), 1, 'omitnan');
                            err = 0.1*std(stimData(iExp).(fnx).X(sel & subsel, :), 0, 1, 'omitnan');
                            ci = quantile(stimData(iExp).(fnx).X(sel & subsel, :), [0.25, 0.75], 1);
                            h(iLine) = plot(ax, t, mu, Color=colorsByPower(iPower+1, :), LineStyle=lineStylesByPower(iPower+1), LineWidth=1.5, DisplayName=label);%sprintf("%s (n=%i)", label, nnz(sel & subsel)));
                            iLine = iLine + 1;
                            patch(ax, [t, flip(t)], [ci(1, :), flip(ci(2, :))], colorsByPower(iPower+1, 1:3), FaceAlpha=0.05, EdgeColor=colorsByPower(iPower+1, 1:3), EdgeAlpha=p.varianceEdgeAlhpa, LineStyle=':', LineWidth=0.5);
                        case {'psmh', 'psmhCumulative', 'psmhCumulativeInhibited', 'psmhCumulativeBreakthrough'}
                            if contains(fn, {'Inhibited', 'Breakthrough'})
                                NMoves = stimData(iExp).psmh.stim.(trialType).N; % trials x t
                                centers = stimData(iExp).psmh.stim.(trialType).edges; centers = 0.5*(centers(1:end-1) + centers(2:end));
                                xline(ax, breakthroughWindow(2), 'k:')
                            end
                            if contains(fn, 'Inhibited')
                                subsel = sum(NMoves(:, isin(centers, [-Inf, breakthroughWindow(2)])), 2) == 0;
                            elseif contains(fn, 'Breakthrough')
                                subsel = sum(NMoves(:, isin(centers, breakthroughWindow)), 2) > 0;
                            else
                                subsel = true(size(sel));
                            end
                            edges = stimData(iExp).psmh.stim.(trialType).edges;
                            N = mean(stimData(iExp).psmh.stim.(trialType).N(sel & subsel, :), 1, 'omitnan');
                            N(~isfinite(N)) = 0;
                            if ismember(fn, {'psmhCumulative', 'psmhCumulativeInhibited', 'psmhCumulativeBreakthrough'})
                                N = cumsum(N);
                            end
                            h(iLine) = histogram(ax, BinEdges=edges, BinCounts=N, DisplayStyle='stairs', EdgeColor=colorsByPower(iPower+1, 1:3), EdgeAlpha=colorsByPower(iPower+1, 4), LineWidth=1.5, DisplayName=label);%sprintf("%s (n=%i)", label, nnz(sel)));
                            iLine = iLine + 1;
                        case {'psmhFirst', 'psmhSecond', 'psmhFirstCumulative', 'psmhSecondCumulative', 'psmhFirstCumulativeInhibited', 'psmhSecondCumulativeInhibited', 'psmhFirstCumulativeBreakthrough', 'psmhSecondCumulativeBreakthrough'}
                            if contains(fn, {'Inhibited', 'Breakthrough'})
                                NMoves = stimData(iExp).psmh.stim.(trialType).N; % trials x t
                                centers = stimData(iExp).psmh.stim.(trialType).edges; centers = 0.5*(centers(1:end-1) + centers(2:end));
                                xline(ax, breakthroughWindow(2), 'k:')
                            end
                            if contains(fn, 'Inhibited')
                                subsel = sum(NMoves(:, isin(centers, [-Inf, breakthroughWindow(2)])), 2) == 0;
                            elseif contains(fn, 'Breakthrough')
                                subsel = sum(NMoves(:, isin(centers, breakthroughWindow)), 2) > 0;
                            else
                                subsel = true(size(sel));
                            end

                            edges = stimData(iExp).psmh.stim.(trialType).edges;
                            N = stimData(iExp).psmh.stim.(trialType).N(sel & subsel, :);
                            N(~isfinite(N)) = 0;
                            N1 = zeros(size(N));
                            N2 = zeros(size(N));
                            assert(all(N<=1, 'all'))
                            for iTrial = 1:size(N, 1)
                                ii = find(N(iTrial, :)>0, 2, 'first');
                                if length(ii) >= 1
                                    N1(iTrial, ii(1)) = 1;
                                end
                                if length(ii) >= 2
                                    N2(iTrial, ii(2)) = 1;
                                end
                            end
                            if contains(fn, 'First')
                                N = mean(N1, 1, 'omitnan');
                            elseif contains(fn, 'Second')
                                N = mean(N2, 1, 'omitnan');
                            end
                            if contains(fn, 'Cumulative')
                                N = cumsum(N);
                            end
                            h(iLine) = histogram(ax, BinEdges=edges, BinCounts=N, DisplayStyle='stairs', EdgeColor=colorsByPower(iPower+1, 1:3), EdgeAlpha=colorsByPower(iPower+1, 4), LineWidth=1.5, DisplayName=label);%sprintf("%s (n=%i)", label, nnz(sel)));
                            iLine = iLine + 1;
                        case {'X', 'XInc', 'XDec', 'XFlt', 'XInhibited', 'XBreakthrough'} % Spike rate
                            t = stimData(iExp).psth.stim.t;
                            switch fn
                                case {'X', 'XInhibited', 'XBreakthrough'}
                                    fnx = 'X';
                                case {'XInc', 'XDec', 'XFlt'}
                                    fnx = char(string(fn) + trialType);
                                    fnx(5) = upper(fnx(5)); % XIncPress
                            end
                            if contains(fn, {'Inhibited', 'Breakthrough'})
                                NMoves = stimData(iExp).psmh.stim.(trialType).N; % trials x t
                                centers = stimData(iExp).psmh.stim.(trialType).edges; centers = 0.5*(centers(1:end-1) + centers(2:end));
                                xline(ax, breakthroughWindow(2), 'k:')
                            end
                            if contains(fn, 'Inhibited')
                                subsel = sum(NMoves(:, isin(centers, [-Inf, breakthroughWindow(2)])), 2) == 0;
                            elseif contains(fn, 'Breakthrough')
                                subsel = sum(NMoves(:, isin(centers, breakthroughWindow)), 2) > 0;
                            else
                                subsel = true(size(sel));
                            end
                            mumu = mean(stimData(iExp).psth.stim.(fnx)(sel & subsel, :, :), [1, 3], 'omitnan');
                            if showIndividualTracesByPower(iPower+1)
                                switch p.showVarianceFor
                                    case "trials"
                                        mu = mean(stimData(iExp).psth.stim.(fnx)(sel & subsel, :, :), 3, 'omitnan'); % avg across units: trials x time x units, then average across 3rd dim -> trials x time
                                    case "units"
                                        mu = mean(permute(stimData(iExp).psth.stim.(fnx)(sel & subsel, :, :), [3, 2, 1]), 3, 'omitnan'); % avg across trials: units x time x trials, then average across 3rd dim -> units x time
                                end
                                if ~isempty(mu)
                                    switch p.showVarianceAs
                                        case "ci"
                                            ci = quantile(mu, [0.25, 0.75], 1);
                                            patch(ax, [t, flip(t)], [ci(1, :), flip(ci(2, :))], colorsByPower(iPower+1, 1:3), FaceAlpha=0.05, EdgeColor=colorsByPower(iPower+1, 1:3), EdgeAlpha=p.varianceEdgeAlhpa, LineStyle=':', LineWidth=0.5);
                                        case "traces"
                                            plot(ax, t, mu', Color=[colorsByPower(iPower+1, 1:3), individualTracesAlpha], LineWidth=0.5);
                                        case "sd"
                                            sd = std(mu, 0, 1, 'omitnan');
                                            patch(ax, [t, flip(t)], [mumu-sd, flip(mumu+sd)], colorsByPower(iPower+1, 1:3), FaceAlpha=0.05, EdgeColor=colorsByPower(iPower+1, 1:3), EdgeAlpha=p.varianceEdgeAlhpa, LineStyle=':', LineWidth=0.5);
                                    end
                                end
                            end
                            h(iLine) = plot(ax, t, mumu, Color=colorsByPower(iPower+1, :), LineWidth=1.5, DisplayName=label);%sprintf("%s (n=%i)", label, nnz(sel)));
                            iLine = iLine + 1;
                            % yline(ax, 0, 'k--')
                    end
                end
                xline(ax, durations, 'k--')
                if iDuration == 1
                    nTrialsPerExp = arrayfun(@(sd) length(sd.trialType), stimData(1:end-1));
                    expToFirstTrialIndex = cumsum([1, nTrialsPerExp(1:end-1)]);
                    if featureDispNames(ifn) == "psmh"
                        switch trialType
                            case "press"
                                ylabel(ax, ["bar-contact", featureUnits(ifn)], Rotation=0)
                            case "lick"
                                ylabel(ax, ["spout-contact", featureUnits(ifn)], Rotation=0)
                        end
                    elseif ismember(string(fn), ["XInc", "XDec", "XFlt"])
                        if iExp < length(stimData)
                            selExp = iExp;
                        else
                            selExp = 1:length(stimData)-1;
                        end
                        fnx = char(string(fn) + trialType);
                        fnx(5) = upper(fnx(5));
                        ylabel(ax, [featureDispNames(ifn), sprintf("n=%i %s", sum(arrayfun(@(sd) size(sd.psth.stim.(fnx), 3), stimData(selExp))), featureUnits(ifn))], Rotation=0)
                    else
                        ylabel(ax, [featureDispNames(ifn), featureUnits(ifn)], Rotation=0)
                    end
                end
                xlim(ax, xl{iDuration})
                ylim(ax, yl{ifn})
                if ifn == length(features)
                    % xticks(ax, xl{iDuration}(1)+1:xl{iDuration}(2)-1)
                    xticks(ax, unique(durations))
                else
                    xticks(ax, [])
                end
                if iDuration == 1
                    yticks(ax, 'auto')
                else
                    yticks(ax, [])
                end
                % switch fn
                %     case {'Jaw', 'HandL', 'HandR', 'JawInhibited', 'HandLInhibited', 'HandRInhibited', 'JawBreakthrough', 'HandLBreakthrough', 'HandRBreakthrough', 'X', 'XInc', 'XDec', 'XFlt'}
                %         lgd.line(iDuration) = legend(h, Location='eastoutside', Orientation='vertical'); lgd.ItemTokenSize = [9, 9];
                %     case {'psmh'}
                %         lgd.hist(iDuration) = legend(h, Location='eastoutside', Orientation='vertical'); lgd.ItemTokenSize = [9, 9];
                % end
            end
        end
        trialTypeDispName = trialType; 
        if trialType == "press"
            trialTypeDispName = "reach";
        end
        title(tl, sprintf('Exp %i - %s - %s (n=%i)', iExp, stimData(iExp).name, trialTypeDispName, nnz(stimData(iExp).trialType==trialType)), Interpreter='none')
        xlabel(tl, 'time to decoder/opto onset (s)')
        fontsize(fig, 10, 'points')
        % print(fig, fullfile(exportPath, sprintf("Exp %i - %s - %s.png", iExp, stimData(iExp).name, trialTypeDispName)), '-dpng')
        linkaxes(AX(:, 1), 'x')
        linkaxes(AX(:, 2), 'x')
    end
end
