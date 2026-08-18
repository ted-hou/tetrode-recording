eu = EphysUnit.load('C:\SERVER\Units\SNr_CoChR_VGATCre\SingleUnit_NonDuplicate_NonDrift_SNr', waveforms=false, spikecounts=false, spikerates=false);


nTrials = arrayfun(@(eu) length(eu.EventTimes.TimeoutOn), eu);
eu(nTrials==0) = [];

[ac.ac, ac.euIndicesFirstInSession, ac.expIndices, ac.uniqueExpNames] = eu.loadArduinoConnection();
eu.alignTimestamps(["WAITFORTOUCH", "TIMEOUT_START"], ac={ac.ac, ac.euIndicesFirstInSession, ac.expIndices, ac.uniqueExpNames}, acRefEventName="TIMEOUT_END", euRefEventName="TimeoutOff");

% Make press/lick trials

for iEu = 1:length(eu)
    selPress = eu(iEu).EventTimes.PressOff - eu(iEu).EventTimes.PressOn > 0.002;
    eu(iEu).Trials.Press = Trial(eu(iEu).EventTimes.TimeoutOn, eu(iEu).EventTimes.Press(selPress), stopMode='first', exclude=eu(iEu).EventTimes.Lick); 
    eu(iEu).Trials.Lick = Trial(eu(iEu).EventTimes.TimeoutOn, eu(iEu).EventTimes.Lick, stopMode='first', exclude=eu(iEu).EventTimes.Press(selPress));
end


%% Make CompleteExperiment3
dlcResultsPath = 'C:\SERVER\DeepLabCut\Results\FourPawsJawTongueSpine_Emma';

exp = CompleteExperiment3(eu, cameras='lr', deeplabcutPath=dlcResultsPath);

exp.alignTimestamps(refEventNameArduino={'REWARD_ON'}, refEventNameEphys={'RewardOn'}, trialDurationTolerance=2);

clear results
results(length(exp)) = struct(name=[], varsL=[], hasTimestampsL=[], varsR=[], hasTimestampsR=[], isValid=[]);
for iExp = 1:length(exp)
    results(iExp).name = exp(iExp).name;
    results(iExp).varsR = string(exp(iExp).vtdR.Properties.VariableNames)';
    results(iExp).hasTimestampsR = ismember('Timestamp', exp(iExp).vtdR.Properties.VariableNames);

    results(iExp).varsL = string(exp(iExp).vtdL.Properties.VariableNames)';
    results(iExp).hasTimestampsL = ismember('Timestamp', exp(iExp).vtdL.Properties.VariableNames);

    results(iExp).isValid = length(results(iExp).varsL) == 23 && length(results(iExp).varsR) == 23;
end
clear iExp
eu = EphysUnit.load('C:\SERVER\Units\SNr_CoChR_VGATCre\SingleUnit_NonDuplicate_NonDrift_SNr', waveforms=false, spikecounts=false, spikerates=false);

%% Cull bad experiments (some are still pending DLC)
assert(isequal({exp.name}, {results.name}));
exp = exp([results.isValid]);
results = results([results.isValid]);
eu = [exp.eu];

% Rename variables
from = ["HandIpsiCam", "HandContraCam", "FootIpsiCam", "FootContraCam", "Jaw", "Tongue", "Spine"];
toL = ["HandL", "HandR", "FootL", "FootR", "Jaw", "Tongue", "Spine"];
toR = ["HandR", "HandL", "FootR", "FootL", "Jaw", "Tongue", "Spine"];
for iExp = 1:length(exp)
    if results(iExp).hasTimestampsL
        vtd = exp(iExp).vtdL;
        for i = 1:length(toL)
            if ismember(sprintf("%s_X", from(i)), results(iExp).varsL)
                vtd = renamevars(vtd, ...
                    [sprintf("%s_X", from(i)), sprintf("%s_Y", from(i)), sprintf("%s_Likelihood", from(i))], ...
                    [sprintf("%s_X", toL(i)), sprintf("%s_Y", toL(i)), sprintf("%s_Likelihood", toL(i))]);
            end
        end
        exp(iExp).vtdL = vtd;
    end
    if results(iExp).hasTimestampsR
        vtd = exp(iExp).vtdR;
        for i = 1:length(toR)
            if ismember(sprintf("%s_X", from(i)), results(iExp).varsR)
                vtd = renamevars(vtd, ...
                    [sprintf("%s_X", from(i)), sprintf("%s_Y", from(i)), sprintf("%s_Likelihood", from(i))], ...
                    [sprintf("%s_X", toR(i)), sprintf("%s_Y", toR(i)), sprintf("%s_Likelihood", toR(i))]);
            end
        end
        exp(iExp).vtdR = vtd;
    end
end

%% Reset results
clearvars -except eu exp results boot c eta p

expIndices = cellfun(@(name) find(strcmpi(name, {exp.name}), 1, 'first'), {eu.ExpName});

%%
for iEu = 1:length(eu)
    try
        selPress = eu(iEu).EventTimes.PressOff - eu(iEu).EventTimes.PressOn > 2e-3;
        eu(iEu).EventTimes.ValidPress = eu(iEu).EventTimes.Press(selPress);
        selLick = eu(iEu).EventTimes.LickOff - eu(iEu).EventTimes.LickOn > 2e-3;
        eu(iEu).EventTimes.ValidLick = eu(iEu).EventTimes.Lick(selLick);  
    
        if isfield(eu(iEu).EventTimes, 'TIMEOUT_START')
            eu(iEu).Trials.Press = Trial(eu(iEu).EventTimes.TIMEOUT_START, eu(iEu).EventTimes.ValidPress, stopMode='first', exclude=eu(iEu).EventTimes.ValidLick);
            eu(iEu).Trials.Lick = Trial(eu(iEu).EventTimes.TIMEOUT_START, eu(iEu).EventTimes.ValidLick, stopMode='first', exclude=eu(iEu).EventTimes.ValidPress);
        else
            eu(iEu).Trials.Press = Trial(eu(iEu).EventTimes.TimeoutOn, eu(iEu).EventTimes.ValidPress, stopMode='first', exclude=eu(iEu).EventTimes.ValidLick);
            eu(iEu).Trials.Lick = Trial(eu(iEu).EventTimes.TimeoutOn, eu(iEu).EventTimes.ValidLick, stopMode='first', exclude=eu(iEu).EventTimes.ValidPress);
        end
    
        eu(iEu).EventTimes.FirstPress = [eu(iEu).Trials.Press.Stop];
        eu(iEu).EventTimes.FirstLick = [eu(iEu).Trials.Lick.Stop];
    catch
        eu(iEu).EventTimes.FirstPress = [];
        eu(iEu).EventTimes.FirstLick = [];
        warning("Could not do unit %i from (exp %i)", iEu, expIndices(iEu))
    end
end


%% Boot movement responses
c.hasPress = false(1, length(eu));
c.hasLick = false(1, length(eu));
for iEu = 1:length(eu)
    c.hasPress(iEu) = nnz(eu(iEu).Trials.Press.duration()>=2) >= 10;
    c.hasLick(iEu) = nnz(eu(iEu).Trials.Lick.duration()>=2) >= 10;
end
p.bootAlpha = 0.01;
p.nboot = 100000;
p.responseWindowPress = [-0.3, 0];
p.responseWindowLick = [-0.3, 0];
assert(isequal(p.responseWindowPress, [-0.3, 0]))
assert(isequal(p.responseWindowLick, [-0.3, 0]))
boot.press = struct('h', NaN(length(eu), 1), 'muDiffCI', NaN(length(eu), 2), 'muDiffObs', NaN(length(eu), 1));
boot.lick = struct('h', NaN(length(eu), 1), 'muDiffCI', NaN(length(eu), 2), 'muDiffObs', NaN(length(eu), 1));
[boot.press.h, boot.press.muDiffCI, boot.press.muDiffObs] = bootstrapMoveResponse( ...
    eu, 'press', nboot=p.nboot, alpha=p.bootAlpha, withReplacement=false, oneSided=false, ...
    responseWindow=p.responseWindowPress);
[boot.lick.h, boot.lick.muDiffCI, boot.lick.muDiffObs] = bootstrapMoveResponse( ...
    eu, 'lick', nboot=p.nboot, alpha=p.bootAlpha, withReplacement=false, oneSided=false, ...
    responseWindow=p.responseWindowLick);
fprintf(1, '\nAll done\n')

%% Get whole session kinematics
p.features = ["HandR", "HandL", "FootR", "FootL", "Tongue", "Jaw", "Spine"];
p.featureStats = ["xPos", "xPos", "xPos", "xPos", "likelihood", "yPos", "yPos"]; % xPos, yPos, likelihood, speed
p.vtdNames = ["vtdR", "vtdL", "vtdR", "vtdL", "both", "both", "both"];
p.minL = [0.5, 0.5, 0.5, 0.5, 0.2, 0.5, 0.5];
p.spikeDataSource = "rate"; % rate, count
p.spikeRes = 0.001;
p.spikeKernelType = 'exponential';
switch p.spikeKernelType
    case 'gaussian'
        p.spikeKernelSigma = 0.075;
        p.spikeKernelWidth = 0.5;
        [~, ~, p.spikeKernel] = eu(1).getSpikeRates('gaussian', p.spikeKernelSigma, p.spikeRes, kernelWidth=p.spikeKernelWidth);
    case 'exponential'
        p.spikeKernelLambda1 = 10;
        p.spikeKernelLambda2 = 100;
        p.spikeKernelWidth = 0.5;
        [~, ~, p.spikeKernel] = eu(1).getSpikeRates('exponential', p.spikeKernelLambda1, p.spikeKernelLambda2, p.spikeRes, kernelWidth=p.spikeKernelWidth);
end

% Smoothing for kinematics
p.kineRes = 0.01;
p.kineKernelType = 'exponential';
switch p.kineKernelType
    case 'gaussian'
        p.kineKernelSigma = 0.075;
        p.kineKernelWidth = 0.5;
        [~, ~, p.kineKernel] = eu(1).getSpikeRates('gaussian', p.kineKernelSigma, p.kineRes, kernelWidth=p.kineKernelWidth);
    case 'exponential'
        p.kineKernelLambda1 = 10;
        p.kineKernelLambda2 = 100;
        p.kineKernelWidth = 0.5;
        [~, ~, p.kineKernel] = eu(1).getSpikeRates('exponential', p.kineKernelLambda1, p.kineKernelLambda2, p.kineRes, kernelWidth=p.kineKernelWidth);
end


tl = tiledlayout(figure(), 2, 1);
ax = nexttile(tl);
title(ax, 'spike rate kernel')
xlabel(ax, 'time (s)')
hold(ax, 'on')
switch p.spikeKernelType
    case 'gaussian'
        plot(ax, p.spikeKernel.t, p.spikeKernel.y, DisplayName=sprintf('\\sigma=%g', p.spikeKernelSigma));
    case 'exponential'
        plot(ax, p.spikeKernel.t, p.spikeKernel.y, DisplayName=sprintf('\\lambda_1=%g, \\lambda_2=%g', p.spikeKernelLambda1, p.spikeKernelLambda2));
end
legend(ax, Interpreter='tex')
ax = nexttile(tl);
title(ax, 'kinematics kernel')
xlabel(ax, 'time (s)')
hold(ax, 'on')
switch p.kineKernelType
    case 'gaussian'
        plot(ax, p.kineKernel.t, p.kineKernel.y, DisplayName=sprintf('\\sigma=%g', p.kineKernelSigma));
    case 'exponential'
        plot(ax, p.kineKernel.t, p.kineKernel.y, DisplayName=sprintf('\\lambda_1=%g, \\lambda_2=%g', p.kineKernelLambda1, p.kineKernelLambda2));
end
legend(ax, Interpreter='tex')
drawnow


% Process kinematics
clear kinematics
kinematics(length(exp)) = struct(HandR=[], HandL=[], FootR=[], FootL=[], Tongue=[]);
for iExp = 1:length(exp)
    for iFeature = 1:length(p.features)
        fn = p.features(iFeature);
        vn = p.vtdNames(iFeature);
        sn = p.featureStats(iFeature);
        if ismember(vn, ["vtdL", "vtdR"])
            vtd = exp(iExp).(vn);
            if isempty(vtd)
                continue
            end
        end
        % Tongue uses bilateral likelihood
        if fn == "Tongue"
            assert(sn == "likelihood");
            iSide = 0;
            t = 0:p.kineRes:max(exp(iExp).vtdL.Timestamp(end), exp(iExp).vtdR.Timestamp(end));
            L = NaN(length(t), 1);
            for side = ["vtdL", "vtdR"]
                iSide = iSide + 1;
                l = exp(iExp).(side).(sprintf("%s_Likelihood", fn));
                l(l<p.minL(iFeature)) = 0;
                l(l>p.minL(iFeature)) = 1;
                L(:, iSide) = interp1(exp(iExp).(side).Timestamp, l, t, 'previous');
            end
            L = sum(L, 2);
            L = single(L > 0);
            clear iSide side l
            % L = smoothdata(L, 'gaussian', p.smoothWindow(iFeature));
            selnan = isnan(L);
            L(selnan) = interp1(t(~selnan), L(~selnan), t(selnan), 'previous');
            L = conv(L, p.kineKernel.y, 'same');
            kinematics(iExp).(fn) = struct(X=single(L(:)), t=single(t(:)));
            clear t
        % Licks from digital events
        elseif fn == "Lick"
            assert(sn == "likelihood");
            L = exp(iExp).eu(1).EventTimes.LickOn;
            t = 0:p.kineRes:exp(iExp).vtdR.Timestamp(end);
            edges = [t - p.kineRes/2, t(end) + p.kineRes/2];
            L = histcounts(L, edges);
            % L = smoothdata(L, 'gaussian', p.smoothWindow(iFeature));
            L = conv(L, p.kineKernel.y, 'same');
            kinematics(iExp).(fn) = struct(X=single(L(:)), t=single(t(:)));
            clear edges t
        % Average displacement from both sides
        elseif vn == "both"
            iSide = 0;
            t = 0:p.kineRes:max(exp(iExp).vtdL.Timestamp(end), exp(iExp).vtdR.Timestamp(end));
            S = NaN(length(t), 1);
            for side = ["vtdL", "vtdR"]
                iSide = iSide + 1;
                x = exp(iExp).(side).(sprintf("%s_X", fn));
                y = exp(iExp).(side).(sprintf("%s_Y", fn));
                l = exp(iExp).(side).(sprintf("%s_Likelihood", fn));
                x(l<p.minL(iFeature)) = NaN;
                y(l<p.minL(iFeature)) = NaN;
                x = (x - mean(x, 'all', 'omitnan')) ./ std(x, 0, 'all', 'omitnan');
                y = (y - mean(y, 'all', 'omitnan')) ./ std(y, 0, 'all', 'omitnan');
                tt = exp(iExp).(side).Timestamp;

                switch sn
                    case "xPos"
                        if vn == "vtdR"
                            s = x;
                        else
                            s = -x;
                        end
                    case "yPos"
                        s = y;
                    case "displacement"
                        s = sqrt(x.^2 + y.^2);
                    otherwise
                        error("Unsupported stat '%s' for feature '%s'", sn, fn)
                end

                % l = smoothdata(l, 'gaussian', 7);
                S(:, iSide) = interp1(tt, s, t, 'previous');
            end
            clear iSide side x y l d
            s = mean(S, 2, 'omitnan');
            % S = smoothdata(S, 'gaussian', p.smoothWindow(iFeature));
            selnan = isnan(s);
            s(selnan) = interp1(t(~selnan), s(~selnan), t(selnan), 'previous');
            s = conv(s, p.kineKernel.y, 'same');
            kinematics(iExp).(fn) = struct(X=single(s(:)), t=single(t(:)));
            clear t
        % Other tracking points use position or speed
        elseif ismember(sprintf("%s_X", fn), vtd.Properties.VariableNames)
            x = vtd.(sprintf("%s_X", fn));
            y = vtd.(sprintf("%s_Y", fn));
            l = vtd.(sprintf("%s_Likelihood", fn));
            t = 0:p.kineRes:vtd.Timestamp(end);
            x(l<p.minL(iFeature)) = NaN;
            y(l<p.minL(iFeature)) = NaN;
            x = (x - mean(x, 'all', 'omitnan')) ./ std(x, 0, 'all', 'omitnan');
            y = (y - mean(y, 'all', 'omitnan')) ./ std(y, 0, 'all', 'omitnan');
            switch sn
                case "xPos"
                    if vn == "vtdR"
                        s = x;
                    else
                        s = -x;
                    end
                case "yPos"
                    s = y;
                case "displacement"
                    s = sqrt(x.^2 + y.^2);
                otherwise
                    error("Unsupported stat '%s' for feature '%s'", sn, fn)
            end
            s = interp1(vtd.Timestamp, s, t, 'previous');
            % s = smoothdata(s, 'gaussian', p.smoothWindow(iFeature));
            selnan = isnan(s);
            s(selnan) = interp1(t(~selnan), s(~selnan), t(selnan), 'previous');
            s = conv(s, p.kineKernel.y, 'same');
            kinematics(iExp).(fn) = struct(X=single(s(:)), t=single(t(:)));
        end
    end
end
clear iExp iFeature vn fn vtd L X Y S selnan


%% Report bootstraped movement response direction


% figure, histogram(boot.press.h)
c.isPressUp = boot.press.h' == 1 & c.hasPress;
c.isPressDown = boot.press.h' == -1 & c.hasPress;
c.isPressResponsive = (c.isPressUp | c.isPressDown);

% figure, histogram(boot.lick.h)
c.isLickUp = boot.lick.h' == 1 & c.hasLick;
c.isLickDown = boot.lick.h' == -1 & c.hasLick;
c.isLickResponsive = (c.isLickUp | c.isLickDown);

fprintf(1, ['%g good units, %g modulated for reach:\n' ...
    '\t%g (%.0f%%) are excited (p<%g);\n' ...
    '\t%g (%.0f%%) are inhibited (p<%g).\n'], ...
    nnz(c.hasPress), nnz(c.isPressResponsive), ...
    nnz(c.isPressUp), 100*nnz(c.isPressUp)/nnz(c.isPressResponsive), p.bootAlpha, ...
    nnz(c.isPressDown), 100*nnz(c.isPressDown)/nnz(c.isPressResponsive), p.bootAlpha);

fprintf(1, ['%g good units, %g modulated for lick:\n' ...
    '\t%g (%.0f%%) are excited (p<%g);\n' ...
    '\t%g (%.0f%%) are inhibited (p<%g).\n'], ...
    nnz(c.hasLick), nnz(c.isLickResponsive), ...
    nnz(c.isLickUp), 100*nnz(c.isLickUp)/nnz(c.isLickResponsive), p.bootAlpha, ...
    nnz(c.isLickDown), 100*nnz(c.isLickDown)/nnz(c.isLickResponsive), p.bootAlpha);

nTotal = nnz(c.isPressResponsive & c.isLickResponsive);
fprintf(1, ['%g good units, %g modulated for reach AND lick:\n' ...
    '\t%g (%.0f%%) are press-excited AND lick-excited;\n' ...
    '\t%g (%.0f%%) are press-inhibited AND lick-inhibited;\n' ...
    '\t%g (%.0f%%) are press-excited AND lick-inhibited;\n' ...
    '\t%g (%.0f%%) are press-inhibited AND lick-excited;\n'], ...
    nnz(c.hasPress & c.hasLick), nTotal, ...
    nnz(c.isPressUp & c.isLickUp), 100*nnz(c.isPressUp & c.isLickUp)/nTotal, ...
    nnz(c.isPressDown & c.isLickDown), 100*nnz(c.isPressDown & c.isLickDown)/nTotal, ...
    nnz(c.isPressUp & c.isLickDown), 100*nnz(c.isPressUp & c.isLickDown)/nTotal, ...
    nnz(c.isPressDown & c.isLickUp), 100*nnz(c.isPressDown & c.isLickUp)/nTotal)   

sel = c.hasPress & c.hasLick;
fprintf('05 Calculate: Of %i: %i (%i%%) showed modulation for BOTH, %i (%i%%) showed modulation for lick only, %i (%i%%) showed modulation for reach only, %i (%i%%) for neither.', ...
    nnz(sel), ...
    nnz(sel & c.isPressResponsive & c.isLickResponsive), round(nnz(sel & c.isPressResponsive & c.isLickResponsive)/nnz(sel)*100), ...
    nnz(sel & c.isLickResponsive & ~c.isPressResponsive), round(nnz(sel & c.isLickResponsive & ~c.isPressResponsive)/nnz(sel)*100), ...
    nnz(sel & c.isPressResponsive & ~c.isLickResponsive), round(nnz(sel & c.isPressResponsive & ~c.isLickResponsive)/nnz(sel)*100), ...
    nnz(sel & ~c.isPressResponsive & ~c.isLickResponsive), round(nnz(sel & ~c.isPressResponsive & ~c.isLickResponsive)/nnz(sel)*100) ...
    )

clear nTotal sel

save('C:\SERVER\Units\meta_SNr_CoChR_VGATCre_ValidVideos.mat', 'boot', 'c', 'eta', 'kinematics', 'p')
eu.save('C:\SERVER\Units\SNr_CoChR_VGATCre\SingleUnit_NonDuplicate_NonDrift_SNr_ValidVideos')

%%
rd = eu.getRasterData('stim', window=[-1, 1], photoelectricBlankDuration=0.5e-3, minTrialDuration=0.75, maxTrialDuration=4);
if ~exist('E:\Figures\SNr_CoChR_VGATCre', 'dir')
    mkdir('E:\Figures\SNr_CoChR_VGATCre')
end

fig = figure(Units='inches', Position=[1, 1, 7, 7]);
ax = axes(fig);
for iEu = 1:length(eu)
    cla(ax)
    try
        EphysUnit.plotRaster(ax, rd(iEu), timeUnit='ms', xlim=[-1000, 4000]);
        legend(ax, 'off')
        print(fig, fullfile('E:\Figures\SNr_CoChR_VGATCre', sprintf('%s.png', eu(iEu).getName())), '-dpng', '-r0')
    end
end


%% Plot PSTH aligned to stim onset
clear artifacts
artifacts(1) = struct(event='StimOn', length=0.5, lengthUnit='ms', direction='right');
artifacts(2) = struct(event='StimOff', length=0.5, lengthUnit='ms', direction='right');
eta.stim = eu.getETA('count', 'stim', [-1, 1], resolution=0.020, alignTo='start', normalize=[-0.5, 0], artifacts=artifacts);
meta.stim = mean(eta.stim.X(:, isin(eta.stim.t, [0.05, 0.2])), 2, 'omitnan');
ax = axes(figure); hold(ax, 'on')
plot(ax, eta.stim.t, eta.stim.X(meta.stim>=0, :), Color=[0.8, 0.2, 0.2, 0.1])
plot(ax, eta.stim.t, eta.stim.X(meta.stim<0, :), Color=[0.2, 0.2, 0.8, 0.1])
plot(ax, eta.stim.t, mean(eta.stim.X, 1, 'omitnan'), Color='k', LineWidth=2)

%% Plot PETH aligned to reach
close all
eta.press = eu.getETA('count', 'press', [-4, 1], resolution=0.2, alignTo='stop', normalize=[-4, -2], minTrialDuration=1);
eta.lick = eu.getETA('count', 'lick', [-4, 1], resolution=0.2, alignTo='stop', normalize=[-4, -2], minTrialDuration=1);
meta.press = mean(eta.press.X(:, isin(eta.press.t, [-0.4, 0.2])), 2, 'omitnan');
meta.lick = mean(eta.lick.X(:, isin(eta.lick.t, [-0.4, 0.2])), 2, 'omitnan');

%%
ax = axes(figure); hold(ax, 'on')
plot(ax, eta.stim.t, eta.stim.X(c.isPressUp, :), Color=[0.8, 0.2, 0.2, 0.1])
plot(ax, eta.stim.t, eta.stim.X(c.isPressDown, :), Color=[0.2, 0.2, 0.8, 0.1])
plot(ax, eta.stim.t, mean(eta.stim.X, 1, 'omitnan'), Color='k', LineWidth=2)
plot(ax, eta.stim.t, mean(eta.stim.X(c.isPressUp, :), 1, 'omitnan'), Color=[0.8, 0.2, 0.2], LineWidth=2)
plot(ax, eta.stim.t, mean(eta.stim.X(c.isPressDown, :), 1, 'omitnan'), Color=[0.2, 0.2, 0.8], LineWidth=2)

%% Plot press PETH
ax = axes(figure); hold(ax, 'on')
plot(ax, eta.press.t, eta.press.X(c.isPressUp, :), Color=[0.8, 0.2, 0.2, 0.1])
plot(ax, eta.press.t, eta.press.X(c.isPressDown, :), Color=[0.2, 0.2, 0.8, 0.1])
plot(ax, eta.press.t, mean(eta.press.X, 1, 'omitnan'), Color='k', LineWidth=2)
xline(ax, 0, 'k-')
title(ax, 'press')
%%
ax = axes(figure); hold(ax, 'on')
plot(ax, eta.press.t, eta.press.X(meta.stim>+0.50, :), Color=[0.8, 0.2, 0.2, 0.1])
plot(ax, eta.press.t, eta.press.X(meta.stim<-0.5, :), Color=[0.2, 0.2, 0.8, 0.1])
plot(ax, eta.press.t, mean(eta.press.X(meta.stim>0, :), 1, 'omitnan'), Color=[0.8, 0.2, 0.2], LineWidth=1)
plot(ax, eta.press.t, mean(eta.press.X(meta.stim<0, :), 1, 'omitnan'), Color=[0.2, 0.2, 0.8], LineWidth=1)
xline(ax, 0, 'k-')
title(ax, 'press')
%%
ax = axes(figure); hold(ax, 'on')
plot(ax, eta.lick.t, eta.lick.X(c.isLickUp, :), Color=[0.8, 0.2, 0.2, 0.1])
plot(ax, eta.lick.t, eta.lick.X(c.isLickDown, :), Color=[0.2, 0.2, 0.8, 0.1])
plot(ax, eta.lick.t, mean(eta.lick.X, 1, 'omitnan'), Color='k', LineWidth=2)
xline(ax, 0, 'k-')
title(ax, 'lick')

ax = axes(figure); hold(ax, 'on')
sel = c.isPressResponsive;
scatter(ax, meta.press(sel), meta.stim(sel))
xline(ax, 0, 'k--')
yline(ax, 0, 'k--')
mdl = fitlm(meta.press(sel), meta.stim(sel));
plot(ax, mdl)

xlabel(ax, 'press')
ylabel(ax, 'stim')


ax = axes(figure); hold(ax, 'on')
sel = c.isLickResponsive;
scatter(ax, meta.lick(sel), meta.stim(sel))
xline(ax, 0, 'k--')
yline(ax, 0, 'k--')
mdl = fitlm(meta.lick(sel), meta.stim(sel));
plot(ax, mdl)

xlabel(ax, 'lick')
ylabel(ax, 'stim')



ax = axes(figure); hold(ax, 'on')
scatter(ax, meta.press, meta.lick)
xline(ax, 0, 'k--')
yline(ax, 0, 'k--')
mdl = fitlm(meta.press, meta.lick);
plot(ax, mdl)

xlabel(ax, 'press')
ylabel(ax, 'lick')