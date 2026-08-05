eu = EphysUnit.load('C:\SERVER\Units\ReachVsLick_1225');
load('C:\SERVER\Units\meta_ReachVsLick_1225_20260728.mat') % 'boot', 'c', 'eta'
load('C:\SERVER\Units\population_dta_rta_25_5_200to800ms_10000boots_20260731.mat')


%% Make CompleteExperiment3
dlcResultsPath = 'C:\SERVER\DeepLabCut\Results\FourPawsTongueJawSpine';

exp = CompleteExperiment3(eu, cameras='lr', deeplabcutPath=dlcResultsPath);

exp.alignTimestamps(refEventNameArduino={'OPTO1_ON', 'OPTO2_ON', 'REWARD_ON'}, refEventNameEphys={'StimOn', 'RewardOn'}, trialDurationTolerance=2);

clear results
results(length(exp)) = struct(name=[], varsL=[], hasTimestampsL=[], varsR=[], hasTimestampsR=[], isValid=[]);
for iExp = 1:length(exp)
    results(iExp).name = exp(iExp).name;
    results(iExp).varsR = string(exp(iExp).vtdR.Properties.VariableNames)';
    results(iExp).hasTimestampsR = ismember('Timestamp', exp(iExp).vtdR.Properties.VariableNames);

    if ~isempty(exp(iExp).vtdL)
        results(iExp).varsL = string(exp(iExp).vtdL.Properties.VariableNames)';
        results(iExp).hasTimestampsL = ismember('Timestamp', exp(iExp).vtdL.Properties.VariableNames);
    else
        results(iExp).hasTimestampsL = false;
    end

    results(iExp).isValid = length(results(iExp).varsL) == 23 && length(results(iExp).varsR) == 23;
end
clear iExp

% Cull bad experiments (some are still pending DLC)
assert(isequal({exp.name}, {results.name}));
exp = exp([results.isValid]);
results = results([results.isValid]);
eu = [exp.eu];

% Rename variables
from = ["HandIpsiCam", "HandContraCam", "LegIpsiCam", "LegContraCam", "Jaw", "Tongue", "Spine"];
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

%% Detect dips in firing
p.features = ["HandR", "HandL", "FootR", "FootL", "Tongue", "Jaw", "Spine"];
p.featureStats = ["xPos", "xPos", "xPos", "xPos", "likelihood", "yPos", "yPos"]; % xPos, yPos, xVel, yVel, likelihood, displacement, speed
p.vtdNames = ["vtdR", "vtdL", "vtdR", "vtdL", "both", "both", "both"];
p.smoothWindow = [20, 20, 20, 20, 5, 20, 20];
p.minL = [0.5, 0.5, 0.5, 0.5, 0.2, 0.5, 0.5];
p.spikeDataSource = "rate"; % rate, count
p.spikeRes = 0.01;
p.spikeKernelType = 'gaussian';
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

ax = axes(figure());
title(ax, 'Spike rate kernel')
xlabel(ax, 'time (s)')
hold(ax, 'on')

switch p.spikeKernelType
    case 'gaussian'
        plot(ax, p.spikeKernel.t, p.spikeKernel.y, DisplayName=sprintf('\\sigma=%g', p.spikeKernelSigma));
    case 'exponential'
        plot(ax, p.spikeKernel.t, p.spikeKernel.y, DisplayName=sprintf('\\lambda_1=%g, \\lambda_2=%g', p.spikeKernelLambda1, p.spikeKernelLambda2));
end
legend(ax, Interpreter='tex')
drawnow

p.xta.res = 1/30;
p.xta.window = [-1, 1];
p.xta.meanWindow = [-0.3, 0.3];

p.xta.dip.samples = [0.2, 0.8]./p.spikeRes; % exceed threshold for 200-800ms
p.xta.dip.nullSamplesPre = 0.2/p.spikeRes; % 200 ms below threshold pre dip
p.xta.dip.nullSamplesPost = 0.2/p.spikeRes; % 200 ms below threshold post dip
p.xta.dip.thresholdQuantile = 0.25; % 0.25
p.xta.dip.thresholdSubQuantile = 0.05;
p.xta.dip.pattern = arrayfun(@(n) [zeros(1, p.xta.dip.nullSamplesPre), ones(1, n), zeros(1, p.xta.dip.nullSamplesPost)] , p.xta.dip.samples(1):p.xta.dip.samples(2), UniformOutput=false);
p.xta.dip.patternOnset = cellfun(@(pat) find(pat, 1, 'first') - 1, p.xta.dip.pattern); % finds the onset

p.xta.rise.samples = [0.2, 0.8]./p.spikeRes; % exceed threshold for 200-800ms
p.xta.rise.nullSamplesPre = 0.2/p.spikeRes; % 200 ms below threshold pre dip
p.xta.rise.nullSamplesPost = 0.2/p.spikeRes; % 200 ms below threshold post dip
p.xta.rise.thresholdQuantile = 1 - p.xta.dip.thresholdQuantile;
p.xta.rise.thresholdSubQuantile = 1 - p.xta.dip.thresholdSubQuantile;
p.xta.rise.pattern = arrayfun(@(n) [zeros(1, p.xta.rise.nullSamplesPre), ones(1, n), zeros(1, p.xta.rise.nullSamplesPost)] , p.xta.rise.samples(1):p.xta.rise.samples(2), UniformOutput=false);
p.xta.rise.patternOnset = cellfun(@(pat) find(pat, 1, 'first') - 1, p.xta.rise.pattern);

p.blank(1).event = "StimOn";
p.blank(1).window = [-1, 1];
p.blank(1).event = "FirstPress";
p.blank(1).window = [-2, 2];
p.blank(1).event = "FirstLick";
p.blank(1).window = [-2, 2];

p.nBoot = 10000;
p.bootAlpha = 0.05;

% Process vtdFeatues
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
            t = 0:p.xta.res:max(exp(iExp).vtdL.Timestamp(end), exp(iExp).vtdR.Timestamp(end));
            L = NaN(length(t), 1);
            for side = ["vtdL", "vtdR"]
                iSide = iSide + 1;
                l = exp(iExp).(side).(sprintf("%s_Likelihood", fn));
                l(l<p.minL(iFeature)) = 0;
                l(l>p.minL(iFeature)) = 1;
                % l = smoothdata(l, 'gaussian', 7);
                L(:, iSide) = interp1(exp(iExp).(side).Timestamp, l, t, 'linear');
            end
            L = sum(L, 2);
            L = single(L > 0);
            clear iSide side l
            L = smoothdata(L, 'gaussian', p.smoothWindow(iFeature));
            kinematics(iExp).(fn) = struct(X=single(L), t=single(t));
            clear t
        % Licks from digital events
        elseif fn == "Lick"
            assert(sn == "likelihood");
            L = exp(iExp).eu(1).EventTimes.LickOn;
            t = 0:p.xta.res:exp(iExp).vtdR.Timestamp(end);
            edges = [t - p.xta.res/2, t(end) + p.xta.res/2];
            L = histcounts(L, edges);
            L = smoothdata(L, 'gaussian', p.smoothWindow(iFeature));
            kinematics(iExp).(fn) = struct(X=single(L), t=single(t));
            clear edges t
        % Average displacement from both sides
        elseif vn == "both"
            iSide = 0;
            t = 0:p.xta.res:max(exp(iExp).vtdL.Timestamp(end), exp(iExp).vtdR.Timestamp(end));
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
                    case "xVel"
                        if vn == "vtdR"
                            s = [NaN; diff(smoothdata(x, 'gaussian', p.smoothWindow(iFeature)))]./[NaN; diff(tt)];
                        else
                            s = -[NaN; diff(smoothdata(x, 'gaussian', p.smoothWindow(iFeature)))]./[NaN; diff(tt)];
                        end
                    case "yVel"
                        s = [NaN; diff(smoothdata(y, 'gaussian', p.smoothWindow(iFeature)))]./[NaN; diff(tt)];
                    case "displacement"
                        s = sqrt(x.^2 + y.^2);
                    case "speed"
                        x = [NaN; diff(smoothdata(x, 'gaussian', p.smoothWindow(iFeature)))]./[NaN; diff(tt)];
                        y = [NaN; diff(smoothdata(y, 'gaussian', p.smoothWindow(iFeature)))]./[NaN; diff(tt)];
                        s = sqrt(x.^2 + y.^2);
                    otherwise
                        error("Unsupported stat '%s' for feature '%s'", sn, fn)
                end

                % l = smoothdata(l, 'gaussian', 7);
                S(:, iSide) = interp1(tt, s, t, 'linear');
            end
            clear iSide side x y l d
            S = mean(S, 2, 'omitnan');
            if ~ismember(sn, ["xVel", "yVel", "speed"])
                S = smoothdata(S, 'gaussian', p.smoothWindow(iFeature));
            end
            kinematics(iExp).(fn) = struct(X=single(S), t=single(t));
            clear t
        % Other tracking points use position or speed
        elseif ismember(sprintf("%s_X", fn), vtd.Properties.VariableNames)
            x = vtd.(sprintf("%s_X", fn));
            y = vtd.(sprintf("%s_Y", fn));
            l = vtd.(sprintf("%s_Likelihood", fn));
            t = vtd.Timestamp;
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
                case "xVel"
                    if vn == "vtdR"
                        s = [NaN; diff(smoothdata(x, 'gaussian', p.smoothWindow(iFeature)))]./[NaN; diff(t)];
                    else
                        s = -[NaN; diff(smoothdata(x, 'gaussian', p.smoothWindow(iFeature)))]./[NaN; diff(t)];
                    end
                case "yVel"
                    s = [NaN; diff(smoothdata(y, 'gaussian', p.smoothWindow(iFeature)))]./[NaN; diff(t)];
                case "displacement"
                    s = sqrt(x.^2 + y.^2);
                case "speed"
                    x = [NaN; diff(smoothdata(x, 'gaussian', p.smoothWindow(iFeature)))]./[NaN; diff(t)];
                    y = [NaN; diff(smoothdata(y, 'gaussian', p.smoothWindow(iFeature)))]./[NaN; diff(t)];
                    s = sqrt(x.^2 + y.^2);
                otherwise
                    error("Unsupported stat '%s' for feature '%s'", sn, fn)
            end
            % selnan = isnan(D);
            % D(selnan) = interp1(vtd.Timestamp(~selnan), D(~selnan), vtd.Timestamp(selnan), 'linear');
            if ~ismember(sn, ["xVel", "yVel", "speed"])
                s = smoothdata(s, 'gaussian', p.smoothWindow(iFeature));
            end
            kinematics(iExp).(fn) = struct(X=single(s), t=single(vtd.Timestamp));
        end
    end
end
clear iExp iFeature vn fn vtd L X Y S selnan

% Process dip/rise-triggered kinematics
tLocal = single(p.xta.window(1):p.xta.res:p.xta.window(2));
clear xta
xta.dip(length(exp)) = struct(iExp=[], params=[], t0=[], duration=[], spikerate=[]);
xta.rise(length(exp)) = struct(iExp=[], params=[], t0=[], duration=[], spikerate=[]);

lineLength = 0;
tTicTotal = tic();
rng(42)
if isempty(gcp('nocreate'))
    parpool('Processes');
end
for iExp = 1:length(exp)
    lastSpike = max(arrayfun(@(eu) eu.SpikeTimes(end), exp(iExp).eu));
    edges = 0:p.spikeRes:lastSpike;
    t = 0.5*(edges(1:end-1) + edges(2:end));
    X = NaN(length(exp(iExp).eu), length(t));
    for iEu = 1:length(exp(iExp).eu)
        switch p.spikeDataSource
            case "count"
                [x, ~] = exp(iExp).eu(iEu).getSpikeCounts(edges);
                x = double(x)./p.spikeRes;
            case "rate"
                switch p.spikeKernelType
                    case 'gaussian'
                        [x, ~] = exp(iExp).eu(iEu).getSpikeRates('gaussian', p.spikeKernelSigma, edges, kernelWidth=p.spikeKernelWidth);
                    case 'exponential'
                        [x, ~] = exp(iExp).eu(iEu).getSpikeRates('exponential', p.spikeKernelLambda1, p.spikeKernelLambda2, edges, kernelWidth=p.spikeKernelWidth);
                end
            otherwise
                error("Unknown p.spikeDataSource=%s", p.spikeDataSource)
        end
        x = (x-mean(x))/std(x, 0);
        X(iEu, :) = x;
    end
    x = mean(X, 1, 'omitnan');
    
    
    t0 = struct(dip=[], rise=[]);
    duration = struct(dip=[], rise=[]);
    nTotal = struct(dip=[], rise=[]);
    threshold = struct(dip=[], rise=[]);
    % Find onset of dips and rises
    for dir = ["dip", "rise"]
        switch dir
            case "dip"
                threshold.(dir) = quantile(x(x<0), p.xta.(dir).thresholdQuantile);
                i0 = arrayfun(@(i) strfind(x<=threshold.(dir), p.xta.(dir).pattern{i}) + p.xta.(dir).patternOnset(i), 1:length(p.xta.(dir).pattern), UniformOutput=false);
            case "rise"
                threshold.(dir) = quantile(x(x>0), p.xta.(dir).thresholdQuantile);
                i0 = arrayfun(@(i) strfind(x>=threshold.(dir), p.xta.(dir).pattern{i}) + p.xta.(dir).patternOnset(i), 1:length(p.xta.(dir).pattern), UniformOutput=false);
        end
        d = cell(size(i0));
        for iPat = 1:length(p.xta.(dir).pattern)
            d{iPat} = repmat(uint8(nnz(p.xta.(dir).pattern{iPat})), size(i0{iPat}));
        end
        i0 = cat(2, i0{:});
        t0.(dir) = t(i0);
        nTotal.(dir) = length(t0.(dir));
        duration.(dir) = cat(2, d{:}); % Multiply by p.spikeRes to get dip duration in seconds

        % Blank out dips around certain behavioral/stim events (opto, etc.)
        for iEvent = 1:length(p.blank)
            tEvent = exp(iExp).eu(1).EventTimes.(p.blank(iEvent).event);
            windows = tEvent(:) + p.blank(iEvent).window;
            for i = 1:length(tEvent)
                sel = t0.(dir)>=windows(i, 1) & t0.(dir)<=windows(i, 2);
                t0.(dir)(sel) = [];
                duration.(dir)(sel) = [];
                clear sel
            end
        end
        clear iEvent tEvent windows i i0
    end


    fprintf(repmat('\b', [1, lineLength]))
    lineLength = fprintf('Session %i/%i; %i(-%i) dips (x<%.2f), %i(-%i) rises (x>%.2f)... %.1fs elapsed...', iExp, length(exp), length(t0.dip), nTotal.dip-length(t0.dip), threshold.dip, length(t0.rise), nTotal.rise-length(t0.rise), threshold.rise, toc(tTicTotal));

    for dir = ["dip", "rise"]
        xta.(dir)(iExp).iExp = iExp;
        xta.(dir)(iExp).params = p;
        xta.(dir)(iExp).t0 = t0.(dir);
        xta.(dir)(iExp).duration = duration.(dir);

        t0Temp = t0.(dir);
        if ~isempty(t0Temp)
            T = tLocal + t0Temp';
            % Spike rates
            xta.(dir)(iExp).spikerate = struct(t=tLocal, X=NaN(length(t0Temp), length(tLocal)));
            for iDip = 1:length(t0Temp)
                xta.(dir)(iExp).spikerate.X(iDip, :) = interp1(t, x, T(iDip, :), 'linear');
            end
            % Subselect dips/rises by quantile
            if ~isnan(p.xta.(dir).thresholdSubQuantile)
                magnitude = mean(xta.(dir)(iExp).spikerate.X(:, xta.(dir)(iExp).spikerate.t > p.xta.meanWindow(1) & xta.(dir)(iExp).spikerate.t < p.xta.meanWindow(2)), 2, 'omitnan');
                switch dir
                    case "dip"
                        sel = magnitude <= quantile(magnitude, p.xta.(dir).thresholdSubQuantile);
                    case "rise"
                        sel = magnitude >= quantile(magnitude, p.xta.(dir).thresholdSubQuantile);
                end
                t0Temp = t0Temp(sel);
                T = T(sel, :);
                xta.(dir)(iExp).spikerate.X = xta.(dir)(iExp).spikerate.X(sel, :);
                xta.(dir)(iExp).t0 = t0Temp;
                lineLength = lineLength + fprintf('\tSubselecting %i/%i %ss...', nnz(sel), length(sel), dir);
            end

            % dip/rise-triggered kinematics
            for fn = p.features
                xta.(dir)(iExp).(fn) = struct(t=tLocal, X=NaN(length(t0Temp), length(tLocal), 'single'));
            end
            for fn = p.features
                if isempty(kinematics(iExp).(fn))
                    xta.(dir)(iExp).(fn) = [];
                    continue
                end
                for iDip = 1:length(t0Temp)
                    xta.(dir)(iExp).(fn).X(iDip, :) = interp1(kinematics(iExp).(fn).t, kinematics(iExp).(fn).X, T(iDip, :), 'linear'); % Consider doing shifts instead of interp1 for faster bootstraping
                end
            end

            % Bootstrap dip-triggered kinematics
            if p.nBoot > 0
                XBoot = NaN(p.nBoot, length(tLocal), length(p.features));
                IT = 1:length(tLocal);
                halfLength = (length(tLocal)-1)/2;
                assert(mod(halfLength, 1) == 0);
                IFEATURE = 1:length(p.features);
                maxT = kinematics(iExp).HandR.t(end);
                kinematicsTemp = kinematics(iExp);
                tTic = tic();
                parfor iBoot = 1:p.nBoot   
                % for iBoot = 1:p.nBoot   
                    t0Boot = rand([length(t0Temp), 1]) * maxT;
                    for iFeature = IFEATURE
                        fn = p.features(iFeature);
                        if isempty(kinematicsTemp.(fn))
                            continue
                        end
                        X = NaN(length(t0Boot), length(tLocal));
                        for iDip = 1:length(t0Boot)
                            [iStart, iStop] = isin(kinematicsTemp.(fn).t, t0Boot(iDip) + p.xta.window, true, true);
                            iMid = round((iStop + iStart)/2);
                            if isempty(iMid)
                                continue
                            elseif iMid-halfLength <= 0
                                iMid = halfLength + 1;
                            elseif iMid+halfLength >= length(kinematicsTemp.(fn).X)
                                iMid = length(kinematicsTemp.(fn).X) - halfLength;
                            end
                            X(iDip, :) = kinematicsTemp.(fn).X(iMid-halfLength:iMid+halfLength);
                        end
                        XBoot(iBoot, IT, iFeature) = mean(X, 1, 'omitnan');
                    end
                end
                for iFeature = 1:length(p.features)
                    fn = p.features(iFeature);
                    if isempty(kinematicsTemp.(fn))
                        continue
                    end
                    xta.(dir)(iExp).(fn).XBoot = XBoot(:, :, iFeature);
                end
                lineLength = lineLength + fprintf('%.1fs;', toc(tTic));
                clear iBoot t0Boot T X XBoot IT IFEATURE maxT kinematicsTemp tTic
            end
        end        
    end
    lineLength = lineLength + fprintf('\n');
end

clear iEu iExp lineLength tLocal tTicTotal x t mu sd threshold iRise t0Temp riseThreshold iRise tRise T magnitude sel riseMagnitude selRises
clear lineLength2 I2 iDip iFeature fn nDipsTotal nRisesTotal


%% Post-hoc do a bootstrap for xta.dip/xta.rise spikerate
lineLength = 0;
tTicTotal = tic();
tLocal = p.xta.window(1):p.xta.res:p.xta.window(2);
if p.nBoot > 0
    for iExp = 1:length(exp)
        fprintf(repmat('\b', [1, lineLength]))
        lineLength = fprintf('Session %i/%i... %.1fs;', iExp, length(exp), toc(tTicTotal));

        % Whole session spike rates
        lastSpike = max(arrayfun(@(eu) eu.SpikeTimes(end), exp(iExp).eu));
        edges = 0:p.spikeRes:lastSpike;
        t = 0.5*(edges(1:end-1) + edges(2:end));
        X = NaN(length(exp(iExp).eu), length(t));
        for iEu = 1:length(exp(iExp).eu)
            switch p.spikeDataSource
                case "count"
                    [x, ~] = exp(iExp).eu(iEu).getSpikeCounts(edges);
                    x = double(x)./p.spikeRes;
                case "rate"
                    switch p.spikeKernelType
                        case 'gaussian'
                            [x, ~] = exp(iExp).eu(iEu).getSpikeRates('gaussian', p.spikeKernelSigma, edges, kernelWidth=p.spikeKernelWidth);
                        case 'exponential'
                            [x, ~] = exp(iExp).eu(iEu).getSpikeRates('exponential', p.spikeKernelLambda1, p.spikeKernelLambda2, edges, kernelWidth=p.spikeKernelWidth);
                    end
                otherwise
                    error("Unknown p.spikeDataSource=%s", p.spikeDataSource)
            end
            x = (x-mean(x))/std(x, 0);
            X(iEu, :) = x;
        end
        x = mean(X, 1, 'omitnan');
    

        maxT = t(end);
        IT = 1:length(tLocal);
        XBoot = NaN(p.nBoot, length(tLocal));

        for dir = ["dip", "rise"]
            if ~isempty(xta.(dir)(iExp).t0)
                t0Temp = xta.(dir)(iExp).t0;
                lineLength = lineLength + fprintf(' %i %ss...', length(t0Temp), dir);
                tTic = tic();
                parfor iBoot = 1:p.nBoot   
                    t0Boot = rand([length(t0Temp), 1]) * maxT;
                    X = NaN(length(t0Boot), length(tLocal));
                    for i0 = 1:length(t0Boot)
                        [iStart, iStop] = isin(t, t0Boot(i0) + p.xta.window, true, true);
                        iMid = round((iStop + iStart)/2);
                        if isempty(iMid)
                            continue
                        elseif iMid-halfLength <= 0
                            iMid = halfLength + 1;
                        elseif iMid+halfLength >= length(x)
                            iMid = length(x) - halfLength;
                        end
                        X(i0, :) = x(iMid-halfLength:iMid+halfLength);
                    end
                    XBoot(iBoot, IT) = mean(X, 1, 'omitnan');
                end
                xta.(dir)(iExp).spikerate.XBoot = XBoot(:, :);
                lineLength = lineLength + fprintf('%.1fs;', toc(tTic));
            end
        end
    end
end
fprintf('\n')
clear tLocal iEu x t mu sd maxT IT XBoot t0Temp tRise tTic iBoot t0Boot T X i0 iRise lineLength tTic tTicTotal


%% Post-hoc calculate the bootstrpped 95%CI of the STD of the displacement traces
p.std.features = ["HandR", "HandL", "Jaw"];
p.std.window = [-0.3, 0.6];

for iExp = 1:length(exp)
    for dir = ["dip", "rise"]
        for fn = p.std.features
            if isempty(xta.(dir)(iExp).(fn))
                continue
            end
            selT = xta.(dir)(iExp).(fn).t >= p.std.window(1) & xta.(dir)(iExp).(fn).t <= p.std.window(2);
            xta.(dir)(iExp).(fn).stats.std = std(mean(xta.(dir)(iExp).(fn).X(:, selT), 1, 'omitnan'), 0, 2, 'omitnan');
            if p.nBoot > 0
                xta.(dir)(iExp).(fn).stats.stdBoot = std(xta.(dir)(iExp).(fn).XBoot(:, selT), 0, 2, 'omitnan');
            else
                xta.(dir)(iExp).(fn).stats.stdBoot = [];
            end
        end
    end
end

% %% Post-hoc boot mi
% p.mi.features = ["HandR", "HandL", "Jaw"];
% p.mi.windowPre = [-0.3, 0];
% p.mi.windowPost = [0, 0.3];
% p.mi.nBoot = 10000;
% 
% for iExp = 1:length(exp)
%     kine = kinematics(iExp);
%     maxT = kine.HandR.t(end);
%     for dir = ["dip", "rise"]
%         miBootTemp = NaN(p.mi.nBoot, length(p.mi.features));
%         nTrials = length(xta.(dir)(iExp).(fn).t0);
%         for fn = p.mi.features
%         end
%     end
% end

%% Save results
exportPath = fullfile("C:\SERVER\Units", sprintf("population_dta_rta_%i_%i_%ito%ims_%iboots_%s.mat", 100*p.xta.dip.thresholdQuantile, 100*p.xta.dip.thresholdSubQuantile, p.spikeRes*1000*p.xta.dip.samples(1), p.spikeRes*1000*p.xta.rise.samples(2), p.nBoot, datetime("now", Format="yyyyMMdd")));
save(exportPath, 'xta', 'kinematics', 'p', '-v7.3')
fprintf("Saved to %s\n", exportPath);


clear d dir duration iPat nTotal s selUnits sn t0 tt




%% Try to get a thingamajig
% Umm. it's a we do the thing.
% Okay
% Here goes
% 1) We calculate ISIs for each neuron in a session
% 2) we then shuffle the ISIs
% 3) reconstruct spike train from shuffled ISI
% 4) redo the spike rates.
% 5) re-detect population dips
% 6) count no. population dips
% 7) repeat 10000 times

p.isi.nBoot = 1000;

% 1) We calculate ISIs for each neuron in a session
lineLength = 0;
xCountBootDip = NaN(p.isi.nBoot, length(exp));
xCountBootRise = NaN(p.isi.nBoot, length(exp));
tTicTotal = tic();
for iExp = 1:length(exp)
    spikeTimes = arrayfun(@(eu) eu.SpikeTimes, exp(iExp).eu, UniformOutput=false);
    lastSpike = max(arrayfun(@(eu) eu.SpikeTimes(end), exp(iExp).eu));
    edges = 0:p.spikeRes:lastSpike;
    t = 0.5*(edges(1:end-1) + edges(2:end));
    for iBoot = 1:p.isi.nBoot
        X = NaN(length(exp(iExp).eu), length(t));
        for iEu = 1:length(exp(iExp).eu)
            % Shuffle spike times
            st0 = spikeTimes{iEu}(1);
            isi = diff(spikeTimes{iEu});
            st = cumsum([st0, isi(randperm(length(isi), length(isi)))]);
            switch p.spikeDataSource
                case "count"
                    [x, ~] = exp(iExp).eu(iEu).getSpikeCounts(edges, spikeTimes=st);
                    x = double(x)./p.spikeRes;
                case "rate"
                    switch p.spikeKernelType
                        case 'gaussian'
                            [x, ~] = exp(iExp).eu(iEu).getSpikeRates('gaussian', p.spikeKernelSigma, edges, kernelWidth=p.spikeKernelWidth, spikeTimes=st);
                        case 'exponential'
                            [x, ~] = exp(iExp).eu(iEu).getSpikeRates('exponential', p.spikeKernelLambda1, p.spikeKernelLambda2, edges, kernelWidth=p.spikeKernelWidth, spikeTimes=st);
                    end
                otherwise
                    error("Unknown p.spikeDataSource=%s", p.spikeDataSource)
            end
            x = (x-mean(x))/std(x, 0);
            X(iEu, :) = x;
        end
        x = mean(X, 1, 'omitnan');

        t0 = struct(dip=[], rise=[]);
        duration = struct(dip=[], rise=[]);
        nTotal = struct(dip=[], rise=[]);
        threshold = struct(dip=[], rise=[]);
        % Find onset of dips and rises
        for dir = ["dip", "rise"]
            switch dir
                case "dip"
                    threshold.(dir) = quantile(x(x<0), p.xta.(dir).thresholdQuantile);
                    i0 = arrayfun(@(i) strfind(x<=threshold.(dir), p.xta.(dir).pattern{i}) + p.xta.(dir).patternOnset(i), 1:length(p.xta.(dir).pattern), UniformOutput=false);
                case "rise"
                    threshold.(dir) = quantile(x(x>0), p.xta.(dir).thresholdQuantile);
                    i0 = arrayfun(@(i) strfind(x>=threshold.(dir), p.xta.(dir).pattern{i}) + p.xta.(dir).patternOnset(i), 1:length(p.xta.(dir).pattern), UniformOutput=false);
            end
            d = cell(size(i0));
            for iPat = 1:length(p.xta.(dir).pattern)
                d{iPat} = repmat(uint8(nnz(p.xta.(dir).pattern{iPat})), size(i0{iPat}));
            end
            i0 = cat(2, i0{:});
            t0.(dir) = t(i0);
            nTotal.(dir) = length(t0.(dir));
            duration.(dir) = cat(2, d{:}); % Multiply by p.spikeRes to get dip duration in seconds

            % Blank out dips around certain behavioral/stim events (opto, etc.)
            for iEvent = 1:length(p.blank)
                tEvent = exp(iExp).eu(1).EventTimes.(p.blank(iEvent).event);
                windows = tEvent(:) + p.blank(iEvent).window;
                for i = 1:length(tEvent)
                    sel = t0.(dir)>=windows(i, 1) & t0.(dir)<=windows(i, 2);
                    t0.(dir)(sel) = [];
                    duration.(dir)(sel) = [];
                end
            end
            switch dir
                case "dip"
                    xCountBootDip(iBoot, iExp) = length(t0.(dir));
                case "rise"
                    xCountBootRise(iBoot, iExp) = length(t0.(dir));
            end
        end

        fprintf(repmat('\b', [1, lineLength]))
        lineLength = fprintf('Session %i/%i (iBoot=%i); %i(-%i) dips (x<%.2f), %i(-%i) rises (x>%.2f)... %.1fs elapsed...', iExp, length(exp), iBoot, length(t0.dip), nTotal.dip-length(t0.dip), threshold.dip, length(t0.rise), nTotal.rise-length(t0.rise), threshold.rise, toc(tTicTotal));

    end
end
% 2) we then shuffle the ISIs
% 3) reconstruct spike train from shuffled ISI
% 4) redo the spike rates.
% 5) re-detect population dips
% 6) count no. population dips
% 7) repeat 10000 times



%% Save results

exportPath = fullfile("C:\SERVER\Units", sprintf("population_dta_rta_%i_%i_%ito%ims_%iboots_%s.mat", 100*p.xta.dip.thresholdQuantile, 100*p.xta.dip.thresholdSubQuantile, p.spikeRes*1000*p.xta.dip.samples(1), p.spikeRes*1000*p.xta.rise.samples(2), p.nBoot, datetime("now", Format="yyyyMMdd")));
save(exportPath, 'xta', 'kinematics', 'p', 'xCountBootRise', 'xCountBootDip', '-v7.3')
fprintf("Saved to %s\n", exportPath);


%% Report nDips, nRises from shuffled spikes
alpha = 0.05;
xCountBoot = struct(dip=xCountBootDip*p.xta.dip.thresholdSubQuantile, rise=xCountBootRise*(1-p.xta.rise.thresholdSubQuantile));
xCountBootCI = struct(dip=[], rise=[]);
h = struct(dip=[], rise=[]);
for dir = ["dip", "rise"]
    xCount.(dir) = arrayfun(@(xta) length(xta.t0), xta.(dir));
    xCountBootCI.(dir) = quantile(xCountBoot.(dir), [alpha/2, 1-alpha/2], 1);
    h.(dir) = xCount.(dir)>=xCountBootCI.(dir)(1, :) & xCount.(dir)<=xCountBootCI.(dir)(2, :);

    fprintf("On average across %i sessions, there are %.3f±%.3f (mean±sd) observed %ss, whereas bootstrapped 95%%CI is [%.3f, %.3f].\n", length(xCount.(dir)), mean(xCount.(dir), 2), std(xCount.(dir), 0, 2), dir, mean(xCountBootCI.(dir)(1, :), 2), mean(xCountBootCI.(dir)(2, :), 2))
end

for dir = ["dip", "rise"]
    fprintf("On average across %i sessions, there are %.0f±%.0f (mean±sd) observed %ss, whereas bootstrapped 95%%CI is [%.0f, %.0f].\n", length(xCount.(dir)), mean(xCount.(dir), 2)*20, std(xCount.(dir), 0, 2)*20, dir, mean(xCountBootCI.(dir)(1, :), 2)*20, mean(xCountBootCI.(dir)(2, :), 2)*20)
end

fig = figure;
tl = tiledlayout(fig, 2, 1);
for dir = ["dip", "rise"]
    ax = nexttile(tl); hold(ax, 'on')
    scatter(ax, 1:length(exp), xCount.(dir), 'ko', DisplayName='observed')
    if any(h.(dir))
        scatter(ax, 1:length(exp), xCount.(dir)(h.(dir)), 'ko', 'filled', DisplayName='observed')
    end
    patch(ax, [1:length(exp), length(exp):-1:1], [xCountBootCI.(dir)(1, :), xCountBootCI.(dir)(2, end:-1:1)], [0.15, 0.15, 0.15], FaceAlpha=0.1, EdgeAlpha=0.1, DisplayName='95%CI')
    xlabel(ax, 'session')
    ylabel(ax, sprintf('no. %ss', dir))
    legend(ax, Location='north')
    title(ax, dir)
end

