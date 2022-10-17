%% Create EphysUnit objects for acute units. This time do not cull ITI 
% spikes, also include lampON/Off events as LIGHT trials. Duplicates are
% not included.

% clear
% 
% whitelist = dir('C:\SERVER\Units\Lite_NonDuplicate\*.mat');
% whitelist = {whitelist.name}';
% whitelist = cellfun(@(x) strsplit(x, '.mat'), whitelist, UniformOutput=false);
% whitelist = cellfun(@(x) x{1}, whitelist, UniformOutput=false);
% 
% ar = AcuteRecording.load();
% for i = 1:length(ar)
%     clear eu
%     try  
%         eu = EphysUnit(ar(i), savepath='C:\SERVER\Units\acute_2cam', whitelist=whitelist, ...
%             cullITI=false, readWaveforms=false);
%     catch ME
%         warning('Error while processing file %g (%s)', i, ar(i).expName);
%     end
% end

%% 1. Load data
%% 1.1. Load acute EU objects (duplicates already removed)
eu = EphysUnit.load('C:\SERVER\Units\acute_2cam'); 

%% Remove multiunit detected by ISI test.
p.ISIThreshold = 0.0015;
for iEu = 1:length(eu)
    st = eu(iEu).SpikeTimes;
    isi = [NaN, diff(st)];
    st(isi == 0) = [];
    isi = [NaN, diff(st)];
    eu(iEu).SpikeTimes = st;
    ISI{iEu} = isi;
end

for iEu = 1:length(eu)
    prcLowISI(iEu) = nnz(ISI{iEu} < p.ISIThreshold) ./ length(ISI{iEu});
end
histogram(prcLowISI, 0:0.01:1)
cat.isMultiUnit = prcLowISI > 0.05;
cat.isSingleUnit = prcLowISI <= 0.05;
eu = eu(cat.isSingleUnit);
clearvars -except p eu

eu = eu';

%%
% 1.2. Load Video Tracking Data (vtd) and ArduinoConnection (ac), and group into experiments
clearvars -except eu
exp = CompleteExperiment(eu);

% 1.3 Align video and ephys timestamps
exp.alignTimestamps();

%% Get first significant arm movemetn befor trials
theta = 2.5;
[velocityKernels, ~, ~] = CompleteExperiment.makeConsineKernels(0, width=0.2, overlap=0.5, direction='both');

TST = cell(1, length(exp));
for iExp = 1:length(exp)
    disp(iExp)
    trials = exp(iExp).eu(1).getTrials('press');
    trueStartTime = NaN(1, length(trials));
    for iTrial = 1:length(trials)
        t = trials(iTrial).Start:1/30:trials(iTrial).Stop;
        F = exp(iExp).getFeatures(timestamps=t, features={'handL', 'handR'}, stats={'xPos', 'yPos'});
        F = CompleteExperiment.convolveFeatures(F, velocityKernels, kernelNames={'_smooth'}, ...
            features={'handL', 'handR'}, ...
            stats={'xPos', 'yPos'}, ...
            mode='replace', normalize='maxabs');
        F.inTrial = [];
        F.t = [];
        F = normalize(F);
        data = table2array(F);
        trueStartIndex = find(all(abs(data) < theta, 2), 1, 'last');
        if ~isempty(trueStartIndex)
            trueStartTime(iTrial) = t(trueStartIndex) - t(end);
        end
    %     lclips = exp(i).getVideoClip(trueStartTime, side='l', numFramesBefore=30);
    %     rclips = exp(i).getVideoClip(trueStartTime, side='r', numFramesBefore=30);
    %     implay(lclips, 30)
    %     implay(rclips, 30)
        % plot(F.t - F.t(end), table2array(F(:, {'handL_xVel', 'handL_yVel', 'handR_xVel', 'handR_yVel'})))
    end
    TST{iExp} = trueStartTime;
end
%%
figure(DefaultAxesFontSize=14)
tst = cat(2, TST{:});
histogram(-tst, 0:0.05:2, FaceColor="auto", Normalization='probability')
xlabel('Movement duration (s)')
ylabel('Probability')
legend(sprintf('%g trials, %g animals', length(tst), length(exp)))

clearvars -except eu exp TST
    
%% 1.4 Get video clips around event times to verify stuff
i = 9;
eventTimes = exp(i).getEventTimestamps('Press');
eventTimes = eventTimes(1);%randi(length(eventTimes)));
lclips = exp(i).getVideoClip(eventTimes, side='l', numFramesBefore=30);
rclips = exp(i).getVideoClip(eventTimes, side='r', numFramesBefore=30);
implay(lclips, 30)
implay(rclips, 30)

%% 2. Build GLM
clearvars -except exp

%% 2.1 Examine data
clearvars -except exp
close all
iExp = 1;
iEu = 1;
pTheta = 0.95;
bodyparts = {'handIpsi', 'footIpsi', 'spine', 'tail', 'nose', 'tongue'};

%% 2.1.1 Plot timecourse of velocity and position data, see if there are correlations (redundancies).
for ibp = 1:length(bodyparts)
    bp = bodyparts{ibp};
    posX = exp(iExp).vtdL.(sprintf('%s_X', bp));
    posY = exp(iExp).vtdL.(sprintf('%s_Y', bp));
    velX = [0; diff(posX)];
    velY = [0; diff(posY)];
    prob = exp(iExp).vtdL.(sprintf('%s_Likelihood', bp));
    isUncertain = prob < pTheta;
    velX(isUncertain) = NaN;
    velY(isUncertain) = NaN;
    posX(isUncertain) = NaN;
    posY(isUncertain) = NaN;


    fig = figure(Units='normalized', OuterPosition=[0.1, 0.2, 0.8, 0.6]);
    suptitle(bp)

    ax = subplot(2, 2, 1); hold on;
    t = exp(iExp).vtdL.Timestamp;
    xname = sprintf('%s_VelX', bp);
    yname = sprintf('%s_VelY', bp);
    plot(ax, t, velX, DisplayName=xname)
    plot(ax, t, velY, DisplayName=yname)
    xlabel(ax, 't', Interpreter='none')
    ylabel(ax, 'vel', Interpreter='none')
    legend(ax, Interpreter='none')
    hold off;
    
    ax = subplot(2, 2, 2);
    scatter(ax, velX, velY)
    xlabel(ax, xname, Interpreter='none')
    ylabel(ax, yname, Interpreter='none')


    ax = subplot(2, 2, 3); hold on;
    t = exp(iExp).vtdL.Timestamp;
    xname = sprintf('%s_X', bp);
    yname = sprintf('%s_Y', bp);
    plot(ax, t, posX, DisplayName=xname)
    plot(ax, t, posY, DisplayName=yname)
    xlabel(ax, 't', Interpreter='none')
    ylabel(ax, 'pos', Interpreter='none')
    legend(ax, Interpreter='none')
    hold off;
    
    ax = subplot(2, 2, 4);
    scatter(ax, posX, posY)
    xlabel(ax, xname, Interpreter='none')
    ylabel(ax, yname, Interpreter='none')
end
clearvars -except exp iExp iEu bodyparts pTheta


%% 2.2 Get Features
% clearvars -except exp
MDL = cell(length(exp), 1);
SRT = cell(length(exp), 1);
FT = cell(length(exp), 1);
TT = cell(length(exp), 1);
NOTNAN = cell(length(exp), 1);

for iExp = 1:length(exp)
    F = exp(iExp).getFeatures(sampleRate=30, trialType={'press'}, stats={'xVel', 'yVel'}, ...
        features={'handL', 'handR', 'footL', 'footR', 'nose', 'spine' 'trialStart', 'pressTrialRamp', 'firstPressRamp'}, ...
        likelihoodThreshold=0.95);
    [velocityKernels, ~, velocityDelays] = CompleteExperiment.makeConsineKernels(4, width=0.1, overlap=0.5, direction='both');
    velocityDelays = round(velocityDelays * 1000);
    F = CompleteExperiment.convolveFeatures(F, velocityKernels, kernelNames=velocityDelays, ...
        features={'handL', 'handR', 'footL', 'footR', 'nose', 'spine'}, ...
        stats={'xVel', 'yVel'}, ...
        mode='replace', normalize='maxabs');

    [eventKernels, ~, eventDelays] = CompleteExperiment.makeConsineKernels(6, width=0.4, overlap=0.75, direction='both');
    eventDelays = round(eventDelays * 1000);
    F = CompleteExperiment.convolveFeatures(F, eventKernels, kernelNames=eventDelays, ...
        features={'trialStart'}, ...
        mode='replace', normalize='none');

    t = F.t;
    inTrial = F.inTrial;
    F.t = [];
    F.inTrial = [];
    Ft = F(inTrial, :);
    tt = t(inTrial);

    Ft.constant = ones(height(Ft), 1);

    % Names of predictors
    % Cue and first move
    % Movement velocities
    % Trial-length ramping signals
    % Short pre-movement ramps
    names = Ft.Properties.VariableNames;
    cuePredictors = names(contains(names, {'trialStart'}))';
    velocityPredictors = names(contains(names, {'Vel'}))';
    timingInvariantRampPredictors = names(contains(names, {'firstPressRamp'}))';
    trialProgressPredictors = names(contains(names, 'TrialRamp'))';


    variantPredictors = { ...
        {'constant'}, ...
        [{'constant'}; cuePredictors], ...
        [{'constant'}; cuePredictors; velocityPredictors], ...
        [{'constant'}; cuePredictors; velocityPredictors; timingInvariantRampPredictors], ...
        [{'constant'}; cuePredictors; velocityPredictors; timingInvariantRampPredictors; trialProgressPredictors]};
    variantNames = {'Constant', '+Cue', '+Velocity', '+TimingInvariantRamp', '+TrialProgress'};
    nVariants = length(variantNames);
    
    % All nan if one nan
    notnan = all(~isnan(Ft(:, variantPredictors{end}).Variables), 2);
    prctNotNan = nnz(notnan) ./ height(Ft); 
    data = Ft.Variables;
    data(~notnan, :) = NaN;
    Ft.Variables = data;

    FT{iExp} = Ft;
    TT{iExp} = tt;

    % 2.3 Fit GLM
    % Grab a unit and start fitting glms!
    fprintf(1, 'Fitting %g modelsets...\n', length(exp(iExp).eu)); tTicAll = tic();
    mdl = cell(length(exp(iExp).eu), nVariants);
    srt = cell(length(exp(iExp).eu), 1);
    warning('off','all')
    for iEu = 1:length(exp(iExp).eu)
        fprintf(1, '\tFitting modelset %g of %g...\n', iEu, length(exp(iExp).eu)); tTic = tic();
        eu = exp(iExp).eu(iEu);
        srTrialAligned = [0, eu.getSpikeRates('gaussian', 0.1, t)]'; 
        srt{iEu} = srTrialAligned(inTrial);
    
        thisF = Ft;
        thisF.SpikeRate = double(srt{iEu});

        for iVariant = 1:nVariants
            mdl{iEu, iVariant} = fitglm(thisF, ResponseVar='SpikeRate', PredictorVars=variantPredictors{iVariant}, Distribution='poisson');
            fprintf(1, '\t\t%s R^2 = %.2f\n', variantNames{iVariant}, mdl{iEu, iVariant}.Rsquared.Ordinary);
        end

        % fprintf(1, '\tDone (%.0f%% not nan) in %.2f sec.\n', prctNotNan*100, toc(tTic));
    end
    warning('on','all')
    warning('query','all')
    fprintf(1, 'Fitted %g units in %.2f seconds.\n', length(exp(iExp).eu), toc(tTicAll));
    MDL{iExp} = cellfun(@compact, mdl, UniformOutput=false);
    SRT{iExp} = srt;
    NOTNAN{iExp} = notnan;
    
%     % Plot R^2 distribution for all units, compare different models.
%     fig = figure();
%     ax = axes(fig);
%     hold(ax, 'on')
%     for iVariant = 1:nVariants
%         N = histcounts(cellfun(@(mdl) mdl.Rsquared.Ordinary, mdl(:, iVariant)), 0:0.1:1);
%         plot(ax, 0.05:0.1:0.95, N, Color=getColor(iVariant, nVariants), LineWidth=2, DisplayName=variantNames{iVariant})
%     end
%     xlabel('R^2')
%     ylabel('Count')
%     legend(ax);
%     title(ax, sprintf('Model R^2 (%g units)', size(mdl, 1)))
%     hold(ax, 'off')
%     print(fig, sprintf('C:\\SERVER\\Figures\\GLM\\%s', exp(iExp).name), '-dpng');
end
mdl = cat(1, MDL{:});
srt = cat(1, SRT{:});
expIndices = zeros(size(mdl, 1), 1);
i = 0;
for iExp = 1:length(exp)
    expIndices(i + 1:i + length(exp(iExp).eu)) = iExp;
    i = i + length(exp(iExp).eu);
end
% clearvars -except exp mdl srt variantNames nVariants FT TT NOTNAN expIndices

% %% Estimate loss
% L = {};
% eu = [exp.eu];
% curExpIndex = -1;
% for iEu = 1:length(eu)
%     iExp = expIndices(iEu);
%     if curExpIndex ~= iEu
%         curExpIndex = iEu;
%         
%         F = exp(iExp).getFeatures(sampleRate=30, stats={'xVel', 'yVel'}, ...
%             features={'handL', 'handR', 'footL', 'footR', 'nose', 'spine' 'trialStart', 'firstPress', 'firstLick', 'pressTrialRamp', 'lickTrialRamp', 'firstPressRamp', 'firstLickRamp'}, ...
%             likelihoodThreshold=0.95);
%         [velocityKernels, ~, velocityDelays] = CompleteExperiment.makeConsineKernels(4, width=0.1, overlap=0.5, direction='both');
%         velocityDelays = round(velocityDelays * 1000);
%         F = CompleteExperiment.convolveFeatures(F, velocityKernels, kernelNames=velocityDelays, ...
%             features={'handL', 'handR', 'footL', 'footR', 'nose', 'spine'}, ...
%             stats={'xVel', 'yVel'}, ...
%             mode='replace', normalize='maxabs');
% 
%         [eventKernels, ~, eventDelays] = CompleteExperiment.makeConsineKernels(6, width=0.4, overlap=0.75, direction='both');
%         eventDelays = round(eventDelays * 1000);
%         F = CompleteExperiment.convolveFeatures(F, eventKernels, kernelNames=eventDelays, ...
%             features={'trialStart'}, ...
%             mode='replace', normalize='none');
%     
%         [eventKernels, ~, eventDelays] = CompleteExperiment.makeConsineKernels(6, width=0.4, overlap=0.75, direction='left');
%         eventDelays = round(eventDelays * 1000);
%         F = CompleteExperiment.convolveFeatures(F, eventKernels, kernelNames=eventDelays, ...
%             features={'firstPress', 'firstLick'}, ...
%             mode='replace', normalize='none');
%     end
% 
%     for iVariant = 1:nVariants
%         ll = loss(mdl(iEu, iVariant), )
%     end
% 
% end


% 2.3 Plot R^2 distribution for all units, compare different models.
close all
edges = 0:0.05:1;
centers = (edges(1:end-1) + edges(2:end))*0.5;

fig = figure(Units='normalized', Position=[0.5, 0, 0.5,0.4], DefaultAxesFontSize=12);

paramNames = {'Ordinary', 'Adjusted', 'AdjGeneralized', 'LLR', 'Deviance'};
np = length(paramNames);

clear ax
np = 1;
for ip = 1:np
%     ax = subplot(np, 3, 3*(ip-1)+1);
%     hold(ax, 'on')
%     for iVariant = 2:nVariants
%         N = histcounts(cellfun(@(mdl) mdl.Rsquared.(paramNames{ip}), mdl(:, iVariant)), edges);
%         plot(ax, centers, N, Color=getColor(iVariant-1, nVariants-1), LineWidth=2, DisplayName=variantNames{iVariant})
%     end
% %     xlabel(sprintf('R^2 %s', paramNames{ip}))
%     xlabel('R^2')
%     ylabel('Count')
%     legend(ax, Location='northeast');
%     hold(ax, 'off')

    ax = subplot(np, 2, 2*(ip-1)+1);
    hold(ax, 'on')
    for iVariant = 2:nVariants
        N = histcounts(cellfun(@(mdl) mdl.Rsquared.(paramNames{ip}), mdl(:, iVariant)), edges, Normalization='probability');
        plot(ax, edges, [0, cumsum(N)], Color=getColor(iVariant-1, nVariants-1), LineWidth=2, DisplayName=variantNames{iVariant})
    end
%     xlabel(sprintf('R^2 %s', paramNames{ip}))
    xlabel('R^2')
    ylabel('Cumulative probability')
    legend(ax, Location='southeast');
    hold(ax, 'off')

    ax = subplot(np, 2, 2*(ip-1)+2);
    hold(ax, 'on')
    R2 = cellfun(@(mdl) mdl.Rsquared.(paramNames{ip}), mdl);
    dR2 = diff(R2, 1, 2);

    x = repmat(1:nVariants-1, [size(dR2, 1), 1]);
    x = x(:);
    y = dR2(:);
    swarmchart(ax, x, y, 3.7, 'filled', 'k')
    
    boxplot(ax, dR2, Symbol='.', OutlierSize=0.000001, Color='k', Whisker=0)
    
    xticks(ax, 1:nVariants-1)
    xticklabels(ax, variantNames(2:end))
    xtickangle(ax, 315)
    ylabel('\DeltaR^2')
    ylim(ax, [0, max(y)+0.1])
    xlim(ax, [0,nVariants])
end


% 2.4 Plot fitted vs observed


%% 2.4.2 Single units, trial averaged, fitted vs observed
nVariants = length(variantNames);
try
    eu = [exp.eu];
catch
    eu = vertcat(exp.eu);
end
for iEu = 84%1:length(mdl) %39 84
    iExp = expIndices(iEu);
    tt = TT{iExp};
    Ft = FT{iExp};

    fig = figure(Units='normalized', OuterPosition=[0, 0, 0.8, 1]);

    ax = subplot(3, 1, 1);
    hold(ax, 'on')
    plot((1:length(srt{iEu}))./30, srt{iEu}, 'k:', LineWidth=2, DisplayName='Observed');
    for iVariant = 2:nVariants
        yHat = predict(mdl{iEu, iVariant}, FT{expIndices(iEu)});
        plot((1:length(srt{iEu}))./30, yHat, Color=getColor(iVariant-1, nVariants-1), LineWidth=1.2, DisplayName=sprintf('%s (R^2=%.2f)', variantNames{iVariant}, mdl{iEu, iVariant}.Rsquared.Ordinary));
    end

    hold(ax, 'off')
    title(ax, eu(iEu).getName('_'), Interpreter='none');
    xlabel(ax, 'Time (s)')
    ylabel(ax, 'Spike rate (sp/s)')
    xlim(ax, [30, 60])
    legend(ax, Location='northwest');

    ax = subplot(3, 1, 3);
    iVariant = nVariants;
    ax.TickLabelInterpreter = 'none';
    hold(ax, 'on')
    x = 2:height(mdl{iEu, iVariant}.Coefficients);
    errorbar(ax, x, mdl{iEu, iVariant}.Coefficients.Estimate(2:end), mdl{iEu, iVariant}.Coefficients.SE(2:end));
    plot(ax, [x(1), x(end)], [0, 0], 'k--', LineWidth=1.5)
    hold(ax, 'off')
    xticks(ax, x);
    xticklabels(ax, mdl{iEu, iVariant}.CoefficientNames(x));
    xlim(ax, [x(1), x(end)])
    ylabel(ax, 'Coefficient +/- SE')
    title(ax, sprintf('Coefficients (%s)', variantNames{iVariant}))
    set(ax, FontSize=9)


    srtHat = NaN(height(Ft), size(mdl, 2));
    for iVariant = 2:size(mdl, 2)
        srtHat(:, iVariant) = predict(mdl{iEu, iVariant}, Ft);
    end

    trials = exp(iExp).eu(1).getTrials('press');

    maxTrialSampleLength = 0;
    for iTrial = 1:length(trials)
        inTrial = tt >= trials(iTrial).Start & tt <= trials(iTrial).Stop;
        maxTrialSampleLength = max(maxTrialSampleLength, nnz(inTrial));
    end

    srTrialAligned = NaN(length(trials), maxTrialSampleLength);
    srTrialAlignedHat = NaN(length(trials), maxTrialSampleLength, size(mdl, 2));
    for iTrial = 1:length(trials)
        inTrial = tt >= trials(iTrial).Start & tt <= trials(iTrial).Stop;
        trialSampleLength = nnz(inTrial);
        srTrialAligned(iTrial, maxTrialSampleLength - trialSampleLength + 1:end) = srt{iEu}(inTrial);
        for iVariant = 2:size(mdl, 2)
            srTrialAlignedHat(iTrial, maxTrialSampleLength - trialSampleLength + 1:end, iVariant) = srtHat(inTrial, iVariant);
        end
    end


    % Trial aligned mean spike rate
    msr = mean(srTrialAligned, 1, 'omitnan');
    msrHat = squeeze(mean(srTrialAlignedHat, 1, 'omitnan'));
    sd = std(srTrialAligned, 0, 1, 'omitnan');
    sdHat = squeeze(std(srTrialAlignedHat, 0, 1, 'omitnan'));
    t = (-maxTrialSampleLength + 1:0)*median(diff(tt));
    
    ax = subplot(3, 1, 2);
    hold(ax, 'on')
    clear h
    h(1) = plot(ax, t, msr, 'k:', LineWidth=3, DisplayName='Observed');
    for iVariant = 2:size(mdl, 2)
        h(iVariant) = plot(ax, t, msrHat(:, iVariant), Color=getColor(iVariant-1, nVariants-1), LineWidth=1.5, DisplayName=sprintf('%s (R^2=%.2f)', variantNames{iVariant}, mdl{iEu, iVariant}.Rsquared.Ordinary));
    end
%         h(end+1) = patch(ax, [t, flip(t)], [msr + sd, flip(msr - sd)], 'k', FaceAlpha=0.1, EdgeAlpha=0, DisplayName='SD');
    for iVariant = 2:size(mdl, 2)
%             patch(ax, [t, flip(t)]', [msrHat(:, iVariant) + sdHat(:, iVariant); flip(msrHat(:, iVariant) - sdHat(:, iVariant))], ...
%                getColor(iVariant-1, nVariants-1), FaceAlpha=0.1, EdgeAlpha=0, DisplayName=sprintf('%s \\pm SD', variantNames{iVariant}));
    end
    legend(ax, h, Location='northwest')
    xlim(ax, [-4, 0])
    xlabel(ax, sprintf('Time relative to %s (s)', 'lever-touch'))
    ylabel(ax, 'Spike rate (sp/s)')
    titleText = sprintf('%s (%s)', eu(iEu).getName('_'), 'lever');
    title(ax, titleText, Interpreter='none')

    if ~isfolder('C:\\SERVER\\Figures\\GLM')
        mkdir('C:\\SERVER\\Figures\\GLM')
    end
    print(fig, sprintf('C:\\SERVER\\Figures\\GLM\\%s', eu(iEu).getName('_')), '-dpng');
%     close(fig)

    clear h ax
end

function c = getColor(i, n)
    c = hsl2rgb([0.88*(i-1)./n, 1, 0.5]);
end