% expDesc = "WT x SNr(AAV-syn-CoChR)";
% eu = EphysUnit.load('C:\SERVER\Units\ReachVsLick_1225');
% load('C:\SERVER\Units\meta_ReachVsLick_1225_20260728.mat') % 'boot', 'c', 'eta'

% Try running script_optrodeSNr_CoChR_VGATCre_20260729 again;
expDesc = "VGAT-Cre x SNr(AAV-flex-CoChR)";
load('C:\SERVER\Units\meta_SNr_CoChR_VGATCre_ValidVideos.mat')
eu = EphysUnit.load('C:\SERVER\Units\SNr_CoChR_VGATCre\SingleUnit_NonDuplicate_NonDrift_SNr_ValidVideos');

expNames = string({eu.ExpName}');
[uniqueExpNames, expToEuIndices, euToExpIndices] = unique(expNames);

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

%% Plot ETAs (reach vs. lick vs. opto; all sessions combined, or per session)
varsSnapshot = who; 

TRIALTYPES = ["press", "lick", "stim", "stim"];
SELTYPES = ["press", "lick", "press", "lick"];

close all
for iExp = 0%:length(uniqueExpNames)
    if iExp == 0
        selUnits = true(1, length(eu));
    else
        selUnits = reshape(euToExpIndices == iExp, size(c.isPressUp));
    end

    fig = figure(Units='inches', Position=[1, 1, 10, 6]);
    tl = tiledlayout(fig, 2, 2); 
    
    for i = 1:length(TRIALTYPES)
        trialType = TRIALTYPES(i);
        selType = SELTYPES(i);
        clear h
        iLine = 1;
        ax = nexttile(tl); hold(ax, 'on')
        switch selType
            case "press"
                selTypeDispName = "reach";
                selUp = c.isPressUp;
                selDown = c.isPressDown;
                selFlat = ~c.isPressResponsive;
            case "lick"
                selTypeDispName = "lick";
                selUp = c.isLickUp;
                selDown = c.isLickDown;
                selFlat = ~c.isLickResponsive;
        end
        switch trialType
            case "press"
                trialTypeDispName = "reach";
                eventDispName = "bar contact";
            case "lick"
                trialTypeDispName = "lick";
                eventDispName = "spout contact";
            case "stim"
                trialTypeDispName = "opto";
                eventDispName = "opto onset";
        end

        plot(ax, eta.(trialType).t, eta.(trialType).X(selUnits & selUp, :), Color=[1, 0, 0, 0.1])
        if any(selUnits & c.isPressDown)
            plot(ax, eta.(trialType).t, eta.(trialType).X(selUnits & selDown, :), Color=[0, 0, 1, 0.1])
        end
        plot(ax, eta.(trialType).t, eta.(trialType).X(selUnits & selFlat, :), Color=[0, 0, 0, 0.1])
        h(iLine) = plot(ax, eta.(trialType).t, mean(eta.(trialType).X(selUnits & selUp, :), 1, 'omitnan'), Color=[0.8, 0.2, 0.2], LineWidth=1, DisplayName=sprintf('%s-inc (n=%i)', selTypeDispName, nnz(selUnits & selUp)));
        iLine = iLine + 1;
        if any(selUnits & c.isPressDown)
            h(iLine) = plot(ax, eta.(trialType).t, mean(eta.(trialType).X(selUnits & selDown, :), 1, 'omitnan'), Color=[0.2, 0.2, 0.8], LineWidth=1, DisplayName=sprintf('%s-dec (n=%i)', selTypeDispName, nnz(selUnits & selDown)));
            iLine = iLine + 1;
        end
        h(iLine) = plot(ax, eta.(trialType).t, mean(eta.(trialType).X(selUnits & selFlat, :), 1, 'omitnan'), Color=[0.2, 0.2, 0.2], LineWidth=1, DisplayName=sprintf('%s-flat (n=%i)', selTypeDispName, nnz(selUnits & selFlat)));
        iLine = iLine + 1;

        xline(ax, 0, 'k:')
        switch trialType
            case {"press", "lick"}
                xlim(ax, [-2, 3])
            case "stim"
                xlim(ax, [-1, 4])
        end
        ylim(ax, [-2, 4])
        lgd = legend(h(isvalid(h)), Location='northeast');
        lgd.ItemTokenSize = [9, 9];
        title(ax, trialTypeDispName)
        xlabel(ax, sprintf('time to %s (s)', eventDispName))
    end
    
    ylabel(tl, 'z-scored spike rate (a.u.)')

    if iExp == 0
        ttl = sprintf("%s - %i sessions", expDesc, length(uniqueExpNames));
    else
        ttl = sprintf("%s - %i %s", expDesc, iExp, eu(expToEuIndices(iExp)).ExpName);
    end
    title(tl, ttl, Interpreter='none');
    fontsize(fig, 9, 'points')
    print(fig, fullfile("C:\Users\AssadLab\Pictures\SNrOpto", sprintf("%s.png", ttl)), '-dpng')
end

clearvars('-except', varsSnapshot{:}, 'finalAnswer');

%% Plot scattered METAs (reach/lick vs. opto)
varsSnapshot = who; 

fig = figure(Units='inches', Position=[1, 1, 10, 5]);
tlp = tiledlayout(fig, 1, 2); 
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
lgd = legend(h, Orientation='horizontal');
lgd.Layout.Tile = 'north';

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
lgd = legend(h, Orientation='horizontal');
lgd.Layout.Tile = 'north';
ax = nexttile(tl(2), 'east'); AX(4) = ax;

ax = nexttile(tl(2), 'east'); hold(ax, 'on'); AX(2) = ax;
edges = -2.5:0.5:6.5;
histogram(ax, meta.stim(c.isLickUp), edges, FaceColor=[.8,.2,.2], EdgeColor=[.8,.2,.2], EdgeAlpha=0.5, Orientation='horizontal')
histogram(ax, meta.stim(c.isLickDown), edges, FaceColor=[.2,.2,.8], EdgeColor=[.2,.2,.8], EdgeAlpha=0.5, Orientation='horizontal')
histogram(ax, meta.stim, edges, DisplayStyle='stairs', EdgeColor=[.2,.2,.2], EdgeAlpha=0.5, Orientation='horizontal')

xlim(AX([1, 3]), [-3, 7])
ylim(AX, [-3, 7])

clearvars('-except', varsSnapshot{:}, 'finalAnswer');
%% Load laserpower data (SLOW!)
V = cell(length(uniqueExpNames), 1);
for iExp = 1:length(uniqueExpNames)
    eu0 = eu(expToEuIndices(iExp));
    stimOn = eu0.EventTimes.LaserModBlueOn;
    stimOff = eu0.EventTimes.LaserModBlueOff;

    tr = TetrodeRecording();
    tr.SelectFiles(NeuropixelPath=fullfile('C:\SERVER', eu0.getAnimalName(), eu0.ExpName))
    tr.LoadNeuropixelIO();
    tr.ParseNeuropixelIO(DigitalChannels={'Sync', 0; 'Lick', 1; 'Press', 2; 'Reward', 3; 'Timeout', 4; 'Mot2Busy', 5; 'CueLeft', 6; 'CueRight', 7});

    assert(strcmpi(tr.AnalogIn.ChannelNames{1}, 'LaserModBlue'))
    assert(tr.AnalogIn.ChannelIndex(1)==1)
    V{iExp} = NaN(length(stimOn), 1);
    for iStim = 1:length(stimOn)
        [a, b] = isin(tr.AnalogIn.Timestamps, [stimOn(iStim), stimOff(iStim)], true, true);
        V{iExp}(iStim) = mean(tr.AnalogIn.Data(1, a:b), 2, 'omitnan');
    end
    clear eu0 tr stimOn stimOff a b iStim
end
clear iExp

%% Parse laser power data 
LP = cell(length(uniqueExpNames), 1);
p.laserPowers = [5, 20]*1e-3;
p.laserPowerVThreshold = 0.6;
for iExp = 1:length(V)
    LP{iExp} = NaN(size(V{iExp}));
    LP{iExp}(V{iExp}<=p.laserPowerVThreshold) = p.laserPowers(1);
    LP{iExp}(V{iExp}>p.laserPowerVThreshold) = p.laserPowers(2);
end

clear stimData
stimData(length(V)) = struct(iExp=[], name=[], power=[], ain=[], stimOn=[], stimOff=[], duration=[]);
for iExp = 1:length(V)
    stimData(iExp).iExp = iExp;
    stimData(iExp).name = eu(expToEuIndices(iExp)).ExpName;
    stimData(iExp).power = LP{iExp};
    stimData(iExp).ain = V{iExp};
    stimData(iExp).stimOn = eu(expToEuIndices(iExp)).EventTimes.LaserModBlueOn(:);
    stimData(iExp).stimOff = eu(expToEuIndices(iExp)).EventTimes.LaserModBlueOff(:);
    stimData(iExp).duration = stimData(iExp).stimOff - stimData(iExp).stimOn;
end

% Trim short pulses
p.pulseDurations = [1, 3];
p.pulseDurationRes = 0.1;
for iExp = 1:length(V)
    d = round(stimData(iExp).duration./p.pulseDurationRes).*p.pulseDurationRes;
    sel = ismember(d, p.pulseDurations);
    fprintf('iExp=%i, removed %i of %i short pulses (mean(removed)=%g, mean(remaining)=%g).\n', iExp, nnz(~sel), length(sel), mean(stimData(iExp).duration(~sel), 'omitnan'), mean(stimData(iExp).duration(sel), 'omitnan'))
    stimData(iExp).power(~sel) = [];
    stimData(iExp).ain(~sel) = [];
    stimData(iExp).stimOn(~sel) = [];
    stimData(iExp).stimOff(~sel) = [];
    stimData(iExp).duration = d;
    stimData(iExp).duration(~sel) = [];
end

fprintf('\n')
% Calculate stim hash
for iExp = 1:length(V)
    [lia, iPower] = ismember(stimData(iExp).power, p.laserPowers);
    assert(all(lia)), clear lia
    [lia, iDuration] = ismember(stimData(iExp).duration, p.pulseDurations);
    assert(all(lia)), clear lia
    hash = iPower*10 + iDuration;
    stimData(iExp).iPower = iPower;
    stimData(iExp).iDuration = iDuration;
    stimData(iExp).hash = hash;
    fprintf('iExp=%i (', iExp)
    [uniqueHashes, ia] = unique(hash);
    for i = 1:length(ia)
        fprintf("[hash=%i, power=%gmW, duration=%gs],\t", uniqueHashes(i), p.laserPowers(iPower(ia(i)))*1e3, p.pulseDurations(iDuration(ia(i))));
    end
    fprintf('\n')
end

%% Check that all pMove thresholds for stim are 0.3;
% for iExp = 1:length(stimData)
%     animalName = strsplit(stimData(iExp).name, '_');
%     animalName = animalName{1};
%     try
%         edit(fullfile("C:\SERVER", animalName, stimData(iExp).name, sprintf('%s.m', stimData(iExp).name)))
%     catch
%         iExp = 5;
%         edit(fullfile("C:\SERVER", 'desmond47', stimData(iExp).name, sprintf('%s.m', stimData(iExp).name)))
%     end
% end
% 
% 'desmond46_20260710', 0.3
% 'desmond46_20260717', 0.3
% 'desmond46_20260722', 0.3
% 'desmond47_20260710', 0.3
% 'desmond47_20260716', 0.3
% 'desmond47_20260724', 0.3
% 'desmond47_20260729', 0.3

%% Get controls for decoder-stim
tLocal = -4:0.01:3+2;

p.decoderDataDelay = 0.66;
p.decoderSampleRate = 50; % 20ms intervals
p.decoderSmoothWindow = 0.1;
for iExp = 1:length(V)
    eu0 = eu(expToEuIndices(iExp));
    stimOn = eu0.EventTimes.LaserModBlueOn;
    stimOff = eu0.EventTimes.LaserModBlueOff;

    [t, P, X] = readDecoderData(eu0, p.decoderDataDelay, p.decoderSampleRate, p.decoderSmoothWindow);

    % ax = axes(figure);
    % hold(ax, 'on')
    % plot(ax, t, (X - mean(X))./std(X, 0), 'k-')
    % plot(ax, t, P, 'b-')
    % ax.InteractionOptions.LimitsDimensions = "x";

    threshold = 0.3;
    tStimCtrl = t(strfind(P >= threshold, [0, 1]) + 1);
    timeoutStart = eu0.EventTimes.TIMEOUT_START;
    trials = Trial(timeoutStart+2, tStimCtrl, 'first', stimOn);
    tStimCtrl = [trials.Stop];
    % [B, t, I] = trials.inTrial(tStimCtrl);
    % [~, ia, ~] = unique(I);
    % tStimCtrl = t(ia);

    for t0 = stimOn(:)'
        tStimCtrl(isin(tStimCtrl, t0 + [-3, 3])) = [];
    end

    stimData(iExp).stimCtrl = tStimCtrl(:);
end

% Calulate a peri-stim movement histogram
edges = -2:0.1:5;
for iExp = 1:length(V)
    eu0 = eu(expToEuIndices(iExp));
    for trialType = ["press", "lick"]
        switch trialType
            case "press"
                tMove = eu0.EventTimes.FirstPress;
            case "lick"
                tMove = eu0.EventTimes.FirstLick;
        end
        T0 = struct(stim=[], ctrl=[]);
        T0.stim = stimData(iExp).stimOn;
        T0.ctrl = stimData(iExp).stimCtrl;
        for cond = ["stim", "ctrl"]
            stimData(iExp).psmh.(cond).(trialType).N = zeros(length(T0.(cond)), length(edges)-1);
            stimData(iExp).psmh.(cond).(trialType).edges = edges;
            for iTrial = 1:length(T0.(cond))
                t0 = T0.(cond)(iTrial);
                N = histcounts(tMove, t0 + edges);
                N(~isfinite(N)) = 0;
                stimData(iExp).psmh.(cond).(trialType).N(iTrial, :) = N;
            end
        end
    end
end
if isfield(stimData, 'psth')
    stimData = rmfield(stimData, 'psth');
end
%% Save stimData:
save('C:\SERVER\Units\meta_SNr_CoChR_VGATCre_ValidVideos_withoutStimData.mat', 'boot', 'c', 'eta', 'meta', 'kinematics', 'p')
save('C:\SERVER\Units\meta_SNr_CoChR_VGATCre_ValidVideos_justStimData.mat', 'stimData', '-v7.3')

%% Can skip all above, just do this to load
% Try running script_optrodeSNr_CoChR_VGATCre_20260729 again;
expDesc = "VGAT-Cre x SNr(AAV-flex-CoChR)";
eu = EphysUnit.load('C:\SERVER\Units\SNr_CoChR_VGATCre\SingleUnit_NonDuplicate_NonDrift_SNr_ValidVideos');

expNames = string({eu.ExpName}');
[uniqueExpNames, expToEuIndices, euToExpIndices] = unique(expNames);

load('C:\SERVER\Units\meta_SNr_CoChR_VGATCre_ValidVideos_withoutStimData.mat'); %'boot', 'c', 'eta', 'meta', 'kinematics', 'p'
load('C:\SERVER\Units\meta_SNr_CoChR_VGATCre_ValidVideos_justStimData.mat'); %'stimData'

%% Calculate psth (spike rates) - we're making large 3D arrays with lots of NaNs, so doing this on the fly is a bit faster than

% Get weak version of bootstrap:
p.useWeakIncDec = false;
meta.press = mean(eta.press.X(:, isin(eta.press.t, p.responseWindowPress)), 2, 'omitnan');
meta.lick = mean(eta.lick.X(:, isin(eta.lick.t, p.responseWindowLick)), 2, 'omitnan');
meta.pressBaseline = mean(eta.press.X(:, isin(eta.press.t, [-4, -2])), 2, 'omitnan');
meta.lickBaseline = mean(eta.lick.X(:, isin(eta.lick.t, [-4, -2])), 2, 'omitnan');

c.isPressUpWeak = meta.press > meta.pressBaseline;
c.isPressDownWeak = meta.press < meta.pressBaseline;
c.isLickUpWeak = meta.lick > meta.lickBaseline;
c.isLickDownWeak = meta.lick < meta.lickBaseline;


p.spikeDataSource = "rate"; % rate, count
p.spikeRes = 0.01; % 0.001
p.spikeKernelType = 'exponential';
switch p.spikeKernelType
    case 'gaussian'
        p.spikeKernelSigma = 0.075;
        p.spikeKernelWidth = 0.5;
        [~, ~, p.spikeKernel] = eu(1).getSpikeRates('gaussian', p.spikeKernelSigma, p.spikeRes, kernelWidth=p.spikeKernelWidth);
    case 'exponential'
        % p.spikeKernelLambda1 = 10;
        % p.spikeKernelLambda2 = 100;
        % p.spikeKernelWidth = 0.5;
        p.spikeKernelLambda1 = 5;
        p.spikeKernelLambda2 = 100;
        p.spikeKernelWidth = 1;
        [~, ~, p.spikeKernel] = eu(1).getSpikeRates('exponential', p.spikeKernelLambda1, p.spikeKernelLambda2, p.spikeRes, kernelWidth=p.spikeKernelWidth);
end
% plot(p.spikeKernel.t, p.spikeKernel.y)
% hold("on")

if isfield(stimData, 'psth')
    stimData = rmfield(stimData, 'psth');
end
tLocal = -4:0.01:3+2;
ll = 0;
for iExp = 1:length(uniqueExpNames)
    euIndicesInExp = find(euToExpIndices(:)'==iExp);
    t = tLocal(1):p.spikeRes:tLocal(end);
    for cond = ["stim", "ctrl"]
        switch cond
            case "stim"
                stimOn = stimData(iExp).stimOn;
                % stimOff = stimData(iExp).stimOff;
                stimOff = stimOn + 1e-6;
            case "ctrl"
                stimOn = stimData(iExp).stimCtrl;
                stimOff = stimOn + 1e-6;
        end
        if p.useWeakIncDec
            isIncPress = c.isPressUpWeak(euIndicesInExp);
            isDecPress = c.isPressDownWeak(euIndicesInExp);
            isIncLick = c.isLickUpWeak(euIndicesInExp);
            isDecLick = c.isLickDownWeak(euIndicesInExp);
        else
            isIncPress = c.isPressUp(euIndicesInExp);
            isDecPress = c.isPressDown(euIndicesInExp);
            isIncLick = c.isLickUp(euIndicesInExp);
            isDecLick = c.isLickDown(euIndicesInExp);
        end
        stimTrials = Trial(stimOn, stimOff, advancedValidation=false);
        assert(length(stimOn) == length(stimTrials))
        assert(p.spikeDataSource=="rate");
        X = NaN(length(stimTrials), length(t)-1, length(euIndicesInExp), 'single');
        for iEuInExp = 1:length(euIndicesInExp)
            [x, ~, ~] = eu(euIndicesInExp(iEuInExp)).getTrialAlignedData('rate', [t(1), t(end)], 'stim', trials=stimTrials, alignTo='start', ...
                resolution=p.spikeRes, includeInvalid=true, kernel=p.spikeKernel, artifacts=p.artifacts);
            xBaseline = eu(euIndicesInExp(iEuInExp)).getTrialAlignedData('rate', [-4, -2], 'stim', trials=eu(euIndicesInExp(iEuInExp)).Trials.Press, alignTo='stop', ...
                resolution=p.spikeRes, includeInvalid=true, kernel=p.spikeKernel, artifacts=p.artifacts);
            X(:, :, iEuInExp) = (x - mean(xBaseline, 'all', 'omitnan'))./std(xBaseline, 0, 'all', 'omitnan');
            fprintf(repmat('\b', [1, ll]))
            ll = fprintf("iExp=%i, %s, iEu=%i\n", iExp, cond, iEuInExp);
        end
        % X: trials x time x units
        stimData(iExp).psth.(cond).t = 0.5*(t(1:end-1) + t(2:end));
        stimData(iExp).psth.(cond).X = X; % mean(X, 3, 'omitnan');
        stimData(iExp).psth.(cond).XIncPress = X(:, :, isIncPress); % mean(X(:, :, isIncPress), 3, 'omitnan');
        stimData(iExp).psth.(cond).XDecPress = X(:, :, isDecPress); % mean(X(:, :, isDecPress), 3, 'omitnan');
        stimData(iExp).psth.(cond).XIncLick = X(:, :, isIncLick); % mean(X(:, :, isIncLick), 3, 'omitnan');
        stimData(iExp).psth.(cond).XDecLick = X(:, :, isDecLick); % mean(X(:, :, isDecLick), 3, 'omitnan');
        % Record number of units per session, so we can do weighted averaging of sessions in the future
        % stimData(iExp).psth.(cond).N = length(euIndicesInExp) * ones(length(stimTrials), 1);
        % stimData(iExp).psth.(cond).NIncPress = nnz(isIncPress) * ones(length(stimTrials), 1);
        % stimData(iExp).psth.(cond).NDecPress = nnz(isDecPress) * ones(length(stimTrials), 1);
        % stimData(iExp).psth.(cond).NIncLick = nnz(isIncLick) * ones(length(stimTrials), 1);
        % stimData(iExp).psth.(cond).NDecLick = nnz(isDecLick) * ones(length(stimTrials), 1);
    end
end


%% Plot opto-aligned movement kinemeatics + spike rates
% close all
xl = {[-2, 3], [-2, 5]};
kineFeatures = ["Jaw", "HandL", "HandR"];
features = ["X", "XInc", "XDec", "move", "Jaw", "HandL", "HandR"];
featureDispNames = ["spike rate", "spike rate (inc)", "spike rate (dec)", "move", "jaw", "l.hand", "r.hand"];
% yl = {[-0.5, 1.5], [-0.5, 1.5], [-0.5, 1.5], [0, 0.2], [0, 12], [0, 8], [0, 8]};
yl = {[-1.5, 3], [-1.5, 3], [-1.5, 3], [0, 0.2], [0, 12], [0, 8], [0, 8]};
featureUnits = ["(a.u.)", "(a.u.)", "(a.u.)", "units", "events/trial", "speed (a.u.)", "speed (a.u.)", "speed (a.u.)"];
p.kinematicDataSource = "spd";
p.showIndividualTracesFor = "units"; % "trials", "units"
p.showIndividualTracesAsCI = true; % true to plot 25%/75% CI as shaded area, false to plot single traces
lineStylesByPower = ["-", "-", "-"]; % ctrl, 5mW, 20mW
colorsByPower = [.2,.2,.2,1; .2,.2,.8,.5; .8,.2,.2,1]; % ctrl, 5mW, 20mW
showIndividualTracesByPower = [true, true, true];
individualTracesAlpha = 0.2;

clear l
l.h = ones(1, length(features));
l.ch = cumsum([1, l.h]);
l.w = cellfun(@diff, xl);
l.cw = cumsum([1, l.w]);
exportPath = "E:\Figures\SNr_VGAT-Cre-CoChR";
if p.useWeakIncDec
    exportPath = fullfile(exportPath, "weakIncDec");
else
    exportPath = fullfile(exportPath, "bootIncDec");
end
if exist(exportPath, 'dir')
    rmdir(exportPath, 's')
end
mkdir(exportPath)
ll = 0;

% Calculate opto-triggered average kinematics
for iExp = 1:length(uniqueExpNames)
    for ifn = 1:length(kineFeatures)
        fn = kineFeatures(ifn);
        % get whole-session continuous kinematics data
        switch p.kinematicDataSource
            case "pos"
                t = kinematics(iExp).(fn).t;
                x = kinematics(iExp).(fn).X(:);
            case "vel"
                t = kinematics(iExp).(fn).t;
                x = diff([NaN; kinematics(iExp).(fn).X(:)]) ./ diff([NaN; kinematics(iExp).(fn).t(:)]);
            case "spd"
                t = kinematics(iExp).(fn).t;
                x = abs(diff([NaN; kinematics(iExp).(fn).X(:)]) ./ diff([NaN; kinematics(iExp).(fn).t(:)]));
            otherwise
                error("Unknown argument p.kinematicDataSource=%s", p.kinematicDataSource)
        end

        % average across trials
        stimData(iExp).(fn).t = tLocal;
        stimData(iExp).(fn).X = NaN(length(stimData(iExp).hash), length(tLocal), 'single');
        for iStim = 1:length(stimData(iExp).hash)
            stimData(iExp).(fn).X(iStim, :) = interp1(t, x, tLocal + stimData(iExp).stimOn(iStim), 'previous');
        end
        stimData(iExp).(fn).XCtrl = NaN(length(stimData(iExp).stimCtrl), length(tLocal), 'single');
        for iStim = 1:length(stimData(iExp).stimCtrl)
            stimData(iExp).(fn).XCtrl(iStim, :) = interp1(t, x, tLocal + stimData(iExp).stimCtrl(iStim), 'previous');
        end
    end

    % Mark down trialtype
    eu0 = eu(expToEuIndices(iExp));
    trialType = repmat("unknown", [length(stimData(iExp).stimOn), 1]);
    trialType(eu0.Trials.Press.inTrial(stimData(iExp).stimOn)) = "press";
    trialType(eu0.Trials.Lick.inTrial(stimData(iExp).stimOn)) = "lick";
    trialTypeCtrl = repmat("unknown", [length(stimData(iExp).stimCtrl), 1]);
    trialTypeCtrl(eu0.Trials.Press.inTrial(stimData(iExp).stimCtrl)) = "press";
    trialTypeCtrl(eu0.Trials.Lick.inTrial(stimData(iExp).stimCtrl)) = "lick";
    stimData(iExp).trialType = categorical(trialType);
    stimData(iExp).trialTypeCtrl = categorical(trialTypeCtrl);
end

% Combine all sessions, store it in stimData(nSessions+1)
stimData = stimData(1:length(expToEuIndices));
for fn = ["power", "ain", "stimOn", "stimOff", "duration", "iPower", "iDuration", "hash", "stimCtrl", "trialType", "trialTypeCtrl"]
    X = {stimData(1:length(expToEuIndices)).(fn)};
    X = cat(1, X{:});
    stimData(length(expToEuIndices)+1).(fn) = X;
end
for fn = ["Jaw", "HandL", "HandR"]
    X = arrayfun(@(stimData) stimData.(fn).X, stimData(1:length(expToEuIndices)), UniformOutput=false);
    X = cat(1, X{:});
    XCtrl = arrayfun(@(stimData) stimData.(fn).XCtrl, stimData(1:length(expToEuIndices)), UniformOutput=false);
    XCtrl = cat(1, XCtrl{:});
    stimData(length(expToEuIndices)+1).(fn).X = X;
    stimData(length(expToEuIndices)+1).(fn).XCtrl = XCtrl;
    stimData(length(expToEuIndices)+1).(fn).t = stimData(1).(fn).t;
end
for trialType = ["press", "lick"]
    for cond = ["stim", "ctrl"]
        N = arrayfun(@(sd) sd.psmh.(cond).(trialType).N, stimData(1:length(stimData)-1), UniformOutput=false);
        stimData(length(expToEuIndices)+1).psmh.(cond).(trialType).N = cat(1, N{:});
        stimData(length(expToEuIndices)+1).psmh.(cond).(trialType).edges = stimData(1).psmh.(cond).(trialType).edges;
    end
end

% Try to merge sessions into one 3d array (trials x timestamps x units; block diagonal on the 1st and 3rd dimension, everything else NaN, so we could still average valid data and use easy trial indexing when splitting by conditions)
% Block diagonal, you think you're better than us don't you (That's prolly not even the right term for it)
for cond = ["stim", "ctrl"]
    nTrials = arrayfun(@(sd) size(sd.psth.(cond).X, 1), stimData(1:length(expToEuIndices)));
    nTimestamps = size(stimData(1).psth.(cond).X, 2);
    nUnits = arrayfun(@(sd) size(sd.psth.(cond).X, 3), stimData(1:length(expToEuIndices)));
    X = NaN(sum(nTrials), nTimestamps, sum(nUnits), 'single');
    iTrial = 0;
    iUnit = 0;
    for iExp = 1:length(expToEuIndices)
        X(iTrial+1:iTrial+nTrials(iExp), :, iUnit+1:iUnit+nUnits(iExp)) = stimData(iExp).psth.(cond).X;
        iTrial = iTrial + nTrials(iExp);
        iUnit = iUnit + nUnits(iExp);
    end
    stimData(length(expToEuIndices)+1).psth.(cond).X = X;
    if p.useWeakIncDec
        stimData(length(expToEuIndices)+1).psth.(cond).XIncPress = X(:, :, c.isPressUpWeak);
        stimData(length(expToEuIndices)+1).psth.(cond).XDecPress = X(:, :, c.isPressDownWeak);
        stimData(length(expToEuIndices)+1).psth.(cond).XIncLick = X(:, :, c.isLickUpWeak);
        stimData(length(expToEuIndices)+1).psth.(cond).XDecLick = X(:, :, c.isLickDownWeak);
    else
        stimData(length(expToEuIndices)+1).psth.(cond).XIncPress = X(:, :, c.isPressUp);
        stimData(length(expToEuIndices)+1).psth.(cond).XDecPress = X(:, :, c.isPressDown);
        stimData(length(expToEuIndices)+1).psth.(cond).XIncLick = X(:, :, c.isLickUp);
        stimData(length(expToEuIndices)+1).psth.(cond).XDecLick = X(:, :, c.isLickDown);
    end
    stimData(length(expToEuIndices)+1).psth.(cond).t = stimData(1).psth.(cond).t;
end

stimData(length(expToEuIndices)+1).iExp = 0;
stimData(length(expToEuIndices)+1).name = 'all sessions';

% Report trial counts
for iExp = 1:length(stimData)
    nTrials = histcounts(stimData(iExp).trialType, ["press", "lick", "unknown"]);
    nTrialsCtrl = histcounts(stimData(iExp).trialTypeCtrl, ["press", "lick", "unknown"]);
    fprintf("iExp=%i (%s):\n\tstim: %i reach, %i lick, %i unknown; ctrl: %i press, %i lick, %i unknown.\n", iExp, stimData(iExp).name, nTrials(1), nTrials(2), nTrials(3), nTrialsCtrl(1), nTrialsCtrl(2), nTrialsCtrl(3))
    if iExp < length(stimData)
        fprintf("\t%i units; reach: %i/%i inc/dec; lick: %i/%i inc/dec.\n", size(stimData(iExp).psth.stim.X, 3), size(stimData(iExp).psth.stim.XIncPress, 3), size(stimData(iExp).psth.stim.XDecPress, 3), size(stimData(iExp).psth.stim.XIncLick, 3), size(stimData(iExp).psth.stim.XDecLick, 3))
    end
end

for iExp = length(uniqueExpNames)+1 % nSessions+1 will plot session average
    x0 = 0;
    for trialType = ["press", "lick"]
        fig = figure(Units='normalized', Position=[x0, 0, 0.5, 1]);
        x0 = x0 + 0.5;
        tl = tiledlayout(fig, sum(l.h), sum(l.w), TileSpacing='compact');
    
        clear ax
        for ifn = 1:length(features)
            fn = char(features(ifn));
            % Plot STA
            for iDuration = 1:2
                ax = nexttile(tl, sum(l.w)*(ifn-1) + l.cw(iDuration), [l.h(ifn), l.w(iDuration)]);
                if ifn == 1
                    title(ax, sprintf('%gs opto', p.pulseDurations(iDuration)))
                end
                clear h
                iLine = 1;
                hold(ax, 'on')
    
                % Ctrk traces
                selCtrl = stimData(iExp).trialTypeCtrl == trialType;
                switch fn
                    case {'Jaw', 'HandL', 'HandR'}
                        mu = mean(stimData(iExp).(fn).XCtrl(selCtrl, :), 1, 'omitnan');
                        err = 0.1*std(stimData(iExp).(fn).XCtrl(selCtrl, :), 0, 1, 'omitnan');
                        ci = quantile(stimData(iExp).(fn).XCtrl(selCtrl, :), [0.25, 0.75], 1);
                        h(iLine) = plot(ax, tLocal, mu, Color=colorsByPower(1, :), LineWidth=1, DisplayName=sprintf("ctrl (n=%i)", nnz(selCtrl)));
                        iLine = iLine + 1;
                        patch(ax, [tLocal, flip(tLocal)], [ci(1, :), flip(ci(2, :))], colorsByPower(1, 1:3), FaceAlpha=0.1, EdgeAlpha=0)
                    case 'move'
                        edges = stimData(iExp).psmh.ctrl.(trialType).edges;
                        N = mean(stimData(iExp).psmh.ctrl.(trialType).N(selCtrl, :), 1, 'omitnan');
                        N(~isfinite(N)) = 0;
                        h(iLine) = histogram(ax, BinEdges=edges, BinCounts=N, DisplayStyle='stairs', EdgeColor=colorsByPower(1, 1:3), EdgeAlpha=0.5, LineWidth=1, DisplayName=sprintf("ctrl (n=%i)", nnz(selCtrl)));
                        iLine = iLine + 1;
                    case {'X', 'XInc', 'XDec'}
                        t = stimData(iExp).psth.ctrl.t;
                        switch fn
                            case 'X'
                                fnx = 'X';
                            case {'XInc', 'XDec'}
                                fnx = char(string(fn) + trialType);
                                fnx(5) = upper(fnx(5)); % XIncPress
                        end
                        if showIndividualTracesByPower(1)
                            switch p.showIndividualTracesFor
                                case "trials"
                                    mu = mean(stimData(iExp).psth.ctrl.(fnx)(selCtrl, :, :), 3, 'omitnan'); % avg across units: trials x time x units, then average across 3rd dim -> trials x time
                                    if ~isempty(mu)
                                        if p.showIndividualTracesAsCI
                                            ci = quantile(mu, [0.25, 0.75], 1);
                                            patch(ax, [t, flip(t)], [ci(1, :), flip(ci(2, :))], colorsByPower(1, 1:3), FaceAlpha=0.05, EdgeAlpha=0);
                                        else
                                            plot(ax, t, mu', Color=[colorsByPower(1, 1:3), individualTracesAlpha], LineWidth=0.5);
                                        end
                                    end
                                case "units"
                                    mu = mean(permute(stimData(iExp).psth.ctrl.(fnx)(selCtrl, :, :), [3, 2, 1]), 3, 'omitnan'); % avg across trials: units x time x trials, then average across 3rd dim -> units x time
                                    if ~isempty(mu)
                                        if p.showIndividualTracesAsCI
                                            ci = quantile(mu, [0.25, 0.75], 1);
                                            patch(ax, [t, flip(t)], [ci(1, :), flip(ci(2, :))], colorsByPower(1, 1:3), FaceAlpha=0.1, EdgeAlpha=0);
                                        else
                                            plot(ax, t, mu', Color=[colorsByPower(1, 1:3), individualTracesAlpha], LineWidth=0.5);
                                        end
                                    end
                            end
                        end

                        mumu = mean(stimData(iExp).psth.ctrl.(fnx)(selCtrl, :, :), [1, 3], 'omitnan');
                        h(iLine) = plot(ax, t, mumu, Color=colorsByPower(1, :), LineWidth=1, DisplayName=sprintf("ctrl (n=%i)", nnz(selCtrl)));
                        iLine = iLine + 1;
                end
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
                        case {'Jaw', 'HandL', 'HandR'}
                            mu = mean(stimData(iExp).(fn).X(sel, :), 1, 'omitnan');
                            err = 0.1*std(stimData(iExp).(fn).X(sel, :), 0, 1, 'omitnan');
                            ci = quantile(stimData(iExp).(fn).X(sel, :), [0.25, 0.75], 1);
                            h(iLine) = plot(ax, tLocal, mu, Color=colorsByPower(iPower+1, :), LineStyle=lineStylesByPower(iPower+1), LineWidth=1, DisplayName=sprintf("%s (n=%i)", label, nnz(sel)));
                            iLine = iLine + 1;
                            % patch(ax, [tLocal, flip(tLocal)], [mu-err, flip(mu+err)], colors(iPower, 1:3), FaceAlpha=0.1, EdgeAlpha=0)
                            patch(ax, [tLocal, flip(tLocal)], [ci(1, :), flip(ci(2, :))], colorsByPower(iPower+1, 1:3), FaceAlpha=0.1, EdgeAlpha=0)
                        case {'move'}
                            edges = stimData(iExp).psmh.stim.(trialType).edges;
                            N = mean(stimData(iExp).psmh.stim.(trialType).N(sel, :), 1, 'omitnan');
                            N(~isfinite(N)) = 0;
                            h(iLine) = histogram(ax, BinEdges=edges, BinCounts=N, DisplayStyle='stairs', EdgeColor=colorsByPower(iPower+1, 1:3), EdgeAlpha=colorsByPower(iPower+1, 4), LineWidth=1, DisplayName=sprintf("%s (n=%i)", label, nnz(sel)));
                        case {'X', 'XInc', 'XDec'} % Spike rate
                            t = stimData(iExp).psth.stim.t;
                            switch fn
                                case 'X'
                                    fnx = 'X';
                                case {'XInc', 'XDec'}
                                    fnx = char(string(fn) + trialType);
                                    fnx(5) = upper(fnx(5)); % XIncPress
                            end
                            if showIndividualTracesByPower(iPower+1)
                                switch p.showIndividualTracesFor
                                    case "trials"
                                        mu = mean(stimData(iExp).psth.stim.(fnx)(sel, :, :), 3, 'omitnan'); % avg across units: trials x time x units, then average across 3rd dim -> trials x time
                                        if ~isempty(mu)
                                            if p.showIndividualTracesAsCI
                                                ci = quantile(mu, [0.25, 0.75], 1);
                                                patch(ax, [t, flip(t)], [ci(1, :), flip(ci(2, :))], colorsByPower(iPower+1, 1:3), FaceAlpha=0.05, EdgeAlpha=0);
                                            else
                                                plot(ax, t, mu', Color=[colorsByPower(iPower+1, 1:3), individualTracesAlpha], LineWidth=0.5);
                                            end
                                        end
                                    case "units"
                                        mu = mean(permute(stimData(iExp).psth.stim.(fnx)(sel, :, :), [3, 2, 1]), 3, 'omitnan'); % avg across trials: units x time x trials, then average across 3rd dim -> units x time
                                        if ~isempty(mu)
                                            if p.showIndividualTracesAsCI
                                                ci = quantile(mu, [0.25, 0.75], 1);
                                                patch(ax, [t, flip(t)], [ci(1, :), flip(ci(2, :))], colorsByPower(iPower+1, 1:3), FaceAlpha=0.1, EdgeAlpha=0);
                                            else
                                                plot(ax, t, mu', Color=[colorsByPower(iPower+1, 1:3), individualTracesAlpha], LineWidth=0.5);
                                            end
                                        end
                                end
                            end

                            mumu = mean(stimData(iExp).psth.stim.(fnx)(sel, :, :), [1, 3], 'omitnan');
                            h(iLine) = plot(ax, t, mumu, Color=colorsByPower(iPower+1, :), LineWidth=1, DisplayName=sprintf("%s (n=%i)", label, nnz(sel)));
                            iLine = iLine + 1;
                            yline(ax, 0, 'k--')
                    end
                end
                xline(ax, durations, 'k--')
                if iDuration == 1
                    nTrialsPerExp = arrayfun(@(sd) length(sd.trialType), stimData(1:end-1));
                    expToFirstTrialIndex = cumsum([1, nTrialsPerExp(1:end-1)]);
                    if featureDispNames(ifn) == "move"
                        switch trialType
                            case "press"
                                ylabel(ax, ["bar-contact", featureUnits(ifn)])
                            case "lick"
                                ylabel(ax, ["spout-contact", featureUnits(ifn)])
                        end
                    elseif ismember(string(fn), ["XInc", "XDec"])
                        if iExp < length(stimData)
                            selExp = iExp;
                        else
                            selExp = 1:length(stimData)-1;
                        end
                        fnx = char(string(fn) + trialType);
                        fnx(5) = upper(fnx(5));
                        ylabel(ax, [featureDispNames(ifn), sprintf("n=%i %s", sum(arrayfun(@(sd) size(sd.psth.stim.(fnx), 3), stimData(selExp))), featureUnits(ifn))])
                    else
                        ylabel(ax, [featureDispNames(ifn), featureUnits(ifn)])
                    end
                end
                xlim(ax, xl{iDuration})
                ylim(ax, yl{ifn})
                if ifn == length(features)
                    xticks(ax, xl{iDuration}(1)+1:xl{iDuration}(2)-1)
                else
                    xticks(ax, [])
                end
                if iDuration == 1
                    yticks(ax, 'auto')
                else
                    yticks(ax, [])
                end
                if ifn == 1
                    lgd = legend(h, Location='northoutside', Orientation='horizontal'); lgd.ItemTokenSize = [9, 9];
                end
            end
        end
        trialTypeDispName = trialType; 
        if trialType == "press"
            trialTypeDispName = "reach";
        end
        title(tl, sprintf('Exp %i - %s - %s (n=%i)', iExp, stimData(iExp).name, trialTypeDispName, nnz(stimData(iExp).trialType==trialType)), Interpreter='none')
        xlabel(tl, 'time to decoder/opto onset (s)')
        fontsize(fig, 9, 'points')
        print(fig, fullfile(exportPath, sprintf("Exp %i - %s - %s.png", iExp, stimData(iExp).name, trialTypeDispName)), '-dpng')
    end
end

%% Plot opt responses as a heatmap
% close all
xl = {[-1, 3], [-1, 3], [-1, 5]}; % ctrl, 5mW 1s, 5mW 3s
l.w = cellfun(@diff, xl);
l.cw = cumsum([1, l.w]);
trialTypes = ["press", "lick"];
trialTypeDispName = ["reach", "lick"];

for iExp = length(uniqueExpNames) + 1
    fig = figure();
    tlp = tiledlayout(fig, 2, 1);
    tl = gobjects(2, 1);
    tl(1) = tiledlayout(tlp, 2, sum(l.w));
    tl(2) = tiledlayout(tlp, 2, sum(l.w));
    tl(1).Layout.Tile = 1; tl(1).Layout.TileSpan = [1, 1];
    tl(2).Layout.Tile = 2; tl(1).Layout.TileSpan = [1, 1];

    for itl = 1:length(trialTypes)
        trialType = trialTypes(itl);
        selCtrl = stimData(iExp).trialTypeCtrl == trialType;
        title(tl(itl), trialTypeDispName(itl))
        AX = gobjects(2, 3);
        for i = 1:2
            for j = 1:3
                AX(i, j) = nexttile(tl(itl), (i-1)*sum(l.w) + l.cw(j), [1, l.w(j)]);
                ax = AX(i, j);
            end
        end
        % AX(2, 1).Visible = false;
        
        % First find order, do not plot

        [uniqueHash, ia] = unique(stimData(iExp).hash);
        for iHash = 1:length(uniqueHash)
            sel = stimData(iExp).hash == uniqueHash(iHash) & stimData(iExp).trialType == trialType;
            iPower = stimData(iExp).iPower(ia(iHash));
            if stimData(iExp).iDuration(ia(iHash)) ~= 2 || iPower ~= 2
                continue
            end
            fnxDec = char("XDec" + trialType);
            fnxDec(5) = upper(fnxDec(5)); % XDecPress
            fnxInc = char("XInc" + trialType);
            fnxInc(5) = upper(fnxInc(5)); % XIncPress
            t = stimData(iExp).psth.stim.t;
            muDec = mean(permute(stimData(iExp).psth.stim.(fnxDec)(sel, :, :), [3, 2, 1]), 3, 'omitnan'); % avg across trials: units x time x trials, then average across 3rd dim -> units x time
            muInc = mean(permute(stimData(iExp).psth.stim.(fnxInc)(sel, :, :), [3, 2, 1]), 3, 'omitnan'); % avg across trials: units x time x trials, then average across 3rd dim -> units x time
            [~, orderDec] = sort(100*mean(muDec(:, isin(t, [0, 3])), 2, 'omitnan') + 1*mean(muDec(:, isin(t, [3, 4])), 2, 'omitnan'), 'ascend');
            [~, orderInc] = sort(100*mean(muInc(:, isin(t, [0, 3])), 2, 'omitnan') + 1*mean(muInc(:, isin(t, [3, 4])), 2, 'omitnan'), 'descend');
        end


        % Plot ctrl vs. 1s vs. 3s
        for iDuration = 0:2 % 0 is ctrl, 1 is 1s, 2 is 3s
            % ctrl
            if iDuration == 0
                for iPowerAirBunnies = 1:2
                    ax = AX(iPowerAirBunnies, 1);
                    label = "ctrl";
                    durations = 0;
                    fnxDec = char("XDec" + trialType);
                    fnxDec(5) = upper(fnxDec(5)); % XDecPress
                    fnxInc = char("XInc" + trialType);
                    fnxInc(5) = upper(fnxInc(5)); % XIncPress
                    muDec = mean(permute(stimData(iExp).psth.ctrl.(fnxDec)(selCtrl, :, :), [3, 2, 1]), 3, 'omitnan'); % avg across trials: units x time x trials, then average across 3rd dim -> units x time
                    muInc = mean(permute(stimData(iExp).psth.ctrl.(fnxInc)(selCtrl, :, :), [3, 2, 1]), 3, 'omitnan'); % avg across trials: units x time x trials, then average across 3rd dim -> units x time
                    mu = vertcat(muDec(orderDec, :), muInc(orderInc, :));
                    imagesc(ax, XData=stimData(iExp).psth.ctrl.t, CData=mu);
                    xline(ax, durations, 'k--')
                    yline(ax, size(muDec, 1)+0.5, 'k--', LineWidth=1)
                    xlim(ax, xl{iDuration+1})
                    ylim(ax, [0.5, size(mu, 1)+0.5])
                    yticks(ax, [1, size(muDec, 1), size(mu, 1)])
                    ax.YAxis.Direction = 'reverse';
                    applyCustomColormap(ax, [-1.5, 1.5], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);
                    title(ax, sprintf("%s", label))
                end
            else
                % stim (1s, 3s)
                [uniqueHash, ia] = unique(stimData(iExp).hash);
                durations = [0]; % For drawing xline at opto onset/offset
                for iHash = 1:length(uniqueHash)
                    sel = stimData(iExp).hash == uniqueHash(iHash) & stimData(iExp).trialType == trialType;
                    iPower = stimData(iExp).iPower(ia(iHash));
                    if stimData(iExp).iDuration(ia(iHash)) ~= iDuration
                        continue
                    end
                    durations = unique([durations, p.pulseDurations(iDuration)]);
                    label = sprintf("%gmw", p.laserPowers(iPower)*1e3);
    
                    ax = AX(iPower, iDuration+1);
                    fnxDec = char("XDec" + trialType);
                    fnxDec(5) = upper(fnxDec(5)); % XDecPress
                    fnxInc = char("XInc" + trialType);
                    fnxInc(5) = upper(fnxInc(5)); % XIncPress
                    muDec = mean(permute(stimData(iExp).psth.stim.(fnxDec)(sel, :, :), [3, 2, 1]), 3, 'omitnan'); % avg across trials: units x time x trials, then average across 3rd dim -> units x time
                    muInc = mean(permute(stimData(iExp).psth.stim.(fnxInc)(sel, :, :), [3, 2, 1]), 3, 'omitnan'); % avg across trials: units x time x trials, then average across 3rd dim -> units x time
                    mu = vertcat(muDec(orderDec, :), muInc(orderInc, :));
                    mu(isnan(mu)) = 0;
                    imagesc(ax, XData=stimData(iExp).psth.stim.t, CData=mu);
                    xline(ax, durations, 'k--')
                    yline(ax, size(muDec, 1)+0.5, 'k--', LineWidth=1)
                    xlim(ax, xl{iDuration+1})
                    ylim(ax, [0.5, size(mu, 1)+0.5])
                    yticks(ax, [1, size(muDec, 1), size(mu, 1)])
                    ax.YAxis.Direction = 'reverse';
                    applyCustomColormap(ax, [-1.5, 1.5], hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);

                    title(ax, sprintf("%s %is", label, durations(end)))
                end
            end
            xticks(ax, durations)
            ylabel(ax, 'units')
        end
        cb = colorbar(ax);
        cb.Layout.Tile = 'east';
        cb.Label.String = 'norm spike rate (a.u.)';
    end
    fontsize(fig, 9, 'points')
end


%% Plot whole session kinematics + spikeRates
% close all

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

p.plotIndividualUnits = false;
p.plotPopulationMeans = true;
p.kinematicDataSource = "pos"; % pos, vel
p.spikeDataSource = "rate"; % rate, count
p.spikeRes = 0.001;
p.spikeKernelType = 'exponential';
p.correctSpikeRateDrift = true; % Subtract a smoothed baseline spike rate
if isfield(p, 'blank')
    p = rmfield(p, 'blank');
end
p.blank(1).event = "StimOn";
p.blank(1).window = [-0.5, 0.5]*1e-3;
p.blank(1).event = "StimOff";
p.blank(1).window = [-0.5, 0.5]*1e-3;
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
trialType = "press"; iExp = 6;
% trialType = "press"; iExp = 5;
clear n ei
isInExp = euToExpIndices == iExp;
n.press.inc = nnz(isInExp & c.isPressUp(:));
n.press.dec = nnz(isInExp & c.isPressDown(:));
n.press.all = nnz(isInExp);
n.lick.inc = nnz(isInExp & c.isLickUp(:));
n.lick.dec = nnz(isInExp & c.isLickDown(:));
n.lick.all = nnz(isInExp);
ei.press.inc = find(isInExp & c.isPressUp(:));
ei.press.dec = find(isInExp & c.isPressDown(:));
ei.press.all = find(isInExp);
ei.lick.inc = find(isInExp & c.isLickUp(:));
ei.lick.dec = find(isInExp & c.isLickDown(:));
ei.lick.all = find(isInExp);
if n.(trialType).inc == 0 || n.(trialType).dec == 0
    error("Cannot plot exp %i for trialtype %s", iExp, trialType)
end
X = struct(dec=[], inc=[]);
x = X;
XBoot = struct(dec=[], inc=[]);
nBoot = 10;
for dir = ["dec", "inc", "all"]
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
    end
    x.(dir) = mean(X.(dir), 1, 'omitnan')';
end
X.dec = X.dec';
X.inc = X.inc';
X.all = X.all';

% Do a PCA on population activity
XMerge = double([X.dec, X.inc]);
XMerge = XMerge - mean(XMerge, 1, 'omitnan');
[coeff, score, ~, ~, explained, mu] = pca(XMerge);

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

features = ["Jaw", "HandR", "HandL", "SpikeRateWithOpto", "DecoderWithStimData"];
featureDispName = ["jaw", "r.hand", "l.hand", "spike rate", "decoder P(Move)"];
nPCs = 1;

l.h = [1, 1, 1, 2, 2];
l.ch = cumsum([1, l.h]);
fig = figure(Units='inches', Position=[1, 1, 10, 6], DefaultAxesFontSize=9);
tl = tiledlayout(fig, sum(l.h), 1, TileSpacing='none', Padding='tight');

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
        case "DecoderWithStimData"
            eu0 = eu(expToEuIndices(iExp));
            [tDec, PDec, XDec] = readDecoderData(eu0, p.decoderDataDelay, p.decoderSampleRate, p.decoderSmoothWindow);
            XDec = (XDec - mean(XDec, 'all', 'omitnan'))./std(XDec, 0, 'all', 'omitnan');

            if p.plotPopulationMeans
                iLine = iLine + 1;
                h(iLine) = plot(ax(i), tDec, XDec, 'k', DisplayName='decoderX', Clipping='on');
                iLine = iLine + 1;
                h(iLine) = plot(ax(i), tDec, PDec, 'r', DisplayName='decoderP', Clipping='on');
                iLine = iLine + 1;
                h(iLine) = plot(ax(i), t, x.all, 'k--', DisplayName=sprintf('all (n=%i)', n.(trialType).all), Clipping='on');
            end
            yline(ax(i), 0, 'k:')
            tOn = stimData(iExp).stimOn;
            tOff = stimData(iExp).stimOff;
            trialType = stimData(iExp).trialType;
            for iStim = 1:length(tOn)
                switch trialType(iStim)
                    case "press"
                        color = [0.8, 0.2, 0.2];
                    case "lick"
                        color = [0.2, 0.2, 0.8];
                    case "unknown"
                        color = [0.2, 0.2, 0.2];
                end
                patch(ax(i), [tOn(iStim), tOff(iStim), tOff(iStim), tOn(iStim)], [-5, -5, 5, 5], [0.2, 0.2, 0.8], FaceAlpha=0.1, EdgeColor=color, EdgeAlpha=0.5)
            end
            tOn = stimData(iExp).stimCtrl;
            tOff = stimData(iExp).stimCtrl + 3;
            trialType = stimData(iExp).trialTypeCtrl;
            for iStim = 1:length(tOn)
                switch trialType(iStim)
                    case "press"
                        color = [0.8, 0.2, 0.2];
                    case "lick"
                        color = [0.2, 0.2, 0.8];
                    case "unknown"
                        color = [0.2, 0.2, 0.2];
                end
                patch(ax(i), [tOn(iStim), tOff(iStim), tOff(iStim), tOn(iStim)], [-5, -5, 5, 5], [0.2, 0.2, 0.2], FaceAlpha=0.1, EdgeColor=color, EdgeAlpha=0.5)
            end
            tTimeoutStart = eu(expToEuIndices(iExp)).EventTimes.TIMEOUT_START;
            for iTrial = 1:length(tTimeoutStart)
                patch(ax(i), tTimeoutStart(iTrial) + [0, 0.1, 0.1, 0], [-5, -5, 5, 5], [0.8, 0.8, 0.2], FaceAlpha=0.1, EdgeColor=[0.8, 0.8, 0.2], EdgeAlpha=0.8)
            end    

            ylim(ax(i), [-3, 3])
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
linkaxes(ax, 'x');
% legend(h);
xlabel(tl, 'time (s)')
box(ax, 'off')

% Thanks AI for drawing pretty buttons! Would've been nicer if they were
% functional but the grad student made it work.
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

function [t, P, X] = readDecoderData(eu, delay, sampleRate, smoothwindow)
    % Parse the decoded data
    % 1. Open the binary file
    fid = fopen(fullfile('C:\SERVER', eu.getAnimalName(), eu.ExpName, 'DecodedData.bin'), 'rb');
    
    % 2. Read all contents as raw unsigned 8-bit bytes
    rawBytes = fread(fid, Inf, '*uint8');
    fclose(fid);
    
    % 3. Calculate total byte-length of one triplet
    % uint32 = 4 bytes, uint16 = 2 bytes, uint8 = 1 byte (Total = 7 bytes per triplet)
    bytesPerTriplet = 7;
    numTriplets = floor(length(rawBytes) / bytesPerTriplet);
    
    % Trim any trailing incomplete bytes
    rawBytes = rawBytes(1 : numTriplets * bytesPerTriplet);
    
    % 4. Reshape bytes into a 7-by-N matrix (column-major layout)
    mat = reshape(rawBytes, bytesPerTriplet, numTriplets);
    
    % 5. Extract and cast fields using byte slices
    I = zeros(1, size(mat, 2), 'uint32');
    X = zeros(1, size(mat, 2), 'uint16');
    P = zeros(1, size(mat, 2), 'uint8');
    for iEvent = 1:size(mat, 2)
        I(iEvent) = typecast(mat(1:4, iEvent), 'uint32'); % First 4 bytes
        X(iEvent) = typecast(mat(5:6, iEvent), 'uint16'); % Next 2 bytes
        P(iEvent)  = mat(7, iEvent);                       % Final 1 byte
    end
    X = single(X)/62235*200;
    P = single(P)/255;
    t = (single(I)-1)./30000;

    t = t - delay;
    P = smoothdata(P, 'gaussian', smoothwindow*sampleRate);
    X = smoothdata(X, 'gaussian', smoothwindow*sampleRate);
end