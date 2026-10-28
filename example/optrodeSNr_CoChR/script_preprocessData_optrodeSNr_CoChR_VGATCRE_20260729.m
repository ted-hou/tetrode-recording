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

pressTrialsWithoutStim = cell(size(eu));
lickTrialsWithoutStim = cell(size(eu));
for iEu = 1:length(eu)
    trials = eu(iEu).Trials.Press;
    [~, ~, I] = trials.inTrial(eu(iEu).EventTimes.LaserModBlueOn);
    trials(I) = [];
    [~, ~, I] = trials.inTrial(eu(iEu).EventTimes.LaserModBlueOff);
    trials(I) = [];
    pressTrialsWithoutStim{iEu} = trials;
    trials = eu(iEu).Trials.Lick;
    [~, ~, I] = trials.inTrial(eu(iEu).EventTimes.LaserModBlueOn);
    trials(I) = [];
    [~, ~, I] = trials.inTrial(eu(iEu).EventTimes.LaserModBlueOff);
    trials(I) = [];
    lickTrialsWithoutStim{iEu} = trials;
end

eta.press = eu.getETA('count', 'press', [-4, 4], resolution=0.1, trials=pressTrialsWithoutStim, alignTo='stop', normalize=[-4, -2], minTrialDuration=1, artifacts=artifacts);
eta.lick = eu.getETA('count', 'lick', [-4, 4], resolution=0.1, trials=lickTrialsWithoutStim, alignTo='stop', normalize=[-4, -2], minTrialDuration=1, artifacts=artifacts);

meta.stim = mean(eta.stim.X(:, isin(eta.stim.t, [0.05, 0.2])), 2, 'omitnan');
meta.press = mean(eta.press.X(:, isin(eta.press.t, [-0.3, 0])), 2, 'omitnan');
meta.lick = mean(eta.lick.X(:, isin(eta.lick.t, [-0.3, 0])), 2, 'omitnan');

%% Plot ETAs (reach vs. lick vs. opto; all sessions combined, or per session)
vs = who; 

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
        h(iLine) = plot(ax, eta.(trialType).t, mean(eta.(trialType).X(selUnits & selUp, :), 1, 'omitnan'), Color=[0.8, 0.2, 0.2], LineWidth=1.5, DisplayName=sprintf('%s-inc (n=%i)', selTypeDispName, nnz(selUnits & selUp)));
        iLine = iLine + 1;
        if any(selUnits & c.isPressDown)
            h(iLine) = plot(ax, eta.(trialType).t, mean(eta.(trialType).X(selUnits & selDown, :), 1, 'omitnan'), Color=[0.2, 0.2, 0.8], LineWidth=1.5, DisplayName=sprintf('%s-dec (n=%i)', selTypeDispName, nnz(selUnits & selDown)));
            iLine = iLine + 1;
        end
        h(iLine) = plot(ax, eta.(trialType).t, mean(eta.(trialType).X(selUnits & selFlat, :), 1, 'omitnan'), Color=[0.2, 0.2, 0.2], LineWidth=1.5, DisplayName=sprintf('%s-flat (n=%i)', selTypeDispName, nnz(selUnits & selFlat)));
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

clearvars('-except', vs{:});

%% Plot scattered METAs (reach/lick vs. opto)
vs = who; 

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

clearvars('-except', vs{:}, 'finalAnswer');
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

%% Get controls for decoder-stim (old method, by redoing decoder data)
% Note: online storage of decoderData has a weird delay ~0.5s
% stimData = stimData(1:length(expToEuIndices));
% tLocal = -4:0.01:3+2;
% if exist('decoderData', 'var') && isfield(decoderData, 'tSham')
%     for iExp = 1:length(stimData)
%         stimData(iExp).stimCtrl = decoderData(iExp).tSham(:);
%     end
% % Fallback method is bad, requires reading MATLAB-recorded decoder pMoves
% % but there's an unknow delay we need to correct by guessing
% else
%     p.decoderDataDelay = 0.66;
%     p.decoderSampleRate = 50; % 20ms intervals
%     p.decoderSmoothWindow = 0.1;
%     for iExp = 1:length(stimData)
%         eu0 = eu(expToEuIndices(iExp));
%         stimOn = eu0.EventTimes.LaserModBlueOn;
%         stimOff = eu0.EventTimes.LaserModBlueOff;
% 
%         [t, P, X] = readDecoderData(eu0, p.decoderDataDelay, p.decoderSampleRate, p.decoderSmoothWindow);
% 
%         % ax = axes(figure);
%         % hold(ax, 'on')
%         % plot(ax, t, (X - mean(X))./std(X, 0), 'k-')
%         % plot(ax, t, P, 'b-')
%         % ax.InteractionOptions.LimitsDimensions = "x";
% 
%         threshold = 0.3;
%         tStimCtrl = t(strfind(P >= threshold, [0, 1]) + 1);
%         timeoutStart = eu0.EventTimes.TIMEOUT_START;
%         trials = Trial(timeoutStart+2, tStimCtrl, 'first', stimOn);
%         tStimCtrl = [trials.Stop];
%         % [B, t, I] = trials.inTrial(tStimCtrl);
%         % [~, ia, ~] = unique(I);
%         % tStimCtrl = t(ia);
% 
%         for t0 = stimOn(:)'
%             tStimCtrl(isin(tStimCtrl, t0 + [-3, 3])) = [];
%         end
% 
%         stimData(iExp).stimCtrl = tStimCtrl(:);
%     end
% end

%% Get controls for decoder-stim (new method: use ArduinoConnection.AnalogOutputEvents)
stimData = stimData(1:length(expToEuIndices));
hasCtrl = [false, true, true, false, false, true, true];
for iExp = 1:length(stimData)
    if hasCtrl(iExp)
        eu0 = eu(expToEuIndices(iExp));
        ac = eu0.loadArduinoConnection();
        stimOn = eu0.EventTimes.LaserModBlueOn;
        stimOff = eu0.EventTimes.LaserModBlueOff;
        assert(size(ac.AnalogOutputEvents, 2) >= 3) % 3rd col is arduino millis(), which we use to align
        aout = ac.AnalogOutputEvents;
        aout(isnan(aout(:, 1)), :) = [];
        aout(:, 5) = eu0.alignTimestamps(aout(:, 3)/1000, acRefEventName="TIMEOUT_END", euRefEventName="TimeoutOff");
        aout = aout(aout(:, 5)>0, :); % Remove calibration/validation events
        % Remove 2nd laser
        aout = aout(aout(:, 1)==1, :);
    
        dPwr = diff([aout(:, 2); Inf]);
        dur = diff([aout(:, 5); Inf]);
        dur = round(dur*10)/10;
    
        selCtrl = (dur==3) & dPwr == 0;
    
        stimData(iExp).stimCtrl = aout(selCtrl, 5);
    else
        stimData(iExp).stimCtrl = [];
    end
end

%% Calulate a peri-stim movement histogram
edges = -2:0.1:5;
for iExp = 1:length(stimData)
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
