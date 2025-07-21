%% Load EphysUnits
% clear files
% files{1} = dir("C:\SERVER\Units\Lite_NonDuplicate_NonDrift\Daisy2_20180514_*.mat");
% files{2} = dir("C:\SERVER\Units\Lite_NonDuplicate_NonDrift\daisy12_20220106_*.mat");
% files{3} = dir("C:\SERVER\Units\Lite_NonDuplicate_NonDrift\daisy14_20220506_*.mat");
% files = cat(1, files{:});
% files = arrayfun(@(f) sprintf('%s\\%s', f.folder, f.name), files, UniformOutput=false);
% eu = EphysUnit.load(files, waveforms=false, spikecounts=false, spikerates=false);

eu = EphysUnit.load('C:\SERVER\Units\Lite_NonDuplicate_NonDrift', waveforms=false, spikecounts=false, spikerates=false);
load('C:\SERVER\Units\meta_Lite_NonDuplicate_NonDrift_20250705.mat')
eu = eu(c.hasPress & c.hasLick);
clearvars -except eu

%% Load euComplete (complete with ITI spikes)
euNames = lower(eu.getName());

files = dir('C:\SERVER\Units\NonLite_PressVsLick\*.mat');
sel = ismember(cellfun(@(n) lower(strrep(n, '.mat', '')), {files.name}, UniformOutput=false), euNames);
files = files(sel);
cd('C:\SERVER\Units\NonLite_PressVsLick\')
euComplete = EphysUnit.load({files.name}, waveforms=false, spikecounts=false, spikerates=false);

% euComplete.save('C:\SERVER\Units\NonLite_PressVsLick_NonDuplicate_NonDrift');

[lia, locb] = ismember(eu.getName(), euComplete.getName());
eu(lia) = euComplete(locb(lia));

clear euComplete files sel lia locb


%% Copy digital events from ArduinoConnection to EphysUnit (using CueOn to correct for clock drift)

clear ac
[ac.ac, ac.euIndicesFirstInSession, ac.expIndices, ac.uniqueExpNames] = eu.loadArduinoConnection();

% Find and correct bad lick labels
for iEu = 1:length(eu)
    eu(iEu).EventTimes.LICK = [];
    eu(iEu).EventTimes.LICK_OFF = [];
end
eu.alignTimestamps(["LICK", "LICK_OFF", "REWARD_ON"], ac={ac.ac, ac.euIndicesFirstInSession, ac.expIndices, ac.uniqueExpNames});
clear nLicks
nLicks.EUvAC = arrayfun(@(eu) [length(eu.EventTimes.Lick), length(eu.EventTimes.LICK)], eu, UniformOutput=false);
nLicks.EUvAC = cat(1, nLicks.EUvAC{:});
nLicks.EUvReward = arrayfun(@(eu) [length(eu.EventTimes.Lick), length(eu.EventTimes.RewardTimes)], eu, UniformOutput=false);
nLicks.EUvReward = cat(1, nLicks.EUvReward{:});
lickIsReward = diff(nLicks.EUvReward, 1, 2) == 0; lickIsReward = lickIsReward(:)';
for iEu = find(lickIsReward)
    eu(iEu).EventTimes.Lick = eu(iEu).EventTimes.LICK;
    eu(iEu).Trials.Lick = Trial(eu(iEu).EventTimes.Cue, eu(iEu).EventTimes.LICK, 'first', eu(iEu).EventTimes.Press);
end

eu.alignTimestamps(["LEVER_RETRACT_START", "LEVER_RETRACTED", "LEVER_RETRACT_END", "LEVER_DEPLOY_START", "LEVER_DEPLOYED", "LEVER_DEPLOY_END", "TUBE_RETRACT_START", "TUBE_RETRACT_END", "TUBE_RETRACTED", "TUBE_DEPLOY_START", "TUBE_DEPLOY_END", "TUBE_DEPLOYED"], ac={ac.ac, ac.euIndicesFirstInSession, ac.expIndices, ac.uniqueExpNames});
for iEu = 1:length(eu)
    try
        if isempty(eu(iEu).EventTimes.LEVER_DEPLOY_START)
            eu(iEu).EventTimes.LEVER_DEPLOY_START = eu(iEu).EventTimes.LEVER_DEPLOYED;
        end
        if isempty(eu(iEu).EventTimes.LEVER_RETRACT_START)
            eu(iEu).EventTimes.LEVER_RETRACT_START = eu(iEu).EventTimes.LEVER_RETRACTED;
        end
        if isempty(eu(iEu).EventTimes.TUBE_DEPLOY_START)
            eu(iEu).EventTimes.TUBE_DEPLOY_START = eu(iEu).EventTimes.TUBE_DEPLOYED;
        end
        if isempty(eu(iEu).EventTimes.TUBE_RETRACT_START)
            eu(iEu).EventTimes.TUBE_RETRACT_START = eu(iEu).EventTimes.TUBE_RETRACTED;
        end
    catch
        warning('Counld not process iEu = %i, "%s"', iEu, eu(iEu).getName())
    end
end

% Sanitize events (timestamps must be positive)
for iEu = 1:length(eu)
    try
        for event = ["LICK", "LICK_OFF", "LEVER_RETRACT_START", "LEVER_RETRACT_END", "LEVER_DEPLOY_START", "LEVER_DEPLOY_END", "TUBE_RETRACT_START", "TUBE_RETRACT_END", "TUBE_DEPLOY_START", "TUBE_DEPLOY_END"]
            e = eu(iEu).EventTimes.(event);
            e(e<0) = [];
            eu(iEu).EventTimes.(event) = e;
        end
    catch
        warning('Counld not process iEu = %i, "%s"', iEu, eu(iEu).getName())
    end
end
clear iEu event

% Make trials
for iEu = 1:length(eu)
    try
        if ~isempty(eu(iEu).EventTimes.LEVER_DEPLOY_END)
            eu(iEu).Trials.LeverDeploy = Trial(eu(iEu).EventTimes.LEVER_DEPLOY_START, eu(iEu).EventTimes.LEVER_DEPLOY_END);
        else
            eu(iEu).Trials.LeverDeploy = Trial(eu(iEu).EventTimes.LEVER_DEPLOY_START, 1 + eu(iEu).EventTimes.LEVER_DEPLOY_START);
        end
        if ~isempty(eu(iEu).EventTimes.LEVER_RETRACT_END)
            eu(iEu).Trials.LeverRetract = Trial(eu(iEu).EventTimes.LEVER_RETRACT_START, eu(iEu).EventTimes.LEVER_RETRACT_END);
        else
            eu(iEu).Trials.LeverRetract = Trial(eu(iEu).EventTimes.LEVER_RETRACT_START, 1 + eu(iEu).EventTimes.LEVER_RETRACT_START);
        end
        if ~isempty(eu(iEu).EventTimes.TUBE_DEPLOY_END)
            eu(iEu).Trials.TubeDeploy = Trial(eu(iEu).EventTimes.TUBE_DEPLOY_START, eu(iEu).EventTimes.TUBE_DEPLOY_END);
        else
            eu(iEu).Trials.TubeDeploy = Trial(eu(iEu).EventTimes.TUBE_DEPLOY_START, 1 + eu(iEu).EventTimes.TUBE_DEPLOY_START);
        end
        if ~isempty(eu(iEu).EventTimes.TUBE_RETRACT_END)
            eu(iEu).Trials.TubeRetract = Trial(eu(iEu).EventTimes.TUBE_RETRACT_START, eu(iEu).EventTimes.TUBE_RETRACT_END);
        else
            eu(iEu).Trials.TubeRetract = Trial(eu(iEu).EventTimes.TUBE_RETRACT_START, 1 + eu(iEu).EventTimes.TUBE_RETRACT_START);
        end
        maxLickInterval = 0.2;
        minLickInterval = 0.05;
        trialsCircLick = eu(iEu).getTrials('circlick');
        trialsCircLick = trialsCircLick(trialsCircLick.duration <= maxLickInterval & trialsCircLick.duration >= minLickInterval);
        eu(iEu).Trials.CircLick = trialsCircLick;
    catch
        warning('Counld not process iEu = %i, "%s"', iEu, eu(iEu).getName())
    end
end
clear trialsCircLick;
clear e iEu lickIsReward minLickInterval maxLickInterval nLicks

%% Clear some RAM
% for iEu = 1:length(eu)
%     eu(iEu).SpikeTimes = [];
% end
% clear iEu

%% Load raw data and save
if ~exist('C:\SERVER\Units\Lite_NonDuplicate_NonDrift\raw', 'dir')
    mkdir('C:\SERVER\Units\Lite_NonDuplicate_NonDrift\raw')
end
for iEu = 1:length(eu)
    % if ismember(eu(iEu).getAnimalName(), {'daisy14', 'daisy15', 'daisy16', 'desmond23', 'desmond24', 'desmond25', 'desmond26', 'desmond27'}) && ~ismember(eu(iEu).ExpName, {'daisy14_20220506', 'daisy16_20220502'})
    clear data t B selT
    try
        [data, t] = eu(iEu).loadRaw();
        selT = false(size(t));
        for trial = ["LeverDeploy", "LeverRetract", "TubeDeploy", "TubeRetract"]
            trials = eu(iEu).Trials.(trial);
            if isempty(trials)
                continue
            end
            B = trials.inTrial(t, window=[-1, 1], windowMode='extend');
            selT = selT | B;
        end

        for trial = ["CircLick"]
            trials = eu(iEu).Trials.(trial);
            if isempty(trials)
                continue
            end
            B = trials.inTrial(t);
            selT = selT | B;
        end

        for trial = ["Press", "Lick"]
            trials = eu(iEu).Trials.(trial);
            if isempty(trials)
                continue
            end
            B = trials.inTrial(t, window=[0, 2], windowMode='extend');
            selT = selT | B;
        end

        data = data(:, selT);
        t = t(:, selT);
        trials = eu(iEu).Trials;
        name = eu(iEu).getName();
        save(sprintf('C:\\SERVER\\Units\\Lite_NonDuplicate_NonDrift\\raw\\raw_%s.mat', eu(iEu).getName()), 'data', 't', 'trials', 'name')
    catch
        warning('Could not process unit %i, %s', iEu, eu(iEu).getName())
    end
    % end
end
clear data t B selT trial trials name

%% Find out which units failed raw data extraction
clear names

names.expected = string(eu.getName()); 
names.expected = names.expected(:);

files = dir("C:\SERVER\Units\Lite_NonDuplicate_NonDrift\raw\raw_*.mat");
names.found = string({files.name}); 
names.found = arrayfun(@(name) strsplit(name, ".mat"), names.found, UniformOutput=false);
names.found = cellfun(@(tokens) tokens{1}, names.found, UniformOutput=false);
names.found = cellfun(@(name) strsplit(name, "_"), names.found, UniformOutput=false);
names.found = cellfun(@(tokens) strjoin(tokens(2:end), "_"), names.found, UniformOutput=false);
names.found = string(names.found(:));

hasRaw = ismember(names.expected, names.found);
clear files

%% Find out which units use intan
isIntan = false(length(eu), 1);
for iEu = 1:length(eu)
    animalName = eu(iEu).getAnimalName();
    expName = eu(iEu).ExpName;
    filesRoot = dir(sprintf("C:\\SERVER\\%s\\%s\\*.rhd", animalName, expName));
    files = dir(sprintf("C:\\SERVER\\%s\\%s\\*\\*.rhd", animalName, expName));
    if ~isempty(files) || ~isempty(filesRoot)
        isIntan(iEu) = true;
    else
        isIntan(iEu) = false;
    end
    % fprintf('%s %i\n', expName, length(files));
end
clear iEu animalName expName files


%% Make behavioral trials and PETH/RD
tEU = eu.alignTimestamps(["REWARD_ON", "LEVER_RELEASED", "LEVER_RETRACT_START", "LEVER_RETRACTED", "LEVER_RETRACT_END", "TUBE_RETRACT_START", "TUBE_RETRACT_END", "LICK", "LICK_OFF"], ac={ac.ac, ac.euIndicesFirstInSession, ac.expIndices, ac.uniqueExpNames});

for iEu = 1:length(eu)
    try
        leverReleaseTimes = eu(iEu).EventTimes.LEVER_RELEASED;
        leverRetractTimes = eu(iEu).EventTimes.LEVER_RETRACT_START;
        if isempty(leverRetractTimes)
            leverRetractTimes = eu(iEu).EventTimes.LEVER_RETRACTED;
        end
        eu(iEu).Trials.Press = Trial(eu(iEu).EventTimes.Cue, eu(iEu).EventTimes.Press, 'first');
        eu(iEu).Trials.PressIncorrect = eu(iEu).Trials.Press(eu(iEu).Trials.Press.duration() < 4 & eu(iEu).Trials.Press.duration() >= 2);
        eu(iEu).Trials.PressCorrect = eu(iEu).Trials.Press(eu(iEu).Trials.Press.duration() >= 4);
        [~, leverRetractTimesIncorrect] = eu(iEu).Trials.PressIncorrect.inTrial(leverRetractTimes, [0, 2], windowMode='stop');
        assert(nnz(leverRetractTimesIncorrect) > 0)
        [~, leverRetractTimesCorrect] = eu(iEu).Trials.PressCorrect.inTrial(leverRetractTimes, [-0.1, 8], windowMode='stop');
        eu(iEu).Trials.RetractReleaseIncorrect = Trial([leverRetractTimesIncorrect, Inf], leverReleaseTimes, 'first');
        eu(iEu).Trials.RetractReleaseCorrect = Trial([leverRetractTimesCorrect, Inf], leverReleaseTimes, 'first');
        eu(iEu).Trials.PressReleaseIncorrect = Trial([eu(iEu).Trials.PressIncorrect.Start, Inf], [eu(iEu).Trials.RetractReleaseIncorrect.Stop], 'first');
        eu(iEu).Trials.PressReleaseCorrect = Trial([eu(iEu).Trials.PressCorrect.Stop, Inf], [eu(iEu).Trials.RetractReleaseCorrect.Stop], 'first');
        eu(iEu).Trials.PressRetractIncorrect = Trial([eu(iEu).Trials.PressIncorrect.Stop], leverRetractTimesIncorrect, 'first');
        eu(iEu).Trials.PressRetractCorrect = Trial([eu(iEu).Trials.PressCorrect.Stop], leverRetractTimesCorrect, 'first');
        % 
        % fprintf(['press=%i\npressIncorrect=%i, pressCorrect=%i;\nleverRetractTimesIncorrect=%i, leverRetractTimesCorrect=%i;\n' ...
        %     'eu(iEu).Trials.RetractReleaseIncorrect=%i, eu(iEu).Trials.RetractReleaseCorrect=%i;\n' ...
        %     'eu(iEu).Trials.PressReleaseIncorrect=%i, eu(iEu).Trials.PressReleaseCorrect=%i\n'], length(eu(iEu).Trials.Press), length(eu(iEu).Trials.PressIncorrect), length(eu(iEu).Trials.PressCorrect), ...
        %     length(leverRetractTimesIncorrect), length(leverRetractTimesCorrect), ...
        %     length(eu(iEu).Trials.RetractReleaseIncorrect), length(eu(iEu).Trials.RetractReleaseCorrect), ...
        %     length(eu(iEu).Trials.PressReleaseIncorrect), length(eu(iEu).Trials.PressReleaseCorrect));
        % disp(1)
    catch
        warning('Counld not process iEu = %i, "%s"', iEu, eu(iEu).getName())
    end
end

% Find and correct bad lick labels
for iEu = 1:length(eu)
    eu(iEu).EventTimes.LICK = [];
    eu(iEu).EventTimes.LICK_OFF = [];
end
eu.alignTimestamps(["LICK", "LICK_OFF"], ac={ac.ac, ac.euIndicesFirstInSession, ac.expIndices, ac.uniqueExpNames});
clear nLicks
nLicks.EUvAC = arrayfun(@(eu) [length(eu.EventTimes.Lick), length(eu.EventTimes.LICK)], eu, UniformOutput=false);
nLicks.EUvAC = cat(1, nLicks.EUvAC{:});
nLicks.EUvReward = arrayfun(@(eu) [length(eu.EventTimes.Lick), length(eu.EventTimes.RewardTimes)], eu, UniformOutput=false);
nLicks.EUvReward = cat(1, nLicks.EUvReward{:});

lickIsReward = diff(nLicks.EUvReward, 1, 2) == 0; lickIsReward = lickIsReward(:)';

for iEu = find(lickIsReward)
    eu(iEu).Trials.Lick = Trial(eu(iEu).EventTimes.Cue, eu(iEu).EventTimes.LICK, 'first', eu(iEu).EventTimes.Press);
end

for iEu = 1:length(eu)
    eu(iEu).Trials.LickIncorrect = eu(iEu).Trials.Lick(eu(iEu).Trials.Lick.duration() < 4 & eu(iEu).Trials.Lick.duration() >= 2);
    eu(iEu).Trials.LickCorrect = eu(iEu).Trials.Lick(eu(iEu).Trials.Lick.duration() >= 4);
end

% Make rd and eta
clear rd
rd.press = eu.getRasterData('press', [-4, 4], minTrialDuration=2, maxTrialDuration=Inf);
rd.releaseCorrect = eu.getRasterData('press_release_correct', [-4, 4], minTrialDuration=0, maxTrialDuration=Inf);
rd.releaseIncorrect = eu.getRasterData('press_release_incorrect', [-4, 4], minTrialDuration=0, maxTrialDuration=Inf);
rd.retractCorrect = eu.getRasterData('press_retract_correct', [-4, 4], minTrialDuration=0, maxTrialDuration=Inf);
rd.retractIncorrect = eu.getRasterData('press_retract_incorrect', [-4, 4], minTrialDuration=0, maxTrialDuration=Inf);
rd.lick = eu.getRasterData('lick', [-4, 4], minTrialDuration=2, maxTrialDuration=Inf);

eta.correctPress = eu.getETA('count', 'press', [-4, 4], minTrialDuration=4, normalize='none', resolution=0.1);
eta.incorrectPress = eu.getETA('count', 'press', [-4, 4], minTrialDuration=2, maxTrialDuration=4, normalize='none', resolution=0.1);

eta.correctRetract = eu.getETA('count', 'press_retract_correct', [-4, 4], alignTo='stop', normalize='none', resolution=0.1);
eta.incorrectRetract = eu.getETA('count', 'press_retract_incorrect', [-4, 4], alignTo='stop', normalize='none', resolution=0.1);

eta.correctRelease = eu.getETA('count', 'press_release_correct', [-4, 4], alignTo='stop', normalize='none', resolution=0.1);
eta.incorrectRelease = eu.getETA('count', 'press_release_incorrect', [-4, 4], alignTo='stop', normalize='none', resolution=0.1);

eta.correctLick = eu.getETA('count', 'lick', [-4, 4], minTrialDuration=4, normalize='none', resolution=0.1);
eta.incorrectLick= eu.getETA('count', 'lick', [-4, 4], minTrialDuration=2, maxTrialDuration=4, normalize='none', resolution=0.1);

for field = ["correctPress", "incorrectPress", "correctLick", "incorrectLick", "correctRelease", "incorrectRelease", "correctRetract", "incorrectRetract"]
    eta.(field).X = eta.(field).X ./ 0.1;
end



%% Plot reach raster, peth, DLC trajectories for reach-dec DLC cells (among 163)
% close all
if ~exist('E:\DATA\Figures\reach_retract_dlc', 'dir')
    mkdir('E:\DATA\Figures\reach_retract_dlc')
end

fig = figure(Units='normalized', InnerPosition=[0 0 1 1]);
tl = tiledlayout(fig, 5, 3, TileIndexing='columnmajor');
ax = gobjects(15, 1);
for iAx = 1:15
    ax(iAx) = nexttile(tl);
end
% for iEu = find(hasRaw(:)' & isIntan(:)')
for iEu = find(hasRaw(:)')
    try
        raw = load(sprintf("C:\\SERVER\\Units\\Lite_NonDuplicate_NonDrift\\raw\\raw_%s.mat", eu(iEu).getName()));
        % 1: Raster
        for iAx = 1:15
            cla(ax(iAx), 'reset');
        end
        EphysUnit.plotRaster(ax(1), rd.press(iEu), xlim=[-4, 4]);
        set(ax(1).Legend, AutoUpdate=false)
        xline(ax(1), 0, 'k-', LineWidth=1.5)
        yline(ax(1), find(rd.press(iEu).duration>4, 1), 'k-', LineWidth=3)
        delete(ax(1).Legend)
        legend(ax(1), ["", "trial start"])
    
        % 2: PETH
        hold(ax(2), 'on')
        h = gobjects(2, 1);
        h(1) = plot(ax(2), eta.incorrectPress.t, eta.incorrectPress.X(iEu, :), Color=hsl2rgb([0.4, 0.4, 0.4]), LineStyle='-', DisplayName=sprintf('Incorrect (n=%i)', eta.incorrectPress.N(iEu)), LineWidth=2);
        h(2) = plot(ax(2), eta.correctPress.t, eta.correctPress.X(iEu, :), 'k-', DisplayName=sprintf('Correct (n=%i)', eta.correctPress.N(iEu)), LineWidth=2);
        legend(h, AutoUpdate=false, Location='southwest')
        xline(ax(2), 0, 'k--')
        yline(ax(2), 0, 'k--')
        xlim(ax(2), [-4, 4])
        % ylim(ax(2), [0, 100])
    
        title(ax(1), 'Reach Raster')
        title(ax(2), 'Reach PETH')

        % 3: Reach incorrect (raw)
        % [xAligned, tAligned] = eu(iEu).getTrialAlignedData(raw.data, raw.t, trials=eu(iEu).Trials.PressIncorrect, alignTo='stop', window=[-1, 1], resolution=1/30000);
        plotRaw(ax(3), eu(iEu), raw, eu(iEu).Trials.PressIncorrect, alignTo='Stop', window=[-0.5, 1.5]);
        title(ax(3), 'PressIncorrect')

        % 4: Reach correct (raw)
        plotRaw(ax(4), eu(iEu), raw, eu(iEu).Trials.PressCorrect, alignTo='Stop', window=[-0.5, 1.5]);
        title(ax(4), 'PressCorrect')

        % 5: Tube servo deploy (raw)
        plotRaw(ax(5), eu(iEu), raw, eu(iEu).Trials.TubeDeploy, alignTo='Start', window=[-0.5, 1.5]);
        title(ax(5), 'TubeDeploy')


        % 6: Raster
        EphysUnit.plotRaster(ax(6), rd.retractCorrect(iEu), xlim=[-4, 4]);
        set(ax(6).Legend, AutoUpdate=false)
        xline(ax(6), 0, 'k-', LineWidth=1.5)
        % yline(ax(6), find(rd.lick(i).duration>4, 1), 'k-', LineWidth=1.5)
        delete(ax(6).Legend)
        legend(ax(6), ["", "bar-contact"])
    
        % 7 PETH
        hold(ax(7), 'on')
        h = gobjects(2, 1);
        h(1) = plot(ax(7), eta.incorrectRetract.t, eta.incorrectRetract.X(iEu, :), Color=hsl2rgb([0.4, 0.4, 0.4]), LineStyle='-', DisplayName=sprintf('Incorrect (n=%i)', eta.incorrectRetract.N(iEu)), LineWidth=2);
        h(2) = plot(ax(7), eta.correctRetract.t, eta.correctRetract.X(iEu, :), 'k-', DisplayName=sprintf('Correct (n=%i)', eta.correctRetract.N(iEu)), LineWidth=2);
        legend(h, AutoUpdate=false, Location='southwest')
        xline(ax(7), 0, 'k--')
        yline(ax(7), 0, 'k--')
        xlim(ax(7), [-4, 4])
        % ylim(ax(7), [0, 100])
    
        % 8: reach retract (raw)
        plotRaw(ax(8), eu(iEu), raw, eu(iEu).Trials.PressRetractCorrect, alignTo='Stop', window=[-0.5, 1.5]);
        title(ax(8), 'PressRetractCorrect')

        % 9: lever servo deploy (raw)
        plotRaw(ax(9), eu(iEu), raw, eu(iEu).Trials.LeverDeploy, alignTo='Start', window=[-0.5, 1.5]);
        title(ax(9), 'LeverDeploy')

        % 10: lever servo retract (raw)
        plotRaw(ax(10), eu(iEu), raw, eu(iEu).Trials.LeverRetract, alignTo='Start', window=[-0.5, 1.5]);
        title(ax(10), 'LeverRetract')
    
        % 11: Raster
        EphysUnit.plotRaster(ax(11), rd.lick(iEu), xlim=[-4, 4]);
        set(ax(11).Legend, AutoUpdate=false)
        xline(ax(11), 0, 'k-', LineWidth=1.5)
        yline(ax(11), find(rd.lick(iEu).duration>4, 1), 'k-', LineWidth=3)
        delete(ax(11).Legend)
        legend(ax(11), ["", "trial start"])
    
        % 12 PETH
        hold(ax(12), 'on')
        h = gobjects(2, 1);
        h(1) = plot(ax(12), eta.incorrectLick.t, eta.incorrectLick.X(iEu, :), Color=hsl2rgb([0.4, 0.4, 0.4]), LineStyle='-', DisplayName=sprintf('Incorrect (n=%i)', eta.incorrectLick.N(iEu)), LineWidth=2);
        h(2) = plot(ax(12), eta.correctLick.t, eta.correctLick.X(iEu, :), 'k-', DisplayName=sprintf('Correct (n=%i)', eta.correctLick.N(iEu)), LineWidth=2);
        legend(h, AutoUpdate=false, Location='southwest')
        xline(ax(12), 0, 'k--')
        yline(ax(12), 0, 'k--')
        xlim(ax(12), [-4, 4])
        % ylim(ax(12), [0, 100])
    
        % 13: Lick incorrect (raw)
        plotRaw(ax(13), eu(iEu), raw, eu(iEu).Trials.LickIncorrect, alignTo='Stop', window=[-0.5, 1.5]);
        title(ax(13), 'LickIncorrect')

        % 14: Lick correct (raw)
        plotRaw(ax(14), eu(iEu), raw, eu(iEu).Trials.LickCorrect, alignTo='Stop', window=[-0.5, 1.5]);
        title(ax(14), 'LickCorrect')

        % 15: Tube servo retract (raw)
        plotRaw(ax(15), eu(iEu), raw, eu(iEu).Trials.TubeRetract, alignTo='Start', window=[-0.5, 1.5]);
        title(ax(15), 'TubeRetract')
    
        title(ax(1), 'Reach Raster')
        title(ax(2), 'Reach PETH')
        title(ax(6), 'Retract Raster (correct trials)')
        title(ax(7), 'Retract PETH')
        title(ax(11), 'Lick Raster')
        title(ax(12), 'Lick PETH')

        title(tl, rd.press(iEu).name, Interpreter='none')
    
        xlim(ax, [-3, 3])
        xlim(ax([3:5, 8:10, 13:15]), [-500, 1500])
    
        ylim(ax([2, 7, 12]), [min([ax(2).YLim, ax(7).YLim, ax(12).YLim]), max([ax(2).YLim, ax(7).YLim, ax(12).YLim])]);
        print(fig, sprintf('E:\\DATA\\Figures\\reach_retract_dlc\\%s.png', eu(iEu).getName()), '-dpng')
    catch
        warning('Error processing unit %i', iEu)
    end
end

%% Plot a few examples
iEu = find(strcmpi(eu.getName(), 'desmond25_20220430_Channel12_Unit1')); % Little lick noise
raw = load(sprintf("C:\\SERVER\\Units\\Lite_NonDuplicate_NonDrift\\raw\\raw_%s.mat", eu(iEu).getName()));
plotRaw(axes(figure), eu(iEu), raw, eu(iEu).Trials.Lick, alignTo='Stop', plotSpikes=true, spacing=400)
xlim([-200, 500])
plotRaw(axes(figure), eu(iEu), raw, eu(iEu).Trials.Lick, alignTo='Stop', plotSpikes=true, spacing=400)
xlim([-25, 100])
save(sprintf("E:\\Data\\%s_3trials.mat", eu(iEu).getName()), 'x', 't', 'st')

%% Big SNR unit with big lick artefact ****** send to Sophie/Prerau
iEu = find(strcmpi(eu.getName(), 'desmond26_20220531_Channel4_Unit1')); % Longer lick noise
raw = load(sprintf("C:\\SERVER\\Units\\Lite_NonDuplicate_NonDrift\\raw\\raw_%s.mat", eu(iEu).getName()));
[x, t, st] = plotRaw(axes(figure), eu(iEu), raw, eu(iEu).Trials.Lick, alignTo='Stop', plotSpikes=true, spacing=400);
xlim([-200, 500])
plotRaw(axes(figure), eu(iEu), raw, eu(iEu).Trials.Lick, alignTo='Stop', plotSpikes=true, spacing=400)
xlim([-25, 100])
save(sprintf("E:\\Data\\%s_3trials.mat", eu(iEu).getName()), 'x', 't', 'st')

%% Big SNR unit with big lick artefact
iEu = find(strcmpi(eu.getName(), 'desmond26_20220531_Channel37_Unit1')); % Longer lick noise
raw = load(sprintf("C:\\SERVER\\Units\\Lite_NonDuplicate_NonDrift\\raw\\raw_%s.mat", eu(iEu).getName()));
[x, t, st] = plotRaw(axes(figure), eu(iEu), raw, eu(iEu).Trials.Lick, alignTo='Stop', plotSpikes=true, spacing=400);
xlim([-200, 500])
plotRaw(axes(figure), eu(iEu), raw, eu(iEu).Trials.Lick, alignTo='Stop', plotSpikes=true, spacing=400)
xlim([-25, 100])
save(sprintf("E:\\Data\\%s_3trials.mat", eu(iEu).getName()), 'x', 't', 'st')

%% Low SNR unit with big lick artefact
iEu = find(strcmpi(eu.getName(), 'daisy14_20220506_Channel29_Unit1')); % Longer lick noise
raw = load(sprintf("C:\\SERVER\\Units\\Lite_NonDuplicate_NonDrift\\raw\\raw_%s.mat", eu(iEu).getName()));
[x, t, st] = plotRaw(axes(figure), eu(iEu), raw, eu(iEu).Trials.Lick, alignTo='Stop', plotSpikes=true, spacing=400);
xlim([-200, 500])
plotRaw(axes(figure), eu(iEu), raw, eu(iEu).Trials.Lick, alignTo='Stop', plotSpikes=true, spacing=400)
xlim([-25, 100])
save(sprintf("E:\\Data\\%s_3trials.mat", eu(iEu).getName()), 'x', 't', 'st')

%% Functions
function [x, t, st] = plotRaw(ax, eu, raw, trials, varargin)
    p = inputParser();
    p.addRequired('ax')
    p.addRequired('eu')
    p.addRequired('raw', @isstruct)
    p.addRequired('trials', @(x) isa(x, 'Trial'))
    p.addParameter('alignTo', 'Start', @(x) ismember(x, {'Start', 'Stop'}))
    p.addParameter('plotSpikes', true, @islogical)
    p.addParameter('spacing', 200, @isnumeric)
    p.addParameter('spacingSigmas', 10, @isnumeric)
    p.addParameter('randomTrials', true, @islogical)
    p.addParameter('window', [-1, 1], @isnumeric)
    p.parse(ax, eu, raw, trials, varargin{:})
    ax = p.Results.ax;
    eu = p.Results.eu;
    raw = p.Results.raw;
    trials = p.Results.trials;
    alignTo = p.Results.alignTo;
    plotSpikes = p.Results.plotSpikes;
    spacing = p.Results.spacing;
    spacingSigmas = p.Results.spacingSigmas;
    randomTrials = p.Results.randomTrials;
    window = p.Results.window;

    if length(trials) < 3
        return
    end
    if randomTrials
        trials = trials(randi(length(trials), [3, 1]));
    else
        trials = trials(1:3);
    end
    [x, t] = eu.getTrialAlignedData(raw.data, raw.t, trials=trials, alignTo=alignTo, window=window, resolution=1/30000);
    
    if ~isnan(spacingSigmas)
        sigma = mad(x, 1, 'all') / 0.67449;
        spacing = round(sigma * spacingSigmas);
    end

    st = eu.SpikeTimes;
    st = cell(3, 1);
    for iTrial = 1:3
        start = trials(iTrial).(alignTo);
        sel = eu.SpikeTimes >= start + window(1) & eu.SpikeTimes <= start + window(2);
        st{iTrial} = eu.SpikeTimes(sel) - start;
    end

    hold(ax, 'on')
    for iTrial = 1:3
        plot(ax, 1000*t, x(iTrial, :) + spacing*(iTrial-2), Color=getColor(iTrial, 3, 0.67))

        if plotSpikes
            scatter(ax, 1000*st{iTrial}, -spacing/2 + spacing*(iTrial-2), 10, MarkerEdgeColor=getColor(iTrial, 3, 0.67))
        end
    end

    xline(ax, 0, 'k--')
    ylim(ax, 2*[-spacing, spacing])
    xlabel(ax, 'Time (ms)')
    ylabel(ax, 'Voltage (uV)')
end

