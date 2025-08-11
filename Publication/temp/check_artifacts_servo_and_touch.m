%% Load EphysUnits
% clear files
% files{1} = dir("C:\SERVER\Units\Lite_NonDuplicate_NonDrift\Daisy2_20180514_*.mat");
% files{2} = dir("C:\SERVER\Units\Lite_NonDuplicate_NonDrift\daisy12_20220106_*.mat");
% files{3} = dir("C:\SERVER\Units\Lite_NonDuplicate_NonDrift\daisy14_20220506_*.mat");
% files = cat(1, files{:});
% files = arrayfun(@(f) sprintf('%s\\%s', f.folder, f.name), files, UniformOutput=false);
% eu = EphysUnit.load(files, waveforms=false, spikecounts=false, spikerates=false);

% eu = EphysUnit.load('C:\SERVER\Units\Lite_NonDuplicate_NonDrift', waveforms=false, spikecounts=false, spikerates=false);
% load('C:\SERVER\Units\meta_Lite_NonDuplicate_NonDrift_20250705.mat')
% eu = eu(c.hasPress & c.hasLick);
% clearvars -except eu
% 
% % Load euComplete (complete with ITI spikes)
% euNames = lower(eu.getName());
% 
% files = dir('C:\SERVER\Units\NonLite_PressVsLick\*.mat');
% sel = ismember(cellfun(@(n) lower(strrep(n, '.mat', '')), {files.name}, UniformOutput=false), euNames);
% files = files(sel);
% cd('C:\SERVER\Units\NonLite_PressVsLick\')
% euComplete = EphysUnit.load({files.name}, waveforms=false, spikecounts=false, spikerates=false);
% 
% euComplete.save('C:\SERVER\Units\NonLite_PressVsLick_NonDuplicate_NonDrift');
% 
% [lia, locb] = ismember(eu.getName(), euComplete.getName());
% eu(lia) = euComplete(locb(lia));
% 
% clearvars -except eu
eu = EphysUnit.load('C:\SERVER\Units\NonLite_PressVsLick_NonDuplicate_NonDrift', waveforms=false, spikecounts=false, spikerates=false);

clear spikeTimesCache
spikeTimesCache(length(eu)) = struct(index=[], name=[], data=[]);

for iEu = 1:length(eu)
    spikeTimesCache(iEu) = struct(index=iEu, name=eu(iEu).getName(), data=eu(iEu).SpikeTimes);
end

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

%% Read raw data and resave (ALREADY DONE, BUT NEED TO REDO FOR BLACKROCK DATA)
% if ~exist('C:\SERVER\Units\Lite_NonDuplicate_NonDrift\raw', 'dir')
%     mkdir('C:\SERVER\Units\Lite_NonDuplicate_NonDrift\raw')
% end
% for iEu = 1:length(eu)
%     % if ismember(eu(iEu).getAnimalName(), {'daisy14', 'daisy15', 'daisy16', 'desmond23', 'desmond24', 'desmond25', 'desmond26', 'desmond27'}) && ~ismember(eu(iEu).ExpName, {'daisy14_20220506', 'daisy16_20220502'})
%     clear data t B selT
%     try
%         [data, t] = eu(iEu).loadRaw();
%         selT = false(size(t));
%         for trial = ["LeverDeploy", "LeverRetract", "TubeDeploy", "TubeRetract"]
%             trials = eu(iEu).Trials.(trial);
%             if isempty(trials)
%                 continue
%             end
%             B = trials.inTrial(t, window=[-1, 1], windowMode='extend');
%             selT = selT | B;
%         end
% 
%         for trial = ["CircLick"]
%             trials = eu(iEu).Trials.(trial);
%             if isempty(trials)
%                 continue
%             end
%             B = trials.inTrial(t);
%             selT = selT | B;
%         end
% 
%         for trial = ["Press", "Lick"]
%             trials = eu(iEu).Trials.(trial);
%             if isempty(trials)
%                 continue
%             end
%             B = trials.inTrial(t, window=[0, 2], windowMode='extend');
%             selT = selT | B;
%         end
% 
%         data = data(:, selT);
%         t = t(:, selT);
%         trials = eu(iEu).Trials;
%         name = eu(iEu).getName();
%         save(sprintf('C:\\SERVER\\Units\\Lite_NonDuplicate_NonDrift\\raw\\raw_%s.mat', eu(iEu).getName()), 'data', 't', 'trials', 'name')
%     catch
%         warning('Could not process unit %i, %s', iEu, eu(iEu).getName())
%     end
%     % end
% end
% clear data t B selT trial trials name

%% Read raw data and resave (JUST FOR BLACKROCK DATA IT IS TRICKY)

% Some sessions had 32 channels, we read all and do common mean ref

% Some sessions had 64 channels (2rigs), we either read 1:32 or 33:64, read
% all and do common mean ref

% Since we're reading all channels anyway, we iterate through sessions
% (rather than iterating through units) and save all requested units. 

% unique(string(eu(~cc.isIntan).getAnimalName))
animalNames.blackrock = ["daisy10", "daisy2", "daisy3", "daisy8", "daisy9", "desmond10", "desmond11", "desmond22"];

animalNames.rig1 = ["desmond10", "desmond11", "desmond12", "daisy4", "desmond14", "desmond16", "desmond18", "daisy7", "desmond21", "desmond22", "daisy9", "daisy11", "daisy12", "daisy13"];
animalNames.rig2 = ["desmond13", "daisy5", "desmond15", "desmond17", "desmond19", "desmond20", "daisy8", "daisy10"];

if ~exist('C:\SERVER\Units\Lite_NonDuplicate_NonDrift\raw', 'dir')
    mkdir('C:\SERVER\Units\Lite_NonDuplicate_NonDrift\raw')
end
selUnits = reshape(find(~cc.isIntan), 1, []);
failedUnits = [];
failedME = {};
for iEu = selUnits
    % if ismember(eu(iEu).getAnimalName(), {'daisy14', 'daisy15', 'daisy16', 'desmond23', 'desmond24', 'desmond25', 'desmond26', 'desmond27'}) && ~ismember(eu(iEu).ExpName, {'daisy14_20220506', 'daisy16_20220502'})
    clear data t B selT
    try
        [data, t] = eu(iEu).loadRaw(chunkDuration=300, subtractMean=true);

        tTic = tic();
        fprintf('Doing trial cropping thing...')
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
        fprintf('Done (%.2fs)\n', toc(tTic));

        tTic = tic();
        fprintf('Saving to file...')
        trials = eu(iEu).Trials;
        name = eu(iEu).getName();
        save(sprintf('C:\\SERVER\\Units\\Lite_NonDuplicate_NonDrift\\raw\\raw_%s.mat', eu(iEu).getName()), 'data', 't', 'trials', 'name')
        fprintf('Done (%.2fs)\n', toc(tTic));
    catch ME
        warning('Could not process unit %i, %s', iEu, eu(iEu).getName())
        failedUnits = [failedUnits, iEu];
        failedME = [failedME, {ME}];
    end
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

cc.hasRaw = reshape(ismember(names.expected, names.found), 1, []);
clear files

% Find out which units use intan
cc.isIntan = false(1, length(eu));
for iEu = 1:length(eu)
    animalName = eu(iEu).getAnimalName();
    expName = eu(iEu).ExpName;
    filesRoot = dir(sprintf("C:\\SERVER\\%s\\%s\\*.rhd", animalName, expName));
    files = dir(sprintf("C:\\SERVER\\%s\\%s\\*\\*.rhd", animalName, expName));
    if ~isempty(files) || ~isempty(filesRoot)
        cc.isIntan(iEu) = true;
    else
        cc.isIntan(iEu) = false;
    end
    % fprintf('%s %i\n', expName, length(files));
end
clear iEu animalName expName files


%% Make behavioral trials and PETH/RD
% tEU = eu.alignTimestamps(["REWARD_ON", "LEVER_RELEASED", "LEVER_RETRACT_START", "LEVER_RETRACTED", "LEVER_RETRACT_END", "TUBE_RETRACT_START", "TUBE_RETRACT_END", "LICK", "LICK_OFF"], ac={ac.ac, ac.euIndicesFirstInSession, ac.expIndices, ac.uniqueExpNames});
% 
% for iEu = 1:length(eu)
%     try
%         leverReleaseTimes = eu(iEu).EventTimes.LEVER_RELEASED;
%         leverRetractTimes = eu(iEu).EventTimes.LEVER_RETRACT_START;
%         if isempty(leverRetractTimes)
%             leverRetractTimes = eu(iEu).EventTimes.LEVER_RETRACTED;
%         end
%         eu(iEu).Trials.Press = Trial(eu(iEu).EventTimes.Cue, eu(iEu).EventTimes.Press, 'first');
%         eu(iEu).Trials.PressIncorrect = eu(iEu).Trials.Press(eu(iEu).Trials.Press.duration() < 4 & eu(iEu).Trials.Press.duration() >= 2);
%         eu(iEu).Trials.PressCorrect = eu(iEu).Trials.Press(eu(iEu).Trials.Press.duration() >= 4);
%         [~, leverRetractTimesIncorrect] = eu(iEu).Trials.PressIncorrect.inTrial(leverRetractTimes, [0, 2], windowMode='stop');
%         assert(nnz(leverRetractTimesIncorrect) > 0)
%         [~, leverRetractTimesCorrect] = eu(iEu).Trials.PressCorrect.inTrial(leverRetractTimes, [-0.1, 8], windowMode='stop');
%         eu(iEu).Trials.RetractReleaseIncorrect = Trial([leverRetractTimesIncorrect, Inf], leverReleaseTimes, 'first');
%         eu(iEu).Trials.RetractReleaseCorrect = Trial([leverRetractTimesCorrect, Inf], leverReleaseTimes, 'first');
%         eu(iEu).Trials.PressReleaseIncorrect = Trial([eu(iEu).Trials.PressIncorrect.Start, Inf], [eu(iEu).Trials.RetractReleaseIncorrect.Stop], 'first');
%         eu(iEu).Trials.PressReleaseCorrect = Trial([eu(iEu).Trials.PressCorrect.Stop, Inf], [eu(iEu).Trials.RetractReleaseCorrect.Stop], 'first');
%         eu(iEu).Trials.PressRetractIncorrect = Trial([eu(iEu).Trials.PressIncorrect.Stop], leverRetractTimesIncorrect, 'first');
%         eu(iEu).Trials.PressRetractCorrect = Trial([eu(iEu).Trials.PressCorrect.Stop], leverRetractTimesCorrect, 'first');
%         % 
%         % fprintf(['press=%i\npressIncorrect=%i, pressCorrect=%i;\nleverRetractTimesIncorrect=%i, leverRetractTimesCorrect=%i;\n' ...
%         %     'eu(iEu).Trials.RetractReleaseIncorrect=%i, eu(iEu).Trials.RetractReleaseCorrect=%i;\n' ...
%         %     'eu(iEu).Trials.PressReleaseIncorrect=%i, eu(iEu).Trials.PressReleaseCorrect=%i\n'], length(eu(iEu).Trials.Press), length(eu(iEu).Trials.PressIncorrect), length(eu(iEu).Trials.PressCorrect), ...
%         %     length(leverRetractTimesIncorrect), length(leverRetractTimesCorrect), ...
%         %     length(eu(iEu).Trials.RetractReleaseIncorrect), length(eu(iEu).Trials.RetractReleaseCorrect), ...
%         %     length(eu(iEu).Trials.PressReleaseIncorrect), length(eu(iEu).Trials.PressReleaseCorrect));
%         % disp(1)
%     catch
%         warning('Counld not process iEu = %i, "%s"', iEu, eu(iEu).getName())
%     end
% end
% 
% % Find and correct bad lick labels
% for iEu = 1:length(eu)
%     eu(iEu).EventTimes.LICK = [];
%     eu(iEu).EventTimes.LICK_OFF = [];
% end
% eu.alignTimestamps(["LICK", "LICK_OFF"], ac={ac.ac, ac.euIndicesFirstInSession, ac.expIndices, ac.uniqueExpNames});
% clear nLicks
% nLicks.EUvAC = arrayfun(@(eu) [length(eu.EventTimes.Lick), length(eu.EventTimes.LICK)], eu, UniformOutput=false);
% nLicks.EUvAC = cat(1, nLicks.EUvAC{:});
% nLicks.EUvReward = arrayfun(@(eu) [length(eu.EventTimes.Lick), length(eu.EventTimes.RewardTimes)], eu, UniformOutput=false);
% nLicks.EUvReward = cat(1, nLicks.EUvReward{:});
% 
% lickIsReward = diff(nLicks.EUvReward, 1, 2) == 0; lickIsReward = lickIsReward(:)';
% 
% for iEu = find(lickIsReward)
%     eu(iEu).Trials.Lick = Trial(eu(iEu).EventTimes.Cue, eu(iEu).EventTimes.LICK, 'first', eu(iEu).EventTimes.Press);
% end
% 
% for iEu = 1:length(eu)
%     eu(iEu).Trials.LickIncorrect = eu(iEu).Trials.Lick(eu(iEu).Trials.Lick.duration() < 4 & eu(iEu).Trials.Lick.duration() >= 2);
%     eu(iEu).Trials.LickCorrect = eu(iEu).Trials.Lick(eu(iEu).Trials.Lick.duration() >= 4);
% end
% 
% % Make rd and eta
% clear rd
% rd.press = eu.getRasterData('press', [-4, 4], minTrialDuration=2, maxTrialDuration=Inf);
% rd.releaseCorrect = eu.getRasterData('press_release_correct', [-4, 4], minTrialDuration=0, maxTrialDuration=Inf);
% rd.releaseIncorrect = eu.getRasterData('press_release_incorrect', [-4, 4], minTrialDuration=0, maxTrialDuration=Inf);
% rd.retractCorrect = eu.getRasterData('press_retract_correct', [-4, 4], minTrialDuration=0, maxTrialDuration=Inf);
% rd.retractIncorrect = eu.getRasterData('press_retract_incorrect', [-4, 4], minTrialDuration=0, maxTrialDuration=Inf);
% rd.lick = eu.getRasterData('lick', [-4, 4], minTrialDuration=2, maxTrialDuration=Inf);
% 
% eta.correctPress = eu.getETA('count', 'press', [-4, 4], minTrialDuration=4, normalize='none', resolution=0.1);
% eta.incorrectPress = eu.getETA('count', 'press', [-4, 4], minTrialDuration=2, maxTrialDuration=4, normalize='none', resolution=0.1);
% 
% eta.correctRetract = eu.getETA('count', 'press_retract_correct', [-4, 4], alignTo='stop', normalize='none', resolution=0.1);
% eta.incorrectRetract = eu.getETA('count', 'press_retract_incorrect', [-4, 4], alignTo='stop', normalize='none', resolution=0.1);
% 
% eta.correctRelease = eu.getETA('count', 'press_release_correct', [-4, 4], alignTo='stop', normalize='none', resolution=0.1);
% eta.incorrectRelease = eu.getETA('count', 'press_release_incorrect', [-4, 4], alignTo='stop', normalize='none', resolution=0.1);
% 
% eta.correctLick = eu.getETA('count', 'lick', [-4, 4], minTrialDuration=4, normalize='none', resolution=0.1);
% eta.incorrectLick= eu.getETA('count', 'lick', [-4, 4], minTrialDuration=2, maxTrialDuration=4, normalize='none', resolution=0.1);
% 
% for field = ["correctPress", "incorrectPress", "correctLick", "incorrectLick", "correctRelease", "incorrectRelease", "correctRetract", "incorrectRetract"]
%     eta.(field).X = eta.(field).X ./ 0.1;
% end



%% Plot reach raster, peth, DLC trajectories for reach-dec DLC cells (among 163)
% % close all
% if ~exist('E:\DATA\Figures\reach_retract_dlc', 'dir')
%     mkdir('E:\DATA\Figures\reach_retract_dlc')
% end
% 
% fig = figure(Units='normalized', InnerPosition=[0 0 1 1]);
% tl = tiledlayout(fig, 5, 3, TileIndexing='columnmajor');
% ax = gobjects(15, 1);
% for iAx = 1:15
%     ax(iAx) = nexttile(tl);
% end
% % for iEu = find(cc.hasRaw(:)' & cc.isIntan(:)')
% for iEu = find(cc.hasRaw(:)')
%     try
%         raw = load(sprintf("C:\\SERVER\\Units\\Lite_NonDuplicate_NonDrift\\raw\\raw_%s.mat", eu(iEu).getName()));
%         % 1: Raster
%         for iAx = 1:15
%             cla(ax(iAx), 'reset');
%         end
%         EphysUnit.plotRaster(ax(1), rd.press(iEu), xlim=[-4, 4]);
%         set(ax(1).Legend, AutoUpdate=false)
%         xline(ax(1), 0, 'k-', LineWidth=1.5)
%         yline(ax(1), find(rd.press(iEu).duration>4, 1), 'k-', LineWidth=3)
%         delete(ax(1).Legend)
%         legend(ax(1), ["", "trial start"])
% 
%         % 2: PETH
%         hold(ax(2), 'on')
%         h = gobjects(2, 1);
%         h(1) = plot(ax(2), eta.incorrectPress.t, eta.incorrectPress.X(iEu, :), Color=hsl2rgb([0.4, 0.4, 0.4]), LineStyle='-', DisplayName=sprintf('Incorrect (n=%i)', eta.incorrectPress.N(iEu)), LineWidth=2);
%         h(2) = plot(ax(2), eta.correctPress.t, eta.correctPress.X(iEu, :), 'k-', DisplayName=sprintf('Correct (n=%i)', eta.correctPress.N(iEu)), LineWidth=2);
%         legend(h, AutoUpdate=false, Location='southwest')
%         xline(ax(2), 0, 'k--')
%         yline(ax(2), 0, 'k--')
%         xlim(ax(2), [-4, 4])
%         % ylim(ax(2), [0, 100])
% 
%         title(ax(1), 'Reach Raster')
%         title(ax(2), 'Reach PETH')
% 
%         % 3: Reach incorrect (raw)
%         % [xAligned, tAligned] = eu(iEu).getTrialAlignedData(raw.data, raw.t, trials=eu(iEu).Trials.PressIncorrect, alignTo='stop', window=[-1, 1], resolution=1/30000);
%         plotRaw(ax(3), eu(iEu), raw, eu(iEu).Trials.PressIncorrect, alignTo='Stop', window=[-0.5, 1.5]);
%         title(ax(3), 'PressIncorrect')
% 
%         % 4: Reach correct (raw)
%         plotRaw(ax(4), eu(iEu), raw, eu(iEu).Trials.PressCorrect, alignTo='Stop', window=[-0.5, 1.5]);
%         title(ax(4), 'PressCorrect')
% 
%         % 5: Tube servo deploy (raw)
%         plotRaw(ax(5), eu(iEu), raw, eu(iEu).Trials.TubeDeploy, alignTo='Start', window=[-0.5, 1.5]);
%         title(ax(5), 'TubeDeploy')
% 
% 
%         % 6: Raster
%         EphysUnit.plotRaster(ax(6), rd.retractCorrect(iEu), xlim=[-4, 4]);
%         set(ax(6).Legend, AutoUpdate=false)
%         xline(ax(6), 0, 'k-', LineWidth=1.5)
%         % yline(ax(6), find(rd.lick(i).duration>4, 1), 'k-', LineWidth=1.5)
%         delete(ax(6).Legend)
%         legend(ax(6), ["", "bar-contact"])
% 
%         % 7 PETH
%         hold(ax(7), 'on')
%         h = gobjects(2, 1);
%         h(1) = plot(ax(7), eta.incorrectRetract.t, eta.incorrectRetract.X(iEu, :), Color=hsl2rgb([0.4, 0.4, 0.4]), LineStyle='-', DisplayName=sprintf('Incorrect (n=%i)', eta.incorrectRetract.N(iEu)), LineWidth=2);
%         h(2) = plot(ax(7), eta.correctRetract.t, eta.correctRetract.X(iEu, :), 'k-', DisplayName=sprintf('Correct (n=%i)', eta.correctRetract.N(iEu)), LineWidth=2);
%         legend(h, AutoUpdate=false, Location='southwest')
%         xline(ax(7), 0, 'k--')
%         yline(ax(7), 0, 'k--')
%         xlim(ax(7), [-4, 4])
%         % ylim(ax(7), [0, 100])
% 
%         % 8: reach retract (raw)
%         plotRaw(ax(8), eu(iEu), raw, eu(iEu).Trials.PressRetractCorrect, alignTo='Stop', window=[-0.5, 1.5]);
%         title(ax(8), 'PressRetractCorrect')
% 
%         % 9: lever servo deploy (raw)
%         plotRaw(ax(9), eu(iEu), raw, eu(iEu).Trials.LeverDeploy, alignTo='Start', window=[-0.5, 1.5]);
%         title(ax(9), 'LeverDeploy')
% 
%         % 10: lever servo retract (raw)
%         plotRaw(ax(10), eu(iEu), raw, eu(iEu).Trials.LeverRetract, alignTo='Start', window=[-0.5, 1.5]);
%         title(ax(10), 'LeverRetract')
% 
%         % 11: Raster
%         EphysUnit.plotRaster(ax(11), rd.lick(iEu), xlim=[-4, 4]);
%         set(ax(11).Legend, AutoUpdate=false)
%         xline(ax(11), 0, 'k-', LineWidth=1.5)
%         yline(ax(11), find(rd.lick(iEu).duration>4, 1), 'k-', LineWidth=3)
%         delete(ax(11).Legend)
%         legend(ax(11), ["", "trial start"])
% 
%         % 12 PETH
%         hold(ax(12), 'on')
%         h = gobjects(2, 1);
%         h(1) = plot(ax(12), eta.incorrectLick.t, eta.incorrectLick.X(iEu, :), Color=hsl2rgb([0.4, 0.4, 0.4]), LineStyle='-', DisplayName=sprintf('Incorrect (n=%i)', eta.incorrectLick.N(iEu)), LineWidth=2);
%         h(2) = plot(ax(12), eta.correctLick.t, eta.correctLick.X(iEu, :), 'k-', DisplayName=sprintf('Correct (n=%i)', eta.correctLick.N(iEu)), LineWidth=2);
%         legend(h, AutoUpdate=false, Location='southwest')
%         xline(ax(12), 0, 'k--')
%         yline(ax(12), 0, 'k--')
%         xlim(ax(12), [-4, 4])
%         % ylim(ax(12), [0, 100])
% 
%         % 13: Lick incorrect (raw)
%         plotRaw(ax(13), eu(iEu), raw, eu(iEu).Trials.LickIncorrect, alignTo='Stop', window=[-0.5, 1.5]);
%         title(ax(13), 'LickIncorrect')
% 
%         % 14: Lick correct (raw)
%         plotRaw(ax(14), eu(iEu), raw, eu(iEu).Trials.LickCorrect, alignTo='Stop', window=[-0.5, 1.5]);
%         title(ax(14), 'LickCorrect')
% 
%         % 15: Tube servo retract (raw)
%         plotRaw(ax(15), eu(iEu), raw, eu(iEu).Trials.TubeRetract, alignTo='Start', window=[-0.5, 1.5]);
%         title(ax(15), 'TubeRetract')
% 
%         title(ax(1), 'Reach Raster')
%         title(ax(2), 'Reach PETH')
%         title(ax(6), 'Retract Raster (correct trials)')
%         title(ax(7), 'Retract PETH')
%         title(ax(11), 'Lick Raster')
%         title(ax(12), 'Lick PETH')
% 
%         title(tl, rd.press(iEu).name, Interpreter='none')
% 
%         xlim(ax, [-3, 3])
%         xlim(ax([3:5, 8:10, 13:15]), [-500, 1500])
% 
%         ylim(ax([2, 7, 12]), [min([ax(2).YLim, ax(7).YLim, ax(12).YLim]), max([ax(2).YLim, ax(7).YLim, ax(12).YLim])]);
%         print(fig, sprintf('E:\\DATA\\Figures\\reach_retract_dlc\\%s.png', eu(iEu).getName()), '-dpng')
%     catch
%         warning('Error processing unit %i', iEu)
%     end
% end

%% Cool units with cool functional responses from John
% coolUnitNames = [
%     % "Daisy2_20180425_Channel17_Unit1", ... 
%     % "Daisy3_20180618_Channel29_Unit1", ... 
%     % "Daisy8_20210708_Channel10_Unit1", ... 
%     "Daisy14_20220506_Channel38_Unit1", ...
%     "Daisy15_20220511_Channel104_Unit1", ...
% ];
% coolUnitNames = [
%     "daisy14_20220506_Channel38_Unit1", ... % Big unit, brief artifact
%     "daisy15_20220511_Channel104_Unit1", ... % Big unit, long artifact
%     "desmond27_20220526_Channel106_Unit1", ... % Randomly chosen medium-SNR unit
%     "desmond25_20220430_Channel124_Unit1", ... % Two units one channel
%     "desmond25_20220430_Channel124_Unit2", ... % Two units one channel 2, electric boogaloo
%     "desmond26_20220531_Channel4_Unit1", ... % Huge unit but missing some spikes due to scaling/shrinkning waveforms
% ];

coolUnitNames = [
    "daisy15_20220601_Channel114_Unit1", ... % Big unit, long artifact
    "daisy14_20220506_Channel38_Unit1", ... % Big unit, brief artifact
    "daisy15_20220511_Channel104_Unit1", ... % Big unit, long artifact
    "desmond27_20220526_Channel106_Unit1", ... % Randomly chosen medium-SNR unit
    "desmond25_20220430_Channel124_Unit1", ... % Two units one channel
    "desmond25_20220430_Channel124_Unit2", ... % Two units one channel 2, electric boogaloo
    "desmond26_20220531_Channel4_Unit1", ... % Huge unit but missing some spikes due to scaling/shrinkning waveforms
];

%% Filter raw and then redo spike detection
close all
fs = 30000;
maxLog = 4;
ksAlpha = 0.01;
ksHitRateThreshold = 0.5;

pTemplateMatching.distanceFactor = 1;
pTemplateMatching.nSigmas = 5;
pTemplateMatching.rateExceed = 0.05;
pTemplateMatching.method = 'euclidean';
pTemplateMatching.maxLog = maxLog;
pTemplateMatching.ksAlpha = ksAlpha;
pTemplateMatching.ksHitRateThreshold = ksHitRateThreshold;
savePath = "C:\SERVER\Figures\lick_artifact_removal\euclidean_run3_just_blackrock";
if ~exist(savePath, 'dir')
    mkdir(savePath);
end
save(sprintf("%s\\pTemplateMatching.mat", savePath), 'pTemplateMatching')

selUnits = find(cc.hasRaw & ~cc.isIntan); useRawCache = false;
% selUnits = find(ismember(eu.getName(), coolUnitNames)); useRawCache = true;
unitNames = eu.getName();
edgesPeriLick = 0:0.01:0.1; nBinsQQPeriLick = min(length(edgesPeriLick) - 1, 10);
n = 11; edgesCircLick = linspace(pi/n, (2-1/n)*pi, n); nBinsQQCircLick = n - 1; 
clear n

clear layout
layout.fig = figure(Unit='normalized', InnerPosition=[0 0 1 1]);
layout.h = [4, 6, 3, 2, 3, 2];
layout.w = [2, 1, 1, 1, 4];
layout.tl = tiledlayout(layout.fig, sum(layout.h), sum(layout.w), TileSpacing='tight', Padding='compact');
layout.ax.waveform = nexttile(layout.tl, [layout.h(1), layout.w(1)]);

% Two QQPlot axis, needs additional schenannnnnnnegans
for iQQ = 1:2
    switch iQQ
        case 1
            nBinsQQ = nBinsQQPeriLick;
            edges = edgesPeriLick;
        case 2
            nBinsQQ = nBinsQQCircLick;
            edges = edgesCircLick;
    end
    layout.qq(iQQ).tl = tiledlayout(layout.tl, nBinsQQ, nBinsQQ, TileSpacing='none', Padding='none');
    layout.qq(iQQ).tl.Layout.Tile = tilenum(layout.tl, 1, sum(layout.w(1:iQQ)) + 1);
    layout.qq(iQQ).tl.Layout.TileSpan = [layout.h(1), layout.w(iQQ+1)];
    layout.qq(iQQ).axDummy = axes(layout.qq(iQQ).tl);
    layout.qq(iQQ).axDummy.Layout.TileSpan = [nBinsQQ, nBinsQQ];
    for i = 1:nBinsQQ
        for j = 1:nBinsQQ
            layout.qq(iQQ).ax(i, j) = axes(layout.qq(iQQ).tl);
            layout.qq(iQQ).ax(i, j).Layout.Tile = tilenum(layout.qq(iQQ).tl, i, j);
            xticks(layout.qq(iQQ).ax(i, j), [])
            yticks(layout.qq(iQQ).ax(i, j), [])
            layout.qq(iQQ).ax(i, j).Visible = false;
        end
    end
    hold(layout.qq(iQQ).ax, 'on')
    switch iQQ
        case 1
            xlim(layout.qq(iQQ).axDummy, round([0, edges(nBinsQQ+1)]*1e3))
            ylim(layout.qq(iQQ).axDummy, round([0, edges(nBinsQQ+1)]*1e3))
            xticks(layout.qq(iQQ).axDummy, [10, 50, 100])
            yticks(layout.qq(iQQ).axDummy, [10, 50, 100])
            xline(layout.qq(iQQ).axDummy, [10], 'k--')
            yline(layout.qq(iQQ).axDummy, [10], 'k--')
            xlabel(layout.qq(iQQ).axDummy, 'time to lick (ms)')
            ylabel(layout.qq(iQQ).axDummy, 'time to lick (ms)')
        case 2
            xlim(layout.qq(iQQ).axDummy, [edges(1), edges(end)])
            ylim(layout.qq(iQQ).axDummy, [edges(1), edges(end)])
            % xticks(layout.qq(iQQ).axDummy, [0, edges(1), pi, edges(end), 2*pi])
            % yticks(layout.qq(iQQ).axDummy, [0, edges(1), pi, edges(end), 2*pi])
            % xticklabels(layout.qq(iQQ).axDummy, ["0", string(arrayfun(@(t) sprintf("%i/%i\\pi", round(t/pi*(nBins+1)), nBins+1), [edges(1), pi, edges(end)])), "2\pi"])
            % yticklabels(layout.qq(iQQ).axDummy, ["0", string(arrayfun(@(t) sprintf("%i/%i\\pi", round(t/pi*(nBins+1)), nBins+1), [edges(1), pi, edges(end)])), "2\pi"])
            xticks(layout.qq(iQQ).axDummy, [edges(1), pi, edges(end)])
            yticks(layout.qq(iQQ).axDummy, [edges(1), pi, edges(end)])
            xticklabels(layout.qq(iQQ).axDummy, [">0", "\pi", "<2\pi"])
            yticklabels(layout.qq(iQQ).axDummy, [">0", "\pi", "<2\pi"])
            xline(layout.qq(iQQ).axDummy, edges(2), 'k--')
            yline(layout.qq(iQQ).axDummy, edges(2), 'k--')
            xlabel(layout.qq(iQQ).axDummy, 'lick phase')
            ylabel(layout.qq(iQQ).axDummy, 'lick phase')
    end
    title(layout.qq(iQQ).axDummy, 'QQ-Plot & -logp_{KS}', FontSize=9)
end
layout.ax.etaCircLick = nexttile(layout.tl, [layout.h(1), layout.w(4)]);
layout.ax.etaLick = nexttile(layout.tl, [layout.h(1), layout.w(5)]);
layout.ax.raw = nexttile(layout.tl, [layout.h(2), sum(layout.w)]);
layout.ax.binnedWaveforms = nexttile(layout.tl, [layout.h(3), sum(layout.w)]);
layout.ax.binnedHistograms = nexttile(layout.tl, [layout.h(4), sum(layout.w)]);
layout.ax.binnedWaveformsByPhase = nexttile(layout.tl, [layout.h(5), sum(layout.w)]);
layout.ax.binnedHistogramsByPhase = nexttile(layout.tl, [layout.h(6), sum(layout.w)]);


if useRawCache
    if ~exist('rawCache', 'var')
        rawCache(length(eu)) = struct(index=[], name=[], data=[], t=[], trials=[]);
    end
else
    clear spikesFiltered spikesRaw
    clear rawCache
    spikesFiltered = struct(index=[], name=[], sampleIndex=[], timestamps=[], waveforms=[], waveformTimestamps=[], isUnit=[]);
    spikesRaw = struct(index=[], name=[], sampleIndex=[], timestamps=[], waveforms=[], waveformTimestamps=[], isUnit=[]);
end
for iUnit = 1:length(selUnits)
    try
        cla(layout.ax.waveform)
        cla(layout.ax.etaLick)
        cla(layout.ax.etaCircLick)
        cla(layout.ax.raw)
        
        clear raw filtered
        iEu = selUnits(iUnit);
        fprintf('Processing unit %i (%i/%i):\n', iEu, iUnit, length(selUnits))
        tTic = tic();
        fprintf('\tLoading raw data...');
        if useRawCache
            if isempty(rawCache(iEu).data)
                raw = load(sprintf("C:\\SERVER\\Units\\Lite_NonDuplicate_NonDrift\\raw\\raw_%s.mat", eu(iEu).getName()));
                raw.index = iEu;
                rawCache(iEu) = raw;
            else
                raw = rawCache(iEu);
            end
        else
            raw = load(sprintf("C:\\SERVER\\Units\\Lite_NonDuplicate_NonDrift\\raw\\raw_%s.mat", eu(iEu).getName()));
            raw.index = iEu;
        end
        fprintf('Done (%.2f s)\n', toc(tTic));
    
        eu(iEu).SpikeTimes = spikeTimesCache([spikeTimesCache.index]==iEu).data;
        oldSpikeTimes = eu(iEu).SpikeTimes;
    
        % Filter the whole continuous data
        tTic = tic();
        fprintf('\tFiltering artifacts...');
        filtered = raw;
        filtered.data = removeArtifact(raw.data, raw.t, method='highpass', highpassCutoff=800, sampleRate=fs);
        stAbs = eu(iEu).SpikeTimes;
        fprintf('Done (%.2f s)\n', toc(tTic));
    
        % Find spikes in raw data
        tTic = tic();
        fprintf('\tDetecting spikes in raw data...');
        if useRawCache && strcmpi(spikesRaw.name, eu(iEu).getName())
            fprintf('Cache found...')
        else
            spikesRaw.index = iEu;
            spikesRaw.name = eu(iEu).getName();
            [spikesRaw.sampleIndex, spikesRaw.timestamps, spikesRaw.waveforms, spikesRaw.waveformTimestamps] = spikeDetect(raw, SampleRate=fs, NumSigmas=2.5, NumSigmasReturn=1.25, NumSigmasReject=20, WaveformWindow=[-0.5, 0.5]);
            spikesRaw.isUnit = ismember(round(spikesRaw.timestamps*fs), round(stAbs*fs));
        end
        fprintf('(%i/%i detected spikes matched timestamps of %i eu spikes)...', nnz(spikesRaw.isUnit), length(spikesRaw.isUnit), length(stAbs))
        spikeTemplateRaw = mean(spikesRaw.waveforms, 1, 'omitnan');
        maxMicroVolts = max(abs(spikeTemplateRaw))*2;
        fprintf('Done (%.2f s)\n', toc(tTic));
    
        % Detect spikes in filtered data
        tTic = tic();
        fprintf('\tDetecting spikes in filtered data...');
        if useRawCache && strcmpi(spikesFiltered.name, eu(iEu).getName())
            fprintf('Cache found...')
        else
            spikesFiltered.index = iEu;
            spikesFiltered.name = eu(iEu).getName();
            [spikesFiltered.sampleIndex, spikesFiltered.timestamps, spikesFiltered.waveforms, spikesFiltered.waveformTimestamps] = spikeDetect(filtered, SampleRate=fs, NumSigmas=2.5, NumSigmasReturn=1.25, NumSigmasReject=20, WaveformWindow=[-0.5, 0.5], MaxMicroVolts=maxMicroVolts);
        end
        fprintf('Done (%.2f s)\n', toc(tTic));

        % % Extract template spikes in filtered data, using existing spiketimes
        % % from eu
        tTic = tic();
        fprintf('\tExtracing waveforms to use as templates from filtered data...');
        [spikeTemplate, ~] = getWaveforms(filtered, [-0.5, 0.5], spikesRaw.sampleIndex(spikesRaw.isUnit), IndexType='SampleIndex');
        [noiseTemplate, tWaveform] = getWaveforms(filtered, [-0.5, 0.5], spikesRaw.sampleIndex(~spikesRaw.isUnit), IndexType='SampleIndex');
        fprintf('Done (%.2f s)\n', toc(tTic));

        spikeTemplate = mean(spikeTemplate, 1, 'omitnan');
        noiseTemplate = mean(noiseTemplate, 1, 'omitnan');
    
        switch pTemplateMatching.method
        % Template matching method 1: Compare to template based on euclidean distance
            case 'euclidean'
                distToSpikeTemplate = sum((spikesFiltered.waveforms - spikeTemplate).^2, 2);
                distToNoiseTemplate = sum((spikesFiltered.waveforms - noiseTemplate).^2, 2);  
                spikesFiltered.isUnit = distToSpikeTemplate < distToNoiseTemplate * pTemplateMatching.distanceFactor; 
                [~, I] = sort(distToSpikeTemplate, 'ascend');
                clear distToSpikeTemplate distToNoiseTemplate
        % Template matching method 1: Compare to template based on corrcoef
            case 'corr'
                warning("Are you sure laddie? Corr does not work as well as euclidean.")
                X = spikesFiltered.waveforms - mean(spikesFiltered.waveforms, 2);
                Y = spikeTemplate - mean(spikeTemplate, 2);
                Z = noiseTemplate - mean(noiseTemplate, 2);
                corrWithSpikeTemplate = (X*Y') ./ (sqrt(sum(X.^2, 2)) * sqrt(sum(Y.^2, 2)));
                corrWithNoiseTemplate = (X*Z') ./ (sqrt(sum(X.^2, 2)) * sqrt(sum(Z.^2, 2)));
                spikesFiltered.isUnit = corrWithSpikeTemplate * pTemplateMatching.distanceFactor > corrWithNoiseTemplate; 
                [~, I] = sort(corrWithSpikeTemplate./corrWithNoiseTemplate, 'descend');
                clear X Y Z corrWithNoiseTemplate corrWithSpikeTemplate
        end

        % Stricter template matching to remove other units/artifacts
        residuals = spikesFiltered.waveforms - spikeTemplate;
        sigma = mad(residuals(spikesFiltered.isUnit, :), 1, 'all') / 0.67449;
        pOutlier = sum(residuals > pTemplateMatching.nSigmas*sigma, 2)./size(residuals, 2);
        spikesFiltered.isUnit = spikesFiltered.isUnit & pOutlier<pTemplateMatching.rateExceed;
        fprintf('\tRemoved %i/%i as outliers.\n', nnz(pOutlier>=pTemplateMatching.rateExceed), length(pOutlier));
    
        eu(iEu).SpikeTimes = spikesFiltered.timestamps(spikesFiltered.isUnit);

        tTic = tic();
        fprintf('\tGenerating plots...');

        hold(layout.ax.waveform, 'on')
        for iWave = 100:100:size(I)
            if mean(spikesFiltered.isUnit(I((iWave-100)+1:iWave))) > 0.5
                plot(layout.ax.waveform, tWaveform, spikesFiltered.waveforms(I(iWave), :), Color=[1 0 0 0.1]) %, Color=[getColor(i/10, ceil(length(I)/10), 0.67), 0.25])
                % plot(layout.ax.waveform, tWaveform, mean(spikesFiltered.waveforms(I(iWave-100+1:iWave), :), 1, 'omitnan'), Color=[1 0 0 0.1]) %, Color=[getColor(i/10, ceil(length(I)/10), 0.67), 0.25])
            else
                plot(layout.ax.waveform, tWaveform, spikesFiltered.waveforms(I(iWave), :), Color=[0.1 0.1 0.1, 0.015]) %Color=[getColor(i/10, ceil(length(I)/10), 0.67, s=0.1, l=0.1), 0.1])
                % plot(layout.ax.waveform, tWaveform, mean(spikesFiltered.waveforms(I(iWave-100+1:iWave), :), 1, 'omitnan'), Color=[0.1 0.1 0.1, 0.025]) %Color=[getColor(i/10, ceil(length(I)/10), 0.67, s=0.1, l=0.1), 0.1])
            end
        end
        plot(layout.ax.waveform, tWaveform, noiseTemplate, LineWidth=2, Color='green')
        plot(layout.ax.waveform, tWaveform, spikeTemplate, LineWidth=2, Color='blue')
        plot(layout.ax.waveform, tWaveform, spikeTemplate + pTemplateMatching.nSigmas*sigma, Color='blue', LineStyle='--')
        plot(layout.ax.waveform, tWaveform, spikeTemplate - pTemplateMatching.nSigmas*sigma, Color='blue', LineStyle='--')
        xlim(layout.ax.waveform, [-0.5, 0.5])
        ylim(layout.ax.waveform, [min(spikeTemplate - pTemplateMatching.nSigmas*sigma), max(spikeTemplate + pTemplateMatching.nSigmas*sigma)])
        xlabel(layout.ax.waveform, 'ms')
        ylabel(layout.ax.waveform, '\muV')
        title(layout.ax.waveform, sprintf('Filtered: %i/%i spikes, %i/%i noise, %i original', nnz(spikesFiltered.isUnit), nnz(spikesRaw.isUnit), nnz(~spikesFiltered.isUnit), nnz(~spikesRaw.isUnit), nnz(stAbs)))

        % Plot raw trace by trial, see if there's false positives
        [xRawAligned, tAligned, stRawAligned, trials] = parseRaw(eu(iEu), raw, eu(iEu).getTrials('circlick_naive', minInterval=0.05, maxInterval=0.20), spikeTimes=oldSpikeTimes, window=[-0, 0.10], alignTo='Start', randomTrials=true, nTrials=15, timestampMode='relative');
        [xFilteredAligned, ~, stFilteredAligned] = parseRaw(eu(iEu), filtered, trials, spikeTimes=spikesFiltered.timestamps(spikesFiltered.isUnit), window=[-0, 0.10], alignTo='Start', randomTrials=false, nTrials='all', timestampMode='relative');
    
        plotRaw(layout.ax.raw, xFilteredAligned, tAligned, stFilteredAligned, plotSpikes=true, spacing=250, xRaw=xRawAligned, stRaw=stRawAligned);
        xlim(layout.ax.raw, [0, 100])
        xticks(layout.ax.raw, [0, 10, 50, 100])
        xline(layout.ax.raw, 0:10:100, 'k--')
    
        % Get the original spike times
        eu(iEu).SpikeTimes = oldSpikeTimes;
        % Calculate ETA
        etaTemp.circLickNaiveRaw = eu(iEu).getETA('count', 'circlick_naive', window=[0, 2*pi], resolution=2*pi/30, normalize='none',  minInterval=0.05, maxInterval=0.20, ...
            lickArtifactLengthType='ms', lickArtifactLength=10, lickArtifactDirection='both', lickOffArtifactLengthType='ms', lickOffArtifactLength=10, lickOffArtifactDirection='both');
        etaTemp.lickRaw = eu(iEu).getETA('count', 'lick', window=[-4, 2], resolution=0.025, normalize='none');
    
        % Get the new filtered spike times
        eu(iEu).SpikeTimes = spikesFiltered.timestamps(spikesFiltered.isUnit);
        etaTemp.circLickNaiveFiltered = eu(iEu).getETA('count', 'circlick_naive', window=[0, 2*pi], resolution=2*pi/30, normalize='none',  minInterval=0.05, maxInterval=0.20, ...
            lickArtifactLengthType='ms', lickArtifactLength=10, lickArtifactDirection='both', lickOffArtifactLengthType='ms', lickOffArtifactLength=10, lickOffArtifactDirection='both');
        etaTemp.lickFiltered = eu(iEu).getETA('count', 'lick', window=[-4, 2], resolution=0.025, normalize='none');
    
        h = gobjects(2, 1);
        hold(layout.ax.etaCircLick, 'on')
        h(1) = plot(layout.ax.etaCircLick, etaTemp.circLickNaiveRaw.t, etaTemp.circLickNaiveRaw.X, 'blue', DisplayName='raw');
        h(2) = plot(layout.ax.etaCircLick, etaTemp.circLickNaiveFiltered.t, etaTemp.circLickNaiveFiltered.X, 'red', DisplayName='filtered');

        xlim(layout.ax.etaCircLick, [0, 2*pi])
        xline(layout.ax.etaCircLick, [0, edgesCircLick, 2*pi], ':', Color=[0.2, 0.2, 0.2, 0.2])
        xticks(layout.ax.etaCircLick, [0, pi, 2*pi])
        xticklabels(layout.ax.etaCircLick, {'0', '\pi', '2\pi'})

        ylabel(layout.ax.etaCircLick, 'spikes/s')
        xlabel(layout.ax.etaCircLick, 'lick phase')
        legend(layout.ax.etaCircLick, h)
        title(layout.ax.etaCircLick, 'PETH (Self-timed lick)')
    
        hold(layout.ax.etaLick, 'on')
        plot(layout.ax.etaLick, etaTemp.lickRaw.t, etaTemp.lickRaw.X./0.025, 'blue', DisplayName='raw')
        plot(layout.ax.etaLick, etaTemp.lickFiltered.t, etaTemp.lickFiltered.X./0.025, 'red', DisplayName='filtered')
        xlabel(layout.ax.etaLick, 'Time to self-timed lick (s)')
        ylabel(layout.ax.etaLick, 'spikes/s')
        legend(layout.ax.etaLick)
        title(layout.ax.etaCircLick, 'PETH (Self-timed lick)')

        ylim(layout.ax.etaLick, 'auto')
        ylim(layout.ax.etaCircLick, 'auto')
        drawnow()
        ylim([layout.ax.etaLick, layout.ax.etaCircLick], ...
            [min(arrayfun(@(ax) ax.YLim(1), [layout.ax.etaLick, layout.ax.etaCircLick])), ...
            max(arrayfun(@(ax) ax.YLim(2), [layout.ax.etaLick, layout.ax.etaCircLick]))]);

        title(layout.tl, eu(iEu).getName(), Interpreter='none')
        
        % Plot waveforms by bin (aligned to lick)
        cla(layout.ax.binnedWaveforms)
        cla(layout.ax.binnedHistograms)
        hold(layout.ax.binnedWaveforms, 'on')
        hold(layout.ax.binnedHistograms, 'on')
        trials = eu(iEu).getTrials('circlick_naive', minInterval=0.05, maxInterval=0.20);
        h = gobjects(2, 1);
        clear binnedPeriLickWaveforms
        nBins = length(edgesPeriLick) - 1;
        binnedPeriLickWaveforms(nBins) = struct(iBin=[], tWave=[], meanWaveform=[], window=[], distToSpikeTemplate=[], distToNoiseTemplate=[], distRatio=[], ratioN=[], ratioEdges=[], nKSTestHits=[]);
        for iBin = 1:nBins
            window = edgesPeriLick(iBin:iBin+1);
            tempTrials = Trial([trials.Start]+window(1), [trials.Start]+window(2), advancedValidation=false);
            [inBin, ~] = tempTrials.inTrial(spikesFiltered.timestamps);
            waveforms = spikesFiltered.waveforms(inBin(:) & spikesFiltered.isUnit(:), :);
            tWave = spikesFiltered.waveformTimestamps*10+mean(window)*1e3;
            mu = mean(waveforms, 1, 'omitnan');

            distToSpikeTemplate = sum((waveforms - spikeTemplate).^2, 2);
            distToNoiseTemplate = sum((waveforms - noiseTemplate).^2, 2);
            distRatio = distToNoiseTemplate ./ distToSpikeTemplate;
            [~, IBinned] = sort(distRatio, 'ascend');
            waveforms = waveforms(IBinned, :);

            waveformBinSize = 10;
            for iWaveform = 1:waveformBinSize:size(waveforms, 1)
                plot(layout.ax.binnedWaveforms, tWave, mean(waveforms(iWaveform:min(iWaveform+waveformBinSize-1, size(waveforms, 1)), :), 1), Color=[getColor(iWaveform/waveformBinSize, length(IBinned)/waveformBinSize, 0.67), 0.25], LineWidth=0.5)
                % plot(layout.ax.binnedWaveforms, tWave, waveforms(iWaveform, :), Color=[getColor(iWaveform/waveformBinSize, length(IBinned)/waveformBinSize, 0.67), 0.25], LineWidth=0.5)
            end
            h(1) = plot(layout.ax.binnedWaveforms, tWave, spikeTemplate, Color='blue', LineWidth=1, LineStyle='--', DisplayName='template');
            h(2) = plot(layout.ax.binnedWaveforms, tWave, mu, Color='red', LineWidth=1, LineStyle='--', DisplayName='unit mean');

            text(layout.ax.binnedWaveforms, tWave(16), min(spikeTemplate - 2*sigma), sprintf('%.1f sp/s', size(waveforms, 1)/length(tempTrials)/0.01), HorizontalAlignment='center', VerticalAlignment='bottom')
            % text(layout.ax.binnedWaveforms, tWave(16), 0, sprintf('%i', size(waveforms, 1)), HorizontalAlignment='center', VerticalAlignment='bottom')

            % Embed a histogram of distanceToTemplate
            [ratioN, ~] = histcounts(log2(distRatio), [0:0.05:maxLog-0.05, Inf], Normalization='probability');
            ratioEdges = (0:0.05:maxLog) + (iBin-1)*maxLog;

            binnedPeriLickWaveforms(iBin) = struct(iBin=iBin, tWave=tWave, meanWaveform=mu, window=window, distToSpikeTemplate=distToSpikeTemplate, distToNoiseTemplate=distToNoiseTemplate, distRatio=distRatio, ratioN=ratioN, ratioEdges=ratioEdges, nKSTestHits=0);
        end
        xticks(layout.ax.binnedWaveforms, 0:10:100)
        xline(layout.ax.binnedWaveforms, 0:10:100, 'k--')
        xlim(layout.ax.binnedWaveforms, [0, 100])
        ylim(layout.ax.binnedWaveforms, [min(spikeTemplate - 2*sigma), max(spikeTemplate + 2*sigma)])
        ylabel(layout.ax.binnedWaveforms, '\muV')
        xlabel(layout.ax.binnedWaveforms, 'Time (ms)')
        legend(layout.ax.binnedWaveforms, h, Orientation='horizontal', Location='northeast')

        for i = 1:nBinsQQPeriLick
            for j = 1:nBinsQQPeriLick
                cla(layout.qq(1).ax(i, j))
            end
        end

        % QQ-plot matrix for waveform distRatio
        for i = 1:nBinsQQPeriLick
            for j = 1:i
                axQQ = layout.qq(1).ax(nBinsQQPeriLick+1-j, i);
                axPKS = layout.qq(1).ax(nBinsQQPeriLick+1-i, j);
                [ks.h, ks.p] = kstest2(log2(binnedPeriLickWaveforms(i).distRatio), log2(binnedPeriLickWaveforms(j).distRatio), Alpha=ksAlpha/(nBinsQQPeriLick*(nBinsQQPeriLick-1)/2));
                Qi = quantile(log2(binnedPeriLickWaveforms(i).distRatio), 0:0.025:1);
                Qj = quantile(log2(binnedPeriLickWaveforms(j).distRatio), 0:0.025:1);
                if ks.h
                    color = 'red';
                    binnedPeriLickWaveforms(i).nKSTestHits = binnedPeriLickWaveforms(i).nKSTestHits + 1;
                    binnedPeriLickWaveforms(j).nKSTestHits = binnedPeriLickWaveforms(j).nKSTestHits + 1;
                else
                    color = 'black';
                end
                plot(axQQ, Qi, Qj, Color=color)
                plot(axPKS, Qi, Qj, Color=color)
                text(axPKS, 0.1, 0.9, sprintf('%i', round(-log10(ks.p))), Color=color, FontSize=6, VerticalAlignment='top', HorizontalAlignment='left', Units='normalized')
            end
        end
        xticks(layout.qq(1).ax, [])
        yticks(layout.qq(1).ax, [])
        set(layout.qq(1).ax, Visible=false)

        % Plot histograms
        for iBin = 1:nBins
            if binnedPeriLickWaveforms(iBin).nKSTestHits / (nBins-1) >= ksHitRateThreshold
                color = 'red';
            else
                color = 'black';
            end
            histogram(layout.ax.binnedHistograms, BinEdges=binnedPeriLickWaveforms(iBin).ratioEdges, BinCounts=binnedPeriLickWaveforms(iBin).ratioN, EdgeColor='none', FaceColor=color);
        end
        xlim(layout.ax.binnedHistograms, [0, maxLog*nBins])
        xticks(layout.ax.binnedHistograms, 0:1:maxLog*nBins)
        xticklabels(layout.ax.binnedHistograms, string(repmat(0:1:(maxLog-1), 1, nBins)))
        xline(layout.ax.binnedHistograms, 0:maxLog:maxLog*nBins, 'k--')
        yticks(layout.ax.binnedHistograms, [])
        ylabel(layout.ax.binnedHistograms, 'p')
        xlabel(layout.ax.binnedHistograms, 'log_2(distToNoise/distToTemplate)')
        layout.ax.binnedHistograms.TickLength = [0, 0];

        % Plot waveforms by bin (aligned to circlick)
        cla(layout.ax.binnedWaveformsByPhase)
        cla(layout.ax.binnedHistogramsByPhase)
        hold(layout.ax.binnedWaveformsByPhase, 'on')
        hold(layout.ax.binnedHistogramsByPhase, 'on')
        trials = eu(iEu).getTrials('circlick_naive', minInterval=0.05, maxInterval=0.20);
        durations = trials.duration();
        % Do lickOff artifact blanking
        if isfield(eu(iEu).EventTimes, 'LickOff')
            lickOff = eu(iEu).EventTimes.LickOff;
        elseif isfield(eu(iEu).EventTimes, 'LICK_OFF')
            lickOff = eu(iEu).EventTimes.LICK_OFF;
        else
            error('Could not find lick off event under either eu.EventTimes.LICK_OFF or eu.EventTimes.LickOff');
        end
        [~, lickOffNonNan, lickOffTrialIndices] = trials.inTrial(lickOff);
        lickOff = NaN(1, length(trials));
        lickOff(lickOffTrialIndices) = lickOffNonNan;
        % Well by difinition there needs to be a lickOff before the next lickOn
        assert(issorted(lickOffTrialIndices, 'ascend') && length(lickOff)==length(trials), "Well by difinition there needs to be a lickOff before the next lickOn")
        assert(size(durations, 1) == 1 && size(lickOffTrialIndices, 1) == 1 && size(lickOff, 1) == 1)
        clear lickOffTrialIndices lickOffNonNan

        h = gobjects(2, 1);
        clear binnedCircLickWaveformsByPhase
        nBins = length(edgesCircLick) - 1;
        binnedCircLickWaveformsByPhase(nBins) = struct( ...
            iBin=[], tWave=[], meanWaveform=[], window=[], windowReal=[], hasArtifact=[],...
            distToSpikeTemplate=[], distToNoiseTemplate=[], distRatio=[], ratioN=[], ratioEdges=[], nKSTestHits=[]);
        for iBin = 1:nBins
            window = edgesCircLick(iBin:iBin+1);
            windowReal = [(window(1)/(2*pi))*durations; (window(2)/(2*pi))*durations];            
            tempTrials = Trial([trials.Start] + windowReal(1, :), [trials.Start] + windowReal(2, :), advancedValidation=false);
            % Do lickOn artifact blanking
            hasArtifactOn = windowReal(1, :) < 10e-3;
            % Do lickOff artifact blanking
            A = [tempTrials.Start];
            B = [tempTrials.Stop];
            C = lickOff - 0.010;
            D = lickOff + 0.010;
            hasArtifactOff = ~(A>=D | B<=C);
            hasArtifactOff = hasArtifactOff & ~isnan(lickOff);
            hasArtifact = hasArtifactOn | hasArtifactOff;
            clear hasArtifactOn hasArtifactOff A B C D

            tempTrials = tempTrials(~hasArtifact);
            if isempty(tempTrials)
                binnedCircLickWaveformsByPhase(iBin) = struct( ...
                    iBin=iBin, tWave=[], meanWaveform=[], window=window, windowReal=windowReal, hasArtifact=hasArtifact, ...
                    distToSpikeTemplate=[], distToNoiseTemplate=[], distRatio=[], ratioN=[], ratioEdges=[], nKSTestHits=[]);
            else
                [inBin, ~] = tempTrials.inTrial(spikesFiltered.timestamps);
                waveforms = spikesFiltered.waveforms(inBin(:) & spikesFiltered.isUnit(:), :);
                tWave = spikesFiltered.waveformTimestamps*diff(window)+mean(window);
                mu = mean(waveforms, 1, 'omitnan');
    
                distToSpikeTemplate = sum((waveforms - spikeTemplate).^2, 2);
                distToNoiseTemplate = sum((waveforms - noiseTemplate).^2, 2);
                distRatio = distToNoiseTemplate ./ distToSpikeTemplate;
                [~, IBinned] = sort(distRatio, 'ascend');
                waveforms = waveforms(IBinned, :);
    
                waveformBinSize = 10;
                for iWaveform = 1:waveformBinSize:size(waveforms, 1)
                    plot(layout.ax.binnedWaveformsByPhase, tWave, mean(waveforms(iWaveform:min(iWaveform+waveformBinSize-1, size(waveforms, 1)), :), 1), Color=[getColor(iWaveform/waveformBinSize, length(IBinned)/waveformBinSize, 0.67), 0.25], LineWidth=0.5)
                    % plot(layout.ax.binnedWaveforms, tWave, waveforms(iWaveform, :), Color=[getColor(iWaveform/waveformBinSize, length(IBinned)/waveformBinSize, 0.67), 0.25], LineWidth=0.5)
                end
                h(1) = plot(layout.ax.binnedWaveformsByPhase, tWave, spikeTemplate, Color='blue', LineWidth=1, LineStyle='--', DisplayName='template');
                h(2) = plot(layout.ax.binnedWaveformsByPhase, tWave, mu, Color='red', LineWidth=1, LineStyle='--', DisplayName='unit mean');
    
                text(layout.ax.binnedWaveformsByPhase, tWave(16), min(spikeTemplate - 2*sigma), sprintf('%.1f sp \\times %i trials', size(waveforms, 1)/length(tempTrials), length(tempTrials)), HorizontalAlignment='center', VerticalAlignment='bottom')
                % text(layout.ax.binnedWaveforms, tWave(16), 0, sprintf('%i', size(waveforms, 1)), HorizontalAlignment='center', VerticalAlignment='bottom')
    
                % Embed a histogram of distanceToTemplate
                [ratioN, ~] = histcounts(log2(distRatio), [0:0.05:maxLog-0.05, Inf], Normalization='probability');
                ratioEdges = (0:0.05:maxLog) + (iBin-1)*maxLog;
                binnedCircLickWaveformsByPhase(iBin) = struct( ...
                    iBin=iBin, tWave=tWave, meanWaveform=mu, window=window, windowReal=windowReal, hasArtifact=hasArtifact, ...
                    distToSpikeTemplate=distToSpikeTemplate, distToNoiseTemplate=distToNoiseTemplate, distRatio=distRatio, ratioN=ratioN, ratioEdges=ratioEdges, nKSTestHits=0);
            end
        end
        xline(layout.ax.binnedWaveformsByPhase, edgesCircLick, 'k--')
        xticks(layout.ax.binnedWaveformsByPhase, [0, edgesCircLick, 2*pi])
        xticklabels(layout.ax.binnedWaveformsByPhase, ["0", string(arrayfun(@(t) sprintf("%i/%i\\pi", round(t/pi*(nBins+1)), nBins+1), edgesCircLick)), "2\pi"])
        xlim(layout.ax.binnedWaveformsByPhase, [edgesCircLick(1), edgesCircLick(end)])
        ylim(layout.ax.binnedWaveformsByPhase, [min(spikeTemplate - 2*sigma), max(spikeTemplate + 2*sigma)])
        ylabel(layout.ax.binnedWaveformsByPhase, '\muV')
        xlabel(layout.ax.binnedWaveformsByPhase, 'lick phase')
        legend(layout.ax.binnedWaveformsByPhase, h, Orientation='horizontal', Location='northeast')

        % QQ-plot matrix for waveform distRatio

        for i = 1:nBinsQQCircLick
            for j = 1:nBinsQQCircLick
                cla(layout.qq(2).ax(i, j))
            end
        end

        for i = 1:nBinsQQCircLick
            for j = 1:i
                if isempty(binnedCircLickWaveformsByPhase(i).distRatio) || isempty(binnedCircLickWaveformsByPhase(j).distRatio)
                    continue;
                end
                axQQ = layout.qq(2).ax(nBinsQQCircLick+1-j, i);
                axPKS = layout.qq(2).ax(nBinsQQCircLick+1-i, j);
                
                [ks.h, ks.p] = kstest2(log2(binnedCircLickWaveformsByPhase(i).distRatio), log2(binnedCircLickWaveformsByPhase(j).distRatio), Alpha=ksAlpha/(nBinsQQCircLick*(nBinsQQCircLick-1)/2));
                Qi = quantile(log2(binnedCircLickWaveformsByPhase(i).distRatio), 0:0.025:1);
                Qj = quantile(log2(binnedCircLickWaveformsByPhase(j).distRatio), 0:0.025:1);
                if ks.h
                    color = 'red';
                    binnedCircLickWaveformsByPhase(i).nKSTestHits = binnedCircLickWaveformsByPhase(i).nKSTestHits + 1;
                    binnedCircLickWaveformsByPhase(j).nKSTestHits = binnedCircLickWaveformsByPhase(j).nKSTestHits + 1;
                else
                    color = 'black';
                end
                plot(axQQ, Qi, Qj, Color=color)
                plot(axPKS, Qi, Qj, Color=color)
                text(axPKS, 0.1, 0.9, sprintf('%i', round(-log10(ks.p))), Color=color, FontSize=6, VerticalAlignment='top', HorizontalAlignment='left', Units='normalized')
            end
        end
        xticks(layout.qq(2).ax, [])
        yticks(layout.qq(2).ax, [])
        set(layout.qq(2).ax, Visible=false)

        % Plot histograms (circlick)
        circLickArtifactRate = 0;
        validCircLickBins = 0;
        for iBin = 1:nBins
            if isempty(binnedCircLickWaveformsByPhase(iBin).distRatio)
                continue
            end
            validCircLickBins = validCircLickBins + 1;
            if binnedCircLickWaveformsByPhase(iBin).nKSTestHits / (nBins-1) >= ksHitRateThreshold
                circLickArtifactRate = circLickArtifactRate + 1;
                color = 'red';
                text(layout.ax.etaCircLick, mean(binnedCircLickWaveformsByPhase(iBin).window), mean(layout.ax.etaCircLick.YLim), '*', VerticalAlignment='bottom', HorizontalAlignment='center', FontSize=14, Color='red')
            else
                color = 'black';
            end
            histogram(layout.ax.binnedHistogramsByPhase, BinEdges=binnedCircLickWaveformsByPhase(iBin).ratioEdges, BinCounts=binnedCircLickWaveformsByPhase(iBin).ratioN, EdgeColor='none', FaceColor=color);
        end
        circLickArtifactRate = circLickArtifactRate/validCircLickBins;
        xlim(layout.ax.binnedHistogramsByPhase, [0, maxLog*nBins])
        xticks(layout.ax.binnedHistogramsByPhase, 0:1:maxLog*nBins)
        xticklabels(layout.ax.binnedHistogramsByPhase, string(repmat(0:1:(maxLog-1), 1, nBins)))
        xline(layout.ax.binnedHistogramsByPhase, 0:maxLog:maxLog*nBins, 'k--')
        yticks(layout.ax.binnedHistogramsByPhase, [])
        ylabel(layout.ax.binnedHistogramsByPhase, 'p')
        xlabel(layout.ax.binnedHistogramsByPhase, 'log_2(distToNoise/distToTemplate)')
        layout.ax.binnedHistogramsByPhase.TickLength = [0, 0];

        fprintf('Done (%.2f s)\n', toc(tTic));
        clear trials iBin window spikeTimesInBin t waveforms IBinned distToSpikeTemplate distToNoiseTemplate distRatio waveformBinSize ks i j Qi Qj color ks

        tTic = tic();
        fprintf('\tSaving plots...');
        print(layout.fig, sprintf('%s\\%s.png', savePath, eu(iEu).getName()), '-dpng', '-r0')
        fprintf('Done (%.2f s)\n', toc(tTic));
    
        index = iEu;
        name = eu(iEu).getName();
        tTic = tic();
        fprintf('\tSaving data...');
        save(sprintf("%s\\%s.mat", savePath, name), 'index', 'name', 'spikesFiltered', 'spikesRaw', 'oldSpikeTimes', 'spikeTemplate', 'noiseTemplate', 'binnedPeriLickWaveforms', 'binnedCircLickWaveformsByPhase', 'circLickArtifactRate', '-v7.3')
        fprintf('Done (%.2f s)\n', toc(tTic));
    catch ME
        warning('Error processing unit %i (%i)', iEu, iUnit);
        warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
    end
end

clear useRawCache iUnit iEu tTic raw filtered stAbs spikeTemplateRaw maxMicroVolts spikeTemplate tWaveform noiseTemplate distToSpikeTemplate distToNoiseTemplate residuals sigma pOutlier I iWave xRawAligned tAligned stRawAligned trials xFilteredAligned stFilteredAligned etaTemp
clear index name tTic
clear trials edgesPeriLick iBin window tempTrials inBin spikeTimesInBin waveforms tWave mu distToSpikeTemplate distToNoiseTemplate distRatio IBinned waveformBinSize iWaveform binnedPeriLickWaveforms nBinsQQPeriLick

%% Load a unit
savePath = "C:\SERVER\Figures\lick_artifact_removal\euclidean_run2";

load(sprintf("%s\\pTemplateMatching.mat", savePath))

selUnits = find(cc.hasRaw & cc.isIntan); useRawCache = false;

clear unit
unit(length(eu)) = struct(index=[], name=[], spikesFiltered=[], oldSpikeTimes=[], spikeTemplate=[], noiseTemplate=[], binnedPeriLickWaveforms=[], binnedCircLickWaveformsByPhase=[], circLickArtifactRate=[]);

lineLength = 0;
tTicTotal = tic();
for iUnit = 1:length(selUnits)
    try
        iEu = selUnits(iUnit);
        tTic = tic();
        S = load(sprintf("%s\\%s.mat", savePath, eu(iEu).getName()), 'index', 'name', 'spikesFiltered', 'oldSpikeTimes', 'spikeTemplate', 'noiseTemplate', 'binnedPeriLickWaveforms', 'binnedCircLickWaveformsByPhase', 'circLickArtifactRate');
        assert(S.index == iEu && strcmpi(S.name, eu(iEu).getName()))
        unit(iEu) = S;        
        fprintf(repmat('\b', 1, lineLength))
        lineLength = fprintf("Loaded unit %i/%i, iEu=%i, name='%s' (%.2fs, %.2fs total)\n", iUnit, length(selUnits), iEu, eu(iEu).getName(), toc(tTic), toc(tTicTotal));
    catch ME
        warning('Error processing unit %i', iUnit);
        warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
    end
end

clear iEu S lineLength tTic tTicTotal iUnit

%% Calculate ETA and bootstrap for cyclic firing
cc.hasLickArtifact = arrayfun(@(unit) isempty(unit.circLickArtifactRate) || unit.circLickArtifactRate > 0, unit);

for iEu = find(~cc.hasLickArtifact)
    eu(iEu).SpikeTimes = unit(iEu).spikesFiltered.timestamps(unit(iEu).spikesFiltered.isUnit);
end
eta.circLickNaiveFiltered = eu.getETA('count', 'circlick_naive', selUnits=~cc.hasLickArtifact, window=[0, 2*pi], resolution=2*pi/30, normalize='none',  minInterval=0.05, maxInterval=0.20, ...
    lickArtifactLengthType='ms', lickArtifactLength=10, lickArtifactDirection='both', lickOffArtifactLengthType='ms', lickOffArtifactLength=10, lickOffArtifactDirection='both');
eta.pressNormFiltered = eu.getETA('count', 'press', window=[-4, 2], normalize=[-4, -2], minTrialDuration=2);
eta.lickBoutNaiveFiltered = eu.getETA('count', 'lickbout_naive', selUnits=~cc.hasLickArtifact, window=[0, 2*pi*4], resolution=2*pi/30, normalize='none',  minInterval=0.05, maxInterval=0.20, ...
    minBoutCycles=2, maxBoutCycles=4, ...
    lickArtifactLengthType='ms', lickArtifactLength=10, lickArtifactDirection='both', lickOffArtifactLengthType='ms', lickOffArtifactLength=10, lickOffArtifactDirection='both');



for iEu = find(~cc.hasLickArtifact)
    eu(iEu).SpikeTimes = unit(iEu).oldSpikeTimes;
end
eta.circLickNaiveRaw = eu.getETA('count', 'circlick_naive', selUnits=~cc.hasLickArtifact, window=[0, 2*pi], resolution=2*pi/30, normalize='none',  minInterval=0.05, maxInterval=0.20, ...
    lickArtifactLengthType='ms', lickArtifactLength=10, lickArtifactDirection='both', lickOffArtifactLengthType='ms', lickOffArtifactLength=10, lickOffArtifactDirection='both');
eta.pressNormRaw = eu.getETA('count', 'press', window=[-4, 2], normalize=[-4, -2], minTrialDuration=2);
eta.lickBoutNaiveRaw = eu.getETA('count', 'lickbout_naive', selUnits=~cc.hasLickArtifact, window=[0, 2*pi*4], resolution=2*pi/30, normalize='none',  minInterval=0.05, maxInterval=0.20, ...
    minBoutCycles=2, maxBoutCycles=4, ...
    lickArtifactLengthType='ms', lickArtifactLength=10, lickArtifactDirection='both', lickOffArtifactLengthType='ms', lickOffArtifactLength=10, lickOffArtifactDirection='both');




%% Bootstrap to find the significance of average Z vector magnitudes (shuffle bins, not trials)
clear bootCLick
[bootCLick.filtered.magH, bootCLick.filtered.magCI, bootCLick.filtered.Z] = bootCircLick(eta.circLickNaiveFiltered, selUnits=~cc.hasLickArtifact, alpha=0.01, nBoot=100000, replace=false, seed=42, interpFirstBin=false, replaceNansWithMean=true);
[bootCLick.raw.magH, bootCLick.raw.magCI, bootCLick.raw.Z] = bootCircLick(eta.circLickNaiveRaw, selUnits=~cc.hasLickArtifact, alpha=0.01, nBoot=100000, replace=false, seed=42, interpFirstBin=false, replaceNansWithMean=true);

cc.isLickFiltered = bootCLick.filtered.magH(:)';
cc.isLickRaw = bootCLick.raw.magH(:)';

fprintf('nCircLick units: raw: %i, filtered: %i, both: %i\n', nnz(cc.isLickRaw), nnz(cc.isLickFiltered), nnz(cc.isLickFiltered & cc.isLickRaw));

%% Plot ETA Heatmap for oscilick units
NAME = ["raw", "filtered"];
SEL = {cc.isLickRaw, cc.isLickFiltered};
ETA = {eta.circLickNaiveRaw, eta.circLickNaiveFiltered};
ETABOUT = {eta.lickBoutNaiveFiltered, eta.lickBoutNaiveRaw};
ETAPRESS = {eta.pressNormFiltered, eta.pressNormRaw};

fig = figure;
tl = tiledlayout(fig, sum(SEL{1}) + sum(SEL{2}), 6, TileIndexing='columnmajor');
ax = gobjects(2, 1);

for iSrc = 1:2
    sel = SEL{iSrc};
    meanZ = bootCLick.(NAME(iSrc)).Z(sel);
    phase = angle(meanZ);
    amp = abs(meanZ);
    phase(phase < 0) = phase(phase < 0) + 2*pi;
    
    phase = phase(:);
    amp = amp(:);
    
    [sortedPhase, I] = sort(phase);
    
    etaTemp.circLickNaiveNorm.(NAME(iSrc)) = ETA{iSrc};
    etaTemp.circLickNaiveNorm.(NAME(iSrc)).X = normalize(ETA{iSrc}.X, 2, 'zscore', 'robust');

    ax = nexttile(tl, [sum(sel), 2]);
    [~, ~] = EphysUnit.plotETA(ax, etaTemp.circLickNaiveNorm.(NAME(iSrc)), sel, order=I, ...
        clim=[-5, 5], xlim=[0, 2*pi], hidecolorbar=true);
    title(ax, sprintf("%s\n(normalized to inter-lick-interval)", NAME(iSrc)));
    xlim(ax, [0, 2*pi])
    xticks(ax, (0:1:2).*pi)
    xticklabels(ax, ["0", "\pi", "2\pi"]);
    xline(ax, (0:1:2).*pi, 'k--')
    xlabel('lick phase')
    ylabel('unit')
end

for iSrc = 1:2
    sel = SEL{iSrc};
    meanZ = bootCLick.(NAME(iSrc)).Z(sel);
    phase = angle(meanZ);
    amp = abs(meanZ);
    phase(phase < 0) = phase(phase < 0) + 2*pi;
    
    phase = phase(:);
    amp = amp(:);
    
    [sortedPhase, I] = sort(phase);
    
    etaTemp.lickBoutNaiveNorm.(NAME(iSrc)) = ETABOUT{iSrc};
    etaTemp.lickBoutNaiveNorm.(NAME(iSrc)).X = normalize(ETABOUT{iSrc}.X, 2, 'zscore', 'robust');

    ax = nexttile(tl, [sum(sel), 4]);
    [~, ~] = EphysUnit.plotETA(ax, etaTemp.lickBoutNaiveNorm.(NAME(iSrc)), sel, order=I, ...
        clim=[-5, 5], xlim=[0, 8*pi], hidecolorbar=false);
    title(ax, sprintf("%s\n(normalized to inter-lick-interval)", NAME(iSrc)));
    xlim(ax, [0, 8*pi])
    xticks(ax, (0:2:8).*pi);
    xticklabels(ax, ["0", arrayfun(@(x) sprintf("%i\\pi", x), 2:2:8)]);
    xline(ax, (0:1:8).*pi, 'k--')
    xlabel('lick phase')
    ylabel('unit')
end

% for iSrc = 1:2
%     sel = SEL{iSrc};
%     meanZ = bootCLick.(NAME(iSrc)).Z(sel);
%     phase = angle(meanZ);
%     amp = abs(meanZ);
%     phase(phase < 0) = phase(phase < 0) + 2*pi;
% 
%     phase = phase(:);
%     amp = amp(:);
% 
%     [sortedPhase, I] = sort(phase);
% 
%     etaTemp.circLickNaiveNormToPressBaseline.(NAME(iSrc)) = ETA{iSrc};
%     etaTemp.circLickNaiveNormToPressBaseline.(NAME(iSrc)).X(sel, :) = (ETA{iSrc}.X(sel, :) - vertcat(ETAPRESS{iSrc}.stats(sel).mean)./0.1) ./ (vertcat(ETAPRESS{iSrc}.stats(sel).sd)./0.1);
% 
%     ax = nexttile(tl, [sum(sel), 2]);
%     [~, ~] = EphysUnit.plotETA(ax, etaTemp.circLickNaiveNormToPressBaseline.(NAME(iSrc)), sel, order=I, ...
%         xlim=[0, 2*pi], clim=[-5, 5], hidecolorbar=true);
%     title(ax, sprintf("%s\n(normalized to [-4, -2] pre-reach)", NAME(iSrc)));
%     xlim(ax, [0, 2*pi])
%     xticks(ax, (0:1:2).*pi)
%     xticklabels(ax, ["0", "\pi", "2\pi"]);
%     xlabel('lick phase')
%     ylabel('unit')
% end

% for iSrc = 1:2
%     sel = SEL{iSrc};
%     meanZ = bootCLick.(NAME(iSrc)).Z(sel);
%     phase = angle(meanZ);
%     amp = abs(meanZ);
%     phase(phase < 0) = phase(phase < 0) + 2*pi;
% 
%     phase = phase(:);
%     amp = amp(:);
% 
%     [sortedPhase, I] = sort(phase);
% 
%     etaTemp.lickBoutNaiveNormToPressBaseline.(NAME(iSrc)) = ETABOUT{iSrc};
%     etaTemp.lickBoutNaiveNormToPressBaseline.(NAME(iSrc)).X(sel, :) = (ETABOUT{iSrc}.X(sel, :) - vertcat(ETAPRESS{iSrc}.stats(sel).mean)./0.1) ./ (vertcat(ETAPRESS{iSrc}.stats(sel).sd)./0.1);
% 
%     ax = nexttile(tl, [sum(sel), 4]);
%     [~, ~] = EphysUnit.plotETA(ax, etaTemp.lickBoutNaiveNormToPressBaseline.(NAME(iSrc)), sel, order=I, ...
%         clim=[-5, 5], xlim=[0, 8*pi], hidecolorbar=false);
%     title(ax, sprintf("%s\n(normalized to [-4, -2] pre-reach)", NAME(iSrc)));
%     xlim(ax, [0, 8*pi])
%     xticks(ax, (0:2:8).*pi);
%     xticklabels(ax, ["0", arrayfun(@(x) sprintf("%i\\pi", x), 2:2:8)]);
%     xlabel('lick phase')
%     ylabel('unit')
% end

%% Functions
function [x, t, st, trials] = parseRaw(eu, raw, trials, varargin)
    p = inputParser();
    p.addRequired('eu')
    p.addRequired('raw', @isstruct)
    p.addRequired('trials', @(x) isa(x, 'Trial'))
    p.addParameter('spikeTimes', [], @isnumeric)
    p.addParameter('alignTo', 'Start', @(x) ismember(x, {'Start', 'Stop'}))
    p.addParameter('window', [-1, 1], @isnumeric)
    p.addParameter('randomTrials', true, @islogical)
    p.addParameter('nTrials', 'all', @(x) isnumeric(x) || strcmpi(x, 'all'))
    p.addParameter('timestampMode', 'relative', @(x) ismember(x, {'absolute', 'relative'}))
    p.parse(eu, raw, trials, varargin{:})
    eu = p.Results.eu;
    raw = p.Results.raw;
    trials = p.Results.trials;
    spikeTimes = p.Results.spikeTimes;
    if isempty(spikeTimes)
        spikeTimes = eu.SpikeTimes;
    end
    alignTo = p.Results.alignTo;
    window = p.Results.window;
    randomTrials = p.Results.randomTrials;
    nTrials = p.Results.nTrials;
    timestampMode = p.Results.timestampMode;

    if ischar(nTrials) && strcmpi(nTrials, 'all')
        nTrials = length(trials);
    end

    if length(trials) < nTrials
        return
    end
    if randomTrials
        trials = trials(randi(length(trials), [nTrials, 1]));
    else
        trials = trials(1:nTrials);
    end

    switch lower(timestampMode)
        case 'absolute'
            [x, ~, ~, t] = eu.getTrialAlignedData(raw.data, raw.t, trials=trials, alignTo=alignTo, window=window, resolution=1/30000);
        case 'relative'
            [x, t] = eu.getTrialAlignedData(raw.data, raw.t, trials=trials, alignTo=alignTo, window=window, resolution=1/30000);
        otherwise
            error('unknown timestampMode %s', timestampMode);
    end

    st = cell(nTrials, 1);
    for iTrial = 1:nTrials
        start = trials(iTrial).(alignTo);
        sel = spikeTimes >= start + window(1) & spikeTimes <= start + window(2);
        switch lower(timestampMode)
            case 'absolute'
                st{iTrial} = spikeTimes(sel);
            case 'relative'
                st{iTrial} = spikeTimes(sel) - start;
        end
    end
    if strcmpi(timestampMode, 'absolute')
        x = reshape(x', 1, []);
        t = reshape(t', 1, []);
        st = cat(2, st{:});
    end
end

function x = removeArtifact(x, t, varargin)
    p = inputParser();
    p.addRequired('x', @isnumeric)
    p.addRequired('t', @isnumeric)
    p.addParameter('method', 'none', @(x) ismember(lower(x), {'none', 'highpass', 'movmean', 'spline'}))
    p.addParameter('highpassCutoff', 400, @isnumeric)
    p.addParameter('sampleRate', 30000, @isnumeric)
    p.parse(x, t, varargin{:})
    x = p.Results.x;
    t = p.Results.t;
    method = p.Results.method;
    highpassCutoff = p.Results.highpassCutoff;
    sampleRate = p.Results.sampleRate;

    xRaw = x;

    switch lower(method)
        case 'highpass'
            % x = bpfft(xRaw', sampleRate, highpassCutoff, sampleRate/2)';
            x = highpass(xRaw', highpassCutoff, sampleRate)';
        case 'movmean'

    end
end

function varargout = plotRaw(ax, x, t, st, varargin)
    p = inputParser();
    p.addRequired('ax')
    p.addRequired('x', @isnumeric)
    p.addRequired('t', @isnumeric)
    p.addRequired('st', @iscell)
    p.addParameter('plotSpikes', true, @islogical)
    p.addParameter('spacing', 200, @isnumeric)
    p.addParameter('spacingSigmas', 10, @isnumeric)
    p.addParameter('artifactWindow', [-0.01, 0.01], @isnumeric)
    p.addParameter('xRaw', [], @isnumeric)
    p.addParameter('stRaw', {}, @iscell)
    p.parse(ax, x, t, st, varargin{:})
    ax = p.Results.ax;
    x = p.Results.x;
    t = p.Results.t;
    st = p.Results.st;
    plotSpikes = p.Results.plotSpikes;
    spacing = p.Results.spacing;
    spacingSigmas = p.Results.spacingSigmas;
    artifactWindow = p.Results.artifactWindow;
    xRaw = p.Results.xRaw;
    stRaw = p.Results.stRaw;

    if ~isnan(spacingSigmas)
        sigma = mad(x(:, t<artifactWindow(1) | t>artifactWindow(2)), 1, 'all') / 0.67449;
        spacing = round(sigma * spacingSigmas);
    end

    hold(ax, 'on')
    nTrials = size(x, 1);
    for iTrial = 1:nTrials
        plot(ax, 1000*t, x(iTrial, :) + spacing*(iTrial-1), Color=getColor(iTrial, nTrials, 0.67))

        if ~isempty(xRaw)
            plot(ax, 1000*t, xRaw(iTrial, :) + spacing*(iTrial-1), Color=[getColor(iTrial, nTrials, 0.67, 0.5, 0.3), 0.25], LineStyle='-')
        end

        if plotSpikes
            if isempty(stRaw)
                scatter(ax, 1000*st{iTrial}, -spacing/2 + spacing*(iTrial-1), 10, MarkerEdgeColor=getColor(iTrial, nTrials, 0.67))
            else
                scatter(ax, 1000*st{iTrial}, -spacing/2 + spacing*(iTrial-1), 80, 'x', MarkerEdgeColor=getColor(iTrial, nTrials, 0.67))
                scatter(ax, 1000*stRaw{iTrial}, -spacing/2 + spacing*(iTrial-1), 20, 'o', MarkerEdgeColor=getColor(iTrial, nTrials, 0.67))
            end
        end
    end

    xline(ax, 0, 'k--')
    ylim(ax, [-spacing, (nTrials)*spacing])
    xlabel(ax, 'time (ms)')
    ylabel(ax, 'voltage (uV)')

    varargout = {x, t, st};
end

% Expand waveform window, fill unavailable data with NaN
function varargout = getWaveforms(raw, waveformWindow, index, varargin)
	p = inputParser;
	addRequired(p, 'raw', @isstruct);
	addRequired(p, 'WaveformWindow', @(x) isnumeric(x) && length(x) == 2);
	addRequired(p, 'Index', @isnumeric);
	addParameter(p, 'IndexType', 'SampleIndex', @ischar);
	addParameter(p, 'SampleRate', 30000, @isnumeric);
	parse(p, raw, waveformWindow, index, varargin{:});
	raw 		    = p.Results.raw;
	waveformWindow 	= p.Results.WaveformWindow;
	index 			= p.Results.Index;
	indexType 		= p.Results.IndexType;
    sampleRate      = p.Results.SampleRate;

	switch indexType
		case 'SampleIndex'
			sampleIndex = index;
			% Only interpolate if sample index in non-integer
			if sum(rem(sampleIndex, 1) == 0) == length(sampleIndex)
			    timestamps = raw.t(sampleIndex);
			else
				timestamps = interp1(1:length(raw.t), raw.t, sampleIndex, 'linear');
			end
		case 'Timestamps'
			% Always interpolate if input index in timestamps. This will always be slower.
			timestamps = index;
			sampleIndex = round(interp1(raw.t, 1:length(raw.t), timestamps, 'linear'));
		otherwise
			error(['Unrecognized index type: ''', indexType, ''', must be ''SampleIndex'' or ''Timestamps''.'])
	end

	sampleRate = sampleRate/1000; % Convert to ms
	t = [flip(0:-1/sampleRate:waveformWindow(1)), 1/sampleRate:1/sampleRate:waveformWindow(2)];
	waveforms = NaN(length(sampleIndex), length(t));
	i = sampleRate*t;
	for iWaveform = 1:length(sampleIndex)
		iQuery = sampleIndex(iWaveform) + i;
		if (sum(rem(iQuery, 1) == 0) == length(iQuery)) && min(iQuery) > 0 && max(iQuery) <= size(raw.data, 2)
			% waveforms(iWaveform, :) = obj.Amplifier.Data(obj.MapChannel_TetrodeToRecorded(channel), iQuery);
			waveforms(iWaveform, :) = raw.data(1, iQuery);
		else
			% waveforms(iWaveform, :) = interp1(1:size(obj.Amplifier.Data, 2), double(obj.Amplifier.Data(obj.MapChannel_TetrodeToRecorded(channel), :)), iQuery, 'pchip', NaN);
			waveforms(iWaveform, :) = interp1(1:size(raw.data, 2), double(raw.data(1, :)), iQuery, 'pchip', NaN);
		end
	end

	% Output
	varargout = {waveforms, t, timestamps, sampleIndex};
end

% Detect spike by simple thresholding
function varargout = spikeDetect(raw, varargin)
	p = inputParser;
    addRequired(p, 'raw', @isstruct)
    addOptional(p, 'sel', [], @(x) islogical(x) & isnumeric(x)) % Selected sampleIndices corresponding to raw struct (can use logical indices)
	addParameter(p, 'SampleRate', 30000, @isnumeric);
	addParameter(p, 'NumSigmas', 4, @isnumeric); % Spike detection threshold = (this*sigma*direction). Sigma is estimated noise standard deviation.
	addParameter(p, 'NumSigmasReturn', 1.25, @isnumeric); % [] to disable. Waveform must return to this*sigma*direction after crossing threshold, helps remove noisy periods with non-zero baseline.
	addParameter(p, 'NumSigmasReject', 40, @isnumeric); % [] to disable. Reject huge waveforms that exceed this many sigmas in either direction
	addParameter(p, 'Direction', 'negative', @ischar);
	addParameter(p, 'WaveformWindow', [-0.5, 0.5], @isnumeric);
    addParameter(p, 'MaxMicroVolts', Inf, @isnumeric);
    addParameter(p, 'MinThresholdMicroVolts', 0, @isnumeric);
    addParameter(p, 'MaxThresholdMicroVolts', Inf, @isnumeric);
    addParameter(p, 'UseClampedThresholdInsteadOfSigma', false, @islogical);
	parse(p, raw, varargin{:});

    raw = p.Results.raw;
    sel = p.Results.sel;
    if isempty(sel)
        sel = true(size(raw.data));
    end
	sampleRate      = p.Results.SampleRate;
	numSigmas 		= p.Results.NumSigmas;
	numSigmasReturn = p.Results.NumSigmasReturn;
	numSigmasReject = p.Results.NumSigmasReject;
	directionMode 	= p.Results.Direction;
	waveformWindow 	= p.Results.WaveformWindow;
	maxMicroVolts 	= p.Results.MaxMicroVolts;
    minThreshold    = p.Results.MinThresholdMicroVolts;
    maxThreshold    = p.Results.MaxThresholdMicroVolts;
    useClampedThresholdInsteadOfSigma = p.Results.UseClampedThresholdInsteadOfSigma;

    cleanedSignal = abs(nonzeros(raw.data(1, sel)));
    cleanedSignal = cleanedSignal(cleanedSignal <= maxMicroVolts);
	sigma = median(cleanedSignal, 'all', 'omitnan')/0.6745;
	threshold = max(minThreshold, numSigmas*sigma);
    threshold = min(maxThreshold, threshold);
    thresholdReturn = numSigmasReturn*sigma;
    thresholdReject = numSigmasReject*sigma;
    if useClampedThresholdInsteadOfSigma
        if ~isempty(thresholdReturn)
            thresholdReturn = threshold * numSigmasReturn / numSigmas;
        end
        if ~isempty(thresholdReject)
            thresholdReject = threshold * numSigmasReject / numSigmas;
        end
    end

    switch lower(directionMode)
	    case 'negative'
		    direction = -1;
	    case 'positive'
		    direction = 1;
	    case 'auto'
            error('Not implemented')
		    % direction = sign(median(obj.Amplifier.Data(iChannel, abs(obj.Amplifier.Data(iChannel, :)) > 1.5*threshold))); % Check if spikes are positive or negative
	    otherwise
		    error(['Unrecognized spike detection mode ''', directionMode, '''.'])
    end
	
	% Find spikes
	[~, sampleIndex] = findpeaks(double(direction*raw.data(1, sel)), 'MinPeakHeight', threshold, 'MinPeakProminence', threshold);

	% Extract waveforms
	[waveforms, t] = getWaveforms(raw, waveformWindow, sampleIndex, IndexType='SampleIndex');

	% Align waveforms to peak
	i = sampleRate*t*1e-3;
	[~, maxIndex] = max(direction*waveforms, [], 2);
	alignmentShift = i(maxIndex);
	[waveforms, t, timestamps, sampleIndex] = getWaveforms(raw, waveformWindow, sampleIndex + alignmentShift, IndexType='SampleIndex');

	% Reject waveforms that do not return to a certain level after crossing threshold
	if ~isempty(numSigmasReturn)
		if direction > 0
			selected = min(waveforms(:, t > 0), [], 2) <= thresholdReturn;
		else
			selected = max(waveforms(:, t > 0), [], 2) >= -thresholdReturn;
		end
		waveforms = waveforms(selected, :);
		timestamps = timestamps(selected);
		sampleIndex = sampleIndex(selected);
	end

	% Reject waveforms that exceed a threshold
	if ~isempty(numSigmasReject)
		selected = max(abs(waveforms), [], 2) < abs(thresholdReject) | max(abs(waveforms), [], 2) < maxMicroVolts;
		waveforms = waveforms(selected, :);
		timestamps = timestamps(selected);
		sampleIndex = sampleIndex(selected);
	end

	% Double data so we can do divisions and stuff
	waveforms = double(waveforms);
    waveformTimestamps = t;

    varargout = {sampleIndex, timestamps, waveforms, waveformTimestamps, direction*threshold, direction*thresholdReturn, thresholdReject};
end

function [magH, magCI, Z] = bootCircLick(eta, varargin)
    parser = inputParser();
    parser.addRequired('eta', @isstruct);
    parser.addParameter('selUnits', [], @(x) islogical(x) || isnumeric(x))
    parser.addParameter('alpha', 0.01, @isnumeric);
    parser.addParameter('nBoot', 10000, @isnumeric);
    parser.addParameter('replace', true, @islogical) % true for bootstrap, false for wda
    parser.addParameter('seed', 42, @isnumeric)
    parser.addParameter('interpFirstBin', false, @islogical)
    parser.addParameter('replaceNansWithMean', true, @islogical)
    parser.parse(eta, varargin{:});
    r = parser.Results;
    X = r.eta.X;
    t = r.eta.t;
    selUnits = r.selUnits;
    alpha = r.alpha;
    nBoot = r.nBoot;
    replace = r.replace;
    seed = r.seed;
    interpFirstBin = r.interpFirstBin;
    replaceNansWithMean = r.replaceNansWithMean;
    rng(seed);

    magH = false(size(X, 1), 1);
    magCI = zeros(size(X, 1), 2);
    Z = NaN(size(X, 1), 1) + 1i*NaN(size(X, 1), 1);

    nUnits = size(X, 1);
    if isempty(selUnits)
        selUnits = 1:nUnits;
    elseif islogical(selUnits)
        selUnits = reshape(find(selUnits), 1, []);
    end
    nBins = length(t);
    tTic = tic();
    fprintf('Bootstrapping %i units...', length(selUnits))
    lineLength = 0;
    for iUnit = selUnits
        fprintf(repmat('\b', 1, lineLength));
        lineLength = fprintf('%i/%i', iUnit, nUnits);
        if replace
            I = randi(nBins, [nBoot, nBins]);
        else
            I = zeros(nBoot, nBins);
            for iBoot = 1:nBoot
                I(iBoot, :) = randperm(nBins);
            end
        end
        x = X(iUnit, :);
        if interpFirstBin
            x(1) = mean(x([2, end])); % First bin has lick artifact usually
        end
        if replaceNansWithMean
            x(isnan(x)) = mean(x, 'omitnan');
            zObs = mean(x.*exp(t*1i));
            zRand = mean(x(I).*exp(t*1i), 2);
        else
            zObs = mean(x.*exp(t*1i));
            zRand = mean(x(I).*exp(t*1i), 2);
        end
        
        magObs = abs(zObs);
        magRand = abs(zRand);
        magCI(iUnit, :) = quantile(magRand, [0, 1 - alpha]);
        magH(iUnit) = magObs > magCI(iUnit, 2);
        Z(iUnit) = zObs;
    end
    fprintf('Done (%.2fs)\n', toc(tTic));
end