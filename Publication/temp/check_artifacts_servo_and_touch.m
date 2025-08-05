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

% Load euComplete (complete with ITI spikes)
euNames = lower(eu.getName());

files = dir('C:\SERVER\Units\NonLite_PressVsLick\*.mat');
sel = ismember(cellfun(@(n) lower(strrep(n, '.mat', '')), {files.name}, UniformOutput=false), euNames);
files = files(sel);
cd('C:\SERVER\Units\NonLite_PressVsLick\')
euComplete = EphysUnit.load({files.name}, waveforms=false, spikecounts=false, spikerates=false);

euComplete.save('C:\SERVER\Units\NonLite_PressVsLick_NonDuplicate_NonDrift');

[lia, locb] = ismember(eu.getName(), euComplete.getName());
eu(lia) = euComplete(locb(lia));

clearvars -except eu

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

% Find out which units use intan
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
% % for iEu = find(hasRaw(:)' & isIntan(:)')
% for iEu = find(hasRaw(:)')
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
    % "daisy14_20220506_Channel38_Unit1", ... % Big unit, brief artifact
    % "daisy15_20220511_Channel104_Unit1", ... % Big unit, long artifact
    % "desmond27_20220526_Channel106_Unit1", ... % Randomly chosen medium-SNR unit
    % "desmond25_20220430_Channel124_Unit1", ... % Two units one channel
    % "desmond25_20220430_Channel124_Unit2", ... % Two units one channel 2, electric boogaloo
    % "desmond26_20220531_Channel4_Unit1", ... % Huge unit but missing some spikes due to scaling/shrinkning waveforms
];

%% Filter raw and then redo spike detection
close all
fs = 30000;
pTemplateMatching.distanceFactor = 1;
pTemplateMatching.nSigmas = 5;
pTemplateMatching.rateExceed = 0.05;
pTemplateMatching.method = 'euclidean';
savePath = "E:\Figures\lick_artifact_removal\euclidean";
if ~exist(savePath, 'dir')
    mkdir(savePath);
end
save(sprintf("%s\\pTemplateMatching.mat", savePath), 'pTemplateMatching')

selUnits = find(hasRaw & isIntan); useRawCache = false;
% selUnits = find(ismember(eu.getName(), coolUnitNames)); useRawCache = true;
unitNames = eu.getName();

clear layout
layout.fig = figure(Unit='normalized', Position=[0 0 1 1]);
layout.h = [4, 10, 4, 1];
layout.w = [2, 5, 2];
layout.tl = tiledlayout(layout.fig, sum(layout.h), sum(layout.w), TileSpacing='tight', Padding='compact');
layout.ax.waveform = nexttile(layout.tl, [layout.h(1), layout.w(1)]);
layout.ax.etaLick = nexttile(layout.tl, [layout.h(1), layout.w(2)]);
layout.ax.etaCircLick = nexttile(layout.tl, [layout.h(1), layout.w(3)]);
layout.ax.raw = nexttile(layout.tl, [layout.h(2), sum(layout.w)]);
layout.ax.binnedWaveforms = nexttile(layout.tl, [layout.h(3), sum(layout.w)]);
layout.ax.binnedHistograms = nexttile(layout.tl, [layout.h(4), sum(layout.w)]);


clear spikesFiltered spikesRaw
if useRawCache
    if ~exist('rawCache', 'var')
        rawCache(length(eu)) = struct(index=[], name=[], data=[], t=[], trials=[]);
    end
else
    clear rawCache
end
spikesFiltered = struct(index=[], name=[], sampleIndex=[], timestamps=[], waveforms=[], waveformTimestamps=[], isUnit=[]);
spikesRaw = struct(index=[], name=[], sampleIndex=[], timestamps=[], waveforms=[], waveformTimestamps=[], isUnit=[]);

for iUnit = 1:length(selUnits)
% for iUnit = 179
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
        spikesRaw.index = iEu;
        spikesRaw.name = eu(iEu).getName();
        [spikesRaw.sampleIndex, spikesRaw.timestamps, spikesRaw.waveforms, spikesRaw.waveformTimestamps] = spikeDetect(raw, SampleRate=fs, NumSigmas=2.5, NumSigmasReturn=1.25, NumSigmasReject=20, WaveformWindow=[-0.5, 0.5]);
        spikesRaw.isUnit = ismember(round(spikesRaw.timestamps*fs), round(stAbs*fs));
        fprintf('(%i/%i detected spikes matched timestamps of %i eu spikes)...', nnz(spikesRaw.isUnit), length(spikesRaw.isUnit), length(stAbs))
        spikeTemplateRaw = mean(spikesRaw.waveforms, 1, 'omitnan');
        maxMicroVolts = max(abs(spikeTemplateRaw))*2;
        fprintf('Done (%.2f s)\n', toc(tTic));
    
        % Detect spikes in filtered data
        tTic = tic();
        fprintf('\tDetecting spikes in filtered data...');
        spikesFiltered.index = iEu;
        spikesFiltered.name = eu(iEu).getName();
        [spikesFiltered.sampleIndex, spikesFiltered.timestamps, spikesFiltered.waveforms, spikesFiltered.waveformTimestamps] = spikeDetect(filtered, SampleRate=fs, NumSigmas=2.5, NumSigmasReturn=1.25, NumSigmasReject=20, WaveformWindow=[-0.5, 0.5], MaxMicroVolts=maxMicroVolts);
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
        [xRawAligned, tAligned, stRawAligned, trials] = parseRaw(eu(iEu), raw, eu(iEu).getTrials('circlick_naive', minInterval=0.05, maxInterval=0.20), spikeTimes=oldSpikeTimes, window=[-0, 0.15], alignTo='Start', randomTrials=true, nTrials=15, timestampMode='relative');
        [xFilteredAligned, ~, stFilteredAligned] = parseRaw(eu(iEu), filtered, trials, spikeTimes=spikesFiltered.timestamps(spikesFiltered.isUnit), window=[-0, 0.15], alignTo='Start', randomTrials=false, nTrials='all', timestampMode='relative');
    
        plotRaw(layout.ax.raw, xFilteredAligned, tAligned, stFilteredAligned, plotSpikes=true, spacing=250, xRaw=xRawAligned, stRaw=stRawAligned);
        xlim(layout.ax.raw, [0, 150])
        xticks(layout.ax.raw, [0, 10, 50, 100, 150])
        xline(layout.ax.raw, 0:10:150, 'k--')
    
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
    
        hold(layout.ax.etaCircLick, 'on')
        plot(layout.ax.etaCircLick, etaTemp.circLickNaiveRaw.t, etaTemp.circLickNaiveRaw.X, 'blue', DisplayName='raw')
        plot(layout.ax.etaCircLick, etaTemp.circLickNaiveFiltered.t, etaTemp.circLickNaiveFiltered.X, 'red', DisplayName='filtered')
        xlim(layout.ax.etaCircLick, [0, 2*pi])
        xticks(layout.ax.etaCircLick, [0, pi, 2*pi])
        xticklabels(layout.ax.etaCircLick, {'0', '\pi', '2\pi'})
        ylabel(layout.ax.etaCircLick, 'spikes/s')
        xlabel(layout.ax.etaCircLick, 'lick phase')
        legend(layout.ax.etaCircLick)
        title(layout.ax.etaCircLick, eu(iEu).getName(), Interpreter='none')
    
        hold(layout.ax.etaLick, 'on')
        plot(layout.ax.etaLick, etaTemp.lickRaw.t, etaTemp.lickRaw.X./0.025, 'blue', DisplayName='raw')
        plot(layout.ax.etaLick, etaTemp.lickFiltered.t, etaTemp.lickFiltered.X./0.025, 'red', DisplayName='filtered')
        xlabel(layout.ax.etaLick, 'Time to self-timed lick (s)')
        ylabel(layout.ax.etaLick, 'spikes/s')
        legend(layout.ax.etaLick)

        ylim(layout.ax.etaLick, 'auto')
        ylim(layout.ax.etaCircLick, 'auto')
        drawnow()
        ylim([layout.ax.etaLick, layout.ax.etaCircLick], ...
            [min(arrayfun(@(ax) ax.YLim(1), [layout.ax.etaLick, layout.ax.etaCircLick])), ...
            max(arrayfun(@(ax) ax.YLim(2), [layout.ax.etaLick, layout.ax.etaCircLick]))]);

        title(layout.tl, eu(iEu).getName(), Interpreter='none')
        fprintf('Done (%.2f s)\n', toc(tTic));
        
        % Plot waveforms by bin (aligned to lick)
        cla(layout.ax.binnedWaveforms)
        cla(layout.ax.binnedHistograms)
        hold(layout.ax.binnedWaveforms, 'on')
        hold(layout.ax.binnedHistograms, 'on')
        trials = eu(iEu).getTrials('circlick_naive', minInterval=0.05, maxInterval=0.20);
        edges = 0:0.01:0.15;
        h = gobjects(2, 1);
        for iBin = 1:length(edges)-1
            window = edges(iBin:iBin+1);
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

            text(layout.ax.binnedWaveforms, tWave(16), min(spikeTemplate - 2*sigma), sprintf('%i', size(waveforms, 1)), HorizontalAlignment='center', VerticalAlignment='bottom')
            % text(layout.ax.binnedWaveforms, tWave(16), 0, sprintf('%i', size(waveforms, 1)), HorizontalAlignment='center', VerticalAlignment='bottom')

            % Embed a histogram of distanceToTemplate
            [ratioN, ~] = histcounts(log2(distRatio), [0:0.05:1.95, Inf]);
            ratioEdges = (0:0.05:2) + (iBin-1)*2;
            histogram(layout.ax.binnedHistograms, BinEdges=ratioEdges, BinCounts=ratioN, EdgeColor='black', FaceColor='black');
        end
        xticks(layout.ax.binnedWaveforms, 0:10:150)
        xline(layout.ax.binnedWaveforms, 0:10:150, 'k--')
        xlim(layout.ax.binnedWaveforms, [0, 150])
        ylim(layout.ax.binnedWaveforms, [min(spikeTemplate - 2*sigma), max(spikeTemplate + 2*sigma)])
        ylabel(layout.ax.binnedWaveforms, '\muV')
        xlabel(layout.ax.binnedWaveforms, 'Time (ms)')
        legend(layout.ax.binnedWaveforms, h, Orientation='horizontal', Location='northeast')

        xlim(layout.ax.binnedHistograms, [0, 2*(length(edges)-1)])
        xticks(layout.ax.binnedHistograms, 0:1:2*(length(edges)-1))
        xticklabels(layout.ax.binnedHistograms, string(repmat([0, 1], 1, length(edges)-1)))
        xline(layout.ax.binnedHistograms, 0:2:2*(length(edges)-1), 'k--')
        ylabel(layout.ax.binnedHistograms, 'SpikeCount')
        xlabel(layout.ax.binnedHistograms, 'log_2(distToNoise/distToTemplate)')

        clear trials edges iBin window iBin spikeTimesInBin t waveforms IBinned distToSpikeTemplate distToNoiseTemplate distRatio waveformBinSize
        
        tTic = tic();
        fprintf('\tSaving plots...');
        print(layout.fig, sprintf('%s\\%s.png', savePath, eu(iEu).getName()), '-dpng', '-r0')
        fprintf('Done (%.2f s)\n', toc(tTic));
    
        index = iEu;
        name = eu(iEu).getName();
        tTic = tic();
        fprintf('\tSaving data...');
        save(sprintf("%s\\%s.mat", savePath, name), 'index', 'name', 'spikesFiltered', 'spikesRaw', 'oldSpikeTimes', '-v7.3')
        fprintf('Done (%.2f s)\n', toc(tTic));
    catch ME
        warning('Error processing unit %i (%i)', iEu, iUnit);
        warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
    end
end

clear useRawCache iUnit iEu tTic raw filtered stAbs spikeTemplateRaw maxMicroVolts spikeTemplate tWaveform noiseTemplate distToSpikeTemplate distToNoiseTemplate residuals sigma pOutlier I iWave xRawAligned tAligned stRawAligned trials xFilteredAligned stFilteredAligned etaTemp
clear index name tTic
clear trials edges iBin window tempTrials inBin spikeTimesInBin waveforms tWave mu distToSpikeTemplate distToNoiseTemplate distRatio IBinned waveformBinSize iWaveform

%% Load a unit
iUnit = 129;
iEu = selUnits(iUnit);
load(sprintf("%s\\%s.mat", savePath, eu(iEu).getName()))


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
    xlabel(ax, 'Time (ms)')
    ylabel(ax, 'Voltage (uV)')

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
