%% Define animal/expname
fdir = uigetdir('C:\SERVER\');
expName = strsplit(fdir, '\');
expName = expName{end};
animalName = strsplit(expName, '_');
animalName = animalName{1};

%% Read auto-sorted(crude) intan data and combine with analog data from blackrock
br = readBlackRock(fdir);
tr = readIntan(fdir);
appendBRtoTR(br, tr);

tr.SpikeClusterAutoReorder([], 'range')
tr.PlotAllChannels('plotMethod', 'mean')

%% Do manual spike sorting
tr.PlotChannel([], 'Reference', 'CueOn', 'Event', 'PressOn', 'Exclude', 'LickOn', 'Event2', '', 'Exclude2', '', 'RasterXLim', [-6, 1], 'ExtendedWindow', [-1, 1], 'PlotStim', true);

%% Save manually sorted units
TetrodeRecording.BatchSave(tr, 'Prefix', 'tr_sorted_', 'DiscardData', false);

%% Load manually sorted results
f = dir(sprintf('C:\\SERVER\\%s\\%s\\SpikeSort\\tr_sorted_%s*.mat', animalName, expName, expName));
load(sprintf('%s\\%s', f.folder, f.name));
clear f

%% Extract stim response and save to file
clear ar bsr probeMap stim sessionInfo
sessionInfo.expName = expName;
sessionInfo.strain = 'A2A-Cre'; % A2A-Cre or D1-Cre;Dlx-Flp;Ai80
sessionInfo.orientation = 'front'; % Whether probe is anterior (front) or posterior (back) to the headstage. Usually this is back, unless recording in anterior SNr.
sessionInfo.ml = 1300; % Center of probe. Negative for left hemisphere.
sessionInfo.dv = -4600; % From surface of brain, not bregma. Assumes Bregma is 200um above brain surface unless otherwise specified for `importProbeMap`.
sessionInfo.ap = -3280+350;
sessionInfo.firstFiber = 'B'; % 'A' or 'B', name of the first patch cable
sessionInfo.firstLight = 0.1; % 

window = [0 0.05];

ar = AcuteRecording(tr, sessionInfo.strain);
stim = ar.extractAllPulses(tr, sessionInfo.firstFiber, sessionInfo.firstLight);
probeMap = ar.importProbeMap(sessionInfo.orientation, sessionInfo.ml, sessionInfo.dv, sessionInfo.ap);
ar.binStimResponse(tr, [], 'Store', true);
ar.summarize(ar.bsr, 'peak', [0, 0.05], 'Store', true);

ar.save();
save(sprintf('C:\\SERVER\\%s\\%s\\AcuteRecording\\sessionInfo_%s.mat', animalName, expName, expName), 'sessionInfo');

%% 
ar = AcuteRecording.load(sprintf('C:\\SERVER\\%s\\%s\\AcuteRecording\\ar_%s.mat', animalName, expName, expName));


% %% Fix probemap ML
% 
% expName = 'daisy15_20220511';
% animalName = strsplit(expName, '_');
% animalName = animalName{1};
% 
% sessionInfo.expName = expName;
% sessionInfo.strain = 'A2A-Cre'; % A2A-Cre or D1-Cre;Dlx-Flp;Ai80
% sessionInfo.orientation = 'back'; % Whether probe is anterior (front) or posterior (back) to the headstage. Usually this is back, unless recording in anterior SNr.
% sessionInfo.ml = 1300; % Center of probe. Negative for left hemisphere.
% sessionInfo.dv = -4500; % From surface of brain, not bregma. Assumes Bregma is 200um above brain surface unless otherwise specified for `importProbeMap`.
% sessionInfo.ap = -3280;
% sessionInfo.firstFiber = 'B'; % 'A' or 'B', name of the first patch cable
% sessionInfo.firstLight = 0.1; % 
% 
% ar = AcuteRecording.load(sprintf('C:\\SERVER\\%s\\%s\\AcuteRecording\\ar_%s.mat', animalName, expName, expName));
% probeMap = ar.importProbeMap(sessionInfo.orientation, sessionInfo.ml, sessionInfo.dv, sessionInfo.ap);
% ar.save();
% save(sprintf('C:\\SERVER\\%s\\%s\\AcuteRecording\\sessionInfo_%s.mat', animalName, expName, expName), 'sessionInfo');
% 
% [bsr, ~] = ar.selectStimResponse('Light', 0.5, 'Duration', 0.01);
% ar.plotMapByStimCondition(bsr, [0.25, 3], 0.25, 'peak', [0 0.05], 0.25, 'UseSignedML', true);
% ar.plotMapByStimCondition(bsr, [0.25, 3], 0.25, 'peak', [0 0.05], 0.25);
%% Load AR file (much faster than loading TR and rebuilding AR)
ar = AcuteRecording.load(sprintf('C:\\SERVER\\%s\\%s\\AcuteRecording\\ar_%s.mat', animalName, expName, expName));

%% One SNr map per Str site (1 figure)
[bsr, ~] = ar.selectStimResponse('Light', 0.5, 'Duration', 0.01);
[statsPeak, conditions] = ar.summarize(bsr, 'peak', window);
statsFirstPeak = ar.summarize(bsr, 'firstPeak', window, 0.25);
statsMean = ar.summarize(bsr, 'mean', window, 0.25);

% Plot and compare mean vs peak response. Check that mean and peak has same sign.
nConditions = length(conditions);
figure;
for i = 1:nConditions
    ax = subplot(nConditions, 1, i);
    hold(ax, 'on')
    plot(ax, statsPeak(:, i), 'DisplayName', sprintf('peak %i-%ims', window(1)*1000, window(2)*1000));
    plot(ax, statsFirstPeak(:, i), 'DisplayName', sprintf('first Peak %i-%ims', window(1)*1000, window(2)*1000));
    plot(ax, statsMean(:, i), 'DisplayName', sprintf('mean %i-%ims', window(1)*1000, window(2)*1000));
    ylabel(ax, '\DeltaActivity')
    xlabel(ax, 'Unit')
    title(ax, conditions(i).label)
    hold(ax, 'off');
    legend(ax);
end
clear ax i nConditions statsPeak statsFirstPeak statsMean

% Plot probe map per stim condition
ar.plotMapByStimCondition(bsr, [0.25, 3], 0.25, 'peak', window, 0.25);

%% One Str map per SNr unit (WARNING: MANY FIGURES)
ar.plotPSTHByStimCondition(bsr, 'CLim', [-.5, .5]);

%% Read multiple files, pool stats and plot in same map.
% plotD1(-3280)
% plotD1(-2930)
plotD1([])

function plotD1(ap)
    if nargin < 1
        ap = [];
    end

    window = [0, 0.05];

    fdir = 'C:\SERVER\Experiment_Galvo_D1Cre;DlxFlp;Ai80\AcuteRecording';
    load(sprintf('%s\\crit.mat', fdir), 'crit');
    load(sprintf('%s\\sessionInfo.mat', fdir), 'sessionInfo');
    ar = AcuteRecording.load(fdir);
    
    if ~isempty(ap)
        sel = [sessionInfo.ap] == ap;
        ar = ar(sel);
        crit = crit(sel);
        sessionInfo = sessionInfo(sel);
    end

    for i = 1:length(ar)
        ar(i).importProbeMap(sessionInfo(i).orientation, sessionInfo(i).ml, sessionInfo(i).dv, sessionInfo(i).ap);
        bsr{i} = ar(i).selectStimResponse('Light', 2, 'Duration', crit(i).duration);
        [stats{i}, conditions{i}] = ar(i).summarize(bsr{i}, 'peak', window);
        titleText = ar(i).plotMapByStimCondition(bsr{i}, [0.25, 1], 0.25, 'peak', window, 0.25);
%         print(sprintf('%s (%.2f AP).png', titleText, ap/1000), '-dpng');
    end
    titleText = ar.plotMapByStimCondition(bsr, [0.25, 1], 0.25, 'peak', window, 0.25);
%     print(sprintf('%s (%.2f AP).png', titleText, ap/1000), '-dpng');
end

%%
function tr = readIntan(fdir)
    if nargin < 1
        fdir = uigetdir('C:\SERVER\');
    end
    expName = strsplit(fdir, '\');
    expName = expName{end};
    tr = TetrodeRecording.BatchLoadSimple(expName, true);
end

function br = readBlackRock(fdir)
    if nargin < 1
        fdir = uigetdir('C:\SERVER\');
    end

    % Read BlackRock files
    nsx = dir(sprintf('%s\\*.ns?', fdir));
    nev = dir(sprintf('%s\\*.nev', fdir));
    nsx = openNSx(sprintf('%s\\%s', nsx.folder, nsx.name));
    nev = openNEV(sprintf('%s\\%s', nev.folder, nev.name), 'nosave', 'nomat');

    % Parse digital data
    digitalChannels = {'Laser', 0; 'Galvo', 1};
    digitalData = flip(dec2bin(nev.Data.SerialDigitalIO.UnparsedData), 2);
    digitalTimestamps = nev.Data.SerialDigitalIO.TimeStampSec;
    digitalSampleIndices = nev.Data.SerialDigitalIO.TimeStamp;
    digitalSamplingRate = nev.MetaTags.TimeRes;
    digitalStartDateTime = nev.MetaTags.DateTimeRaw;


    for iChannel = 1:size(digitalChannels, 1)
        iBit = digitalChannels{iChannel, 2} + 1;
        channelName = digitalChannels{iChannel, 1};
        [digitalEvents.([channelName, 'On']), digitalEvents.([channelName, 'Off'])] = TetrodeRecording.FindEdges(transpose(digitalData(:, iBit)), digitalTimestamps);
    end
    clear iBit channelName

    % Parse analog data
    analogChannels = {'Laser', 1; 'Galvo', 2};
    analogConversionFactor = double([nsx.ElectrodesInfo.MaxAnalogValue])./double([nsx.ElectrodesInfo.MaxDigiValue]);
    analogData = double(nsx.Data) .* analogConversionFactor';
    analogSamplingRate = nsx.MetaTags.SamplingFreq;
    analogTimestamps = (0:nsx.MetaTags.DataPoints - 1) / analogSamplingRate;
    
    br.digitalEvents = digitalEvents;
    br.digitalChannels = digitalChannels;
    br.analogData = analogData;
    br.analogChannels = analogChannels;
    br.analogTimestamps = analogTimestamps;
end

% Append blackrock/galvo digital events to TR object (intan recording),
% align to intan file timestamps.
function appendBRtoTR(br, tr)
    brTimeOffset = tr.DigitalEvents.StimOn - br.digitalEvents.LaserOn;
    tr.DigitalEvents.LaserOn = br.digitalEvents.LaserOn + brTimeOffset;
    tr.DigitalEvents.LaserOff = br.digitalEvents.LaserOff + brTimeOffset;
    tr.DigitalEvents.GalvoOn = br.digitalEvents.GalvoOn + interp1(br.digitalEvents.LaserOn, brTimeOffset, br.digitalEvents.GalvoOn, 'linear', 'extrap');
    tr.DigitalEvents.GalvoOff = br.digitalEvents.GalvoOff + interp1(br.digitalEvents.LaserOn, brTimeOffset, br.digitalEvents.GalvoOff, 'linear', 'extrap');
    tr.AnalogIn.Channels = br.analogChannels;
    tr.AnalogIn.Timestamps = br.analogTimestamps + interp1(br.digitalEvents.LaserOn, brTimeOffset, br.analogTimestamps, 'linear', 'extrap');
    tr.AnalogIn.Data = br.analogData;
end

