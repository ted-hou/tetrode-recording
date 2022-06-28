fdir = 'C:\SERVER\daisy14\daisy14_20220506';
br = readBlackRock(fdir);
tr = readIntan(fdir);
appendBRtoTR(br, tr);

tr.SpikeClusterAutoReorder([], 'range')
tr.PlotAllChannels('plotMethod', 'mean')

%%
tr.PlotChannel([], 'Reference', 'CueOn', 'Event', 'PressOn', 'Exclude', 'LickOn', 'Event2', '', 'Exclude2', '', 'RasterXLim', [-6, 1], 'ExtendedWindow', [-1, 1], 'PlotStim', true);

%%
TetrodeRecording.BatchSave(tr, 'Prefix', 'tr_sorted_', 'DiscardData', false);

%%
load('C:\SERVER\daisy14\daisy14_20220506\SpikeSort\tr_sorted_daisy14_20220506_220506_210458.mat')

%%
clear ar bsr probeMap stim

ar = AcuteRecording(tr, 'D1-Cre;Dlx-Flp;Ai80');
stim = ar.extractAllPulses('A', 0.5);
probeMap = ar.importProbeMap('back', 1300, -4300, -3280);
[bsrFull, m, s] = ar.binStimResponse([]);
% ar.plotPSTHByStimCondition(bsrFull, 'CLim', [-.5, .5]);

%% Calculate main effect, use 1d stat (mean or peak)
window = [0 0.05];
[bsr, I] = ar.selectStimResponse(bsrFull, 'Light', 0.5, 'Duration', 0.01);
[statsPeak, bsr] = AcuteRecording.summarizeStimResponse(bsr, 'Window', window, 'Method', 'peak', 'Normalized', true);
statsMean = ar.summarizeStimResponse(bsr, 'Window', window, 'Method', 'mean', 'Normalized', true);
% ar.plotPSTHByStimCondition(bsr, 'CLim', [-.5, .5]);

%% Plot and compare mean vs peak response. Check that mean and peak has same sign.
figure, ax = axes();
hold(ax, 'on')
plot(ax, statsPeak, 'DisplayName', sprintf('peak %i-%ims', window(1)*1000, window(2)*1000));
plot(ax, statsMean, 'DisplayName', sprintf('mean %i-%ims', window(1)*1000, window(2)*1000));
ylabel(ax, 'Activity change (a.u.)')
xlabel(ax, 'Unit')
hold(ax, 'off');
legend(ax);
clear ax

%% Plot heatmap per stim condition
coords = ar.getProbeCoords([bsr.channel]);
AcuteRecording.plotMap(coords, statsPeak);

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

