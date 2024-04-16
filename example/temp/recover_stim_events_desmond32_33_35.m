expNames = { ...
    'desmond32_20231128'; ...
    'desmond35_20231130'; ...
    };

warning('This fix only works for desmond32_20231128 and desmond35_20231130.')

for iExp = 1:length(expNames)
    tr(iExp) = TetrodeRecording.BatchLoadSimple(expNames{iExp}, false, 'tr_');
    % Load tce(iExp) (TwoColorExperiment)
    file = dir(sprintf('C:\\SERVER\\%s\\%s\\twocolor_%s*.mat', tr(iExp).GetAnimalName(), expNames{iExp}, expNames{iExp}));
    assert(length(file) == 1)
    assert(file.folder(end) ~= '\')
    file = sprintf('%s\\%s', file.folder, file.name);
    T = load(file);
    tce(iExp) = T.obj;
    clear file
    
    % 1. TetrodeRecording.tRef = Threshold Intan TetrodeRecording.AnalogIn (laserModBlue/laserModOrange) into digital events
    % Timestamp naming convention:
    %   tXXX: ephys time (s)
    %   millisXXX: arduino time (ms, int)
    %   datetimeXXX: windows time (MATLAB datetime, might be inaccurate for received serial messages due to buffering, should be accurate for MATLAB output: laser modulation, mirror controls)
    assert(strcmpi(tr(iExp).AnalogIn.ChannelNames{1}, 'LaserModBlue'))
    assert(strcmpi(tr(iExp).AnalogIn.ChannelNames{2}, 'LaserModOrange'))
    
    nPulsestce(iExp) = sum(arrayfun(@(x) x.nPulses, [tce(iExp).Log.params]));
    
    eventIndexStimOn = find(strcmp('OPTOGEN_STIM_1_ON', tce(iExp).LaserArduino.EventMarkerNames));
    eventIndexStimOff = find(strcmp('OPTOGEN_STIM_1_OFF', tce(iExp).LaserArduino.EventMarkerNames));
    stimOn.millis = tce(iExp).LaserArduino.EventMarkersUntrimmed(tce(iExp).LaserArduino.EventMarkersUntrimmed(:, 1) == eventIndexStimOn, 2);
    stimOff.millis = tce(iExp).LaserArduino.EventMarkersUntrimmed(tce(iExp).LaserArduino.EventMarkersUntrimmed(:, 1) == eventIndexStimOff, 2);
    stimOn.datetime = datetime(tce(iExp).LaserArduino.EventMarkersUntrimmed(tce(iExp).LaserArduino.EventMarkersUntrimmed(:, 1) == eventIndexStimOn, 3), ...
        ConvertFrom='datenum');
    stimOff.datetime = datetime(tce(iExp).LaserArduino.EventMarkersUntrimmed(tce(iExp).LaserArduino.EventMarkersUntrimmed(:, 1) == eventIndexStimOff, 3), ...
        ConvertFrom='datenum');
    
    % Windows (input) timestamps are noisy, this can be shown by comparing arduino-recorded stim durations vs. windows recorded stim durations (logging datetime('now') upon receiving arduino messages, i.e., BytesAvailableFcn)
    stimDuration.millis = stimOff.millis - stimOn.millis;
    stimDuration.datetime = milliseconds(stimOff.datetime - stimOn.datetime);
    histogram(stimDuration.millis - stimDuration.datetime, -50:5:50)
    
    % Now we check windows output timestamps agains Intan analog read (laser modulation signal)
    % We expect a 8100ms delay between laserOn and trainOn (both are initiated
    % by MATLAB, with a pause(8.0) and pause(0.1) in between)
    blueLaserOn.t = tr(iExp).AnalogIn.Timestamps(strfind(tr(iExp).AnalogIn.Data(1, :) > 0.1, [0, 1])); nnz(blueLaserOn.t)
    isBlueTrain = [tce(iExp).Log.wavelength]==473;
    try
        assert(nnz(blueLaserOn.t) == nnz(isBlueTrain))
    catch
        df1 = diff(blueLaserOn.t);
        df2 = seconds(diff([tce(iExp).Log(isBlueTrain).laserOnTime]));
        if length(df1) == length(df2) + 1
            dfa = df1(1:end-1) - df2;
            dfb = df1(2:end) - df2;
            if mean(dfa) < mean(dfb)
                error()
            else
                warning('AnalogIn/blueLaserOn.t has %i events, tce(iExp).Log has %i events, we ignore the first one.', nnz(blueLaserOn.t), nnz(isBlueTrain))
                blueLaserOn.t = blueLaserOn.t(2:end);
            end
            assert(nnz(blueLaserOn.t) == nnz(isBlueTrain))
        else
            error()
        end
    end

    blueLaserOn.datetime = [tce(iExp).Log(isBlueTrain).trainOnTime];
    
    t0.t = blueLaserOn.t(1);
    t0.datetime = blueLaserOn.datetime(1);
    
    iEnd = find(isnan(tce(iExp).LaserArduino.AnalogOutputEvents(:, 1)), 1, 'first');
    % The first two events are turning both lasers off, followed by the
    % tce(iExp).Log stim trains
    if (iEnd-1)/2 == length(tce(iExp).Log) + 1 && all(tce(iExp).LaserArduino.AnalogOutputEvents(1, 1:2) == [1, 0]) && all(tce(iExp).LaserArduino.AnalogOutputEvents(2, 1:2) == [2, 0])
        iStart = find(tce(iExp).LaserArduino.AnalogOutputEvents(:, 1) == 1 & tce(iExp).LaserArduino.AnalogOutputEvents(:, 2) > 0, 1, 'first');
        assert(tce(iExp).LaserArduino.AnalogOutputEvents(iStart+1, 2) == 0)
    % The first batch of events are laser calibration, but the last event
    % corresponds to the last train in tce(iExp).Log
    else
        iStart = iEnd - length(tce(iExp).Log)*2;
        assert(tce(iExp).LaserArduino.AnalogOutputEvents(iStart, 1) == 1 && tce(iExp).LaserArduino.AnalogOutputEvents(iStart, 2) > 0 && tce(iExp).LaserArduino.AnalogOutputEvents(iStart+1, 2) == 0)
    end
    t0.millis = tce(iExp).LaserArduino.AnalogOutputEvents(iStart, 3);
    stimOn.t = t0.t + (double(stimOn.millis) - double(t0.millis)).*1e-3; % TADA! This is it
    stimOff.t = t0.t + (double(stimOff.millis) - double(t0.millis)).*1e-3;

    tr(iExp).DigitalEvents.StimOn = stimOn.t(stimOn.t > 0);
    tr(iExp).DigitalEvents.StimOff = stimOff.t(stimOff.t > 0);
    assert(length(tr(iExp).DigitalEvents.StimOn) == length(tr(iExp).DigitalEvents.StimOff))

    TetrodeRecording.BatchSave(tr(iExp), Prefix='tr_fixed_', DiscardData=false);
end

%% Categorize stim pulses into blue and orange.
for iExp = 1:length(tr)
    % Inter-train intervals can be used as a signature to see if MATLAB
    % based stimOn (tce.Log.trainOnTime)
    % agrees with arduino-based stimOn (DigitalEvents.StimOn).
    iti1 = diff(tr(iExp).DigitalEvents.StimOn(1:10:end));
    iti2 = seconds(diff([tce(iExp).Log.trainOnTime]))';

    % Cull extra arduino-based stim events
    if length(iti1) == length(iti2)
        assert(all(abs(iti1 - iti2) < 0.005))
    else
        % desmond35_20231130 log: first 10 pulses did not work because
        % arduino was stuck, arduino (iti1) should have 10 extra pulses at
        % the beginning, we will discard those
        assert(strcmpi(tr(iExp).GetExpName(includeSuffix=false), 'desmond35_20231130'))
        iti1(1) = [];
        assert(length(iti1) == length(iti2))
        assert(all(abs(iti1 - iti2) < 0.005))

        tr(iExp).DigitalEvents.StimOn(1:10) = [];
        tr(iExp).DigitalEvents.StimOff(1:10) = [];
    end

    nPulses = length(tr(iExp).DigitalEvents.StimOn);
    assert(nPulses == sum(arrayfun(@(x) x.params.nPulses, tce(iExp).Log)));
    pulseWavelength = NaN(length(nPulses), 1);
    iPulse = 0;
    for iTrain = 1:length(tce(iExp).Log) - 3
        pulseWavelength(iPulse + 1:iPulse + tce(iExp).Log(iTrain).params.nPulses) = tce(iExp).Log(iTrain).wavelength;
        iPulse = iPulse + tce(iExp).Log(iTrain).params.nPulses;
    end
    for iTrain = length(tce(iExp).Log) - 2:length(tce(iExp).Log)
        pulseWavelength(iPulse + 1:iPulse + tce(iExp).Log(iTrain).params.nPulses) = 0;
        iPulse = iPulse + tce(iExp).Log(iTrain).params.nPulses;
    end
    tr(iExp).DigitalEvents.StimOnBlue = tr(iExp).DigitalEvents.StimOn(pulseWavelength==473);
    tr(iExp).DigitalEvents.StimOffBlue = tr(iExp).DigitalEvents.StimOff(pulseWavelength==473);
    tr(iExp).DigitalEvents.StimOnOrange = tr(iExp).DigitalEvents.StimOn(pulseWavelength==593);
    tr(iExp).DigitalEvents.StimOffOrange = tr(iExp).DigitalEvents.StimOff(pulseWavelength==593);
    tr(iExp).DigitalEvents.StimOnShutterControl = tr(iExp).DigitalEvents.StimOn(pulseWavelength==0);
    tr(iExp).DigitalEvents.StimOffShutterControl = tr(iExp).DigitalEvents.StimOff(pulseWavelength==0);

    assert(ismember(tr(iExp).GetExpName(includeSuffix=false), {'desmond32_20231128', 'desmond35_20231130'}))
    

    fprintf('Exp %i (%s) passed test.\n', iExp, tr(iExp).GetExpName(includeSuffix=false))
    TetrodeRecording.BatchSave(tr(iExp), Prefix='tr_fixed_', DiscardData=false);
end