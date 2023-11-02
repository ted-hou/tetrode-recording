files = { ...
    '\\research.files.med.harvard.edu\neurobio\NEUROBIOLOGY SHARED\Assad Lab\Lingfeng\Data\daisy17\daisy17_20230801\SpikeSort\tr_sorted_daisy17_20230801_230801_154516.mat', ...
    '\\research.files.med.harvard.edu\neurobio\NEUROBIOLOGY SHARED\Assad Lab\Lingfeng\Data\daisy17\daisy17_20230804\SpikeSort\tr_sorted_daisy17_20230804_230804_153713.mat', ...
    '\\research.files.med.harvard.edu\neurobio\NEUROBIOLOGY SHARED\Assad Lab\Lingfeng\Data\daisy18\daisy18_20230731\SpikeSort\tr_sorted_daisy18_20230731_230731_161918.mat', ...
    '\\research.files.med.harvard.edu\neurobio\NEUROBIOLOGY SHARED\Assad Lab\Lingfeng\Data\daisy18\daisy18_20230728\SpikeSort\tr_sorted_daisy18_20230728_230728_151933.mat', ...
    '\\research.files.med.harvard.edu\neurobio\NEUROBIOLOGY SHARED\Assad Lab\Lingfeng\Data\daisy18\daisy18_20230727\SpikeSort\tr_sorted_daisy18_20230727_230727_180050.mat', ...
    '\\research.files.med.harvard.edu\neurobio\NEUROBIOLOGY SHARED\Assad Lab\Lingfeng\Data\daisy18\daisy18_20230802\SpikeSort\tr_sorted_daisy18_20230802_230802_140526.mat', ...
    '\\research.files.med.harvard.edu\neurobio\NEUROBIOLOGY SHARED\Assad Lab\Lingfeng\Data\desmond31\desmond31_20230801\SpikeSort\tr_sorted_desmond31_20230801_230801_171746.mat', ...
    '\\research.files.med.harvard.edu\neurobio\NEUROBIOLOGY SHARED\Assad Lab\Lingfeng\Data\desmond31\desmond31_20230804\SpikeSort\tr_sorted_desmond31_20230804_230804_173252.mat' ...
    };

for iTr = 1:length(files)
    S = load(files{iTr});
    tr = S.tr;
    channels = [tr.Spikes.Channel];
    for iChn = channels
        refChn = max(tr.Spikes(iChn).Cluster.Classes);
        for iUnit = 1:refChn - 1
            hFigure = tr.PlotChannel(iChn, 'Reference', 'TRIAL_START', 'Event', 'PressOn', 'Exclude', '', 'Event2', '', 'Exclude2', '', ...
                'RasterXLim', [-4, 1], 'ExtendedWindow', [-0, 1], 'WaveformYLim', [-200, 200], 'PlotStim', false, 'Bins', 3, ...
                'Clusters', iUnit, 'ReferenceCluster', refChn, 'PrintMode', true, 'FontSize', 11);
            tr.GUISavePlot([], [], hFigure, Filename=['C:\SERVER\Figures\Spontaneous\', tr.GetExpName(), '_Chn', num2str(iChn), 'Unit', num2str(iUnit)], ...
                Reformat='RasterAndPETHAndWaveform')
            close all
        end
    end
    clear S tr
end

%% Create EphysUnit objects from TetrodeRecording.
for iTr = 1:length(files)
    clear S tr ar eu
    try
        S = load(files{iTr});
        tr = S.tr;
        ar = AcuteRecording(tr, 'N/A');
        [~, ~, ~, tRef, tMove] = ar.binMoveResponse(tr, 'press_spontaneous', Store=true);
        eu = EphysUnit(ar, savepath='C:\SERVER\Units\acute_spontaneous_reach', ...
            cullITI=false, readWaveforms=false, isSpontaneousPress=true, tr=tr, extendedWindow=[0, 0]);
    catch ME
        warning('Error while processing file %g (%s)', iTr, files{iTr});
    end
end

clear S tr ar eu iTr

%% 1. Load data
%% 1.1. Load acute EU objects (duplicates already removed)
eu = EphysUnit.load('C:\SERVER\Units\acute_spontaneous_reach'); 

% Remove multiunit detected by ISI test.
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

%
% 1.2. Load Video Tracking Data (vtd) and ArduinoConnection (ac), and group into experiments
clearvars -except eu
exp = CompleteExperiment2(eu);

% 1.3 Align video and ephys timestamps
exp.alignTimestamps();

%% Average movement trajectories across lever-touch "trials"
close all
for iExp = 1:length(exp)
    touchTimes = [exp(iExp).eu(1).Trials.Press.Stop];
    F = cell(length(touchTimes), 1);
    for iTrial = 1:length(touchTimes)
        F{iTrial} = exp(iExp).getFeatures(timestamps=touchTimes(iTrial)-4:1/30:touchTimes(iTrial), features={'handL', 'handR', 'nose', 'jaw'}, stats={'xPos', 'yPos'});
    end
    
    data = cellfun(@(f) f.handL_xPos, F, UniformOutput=false);
    handL_x = cat(2, data{:});
    
    data = cellfun(@(f) f.handL_yPos, F, UniformOutput=false);
    handL_y = cat(2, data{:});
    
    data = cellfun(@(f) f.handR_xPos, F, UniformOutput=false);
    handR_x = cat(2, data{:});
    
    data = cellfun(@(f) f.handR_yPos, F, UniformOutput=false);
    handR_y = cat(2, data{:});
    
    data = cellfun(@(f) f.nose_xPos, F, UniformOutput=false);
    nose_x = cat(2, data{:});
    
    data = cellfun(@(f) f.nose_yPos, F, UniformOutput=false);
    nose_y = cat(2, data{:});
    
    data = cellfun(@(f) f.jaw_xPos, F, UniformOutput=false);
    jaw_x = cat(2, data{:});
    
    data = cellfun(@(f) f.jaw_yPos, F, UniformOutput=false);
    jaw_y = cat(2, data{:});
    
    t = -4:1/30:0;
    
    ax = axes(figure(iExp));
    hold(ax, 'on')
    plot(t, mean(handL_x, 2, 'omitnan'), Color=getColor(1, 4, 0.7), DisplayName='handL_x')
    plot(t, mean(handL_y, 2, 'omitnan'), Color=getColor(1, 4, 0.7), DisplayName='handL_y')
    plot(t, mean(handR_x, 2, 'omitnan'), Color=getColor(2, 4, 0.7), DisplayName='handR_x')
    plot(t, mean(handR_y, 2, 'omitnan'), Color=getColor(2, 4, 0.7), DisplayName='handR_y')
    plot(t, mean(nose_x, 2, 'omitnan'), Color=getColor(3, 4, 0.7), DisplayName='nose_x')
    plot(t, mean(nose_y, 2, 'omitnan'), Color=getColor(3, 4, 0.7), DisplayName='nose_y')
    plot(t, mean(jaw_x, 2, 'omitnan'), Color=getColor(4, 4, 0.7), DisplayName='jaw_x')
    plot(t, mean(jaw_y, 2, 'omitnan'), Color=getColor(4, 4, 0.7), DisplayName='jaw_y')
    legend(ax, Location='northoutside', Orientation='horizontal', Interpreter='none')
end
% [clip, t] = exp(iExp).getVideoClip(touchTimes(1:10:end), side='l', bodyParts={'handIpsiCam', 'handContraCam', 'nose', 'jaw'});


%% Single trials
close all
[velocityKernels, ~, ~] = CompleteExperiment.makeConsineKernels(0, width=0.1, overlap=0.5, direction='both');

for iExp = 1:length(exp)
    touchTimes = [exp(iExp).eu(1).Trials.Press.Stop];
    F = cell(length(touchTimes), 1);
    for iTrial = 1:length(touchTimes)
        F{iTrial} = exp(iExp).getFeatures(timestamps=touchTimes(iTrial)-4:1/30:touchTimes(iTrial), features={'handL', 'handR'}, stats={'xPos', 'yPos'}, ...
            useGlobalNormalization=false);
        F{iTrial} = CompleteExperiment.convolveFeatures(F{iTrial}, velocityKernels, kernelNames={'_smooth'}, ...
            features={'handL', 'handR'}, ...
            stats={'xPos', 'yPos'}, ...
            mode='replace', normalize='zscore');
    end
    
    data = cellfun(@(f) f.handL_xPos_smooth, F, UniformOutput=false);
    handL_x = cat(2, data{:});
    
    data = cellfun(@(f) f.handL_yPos_smooth, F, UniformOutput=false);
    handL_y = cat(2, data{:});
    
    data = cellfun(@(f) f.handR_xPos_smooth, F, UniformOutput=false);
    handR_x = cat(2, data{:});
    
    data = cellfun(@(f) f.handR_yPos_smooth, F, UniformOutput=false);
    handR_y = cat(2, data{:});
    
    t = -4:1/30:0;
    
    ax = axes(figure(iExp));
    hold(ax, 'on')
    plot(handL_x, DisplayName='handL_x')
    plot(handL_y, DisplayName='handL_y')
    plot(handR_x, DisplayName='handR_x')
    plot(handR_y, DisplayName='handR_y')
%     legend(ax, Location='northoutside', Orientation='horizontal', Interpreter='none')
end

%% Single trials
close all
theta = 1;
% thetaPrct = 0.1;
[velocityKernels, ~, ~] = CompleteExperiment.makeConsineKernels(0, width=0.1, overlap=0.5, direction='both');

TST = cell(1, length(exp));
for iExp = 1:length(exp)
    disp(iExp)
    trials = exp(iExp).eu(1).getTrials('press');
    trueStartTime = NaN(1, length(trials));
    for iTrial = 1:length(trials)
        t = trials(iTrial).Stop-4:1/30:trials(iTrial).Stop;
        F = exp(iExp).getFeatures(timestamps=touchTimes(iTrial)-4:1/30:touchTimes(iTrial), features={'handL', 'handR'}, stats={'xPos', 'yPos'}, ...
            useGlobalNormalization=false);
        F = CompleteExperiment.convolveFeatures(F, velocityKernels, kernelNames={'_smooth'}, ...
            features={'handL', 'handR'}, ...
            stats={'xPos', 'yPos'}, ...
            mode='replace', normalize='zscore');
        F.inTrial = [];
        F.t = [];
        F = normalize(F);
        data = flip(table2array(F), 1);
        isAbove = abs(data) >= theta;
%         isAbove = abs(data) >= abs(data(1, :)) * thetaPrct;
        isAboveConseq = any(isAbove & [0, 0, 0, 0; diff(isAbove)] == 0, 2);
        t = flip(t);
        t = t - t(1);
        trueStartIndex = find(~isAboveConseq, 1, 'first') - 1;
        if ~isempty(trueStartIndex)
            trueStartTime(iTrial) = t(max(trueStartIndex, 1));
        end
%         close all
%         plot(t, data, 'r')
%         hold('on')
%         plot(t, isAboveConseq, 'g')
%         plot([trueStartTime(iTrial), trueStartTime(iTrial)], [-4, 4], 'k:')
    %     lclips = exp(i).getVideoClip(trueStartTime, side='l', numFramesBefore=30);
    %     rclips = exp(i).getVideoClip(trueStartTime, side='r', numFramesBefore=30);
    %     implay(lclips, 30)
    %     implay(rclips, 30)
        % plot(F.t - F.t(end), table2array(F(:, {'handL_xVel', 'handL_yVel', 'handR_xVel', 'handR_yVel'})))
    end
    TST{iExp} = trueStartTime;
end


figure(1)
for iExp = 1:length(exp)
    ax = subplot(2, 4, iExp);
    histogram(ax, TST{iExp})
    xlim(ax, [-4, 0])
end

ax = axes(figure(2));
histogram(ax, cat(2, TST{:}));
xlim(ax, [-4, 0])


%% Make PETH heatmaps with corrected movement initiation times
clear eta
ETA = cell(length(exp), 1);
for iExp = 1:length(exp)
    ETA{iExp} = exp(iExp).eu.getETA('count', 'press', window=[-4, 0], alignTo='stop', includeInvalid=true, correction=TST{iExp}, ...
        normalize=[-4, -2]);
    
%     EphysUnit.plotETA(ETA{iExp}, xlim=[-4,0], clim=[-2, 2], sortWindow=[-2, 0], signWindow=[-0.5, 0], sortThreshold=0.6, negativeSortThreshold=0.3);
end


eta.X = cellfun(@(e) e.X, ETA, UniformOutput=false);
eta.X = cat(1, eta.X{:});
eta.t = ETA{1}.t;
eta.N = repmat(size(eta.X, 1), size(eta.X, 1), 1);

EphysUnit.plotETA(eta, xlim=[-4,0], clim=[-1.5, 1.5], sortWindow=[-2, 0], signWindow=[-0.5, 0], sortThreshold=0.5, negativeSortThreshold=0.25);

%% Calculate InterPressIntervals
close all

ipi = cell(length(exp), 1);
for iExp = 1:length(exp)
    ipi{iExp} = diff([exp(iExp).eu(1).Trials.Press.Stop]);
end

ipiAggr = cat(2, ipi{:});

edges = 0:2:50;

fig = figure(Units='normalized', OuterPosition=[0, 0, 1, 0.67]);
for iExp = 1:length(exp)
    ax = subplot(2, 4, iExp);
    histogram(ax, ipi{iExp}, edges, Normalization='pdf')
end

fig = figure(Units='normalized', OuterPosition=[0, 0.67, 1, 0.33]);
ax = axes(fig);
histogram(ax, ipiAggr, edges, Normalization='pdf')