%% Spike detection
% pause(30*60);

clear, clc
folders = { ...
    % 'C:\SERVER\daisy29\daisy29_20251023', ...
    % 'C:\SERVER\daisy29\daisy29_20251024', ...
    % 'C:\SERVER\daisy29\daisy29_20251025', ...
    % 'C:\SERVER\daisy29\daisy29_20251027', ...
    % 'C:\SERVER\desmond41\desmond41_20251028', ...
    % 'C:\SERVER\daisy29\daisy29_20251028', ...
    % 'C:\SERVER\desmond41\desmond41_20251029', ...
    % 'C:\SERVER\daisy30\daisy30_20251029', ...
    % 'C:\SERVER\daisy30\daisy30_20251030', ...
    % 'C:\SERVER\desmond42\desmond42_20251030', ...
    % 'C:\SERVER\daisy29\daisy29_20251031', ...
    % 'C:\SERVER\desmond42\desmond42_20251031', ...
    % 'C:\SERVER\daisy29\daisy29_20251103', ...
    % 'C:\SERVER\desmond41\desmond41_20251104', ...
    % 'C:\SERVER\desmond41\desmond41_20251117', ...
    % 'C:\SERVER\daisy29\daisy29_20251118', ...
    % 'C:\SERVER\daisy30\daisy30_20251119', ...
    'C:\SERVER\daisy29\daisy29_20251120', ...
    'C:\SERVER\daisy30\daisy30_20251121', ...
    'C:\SERVER\desmond42\desmond42_20251121', ...
    'C:\SERVER\desmond41\desmond41_20251124', ...
    };

for iSession = 1:length(folders)
    try
        tr = TetrodeRecording;
        tr.SelectFiles(NeuropixelPath=folders{iSession});
        tr.ReadFiles(Duration=240, NumSigmas=4, NumSigmasReturn=1.5, NumSigmasReject=40, WaveformWindow=[-0.5, 1])
        tr.SaveNeuropixelIO()
        
        % Read detected spikes and NIDQ digital/analog channels
        % tr.LoadNeuropixelIO();
        tr.ParseNeuropixelIO(DigitalChannels={'Sync', 0; 'Lick', 1; 'Press', 2; 'Reward', 3; 'Timeout', 4; 'Mot2Busy', 5; 'CueLeft', 6; 'CueRight', 7});
        
        % Spike sort
        tr.LoadSpikes(1:384);

        % % Remove spikes after train 80, i.e. only keep 10ms pulses
        % [tce, stimOn, stimOff, ~] = tr.LoadTwoColorExperiment();
        % selPulses = abs((stimOff - stimOn) - 0.010) < 0.001;
        % iLastPulse = find(selPulses, 1, 'last');
        % for iChn = 1:384
        %     toDiscard = tr.Spikes(iChn).Timestamps > stimOff(iLastPulse) + 10;
        %     tr.Spikes(iChn).SampleIndex(toDiscard) = [];
        %     tr.Spikes(iChn).Timestamps(toDiscard) = [];
        %     tr.Spikes(iChn).Waveforms(toDiscard, :) = [];
        % end
        % clear tce stimOn stimOff selPulses iLastPulse iChn toDiscard

        tr.IterativeArtifactRemoval(1:384, MinSpikeRate=0.5, KIterative=4, KFinal=2, MaxIters=5, ...
            DimensionIterative=3, DimensionFinal=10, FeatureMethod='PCA', ClusterMethod='kmeans', ...
            WaveformWindow=[-0.5, 0.5]);
        tr.SaveSpikes(Path='Spikes_AutoSortedIterative');
    catch ME
        warning('Could not process folder: %s', folders{iSession})
        warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
    end
end

% %% Remove spikes after train 80 (i.e. only keep 10ms pulses) and redo spike sorting, this helps remove drift
% folders = { ...
%     'C:\SERVER\daisy29\daisy29_20251023', ...
%     'C:\SERVER\daisy29\daisy29_20251024', ...
%     'C:\SERVER\daisy29\daisy29_20251027', ...
%     };
% 
% chunkSize = 32; % NumChannelsPerChunk
% for iSession = 1:length(folders)
%     try
%         tr = TetrodeRecording();
%         tr.SelectFiles(NeuropixelPath=folders{iSession});
%         tr.LoadNeuropixelIO();
%         tr.ParseNeuropixelIO(DigitalChannels={'Sync', 0; 'Lick', 1; 'Press', 2; 'Reward', 3; 'Timeout', 4; 'Mot2Busy', 5; 'CueLeft', 6; 'CueRight', 7});
% 
%         % IterativeArtifactRemoval (OnionPeeling): Load spikes
%         for iChunk = 1:(384/chunkSize)
%             channels = (iChunk-1)*chunkSize + 1 : iChunk*chunkSize;
%             tr.LoadSpikes(channels, Path='Spikes');
%             channels = [tr.Spikes.Channel];
% 
%             if isempty(channels)
%                 tr.Spikes = [];
%                 continue
%             end
% 
%             % Remove spikes after train 80, i.e. only keep 10ms pulses
%             [tce, stimOn, stimOff, ~] = tr.LoadTwoColorExperiment();
%             selPulses = abs((stimOff - stimOn) - 0.010) < 0.001;
%             iLastPulse = find(selPulses, 1, 'last');
%             for iChn = channels(:)'
%                 toDiscard = tr.Spikes(iChn).Timestamps > stimOff(iLastPulse) + 10;
%                 tr.Spikes(iChn).SampleIndex(toDiscard) = [];
%                 tr.Spikes(iChn).Timestamps(toDiscard) = [];
%                 tr.Spikes(iChn).Waveforms(toDiscard, :) = [];
%             end
%             clear tce stimOn stimOff selPulses iLastPulse iChn toDiscard
% 
%             % IterativeArtifactRemoval (OnionPeeling): kmeans, pca
%             tr.IterativeArtifactRemoval(channels, MinSpikeRate=0.5, KIterative=4, KFinal=2, MaxIters=5, ...
%                 DimensionIterative=3, DimensionFinal=10, FeatureMethod='PCA', ClusterMethod='kmeans', ...
%                 WaveformWindow=[-0.5, 0.5]);
% 
%             tr.SaveSpikes(Channels=channels, Path='Spikes_AutoSortedIterative')
%             tr.Spikes = [];
%         end
%     catch ME
%         warning('Could not process folder: %s', folders{iSession})
%         warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
%     end
% end


%% Load sorted data on a different PC

clear, clc
tr = TetrodeRecording();
tr.SelectFiles(NeuropixelPath='C:\SERVER\desmond42\desmond42_20251121')
tr.LoadNeuropixelIO();
tr.ParseNeuropixelIO(DigitalChannels={'Sync', 0; 'Lick', 1; 'Press', 2; 'Reward', 3; 'Timeout', 4; 'Mot2Busy', 5; 'CueLeft', 6; 'CueRight', 7});

channels = 1:128;
tr.LoadSpikes(channels, Path='Spikes_AutoSortedIterative');
tr.PlotAllChannels(Channels=channels, plotMethod='mean')

%%
tr.SaveSpikes(Channels=channels, Path='Spikes_Sorted')
channels = 129:256;
tr.Spikes = [];
tr.LoadSpikes(channels, Path='Spikes_AutoSortedIterative');
tr.PlotAllChannels(Channels=channels, plotMethod='mean')
%%
tr.SaveSpikes(Channels=channels, Path='Spikes_Sorted')
channels = 257:384;
tr.Spikes = [];
tr.LoadSpikes(channels, Path='Spikes_AutoSortedIterative');
tr.PlotAllChannels(Channels=channels, plotMethod='mean')
%%
tr.SaveSpikes(Channels=channels, Path='Spikes_Sorted')

%% Convert to EphysUnits
folders = { ...
    ... 'C:\SERVER\daisy29\daisy29_20251023', ... SORTED, EU
    ... 'C:\SERVER\desmond41\desmond41_20251028', ... SORTED, EU
    ... 'C:\SERVER\desmond42\desmond42_20251030', ... SORTED, EU
    ... 'C:\SERVER\daisy29\daisy29_20251031', ... SORTED, EU
    ... 'C:\SERVER\desmond41\desmond41_20251104', ... SORTED, EU
    ... 'C:\SERVER\desmond41\desmond41_20251117', ... SORTED, EU
    ... 'C:\SERVER\daisy29\daisy29_20251118', ... SORTED, EU
    ... 'C:\SERVER\daisy30\daisy30_20251119', ... SORTED, EU
    ... 'C:\SERVER\daisy29\daisy29_20251024', ...
    ... 'C:\SERVER\daisy29\daisy29_20251025', ...
    ... 'C:\SERVER\daisy29\daisy29_20251027', ...
    ... 'C:\SERVER\daisy29\daisy29_20251028', ...
    ... 'C:\SERVER\desmond41\desmond41_20251029', ...
    ... 'C:\SERVER\daisy30\daisy30_20251029', ...
    ... 'C:\SERVER\daisy30\daisy30_20251030', ...
    ... 'C:\SERVER\desmond42\desmond42_20251031', ...
    ... 'C:\SERVER\daisy29\daisy29_20251103', ...
    ... 'C:\SERVER\daisy29\daisy29_20251120', ...
    ... 'C:\SERVER\daisy30\daisy30_20251121', ...
    ... 'C:\SERVER\desmond42\desmond42_20251121', ...
    ... 'C:\SERVER\desmond41\desmond41_20251124', ...
    };

chunkSize = 32; % NumChannelsPerChunk
for iSession = 1:length(folders)
    try
        tr = TetrodeRecording();
        tr.SelectFiles(NeuropixelPath=folders{iSession});
        tr.LoadNeuropixelIO();
        tr.ParseNeuropixelIO(DigitalChannels={'Sync', 0; 'Lick', 1; 'Press', 2; 'Reward', 3; 'Timeout', 4; 'Mot2Busy', 5; 'CueLeft', 6; 'CueRight', 7});

        % IterativeArtifactRemoval (OnionPeeling): Load spikes
        for iChunk = 1:(384/chunkSize)
            channels = (iChunk-1)*chunkSize + 1 : iChunk*chunkSize;
            tr.LoadSpikes(channels, Path='Spikes_Sorted');

            if isempty(tr.Spikes) || isempty([tr.Spikes.Channel])
                tr.Spikes = [];
                continue
            end

            channels = [tr.Spikes.Channel];

            ar = AcuteRecording(tr, 'N/A');
            ar.binMoveResponse(tr, 'none', Window=[-1, 0], Store=true);
            eu = EphysUnit(ar, readWaveforms=false, cullITI=false, savepath='C:\SERVER\Units\TwoColor_Striatonigral', tr=tr);

            tr.Spikes = [];
            clear eu
        end
    catch ME
        warning('Could not process folder: %s', folders{iSession})
        warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
    end
end

%%
eu = EphysUnit.load('C:\SERVER\Units\TwoColor_Striatonigral', waveforms=false, spikecounts=false, spikerates=false);

%% Remove multiunits, fast (ISS test)
eu = eu.removeMultiUnits(cullZeros=true);

% Remove drift, low spike rate units, fast
clear c
% Remove drift
c.isDrifting = detectDriftingUnits(eu, smoothWindow=300, tolerance=0.05, spikeRateThreshold=5, includeITI=true);

% Filter by spike rate
msr = arrayfun(@(eu) eu.SpikeRateStats.median, eu);
p.minSpikeRate = 15;
c.isSNr = msr >= p.minSpikeRate;

eu = eu(c.isSNr & ~c.isDrifting);

% Remove duplicates (slow, pairwise comparisons)
[eu, isDuplicate] = eu.removeDuplicates(0.7);

eu.save('C:\SERVER\Units\TwoColor_Striatonigral\SingleUnit_NonDuplicate_NonDrift_SNr')

%% ETA Stim (by ISI)
p.isiBaselineWindow = [-0.2, 0];
p.stimBluePowers = [25]*1e-6; 
p.stimRedPowers = [25]*1e-6;
p.stimBlueDurations = [10]*1e-3;
p.stimRedDurations = [10]*1e-3;

p.isiWindow = [-0.4, 0.4];
p.isiRes = 1e-3;
p.xlim.stim = [-0.1, 0.3];
p.xlim.move = [-4, 2];
p.path = 'C:\SERVER\Figures\TwoColor_Striatonigral\TestBatch2';
p.rasterSzStim = 1;
p.rasterSzMove = 1;

close all
XBlue = cell(length(eu), 1);
XRed = cell(length(eu), 1);
for iEu = 1:length(eu)
    groupsBlue = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'duration'}, selectBy=struct(power=p.stimBluePowers, duration=p.stimBlueDurations, location=[], wavelength=[470, 473]));
    groupsRed = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'duration'}, selectBy=struct(power=p.stimRedPowers, duration=p.stimRedDurations, location=[], wavelength=[593, 635]));

    % rd = eu(iEu).getRasterData('stimtwocolor', p.isiWindow, trials=[groupsRed.trials], alignTo='start', shutterDelay=0, sort=false, photoelectricBlankDuration=0.5e-3);
    % EphysUnit.plotRaster(rd)

    [isi, t] = eu(iEu).getMeanPEISI('stimtwocolor', [groupsBlue.trials], window=p.isiWindow, resolution=p.isiRes, photoelectricBlankDuration=0);
    selBaseline = t>p.isiBaselineWindow(1) & t<p.isiBaselineWindow(2);
    normSR = (1./isi - mean(1./isi(:, selBaseline), 'omitnan')) ./ std(1./isi(:, selBaseline), 0, 2, 'omitnan');
    XBlue{iEu} = normSR;

    [isi, t] = eu(iEu).getMeanPEISI('stimtwocolor', [groupsRed.trials], window=p.isiWindow, resolution=p.isiRes, ...
        photoelectricBlankDuration=1.5e-3);
    selBaseline = t>p.isiBaselineWindow(1) & t<p.isiBaselineWindow(2);
    normSR = (1./isi - mean(1./isi(:, selBaseline), 'omitnan')) ./ std(1./isi(:, selBaseline), 0, 2, 'omitnan');
    XRed{iEu} = normSR;

    % fprintf('nan=%i, nan=%i\n', nnz(isnan(XBlue{iEu})), nnz(isnan(XRed{iEu})))
end

eta.stimBlueISI = struct(X=cat(1, XBlue{:}), t=t, N=[], D=[], stats=[]);
eta.stimRedISI = struct(X=cat(1, XRed{:}), t=t, N=[], D=[], stats=[]);

clear XBlue XRed iEu groupsBlue groupsRed isi t selBaseline normSR

%% Calculate META
p.metaWindowStim = [0.005, 0.019];
p.posRespThresholdStim = 3;
p.negRespThresholdStim = -2;

t = eta.stimBlue.t;
meta.stimBlueISI = mean(eta.stimBlueISI.X(:, t>=p.metaWindowStim(1) & t<=p.metaWindowStim(2)), 2, 'omitnan');
t = eta.stimRed.t;
meta.stimRedISI = mean(eta.stimRedISI.X(:, t>=p.metaWindowStim(1) & t<=p.metaWindowStim(2)), 2, 'omitnan');
clear t
% 
% c.isPressUp =         meta.press >= p.posRespThreshold;
% c.isPressDown =       meta.press <= p.negRespThreshold;
% c.isPressResponsive = c.isPressUp | c.isPressDown;
% c.isLickUp =          meta.lick >= p.posRespThreshold;
% c.isLickDown =        meta.lick <= p.negRespThreshold;
% c.isLickResponsive =  c.isLickUp | c.isLickDown;
% c.isLickUnresponsiveButUp = ~c.isLickResponsive & meta.lick > 0;
% c.isLickUnresponsiveButDown = ~c.isLickResponsive & meta.lick < 0;
% c.isPressUnresponsiveButUp = ~c.isPressResponsive & meta.press > 0;
% c.isPressUnresponsiveButDown = ~c.isPressResponsive & meta.press < 0;

c.isStimBlueUp = meta.stimBlueISI >= p.posRespThresholdStim;
c.isStimBlueDown = meta.stimBlueISI <= p.negRespThresholdStim;
c.isStimRedUp = meta.stimRedISI >= p.posRespThresholdStim;
c.isStimRedDown = meta.stimRedISI <= p.negRespThresholdStim;

c.isStimBlueUpRedUpThereforeChrimsonMaybe = c.isStimBlueUp & c.isStimRedUp;
c.isStimBlueUpRedNotUpThereforeCoChrMaybe = c.isStimBlueUp & ~c.isStimRedUp;
c.isStimBlueNotUpRedUpThereforeChrimsonMaybe = ~c.isStimBlueUp & c.isStimRedUp;
c.isStimBlueNotUpRedNotUp = ~c.isStimBlueUp & ~c.isStimRedUp;

fprintf('ChrimsonR %i, CoChR %i\n', nnz(c.isStimRedUp), nnz(c.isStimBlueUpRedNotUpThereforeCoChrMaybe))

%% Plot opto heatmap (ISI)
theta = 0;

selUnitsStim = true(size(eu));c.isStimBlueUp | c.isStimRedUp | c.isStimBlueDown | c.isStimRedDown;
% close all
fig = figure(Units='inches', Position=[1, 1, 4.5, 5]);
ETASORT = {eta.stimBlueISI, eta.stimRedISI};
SORTWINDOW = {[0, 100]*1e-3, [0, 100]*1e-3};
ETA = {eta.stimBlueISI, eta.stimRedISI};
NAME = [sprintf("470nm\n%suW, %sms", string(p.stimBluePowers*1e6).join(', '), string(p.stimBlueDurations*1e3).join(', ')), ...
    sprintf("593nm\n%suW, %sms", string(p.stimRedPowers*1e6).join(', '), string(p.stimRedDurations*1e3).join(', '))];
NAMECOLOR = ["blue", "red"];
ZEROLABEL = ["stim", "stim"];
XLIM = {[-200, 200]*1e-3, [-200, 200]*1e-3};
XTICKS = {[0, 10, 50, 100]*1e-3, [0, 10, 50, 100]*1e-3};
CLIM = {[-5, 5], [-5, 5]};
CBTILE = {'', 'east'};
w = cellfun(@(xl) round(10*diff(xl)), XLIM);
cw = [0, cumsum(w)];


% Combine ETA, PCA, and sort along 1st dimension
etaCombined = struct(X=[], t=[]);
etaCombined.X = cellfun(@(eta) eta.X, ETASORT, UniformOutput=false);
etaCombined.X = cat(2, etaCombined.X{:});
etaCombined.t = cellfun(@(eta) eta.t, ETASORT, UniformOutput=false);
etaCombined.t = cat(2, etaCombined.t{:});
etaCombined.epoch = arrayfun(@(i) i*ones(1, length(ETASORT{i}.t)), 1:length(ETASORT), UniformOutput=false);
etaCombined.epoch = cat(2, etaCombined.epoch{:});
etaCombined.X(etaCombined.X>3) = 3;
etaCombined.X(etaCombined.X<-1.5) = -1.5;
etaCombined.X = etaCombined.X(selUnitsStim, :);

% For sorting, make templates to dot-product with
clear template
template(length(ETASORT)) = struct(t=[], x=[]);
for iETA = 1:length(ETASORT)
    template(iETA).t = etaCombined.t;
    template(iETA).x = zeros(1, length(etaCombined.t));
    template(iETA).x(1, isin(etaCombined.t, SORTWINDOW{iETA}) & etaCombined.epoch==iETA) = 1;
end

score = zeros(size(etaCombined.X, 1), length(ETASORT));
etaCombined.X(isnan(etaCombined.X)) = 0;
for iETA = 1:length(ETASORT)
    score(:, iETA) = etaCombined.X * template(iETA).x';
end
groupVar = arrayfun(@(i) bitshift(int16(score(:, i)>theta), length(ETASORT)-i), 1:size(score, 2), UniformOutput=false);
groupVar = sum(horzcat(groupVar{:}), 2);

% First, sort by number of negative modulations
numNeg = sum(score>theta, 2);
[uniqueGroupVars, ia] = unique(groupVar);
[~, I] = sort(numNeg(ia), 'ascend');
groupVar = changem(groupVar, 0:length(uniqueGroupVars)-1, uniqueGroupVars(I));

% Then, put all small groups (excluding single neg ones) at the bottom
[uniqueGroupVars, ia] = unique(groupVar);
assert(length(uniqueGroupVars) == max(groupVar)+1);
groupSize = histcounts(groupVar, 0:length(uniqueGroupVars));

numUnitsInSameGroup = arrayfun(@(gv) nnz(groupVar==gv), groupVar);
isRare = numUnitsInSameGroup < 3;
isSingleNeg = numNeg==1;
groupVar(isRare & ~isSingleNeg) = max(groupVar)+1;
% Tighten up the groupvars
uniqueGroupVars = unique(groupVar);
groupVar = changem(groupVar, 0:length(uniqueGroupVars)-1, uniqueGroupVars);
groupSize = histcounts(groupVar, 0:length(uniqueGroupVars));
groupSizeCum = cumsum(groupSize);


% groupVar = zeros(length(eu), 1);
% groupVar(c.isStimRedUp) = groupVar(c.isStimRedUp) + 1;
% groupVar(c.isStimBlueUp) = groupVar(c.isStimBlueUp) + 2;
% groupVar = groupVar(selUnitsStim);
% uniqueGroupVars = unique(groupVar);
% groupVar = changem(groupVar, 0:length(uniqueGroupVars)-1, uniqueGroupVars);
% groupSize = histcounts(groupVar, 0:length(uniqueGroupVars));
% groupSizeCum = cumsum(groupSize);
[~, sortOrder] = sort(double(groupVar)*10 + score(:, 1)./max(abs(score(:, 1))), 'ascend');


tl = tiledlayout(fig, 1, sum(w), TileSpacing='tight', Padding='tight');
ax = gobjects(1, length(ETA));
for iAx = 1:length(ETA)
    hidecb = isempty(CBTILE{iAx});
    ax(iAx) = nexttile(tl, 1 + cw(iAx), [1, w(iAx)]);
    EphysUnit.plotETA(ax(iAx), ETA{iAx}, selUnitsStim, xlim=XLIM{iAx}, clim=[-1.5, 1.5], order=sortOrder, hidecolorbar=hidecb);
    % applyCustomColormap(ax(iAx), [-1.5, 1.5], hlim=[0.375, 0, 0, -0.375], llim=[0.125, 0.5, 0.5, 0.25], hpwr=.5, lpwr=1, h0=0.33);
    applyCustomColormap(ax(iAx), CLIM{iAx}, hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);
    if ~hidecb
        ax(iAx).Colorbar.Layout.Tile = CBTILE{iAx};
    end
    if iAx > 1
        yticks(ax(iAx), [])
    else
        yticks(ax(iAx), groupSizeCum(1:end)+0.5)
        yticklabels(ax(iAx), string(groupSizeCum(1:end)))
    end
    title(ax(iAx), strsplit(NAME(iAx), "\\n"), Color=NAMECOLOR(iAx))
    xlabel(ax(iAx), "")
    ylabel(ax(iAx), "")
    xticks(ax(iAx), XTICKS{iAx})
    xticklabels(ax(iAx), string(1e3*XTICKS{iAx}))
    % xticks(ax(iAx), [0, 0.1])
    % xticklabels(ax(iAx), [string(XLIM{iAx}(1)), ZEROLABEL(iAx), string(XLIM{iAx}(2))])
    xtickangle(ax(iAx), 0)
    xline(ax(iAx), 0, 'k-')
    if iAx <= 2
        xline(ax(iAx), 0.01, 'k-')
    end
    yline(ax(iAx), groupSizeCum(1:end-1)+0.5, 'k:', LineWidth=1.5)
    ax(iAx).YAxis.TickLength = [0, 0];
end
ax(2).Colorbar.Label.String = 'opto response (a.u.)';

xlabel(tl, "time (ms)")
ylabel(ax(1), "unit")
fontsize(fig, 14, 'points')

copygraphics(fig, ContentType='vector', BackgroundColor='none')
    
%% ETA Stim (binned spike counts)
clear sortOrder groupSizeCum
for power = [2000, 500, 100, 50, 25]    
    p.etaBaselineWindow = [-0.5, -0.1];
    p.stimBluePowers = [power]*1e-6; 
    p.stimRedPowers = [power]*1e-6;
    p.stimBlueDurations = [10]*1e-3;
    p.stimRedDurations = [10]*1e-3;
    
    p.isiWindow = [-0.5, 0.5];
    p.etaRes = 0.025;
    p.xlim.stim = [-0.1, 0.3];
    p.xlim.move = [-4, 2];
    p.path = 'C:\SERVER\Figures\TwoColor_Striatonigral\TestBatch2';
    p.rasterSzStim = 1;
    p.rasterSzMove = 1;
    
    XBlue = cell(length(eu), 1);
    XRed = cell(length(eu), 1);
    linelength = 0;
    for iEu = 1:length(eu)
        groupsBlue = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'duration'}, selectBy=struct(power=p.stimBluePowers, duration=p.stimBlueDurations, location=[], wavelength=[470, 473]));
        groupsRed = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'duration'}, selectBy=struct(power=p.stimRedPowers, duration=p.stimRedDurations, location=[], wavelength=[593, 635]));
    
        % rd = eu(iEu).getRasterData('stimtwocolor', p.isiWindow, trials=[groupsRed.trials], alignTo='start', shutterDelay=0, sort=false, photoelectricBlankDuration=0.5e-3);
        % EphysUnit.plotRaster(rd)
    
        % etaTemp = eu(iEu).getETA('count', 'stim', window=p.isiWindow, resolution=p.etaRes, trials=[groupsBlue.trials], alignTo='start', normalize=p.etaBaselineWindow);
        etaTemp = eu(iEu).getETA('count', 'stim', window=p.isiWindow, resolution=p.etaRes, trials=[groupsBlue.trials], alignTo='start', normalize='none');...p.etaBaselineWindow);
        XBlue{iEu} = (etaTemp.X - mean(etaTemp.X(isin(etaTemp.t, p.etaBaselineWindow))))./p.etaRes;
    
        % etaTemp = eu(iEu).getETA('count', 'stim', window=p.isiWindow, resolution=p.etaRes, trials=[groupsRed.trials], alignTo='start', normalize=p.etaBaselineWindow);
        etaTemp = eu(iEu).getETA('count', 'stim', window=p.isiWindow, resolution=p.etaRes, trials=[groupsRed.trials], alignTo='start', normalize='none');...p.etaBaselineWindow);
        XRed{iEu} = (etaTemp.X - mean(etaTemp.X(isin(etaTemp.t, p.etaBaselineWindow))))./p.etaRes;
        fprintf(repmat('\b', 1, linelength))
        linelength = fprintf('%i of %i...\n', iEu, length(eu));
        % plot(ax, etaTemp.t, [XRed{iEu}', XBlue{iEu}']./p.etaRes)
        % ylim(ax, [-20, 20])
        % cla(ax)
    end
    
    eta.stimBlue = struct(X=cat(1, XBlue{:}), t=etaTemp.t, N=[], D=[], stats=[]);
    eta.stimRed = struct(X=cat(1, XRed{:}), t=etaTemp.t, N=[], D=[], stats=[]);
    
    clear XBlue XRed iEu groupsBlue groupsRed isi t selBaseline normSR etaTemp linelength
    
    % Plot opto heatmap (ETA)
    theta = 0;
    
    selUnitsStim = true(size(eu));c.isStimBlueUp | c.isStimRedUp | c.isStimBlueDown | c.isStimRedDown;
    fig = figure(Units='inches', Position=[1, 1, 4.5, 5]);
    ETASORT = {eta.stimBlue, eta.stimRed};
    SORTWINDOW = {[0, 100]*1e-3, [0, 100]*1e-3};
    ETA = {eta.stimBlue, eta.stimRed};
    NAME = [sprintf("470nm\n%suW, %sms", string(p.stimBluePowers*1e6).join(', '), string(p.stimBlueDurations*1e3).join(', ')), ...
        sprintf("593nm\n%suW, %sms", string(p.stimRedPowers*1e6).join(', '), string(p.stimRedDurations*1e3).join(', '))];
    NAMECOLOR = ["blue", "red"];
    ZEROLABEL = ["stim", "stim"];
    XLIM = {[-500, 200]*1e-3, [-500, 200]*1e-3};
    XTICKS = {[0, 10, 50, 100]*1e-3, [0, 10, 50, 100]*1e-3};
    CLIM = {[-20, 20], [-20, 20]};
    CBTILE = {'', 'east'};
    w = cellfun(@(xl) round(10*diff(xl)), XLIM);
    cw = [0, cumsum(w)];
    
    
    % Combine ETA, PCA, and sort along 1st dimension
    if ~exist('sortOrder', 'var')
        etaCombined = struct(X=[], t=[]);
        etaCombined.X = cellfun(@(eta) eta.X, ETASORT, UniformOutput=false);
        etaCombined.X = cat(2, etaCombined.X{:});
        etaCombined.t = cellfun(@(eta) eta.t, ETASORT, UniformOutput=false);
        etaCombined.t = cat(2, etaCombined.t{:});
        etaCombined.epoch = arrayfun(@(i) i*ones(1, length(ETASORT{i}.t)), 1:length(ETASORT), UniformOutput=false);
        etaCombined.epoch = cat(2, etaCombined.epoch{:});
        etaCombined.X = etaCombined.X(selUnitsStim, :);
        
        % For sorting, make templates to dot-product with
        clear template
        template(length(ETASORT)) = struct(t=[], x=[]);
        for iETA = 1:length(ETASORT)
            template(iETA).t = etaCombined.t;
            template(iETA).x = zeros(1, length(etaCombined.t));
            template(iETA).x(1, isin(etaCombined.t, SORTWINDOW{iETA}) & etaCombined.epoch==iETA) = 1;
        end
        
        score = zeros(size(etaCombined.X, 1), length(ETASORT));
        etaCombined.X(isnan(etaCombined.X)) = 0;
        for iETA = 1:length(ETASORT)
            score(:, iETA) = etaCombined.X * template(iETA).x';
        end
        groupVar = arrayfun(@(i) bitshift(int16(score(:, i)>theta), length(ETASORT)-i), 1:size(score, 2), UniformOutput=false);
        groupVar = sum(horzcat(groupVar{:}), 2);
        
        % First, sort by number of negative modulations
        numNeg = sum(score>theta, 2);
        [uniqueGroupVars, ia] = unique(groupVar);
        [~, I] = sort(numNeg(ia), 'ascend');
        groupVar = changem(groupVar, 0:length(uniqueGroupVars)-1, uniqueGroupVars(I));
        
        % Then, put all small groups (excluding single neg ones) at the bottom
        [uniqueGroupVars, ia] = unique(groupVar);
        assert(length(uniqueGroupVars) == max(groupVar)+1);
        groupSize = histcounts(groupVar, 0:length(uniqueGroupVars));
        
        numUnitsInSameGroup = arrayfun(@(gv) nnz(groupVar==gv), groupVar);
        isRare = numUnitsInSameGroup < 3;
        isSingleNeg = numNeg==1;
        groupVar(isRare & ~isSingleNeg) = max(groupVar)+1;
        % Tighten up the groupvars
        uniqueGroupVars = unique(groupVar);
        groupVar = changem(groupVar, 0:length(uniqueGroupVars)-1, uniqueGroupVars);
        groupSize = histcounts(groupVar, 0:length(uniqueGroupVars));
        groupSizeCum = cumsum(groupSize);
        
        
        % groupVar = zeros(length(eu), 1);
        % groupVar(c.isStimRedUp) = groupVar(c.isStimRedUp) + 1;
        % groupVar(c.isStimBlueUp) = groupVar(c.isStimBlueUp) + 2;
        % groupVar = groupVar(selUnitsStim);
        % uniqueGroupVars = unique(groupVar);
        % groupVar = changem(groupVar, 0:length(uniqueGroupVars)-1, uniqueGroupVars);
        % groupSize = histcounts(groupVar, 0:length(uniqueGroupVars));
        % groupSizeCum = cumsum(groupSize);
        [~, sortOrder] = sort(double(groupVar)*10 + score(:, 1)./max(abs(score(:, 1))), 'ascend');
    end
    
    
    tl = tiledlayout(fig, 1, sum(w), TileSpacing='compact', Padding='tight');
    ax = gobjects(1, length(ETA));
    for iAx = 1:length(ETA)
        hidecb = isempty(CBTILE{iAx});
        ax(iAx) = nexttile(tl, 1 + cw(iAx), [1, w(iAx)]);
        EphysUnit.plotETA(ax(iAx), ETA{iAx}, selUnitsStim, xlim=XLIM{iAx}, clim=[-1.5, 1.5], order=sortOrder, hidecolorbar=hidecb);
        % applyCustomColormap(ax(iAx), [-1.5, 1.5], hlim=[0.375, 0, 0, -0.375], llim=[0.125, 0.5, 0.5, 0.25], hpwr=.5, lpwr=1, h0=0.33);
        applyCustomColormap(ax(iAx), CLIM{iAx}, hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);
        if ~hidecb
            ax(iAx).Colorbar.Layout.Tile = CBTILE{iAx};
        end
        if iAx > 1
            yticks(ax(iAx), [])
        else
            yticks(ax(iAx), groupSizeCum(1:end)+0.5)
            yticklabels(ax(iAx), string(groupSizeCum(1:end)))
        end
        title(ax(iAx), strsplit(NAME(iAx), "\\n"), Color=NAMECOLOR(iAx))
        xlabel(ax(iAx), "")
        ylabel(ax(iAx), "")
        xticks(ax(iAx), XTICKS{iAx})
        xticklabels(ax(iAx), string(1e3*XTICKS{iAx}))
        % xticks(ax(iAx), [0, 0.1])
        % xticklabels(ax(iAx), [string(XLIM{iAx}(1)), ZEROLABEL(iAx), string(XLIM{iAx}(2))])
        xtickangle(ax(iAx), 0)
        xline(ax(iAx), 0, 'k-')
        if iAx <= 2
            xline(ax(iAx), 0.01, 'k-')
        end
        yline(ax(iAx), groupSizeCum(1:end-1)+0.5, 'k:', LineWidth=1.5)
        ax(iAx).YAxis.TickLength = [0, 0];
    end
    ax(2).Colorbar.Label.String = 'opto response (a.u.)';
    
    xlabel(tl, "time (ms)")
    ylabel(ax(1), "unit")
    fontsize(fig, 14, 'points')
    
    % copygraphics(fig, ContentType='vector', BackgroundColor='none')
end


%% ETA Stim (binned spike counts)
clear sortOrder groupSizeCum
powers = [2000, 500, 100, 50, 25]*1e-6;
locations = [-2400, -2600, -2800, -3000, -600, -400, -200, 0];
durations = [10]*1e-3;

expNames = string({eu.ExpName});
[uniqueExpNames, ia, ic] = unique(expNames);



for iPower = 1:length(POWERS)
    power = POWERS(iPower);
    for iLocation = 1:length(LOCATIONS)
        location = LOCATIONS(iLocation);
        p.etaBaselineWindow = [-0.5, -0.1];
        p.stimBluePowers = [power]*1e-6; 
        p.stimRedPowers = [power]*1e-6;
        p.stimBlueDurations = [10]*1e-3;
        p.stimRedDurations = [10]*1e-3;
        
        p.isiWindow = [-0.5, 0.5];
        p.etaRes = 0.025;
        p.xlim.stim = [-0.1, 0.3];
        p.xlim.move = [-4, 2];
        p.path = 'C:\SERVER\Figures\TwoColor_Striatonigral\TestBatch2';
        p.rasterSzStim = 1;
        p.rasterSzMove = 1;
        
        XBlue = cell(length(eu), 1);
        XRed = cell(length(eu), 1);
        linelength = 0;
        for iEu = 1:length(eu)
            groupsBlue = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'duration'}, selectBy=struct(power=p.stimBluePowers, duration=p.stimBlueDurations, location=location, wavelength=[470, 473]));
            groupsRed = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'duration'}, selectBy=struct(power=p.stimRedPowers, duration=p.stimRedDurations, location=location, wavelength=[593, 635]));
        
            % rd = eu(iEu).getRasterData('stimtwocolor', p.isiWindow, trials=[groupsRed.trials], alignTo='start', shutterDelay=0, sort=false, photoelectricBlankDuration=0.5e-3);
            % EphysUnit.plotRaster(rd)
        
            % etaTemp = eu(iEu).getETA('count', 'stim', window=p.isiWindow, resolution=p.etaRes, trials=[groupsBlue.trials], alignTo='start', normalize=p.etaBaselineWindow);
            etaTemp = eu(iEu).getETA('count', 'stim', window=p.isiWindow, resolution=p.etaRes, trials=[groupsBlue.trials], alignTo='start', normalize='none');...p.etaBaselineWindow);
            XBlue{iEu} = (etaTemp.X - mean(etaTemp.X(isin(etaTemp.t, p.etaBaselineWindow))))./p.etaRes;
        
            % etaTemp = eu(iEu).getETA('count', 'stim', window=p.isiWindow, resolution=p.etaRes, trials=[groupsRed.trials], alignTo='start', normalize=p.etaBaselineWindow);
            etaTemp = eu(iEu).getETA('count', 'stim', window=p.isiWindow, resolution=p.etaRes, trials=[groupsRed.trials], alignTo='start', normalize='none');...p.etaBaselineWindow);
            XRed{iEu} = (etaTemp.X - mean(etaTemp.X(isin(etaTemp.t, p.etaBaselineWindow))))./p.etaRes;
            fprintf(repmat('\b', 1, linelength))
            linelength = fprintf('%i of %i...\n', iEu, length(eu));
            % plot(ax, etaTemp.t, [XRed{iEu}', XBlue{iEu}']./p.etaRes)
            % ylim(ax, [-20, 20])
            % cla(ax)
        end
        
        eta.stimBlue = struct(X=cat(1, XBlue{:}), t=etaTemp.t, N=[], D=[], stats=[]);
        eta.stimRed = struct(X=cat(1, XRed{:}), t=etaTemp.t, N=[], D=[], stats=[]);
    end
end