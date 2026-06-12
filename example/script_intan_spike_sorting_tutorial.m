
tr = TetrodeRecording;


tr.SelectFiles();
tr.ReadFiles(Duration=240, NumSigmas=4, NumSigmasReturn=1.5, NumSigmasReject=40, WaveformWindow=[-0.5, 1])

%%


%%
for iSession = 1:length(folders)
    try

        % Read detected spikes and NIDQ digital/analog channels
        % tr.LoadNeuropixelIO();

        % Spike sort
        tr.LoadSpikes(1:384);

        tr.IterativeArtifactRemoval(1:384, MinSpikeRate=0.5, KIterative=4, KFinal=2, MaxIters=5, ...
            DimensionIterative=3, DimensionFinal=10, FeatureMethod='PCA', ClusterMethod='kmeans', ...
            WaveformWindow=[-0.5, 0.5]);
        tr.SaveSpikes(Path='Spikes_AutoSortedIterative');
    catch ME
        warning('Could not process folder: %s', folders{iSession})
        warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
    end
end


%% Load sorted data on a different PC

clear, clc
tr = TetrodeRecording();
tr.SelectFiles(NeuropixelPath='C:\SERVER\daisy34\daisy34_20260320')
tr.LoadNeuropixelIO();
tr.ParseNeuropixelIO(DigitalChannels={'Sync', 0; 'Lick', 1; 'Press', 2; 'Reward', 3; 'Timeout', 4; 'Mot2Busy', 5; 'CueLeft', 6; 'CueRight', 7});

channels = 1:128;
tr.LoadSpikes(channels, Path='Spikes_AutoSortedIterative');
% tr.LoadSpikes(channels, Path='Spikes_AutoSortedIterative', ExpName='daisy29_20251120');
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
    ... 'C:\SERVER\daisy33\daisy33_20260226', ...
    ... 'C:\SERVER\daisy33\daisy33_20260227', ...
    ... 'C:\SERVER\daisy33\daisy33_20260302', ...
    ... 'C:\SERVER\daisy33\daisy33_20260303', ... 16 channels only
    ... 'C:\SERVER\daisy34\daisy34_20260303', ...
    ... 'C:\SERVER\daisy33\daisy33_20260304', ...
    ... 'C:\SERVER\daisy34\daisy34_20260304', ...
    ... 'C:\SERVER\daisy34\daisy34_20260305', ...
    ... 'C:\SERVER\daisy34\daisy34_20260306', ...
    ... 'C:\SERVER\daisy34\daisy34_20260309', ...
    ... 'C:\SERVER\desmond43\desmond43_20260309'...
    ... 'C:\SERVER\desmond43\desmond43_20260310'...
    ... 'C:\SERVER\daisy34\daisy34_20260325_incompleteStim' ...
    'C:\SERVER\daisy33\daisy33_20260320', ... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch2)
    'C:\SERVER\daisy34\daisy34_20260320', ... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch2)
    'C:\SERVER\daisy33\daisy33_20260323', ... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch2)
    ... 'C:\SERVER\daisy34\daisy34_20260323', ... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch1)
    'C:\SERVER\daisy33\daisy33_20260325', ... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch2)
    ... 'C:\SERVER\daisy33\daisy33_20260326', ... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch1)
    ... 'C:\SERVER\daisy34\daisy34_20260326', ... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch1)
    ... 'C:\SERVER\daisy33\daisy33_20260327', ... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch1)
    ... 'C:\SERVER\daisy34\daisy34_20260327', ... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch1)
    ... 'C:\SERVER\desmond43\desmond43_20260311'... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch1)
    ... 'C:\SERVER\desmond43\desmond43_20260312'... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch1)
    ... 'C:\SERVER\desmond43\desmond43_20260313'... HIGHERPOWER, SORTED, CONVERTED TO EU (Batch1)
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
            eu = EphysUnit(ar, readWaveforms=false, cullITI=false, savepath='C:\SERVER\Units\TwoColor_Striatonigral\Batch2', tr=tr);

            tr.Spikes = [];
            clear eu
        end
    catch ME
        warning('Could not process folder: %s', folders{iSession})
        warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
    end
end

%%
clear
eu = EphysUnit.load('C:\SERVER\Units\TwoColor_Striatonigral\Batch2', waveforms=false, spikecounts=false, spikerates=false);
% Remove multiunits, fast (ISS test)
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

%%

clear clc
eu = EphysUnit.load('C:\SERVER\Units\TwoColor_Striatonigral\SingleUnit_NonDuplicate_NonDrift_SNr', waveforms=false, spikecounts=false, spikerates=false, ...
    animalNames={'daisy33', 'daisy34', 'desmond43'});

% Do psth
rd.stim = eu.getRasterData('stimtwocolor', window=[-0.1, 0.4], durErr=1e-3, shutterDelay=0, photoelectricBlankDuration=0.5e-3, photoelectricOffsetBlankWindow=[10e-3, 10.5e-3]);

%% ETA Stim
close all

p.path = 'E:\DATA\Figures\TwoColor_SNr_Striatonigral\ErinsBatch\WavelengthPowerLocation\PSTH';

p.isiWindow = [-0.5, 0.5];
p.isiRes = 1e-3;
p.isiBaselineWindow = [-0.2, 0];
p.xlim.stim = [-0.05, 0.1];
p.xlim.move = [-4, 2];
p.rasterSzStim = 3;
p.rasterSzMove = 1;
p.mode = "psth"; % "isi", "psth"
p.clim = [-2, 2];
% p.mode = "isi"; % "isi", "psth"
% p.clim= [-8, 8];

p.artifacts = struct(event=[], length=[], lengthUnit=[], direction=[]);
p.artifacts(1) = struct(event='StimOn', length=0.5, lengthUnit='ms', direction='right');
p.artifacts(2) = struct(event='StimOff', length=0.5, lengthUnit='ms', direction='right');

% Combined PEISI and Stim Raster
% Stim Rasters

fig = figure(Units='inches', Position=[1, 1, 5.5, 10]);
clear layout

layout.w = [3]; % Stim
layout.h = [1, 3, 7]; % Raster, peth

layout.tl = tiledlayout(fig, sum(layout.h), sum(layout.w), TileSpacing='tight', Padding='tight', TileIndexing='columnmajor');
layout.ax = gobjects(length(layout.h), length(layout.w));
layout.ax(1, 1) = nexttile(layout.tl, [sum(layout.h(1:2)), layout.w(1)]);
layout.ax(2, 1) = nexttile(layout.tl, [layout.h(3), layout.w(1)]);
% layout.ax(1, 2) = nexttile(layout.tl, [layout.h(1), layout.w(2)]);
% layout.ax(2, 2) = nexttile(layout.tl, [layout.h(2), layout.w(2)]);
% layout.ax(3, 2) = nexttile(layout.tl, [layout.h(3), layout.w(2)]);

if ~exist(p.path, 'dir')
    mkdir(p.path)
end
for iEu = 1:length(eu)
    % Make Raster
    ax = layout.ax(1, 1);
    cla(ax)
    EphysUnit.plotRaster(ax, rd.stim(iEu), xlim=p.xlim.stim*1e3, sz=p.rasterSzStim, timeUnit='ms');
    % delete(ax.Legend)
    % ax.Legend.Location = 'northeast';
    % ax.Legend.FontSize = 6;
    title(ax, eu(iEu).getName(), Interpreter="none")

    % Make ETA
    ax = layout.ax(2, 1);
    cla(ax)
    groups = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'location'});
    sepPower = find([groups.location] == min([groups.location]));
    sepWavelength = find([groups.wavelength] == groups(1).wavelength, 1, 'last');
    switch p.mode
        case "isi"
            isi = NaN(length(groups), length(p.isiWindow(1):p.isiRes:p.isiWindow(2)));
            deltaSR = isi;
            for iGrpMin = 1:length(groups)
                [isi(iGrpMin, :), t] = eu(iEu).getMeanPEISI('stimtwocolor', groups(iGrpMin).trials, window=p.isiWindow, resolution=p.isiRes, ...
                    photoelectricBlankDuration=0.5e-3, photoelectricOffsetBlankWindow=[10e-3, 10.5e-3]);
                baseline = 1./isi(iGrpMin, isin(t, p.isiBaselineWindow));
                deltaSR(iGrpMin, :) = (1./isi(iGrpMin, :) - mean(baseline, 'all', 'omitnan')) ./ std(baseline, 0, 'all', 'omitnan');
                clear baseline
            end
        case "psth"
            clear deltaSR
            for iGrpMin = 1:length(groups)
                etaTemp = eu(iEu).getETA('count', 'stimtwocolor', [-0.5, 0.5], normalize=[-0.5, 0], resolution=0.010, alignTo='start', ...
                    trials = groups(iGrpMin).trials, artifacts=p.artifacts);
                t = etaTemp.t;
                deltaSR(iGrpMin, :) = etaTemp.X;
            end
    end
    imagesc(ax, 1e3*t, [], deltaSR)
    xline(ax, 0)
    yline(ax, sepPower + 0.5, 'k--', LineWidth=1)
    yline(ax, sepWavelength + 0.5, 'k-', LineWidth=1.5)
    xlim(ax, p.xlim.stim*1e3)
    ax.YAxisLocation = 'right';
    applyCustomColormap(ax, p.clim, hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);
    h = colorbar(ax, 'westoutside');
    h.Label.String = '\DeltaSR (a.u.)';
    sepPower = [0, sepPower];
    yticks(ax, 0.5*(sepPower(1:end-1) + sepPower(2:end)));
    yticklabels(ax, arrayfun(@(grp) sprintf("%gmW %inm", grp.power*1e3, grp.wavelength), groups(sepPower(2:end))));
    ax.YAxis.TickLength = [0, 0];
    % yticks(ax, 1:length(groups));
    % yticklabels(ax, {groups.label})
    xlabel(ax, 'Time from opto onset (ms)')
    switch p.mode
        case "isi"
            title(ax, 'PSTH (from ISI)')
        case "psth"
            title(ax, 'PSTH')
    end

    fontsize(fig, 9, 'points');
    print(fig, sprintf('%s\\%s.png', p.path, eu(iEu).getName()), '-dpng', '-r0')
end

clear fig layout iEu ax groups sepPower sepWavelength deltaSR iGrpMin etaTemp t h isi

%% Calculate data for population heatmap
switch p.mode
    case "isi"
        t = p.isiWindow(1):p.isiRes:p.isiWindow(2);
    case "psth"
        t = -0.5:0.01:0.5;
        t = 0.5*(t(1:end-1) + t(2:end));
end
X = NaN(length(eu), 80, length(t));
groupsGlobal = eu(1).groupTwoColorStimTrials({'wavelength', 'power', 'location'});
groupsGlobal = rmfield(groupsGlobal, {'trials', 'pulseIndices'});
assert(length(groupsGlobal) == 80)

ll = 0;
tTicTotal = tic();
for iEu = 1:length(eu)
    fprintf(repmat('\b', [1, ll]))
    ll = fprintf("iEu=%i/%i...(%.1fs)\n", iEu, length(eu), toc(tTicTotal));
    groups = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'location'});
    switch p.mode
        case "isi"
            isi = NaN(length(groups), length(p.isiWindow(1):p.isiRes:p.isiWindow(2)));
            deltaSR = NaN(length(groupsGlobal), length(p.isiWindow(1):p.isiRes:p.isiWindow(2)));
            for iGrpMin = 1:length(groups)
                iGrpGlobal = find(string({groupsGlobal.label}) == string(groups(iGrpMin).label));
                assert(~isempty(iGrpGlobal))
                [isi(iGrpMin, :), t] = eu(iEu).getMeanPEISI('stimtwocolor', groups(iGrpMin).trials, window=p.isiWindow, resolution=p.isiRes, ...
                    photoelectricBlankDuration=0.5e-3, photoelectricOffsetBlankWindow=[10e-3, 10.5e-3]);
                baseline = 1./isi(iGrpMin, isin(t, p.isiBaselineWindow));
                deltaSR(iGrpGlobal, :) = (1./isi(iGrpMin, :) - mean(baseline, 'all', 'omitnan')) ./ std(baseline, 0, 'all', 'omitnan');
                clear baseline
            end
        case "psth"
            deltaSR = NaN(length(groupsGlobal), size(X, 3));
            for iGrpMin = 1:length(groups)
                iGrpGlobal = find(string({groupsGlobal.label}) == string(groups(iGrpMin).label));
                assert(~isempty(iGrpGlobal))
                etaTemp = eu(iEu).getETA('count', 'stimtwocolor', [-0.5, 0.5], normalize=[-0.5, 0], resolution=0.010, alignTo='start', ...
                    trials = groups(iGrpMin).trials, artifacts=p.artifacts);
                t = etaTemp.t;
                deltaSR(iGrpGlobal, :) = etaTemp.X;
            end
    end
    X(iEu, :, :) = deltaSR;
end

psth = struct(X=X, t=t, groups=groupsGlobal);

clear ll iEu isi deltaSR iGrpMin etaTemp tTicTotal groups iGrpGlobal X t groupsGlobal

%% Plot population heatmap
close all
p.metaWindow = [0.01, 0.03];
p.metaThreshold = 0.2;
p.clim = [-1.5, 1.5];

metaX = transpose(mean(psth.X(:, :, isin(psth.t, p.metaWindow)), 3, 'omitnan')); % Location x Power x Wavelength
metaX(isnan(metaX)) = 0;
metaXHiPower = squeeze(mean(reshape(metaX([33:40, 73:80], :, :), 8, 2, length(eu)), 1, 'omitnan'));
B = sign(metaXHiPower) .* (abs(metaXHiPower)>p.metaThreshold);
B(isnan(B)) = 0;
B = B + 1;
hash = B(1, :)*3 + B(2, :);

[~, order] = sort(hash);

fig = figure(Units='inches', Position=[0, 1, 14, 4]);
ax = axes(fig);

groups = psth.groups;
sepPower = find([groups.location] == min([groups.location]));
sepWavelength = find([groups.wavelength] == groups(1).wavelength, 1, 'last');
sepHash = arrayfun(@(i) find(hash(order)==i, 1, 'last'), min(hash):max(hash), UniformOutput=false);
sepHash = [sepHash{:}];

imagesc(ax, metaX(:, order))
xline(ax, sepHash + 0.5, 'k-', LineWidth=1)
yline(ax, sepPower + 0.5, 'k--', LineWidth=1)
yline(ax, sepWavelength + 0.5, 'k-', LineWidth=1.5)
ax.YAxisLocation = 'right';
applyCustomColormap(ax, p.clim, hlim=[0.375, 0, 0, -0.375], llim=[0.2, 1, 1, 0.3], hpwr=.3, lpwr=0.33, h0=0.33);
h = colorbar(ax, 'westoutside');
h.Label.String = 'normalized stim response (a.u.)';
sepPower = [0, sepPower];
yticks(ax, 0.5*(sepPower(1:end-1) + sepPower(2:end)));
yticklabels(ax, arrayfun(@(grp) sprintf("%gmW %inm", grp.power*1e3, grp.wavelength), groups(sepPower(2:end))));
ax.YAxis.TickLength = [0, 0];
xticks(ax, sepHash+0.5)
xticklabels(ax, sepHash);

% Binary labesl
uniqueB = B(:, order) - 1;
uniqueB = uniqueB(:, sepHash);
for iCat = 1:length(sepHash)
    if iCat == 1
        center = (sepHash(1) + 0) / 2;
    else
        center = mean(sepHash(iCat-1:iCat));
    end
    text(ax, center, 20.5, string(uniqueB(1, iCat)), HorizontalAlignment='center', VerticalAlignment='middle', FontWeight='bold', FontSize=9);
    text(ax, center, 60.5, string(uniqueB(2, iCat)), HorizontalAlignment='center', VerticalAlignment='middle', FontWeight='bold', FontSize=9);
end
clear iCat center

xlabel(ax, 'unit')
ylabel(ax, 'wavelength x power x position', Rotation=-90)
switch p.mode
    case "isi"
        title(ax, 'Stim Response (from ISI)')
    case "psth"
        title(ax, sprintf('SNr response to striatal stimulation %i-%i ms', 1e3*p.metaWindow(1), 1e3*p.metaWindow(2)))
end

fontsize(fig, 9, 'points');
copygraphics(fig, ContentType='vector', BackgroundColor='none')

B = B-1;
meta = struct(X=metaX, XHiPower=metaXHiPower, B=B, hash=hash, order=order);

clear metaX metaXHiPower B uniqueB hash order fig ax groups sepPower sepWavelength sepHash h

%% Plot receptive fields (per power x wavelength, averaged across units)
close all

path = 'E:\DATA\Figures\TwoColor_SNr_Striatonigral\ErinsBatch\RF';
if ~isfolder(path)
    mkdir(path)
end
if ~isfolder(fullfile(path, 'avgAcrossPwr'))
    mkdir(fullfile(path, 'avgAcrossPwr'))
end

selUnits = find(any(meta.B ~= 0, 1));
xAligned = NaN(8, length(selUnits), 2, 5); % Location x unit x wavelength x power
xRaw = xAligned;
for iPwr = 1:5
    for iWavelength = 1:2
        selConds = 8*(iPwr-1) + 1 : 8*iPwr;
        selConds = selConds + (iWavelength-1)*40;
        x = abs(meta.X(selConds, selUnits));
        for iUnit = 1:size(x, 2)
            [~, iPeak] = max(x(:, iUnit));
            xAligned(:, iUnit, iWavelength, iPwr) = circshift(x(:, iUnit), 4-iPeak, 1);
            xRaw(:, iUnit, iWavelength, iPwr) = meta.X(selConds, iUnit);
        end
    end
end

% Plot one figure, average across units
fig = figure(Units='inches', Position=[1, 1, 7, 3]);
ax = axes(fig);
h = gobjects(5, 2);
hold(ax, 'on')
hues = [0.61, 0.14];
xMed = mean(xAligned(:, :, :, 5), 2, 'omitnan');
for iPwr = 1:5
    for iWavelength = 1:2
        color = hsl2rgb([hues(iWavelength), 1, 1-(iPwr-1)/8]);
        x = mean(xAligned(:, :, iWavelength, iPwr), 2);
        iGrpMin = 1+8*(iPwr-1)+40*(iWavelength-1);
        name = sprintf("%inm %gmW", psth.groups(iGrpMin).wavelength, psth.groups(iGrpMin).power*1e3);
        h(iPwr, iWavelength) = plot(ax, 1:8, (x-min(xMed, [], 'all'))./(max(xMed, [], 'all') - min(xMed, [], 'all')), Color=color, Marker='o', MarkerSize=iPwr*1.5, DisplayName=name);
    end
end
hold(ax, 'off')
xlabel(ax, 'stim location')
ylabel(ax, 'stim response (a.u.)')
legend(ax, h, Location='eastoutside')
title(ax, sprintf("%i units (average)", length(selUnits)), Interpreter='none')
print(fig, fullfile(path, sprintf('%i units (average).png', length(selUnits))), '-dpng', '-r0')

% % Plot one figure per unit
% doNorm = false;
% for iUnit = 1:size(xAligned, 2)
%     cla(ax)
%     hold(ax, 'on')
%     xMin = min(xRaw(:, iUnit, :, :), [], 'all');
%     xMax = max(xRaw(:, iUnit, :, :), [], 'all');
%     for iPwr = 1:5
%         for iWavelength = 1:2
%             color = hsl2rgb([hues(iWavelength), 1, 1-(iPwr-1)/8]);
%             x = xRaw(:, iUnit, iWavelength, iPwr);
%             iGrp = 1+8*(iPwr-1)+40*(iWavelength-1);
%             name = sprintf("%inm %gmW", psth.groups(iGrp).wavelength, psth.groups(iGrp).power*1e3);
%             if doNorm
%                 plot(ax, 1:4, (x(1:4)-xMin)./(xMax-xMin), Color=color, Marker='o', MarkerSize=iPwr*1.5);
%                 plot(ax, 5:8, (x(5:8)-xMin)./(xMax-xMin), Color=color, Marker='o', MarkerSize=iPwr*1.5);
%             else
%                 plot(ax, 1:4, x(1:4), Color=color, Marker='o', MarkerSize=iPwr*1.5);
%                 h(iPwr, iWavelength) = plot(ax, 5:8, x(5:8), Color=color, Marker='o', MarkerSize=iPwr*1.5, DisplayName=name);
%             end
%         end
%     end
%     yline(ax, 0, 'k--')
%     xticks(ax, [1, 4, 5, 8])
%     xticklabels(ax, ["DLS", "VLS", "DMS", "VMS"])
%     xlabel(ax, 'stim location')
%     ylabel(ax, 'stim response (a.u.)')
%     legend(ax, h, Location='eastoutside')
%     hold(ax, 'off')
%     title(ax, eu(selUnits(iUnit)).getName(), Interpreter='none')
%     ylim(ax, [-1, 1])
%     xlim(ax, [0, 9])
%     print(fig, fullfile(path, sprintf('%s.png', eu(selUnits(iUnit)).getName())), '-dpng', '-r0')
% end

% Plot one figure per unit, average scross powers
h = gobjects(1, 2);
pwrRange = [3, 5];
for iUnit = 1:size(xAligned, 2)
    cla(ax)
    hold(ax, 'on')
    for iWavelength = 1:2
        color = hsl2rgb([hues(iWavelength), 1, 0.5]);
        x = mean(xRaw(:, iUnit, iWavelength, pwrRange(1):pwrRange(2)), 4, 'omitnan');
        iGrpMin = 1+8*(pwrRange(1)-1)+40*(iWavelength-1);
        iGrpMax = 1+8*(pwrRange(2)-1)+40*(iWavelength-1);
        name = sprintf("%inm %g-%gmW", psth.groups(iGrpMin).wavelength, psth.groups(iGrpMin).power*1e3, psth.groups(iGrpMax).power*1e3);
        plot(ax, 1:4, x(1:4), Color=color, Marker='o', MarkerSize=iPwr*1.5);
        h(iWavelength) = plot(ax, 5:8, x(5:8), Color=color, Marker='o', MarkerSize=5*1.5, DisplayName=name);
    end
    yline(ax, 0, 'k--')
    xticks(ax, [1, 4, 5, 8])
    xticklabels(ax, ["DLS", "VLS", "DMS", "VMS"])
    xlabel(ax, 'stim location')
    ylabel(ax, 'stim response (a.u.)')
    legend(ax, h, Location='eastoutside')
    hold(ax, 'off')
    title(ax, eu(selUnits(iUnit)).getName(), Interpreter='none')
    ylim(ax, [-1, 1])
    xlim(ax, [0, 9])
    print(fig, fullfile(path, 'avgAcrossPwr', sprintf('%s.png', eu(selUnits(iUnit)).getName())), '-dpng', '-r0')
end

clear selUnits iPwr iWavelength selConds x iUnit iPeak xMax xMin hues h name iGrpMin color xAligned xMed xRaw fig ax doNorm path iGrpMin iGrpMax
