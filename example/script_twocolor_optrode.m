%% Load TetrodeRecording objects and convert to EphysUnits
sessionInfo(1) = struct('name', 'daisy19_20240410',     'ml', +2300, 'ap', 600, 'dv', -3000, 'duraOffset', -950, 'facing', 'left', 'model', '128DN');
sessionInfo(2) = struct('name', 'daisy19_20240411',     'ml', +2300, 'ap', 600, 'dv', -3000, 'duraOffset', -950, 'facing', 'left', 'model', '128DN');
sessionInfo(3) = struct('name', 'daisy20_20240416',     'ml', -2300, 'ap', 600, 'dv', -3000, 'duraOffset', -950, 'facing', 'left', 'model', '128DN');
sessionInfo(4) = struct('name', 'desmond36_20240419',   'ml', +2300, 'ap', 600, 'dv', -3000, 'duraOffset', -950, 'facing', 'left', 'model', '128DN');
sessionInfo(5) = struct('name', 'desmond36_20240425',   'ml', +2300, 'ap', 600, 'dv', -3000, 'duraOffset', -950, 'facing', 'left', 'model', '128DN');
sessionInfo(6) = struct('name', 'desmond37_20240415',   'ml', -2300, 'ap', 600, 'dv', -3000, 'duraOffset', -950, 'facing', 'left', 'model', '128DN');

for iExp = 2:length(sessionInfo)
    try
        files = dir(sprintf('Y:\\tr_sorted_%s*.mat', sessionInfo(iExp).name));
        files = arrayfun(@(x) [x.folder, x.name], files, UniformOutput=false);
        tr = TetrodeRecording.BatchLoadSimple(files);
        ar = AcuteRecording(tr, 'N/A');
        ar.binMoveResponse(tr, 'none', Window=[-1, 0], Store=true);
        eu = EphysUnit(ar, readWaveforms=true, cullITI=false, savepath='Y:\Units\TwoColor_Optrode', tr=tr, sessionInfo=sessionInfo(iExp));
        clear tr ar eu
    catch
        warning('Failed to process session %i', iExp)
    end
end
clear iExp
%% Load EU
eu = EphysUnit.load('Y:\Units\TwoColor_Optrode\NonDuplicate_SingleUnit');

%% Make rasters and save to disk
rd = eu.getRasterData('stimtwocolor', window=[-0.1, 0.4], durErr=1e-2, shutterDelay=0);
%%
fig = figure(Units='inches', Position=[0, 0, 6, 8]);
ax = axes(fig);
if ~exist('Y:\Figures\TwoColor_Optrode', 'dir')
    mkdir('Y:\Figures\TwoColor_Optrode')
end
for iEu = 1:length(eu)
    cla(ax)
    EphysUnit.plotRaster(ax, rd(iEu), xlim=[-0.1, 0.3], sz=2);
    print(fig, sprintf('Y:\\Figures\\TwoColor_Optrode\\%s.png', eu(iEu).getName()), '-dpng');
end

%% Make PETH and save to disk

%% Calculate peristim ISI, average across trials, group by condition. We use this in lieu of PSTHs


%% Plot 


%  tr = TetrodeRecording.BatchLoadSimple();
%  % Try daisy19_20240411, Chn15, Cluster 2
% %% Load two color experiment and validate stim pulse timing
% file = dir(sprintf('%s\\..\\%s.mat', tr.Path, tr.GetExpName(includeSuffix=false)));
% assert(~isempty(file))
% 
% tce = load(sprintf('%s\\%s', file.folder, file.name));
% tce = tce.obj;
% assert(isa(tce, 'TwoColorExperiment'))
% clear file
% 
% tOn = tr.DigitalEvents.StimOn;
% tOff = tr.DigitalEvents.StimOff;
% 
% nPulsesPerTrain = arrayfun(@(log) log.params.nPulses, tce.Log);
% pulseWidth = arrayfun(@(log) repmat(log.params.pulseWidth, [1, log.params.nPulses]), tce.Log, UniformOutput=false);
% % pulseWidth = cat(2, pulseWidth{:});
% 
% % Verify congruence between Intan shutter events and TwoColorExperiment.Log
% assert(sum(nPulsesPerTrain) == nnz(tOn), 'Intan recorded %i shutterOn events while TwoColorExperiment.Log has %i.', nnz(tOn), sum(nPulsesPerTrain));
% assert(sum(nPulsesPerTrain) == nnz(tOff), 'Intan recorded %i shutterOn events while TwoColorExperiment.Log has %i.', nnz(tOn), sum(nPulsesPerTrain));
% nPulses = nnz(tOn);
% fprintf('Intan & TwoColorExperiment.Log both recorded %i shutter pulses.\n', nnz(tOn))
% 
% % Pulse durations should match between Intan and TCE as well.
% pulseWidthErrorMargin = 1e-3;
% assert(all(abs(tOff - tOn - pulseWidth) < pulseWidthErrorMargin), 'Only %i/%i pulseWidths agree (df<%gs).', nnz(abs(tOff - tOn - pulseWidth) < pulseWidthErrorMargin), nPulses, pulseWidthErrorMargin)
% fprintf('All %i pulseWidths agree (df<%gs).\n', nPulses, pulseWidthErrorMargin)
% 
% % Generate pulse->train map
% iPulse = 0;
% pulseToTrain = NaN(1, nPulses);
% for iTrain = 1:length(tce.Log)
%     pulseToTrain(iPulse + 1:iPulse + nPulsesPerTrain(iTrain)) = iTrain;
%     iPulse = iPulse + nPulsesPerTrain(iTrain);
% end
% 
% clear iPulse iTrain pulseWidth
% 
% %% Make one raster
% iChn = 15;
% iUnit = [3];
% % Stim condition has
% % 1000000*wavelength (1, 2) + 100000*location (1, 2, 3, 4) + 1000*duration (1, 2, 4, 6, 10, 20, 40) + 1*nominalpower (1, 2, 4, 20, 80, 320)
% 
% trainHash = tce.getStimHash();
% pulseHash = trainHash(pulseToTrain);
% 
% [~, trainOrder] = sort(trainHash, 'ascend');
% [~, pulseOrder] = sort(pulseHash, 'ascend');
% 
% % Get pulse on/pulse off times, sorted by hash
% % edges = [tOn(pulseOrder); tOff(pulseOrder)];
% 
% % assert(sizes(edges, 2) == nPulses)
% % assert(size(edges, 1) == 2)
% tOnSorted = tOn(pulseOrder);
% tOffSorted = tOff(pulseOrder);
% 
% st = tr.Spikes(iChn).Timestamps(ismember(tr.Spikes(iChn).Cluster.Classes, iUnit));
% stAligned = st;
% I = NaN(size(st));
% for iPulse = 1:nPulses
%     sel = st >= tOnSorted(iPulse) - 0.25 & st <= tOnSorted(iPulse) + 0.5;
%     I(sel) = iPulse;
%     stAligned(sel) = st(sel) - tOnSorted(iPulse);
% end
% 
% close all
% fig = figure();
% ax = axes(fig);
% hold(ax, 'on')
% scatter(ax, stAligned, I, 5, 'k', 'filled');
% ylim(ax, [0, nPulses + 1])
% ax.YAxis.Direction = 'reverse';
% 
% for iTrain = 1:length(tce.Log)
%     selPulse = pulseToTrain == iTrain;
%     iPulseStart = strfind(selPulse, [0, 1]) + 1;
%     iPulseEnd = strfind(selPulse, [1, 0]);
%     if isempty(iPulseStart)
%         iPulseStart = 1;
%     end
%     if isempty(iPulseEnd)
%         iPulseEnd = nPulses;
%     end
%     iPulseStart = find(pulseOrder == iPulseStart);
%     iPulseEnd = find(pulseOrder == iPulseEnd);
%     pulseWidth = tce.Log(iTrain).params.pulseWidth;
%     switch tce.Log(iTrain).wavelength
%         case 473
%             color = [0.2, 0.2, 0.8];
%         case 593
%             color = [0.8, 0.2, 0.2];
%     end
%     alpha = 0.2 + 0.6*((tce.Log(iTrain).params.iPower - 1)./(length(tce.Params.targetPowers) - 1));
%     % assert(tce.Log(iTrain).targetPower == 25e-6)
%     patch(ax, [0, pulseWidth, pulseWidth, 0], [iPulseStart, iPulseStart, iPulseEnd, iPulseEnd], color, ...
%         FaceAlpha=alpha, EdgeColor='none')
% end
% hold(ax, 'off')