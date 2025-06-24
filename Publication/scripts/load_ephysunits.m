%% 1.2Alt Or just load lite version, without non-SNr cells, without waveforms, spikecounts or spikerates.
if exist('E:\Data', 'dir') 
    eu = EphysUnit.load('E:\Data\Units\SNr_nonDuplicate_nonDrift_withITI', waveforms=false, spikecounts=false, spikerates=false);
else
    eu = EphysUnit.load('C:\SERVER\Units\Lite_NonDuplicate_NonDrift', waveforms=false, spikecounts=false, spikerates=false);
end

%% Load euComplete (complete with ITI spikes)
if ~exist('E:\Data', 'dir') 
    euNames = lower(eu.getName());
    
    files = dir('C:\SERVER\Units\NonLite_PressVsLick\*.mat');
    sel = ismember(cellfun(@(n) lower(strrep(n, '.mat', '')), {files.name}, UniformOutput=false), euNames);
    files = files(sel);
    cd('C:\SERVER\Units\NonLite_PressVsLick\')
    euComplete = EphysUnit.load({files.name}, waveforms=false, spikecounts=false, spikerates=false);
    
    % euComplete.save('C:\SERVER\Units\NonLite_PressVsLick_NonDuplicate_NonDrift');
    
    [lia, locb] = ismember(eu.getName(), euComplete.getName());
    eu(lia) = euComplete(locb(lia));
    
    clear euComplete
end

%% Load metadata
if exist('E:\Data', 'dir')
    load('E:\Data\Units\meta_Lite_NonDuplicate_NonDrift.mat')
else
    load('C:\SERVER\Units\meta_Lite_NonDuplicate_NonDrift.mat')
end
% save('C:\SERVER\Units\meta_Lite_NonDuplicate_NonDrift.mat')
