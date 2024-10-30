load_ephysunits
euNames = lower(eu.getName());

%%
sel = c.hasPress & c.hasLick;
expNames = unique(lower({eu(sel).ExpName}));
%%
for iExp = 4:length(expNames)
    try
        clear tr ar bmr eu
        tr = TetrodeRecording.BatchLoadSimple(expNames{iExp}, true);

        % Make AR
        ar = AcuteRecording(tr, 'WT');
        bmr = ar.binMoveResponse(tr, 'Press', 'Window', [-1, 0], 'BaselineWindow', [-1, 0], 'Store', true);
        
        % Make EphysUnits from AR
        eu = EphysUnit(ar, tr=tr, savepath='C:\SERVER\Units\NonLite_PressVsLick', cullITI=false, readWaveforms=false);
        
    catch ME
        warning('Error while processing file %g (%s)', iExp, expNames{iExp});
    end
end



% files = dir('C:\SERVER\Units\Old\*.mat');
% sel = ismember(cellfun(@(n) lower(strrep(n, '.mat', '')), {files.name}, UniformOutput=false), euNames);
% files = files(sel);
% cd('C:\SERVER\Units\Old\')
% euOld = EphysUnit.load({files.name});
% 
% %% Combined old and new units (old units have all spikes, new ones have ITI culled)
% oldNames = euOld.getName();
% newNames = eu.getName();
% 
% [lia, locb] = ismember(newNames, oldNames);
% sel = c.hasPress & c.hasLick & c.isLick;
% 
% euMixed = eu;
% euMixed(lia) = euOld(locb(lia));
% 
% %%

%% Load euComplete (complete with ITI spikes)

files = dir('C:\SERVER\Units\NonLite_PressVsLick\*.mat');
sel = ismember(cellfun(@(n) lower(strrep(n, '.mat', '')), {files.name}, UniformOutput=false), euNames);
files = files(sel);
cd('C:\SERVER\Units\NonLite_PressVsLick\')
euComplete = EphysUnit.load({files.name});

% euComplete.save('C:\SERVER\Units\NonLite_PressVsLick_NonDuplicate_NonDrift');

%% Make ETA
euMixed = eu;

[lia, locb] = ismember(euComplete.getName(), eu.getName());
euMixed(lia) = euComplete(locb(lia));

%%
eta.lickBoutEnd = euMixed.getETA('count', 'lickboutend', window=[0, 2], resolution=[2*pi/12, 0.01], normalize='none', ...
    maxInterval=0.2, minInterval=0.05, minBoutCycles=8);
eta.lickBoutEndNorm = euMixed.getETA('count', 'lickboutend', window=[0, 2], resolution=[2*pi/12, 0.01], normalize=[1, 2], ...
    maxInterval=0.2, minInterval=0.05, minBoutCycles=8);