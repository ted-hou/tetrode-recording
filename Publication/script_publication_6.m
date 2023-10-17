%%% Figure 6. Lick vs. Reach ephys + Osci lick and somatotopy

%% Load all units
load_ephysunits;
boot_response_dir;


%% plot params
p.fontSize = 8;
p.width = 8;
p.height = 10;
p.heightFirstRow = 3;


%% 6a. Double rasters (reach vs lick) 4 example untis
unitNames = { ...
    'desmond24_20220510_Channel44_Unit1'; ...
    'Daisy2_20180420_Channel14_Unit1'; ...
    'daisy13_20220106_Electrode97_Unit1'; ...
    'daisy8_20210709_Channel7_Unit1'; ...
    };
files = cellfun(@(name) sprintf('C:\\SERVER\\Units\\Lite_NonDuplicate\\%s.mat', name), unitNames, UniformOutput=false);
euEg_lickVsPress = EphysUnit.load(files);

close all
for iEu = 1:length(unitNames)
    thisRdPress = euEg_lickVsPress(iEu).getRasterData('press', window=[0, 2], sort=true);
    thisRdLick = euEg_lickVsPress(iEu).getRasterData('lick', window=[0, 2], sort=true);
    ax = EphysUnit.plotMultiRaster([thisRdPress, thisRdLick], label=["Reach", "Lick"], ...
        xlim=[-3, 2], iti=false, sz=0.75, maxTrials=80, ...
        figUnits='inches', figPosition=[0, 0, p.width/length(unitNames), p.heightFirstRow], figMargins=1.2*[0.16, 0.09]);
    l = findobj(Type='Legend', Parent=ax(1).Parent);
    if iEu == 1
        delete(l(1))
        l(2).Orientation = 'horizontal';
        l(2).Position = [0.052677548682703,0.947314293123628,0.906249988824129,0.0442708323244];
    else
        delete(l)
        ylabel(ax, '')
    end
    xlabel(ax(1), '')
    xlabel(ax(2), 'Time to contact (s)')
    l = findobj(Type='Legend', Parent=ax(1).Parent);
    fontsize([ax; l], p.fontSize, 'points')
    fontname([ax; l], 'Arial')
end
clear thisRd ax iEu

%% 6b. PETH Reach vs. Lick vs. Osci Lick
