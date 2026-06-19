%% Set root
ROOTPATH = "C:\SERVER";
% ROOTPATH = 'E:\DATA';

%% Load eu, exp
folders = [ ...
    % "C:\SERVER\daisy26\daisy26_20250424", ...
    % "C:\SERVER\daisy26\daisy26_20250425", ...
    % "C:\SERVER\desmond38\desmond38_20250401", ...
    % "C:\SERVER\desmond38\desmond38_20250402", ...
    % "C:\SERVER\desmond38\desmond38_20250407", ...
    % "C:\SERVER\desmond38\desmond38_20250417", ...
    % "C:\SERVER\desmond39\desmond39_20250416", ...
    % "C:\SERVER\desmond39\desmond39_20250423", ...
    "C:\SERVER\daisy27\daisy27_20250624", ...
    "C:\SERVER\daisy27\daisy27_20250626", ...
    "C:\SERVER\daisy27\daisy27_20250707", ...
    "C:\SERVER\daisy27\daisy27_20250715", ...
    "C:\SERVER\daisy27\daisy27_20250717", ...
    "C:\SERVER\daisy27\daisy27_20250721", ...
    "C:\SERVER\daisy27\daisy27_20250724", ...
    "C:\SERVER\daisy28\daisy28_20250701", ...
    "C:\SERVER\daisy28\daisy28_20250702", ...
    "C:\SERVER\daisy28\daisy28_20250716", ...
    "C:\SERVER\daisy28\daisy28_20250718", ...
    "C:\SERVER\daisy28\daisy28_20250723", ...
    "C:\SERVER\daisy28\daisy28_20250725", ...
    "C:\SERVER\daisy28\daisy28_20250729", ...
    "C:\SERVER\daisy28\daisy28_20250714", ...
    "C:\SERVER\daisy29\daisy29_20251023", ...
    "C:\SERVER\daisy29\daisy29_20251024", ...
    "C:\SERVER\daisy29\daisy29_20251025", ...
    "C:\SERVER\daisy29\daisy29_20251027", ...
    "C:\SERVER\daisy29\daisy29_20251028", ...
    "C:\SERVER\daisy29\daisy29_20251031", ...
    "C:\SERVER\daisy29\daisy29_20251118", ...
    "C:\SERVER\daisy29\daisy29_20251120", ...
    "C:\SERVER\daisy29\daisy29_20251103", ...
    "C:\SERVER\daisy30\daisy30_20251029", ...
    "C:\SERVER\daisy30\daisy30_20251030", ...
    "C:\SERVER\daisy30\daisy30_20251119", ...
    "C:\SERVER\daisy30\daisy30_20251121", ...
    "C:\SERVER\desmond41\desmond41_20251028", ...
    "C:\SERVER\desmond41\desmond41_20251029", ...
    "C:\SERVER\desmond42\desmond42_20251030", ...
    "C:\SERVER\desmond42\desmond42_20251031", ...
    "C:\SERVER\desmond41\desmond41_20251104", ...
    "C:\SERVER\desmond41\desmond41_20251117", ...
    "C:\SERVER\desmond42\desmond42_20251121", ...
    "C:\SERVER\desmond41\desmond41_20251124", ...
    ];

euFolders = [ ...
    % "C:\SERVER\Units\TwoColor_SNr_SCRetro\SingleUnit_NonDuplicate_NonDrift_SNr", ... daisy26, desmond38, desmond39
    "C:\SERVER\Units\TwoColor_SNr_SCRetro\ReverseInjection\SingleUnit_NonDuplicate_NonDrift_SNr_withTrials", ... daisy27, 28
    "C:\SERVER\Units\TwoColor_Striatonigral\SingleUnit_NonDuplicate_NonDrift_SNr", ... daisy29, 30, desmond41, 42
];

dlcResultsPath = 'C:\SERVER\DeepLabCut\Results\FourPawsTongueJawSpine';

% Load ephysunits
eu = cell(length(euFolders), 1);
for iDir = 1:length(euFolders)
    eu{iDir} = EphysUnit.load(euFolders(iDir));
end

eu = cat(2, eu{:});


% Make CompleteExperiment3

exp = CompleteExperiment3(eu, cameras='lr', deeplabcutPath=dlcResultsPath);

exp.alignTimestamps(refEventNameArduino={'OPTO1_ON', 'OPTO2_ON', 'REWARD_ON'}, refEventNameEphys={'StimOn', 'RewardOn'}, trialDurationTolerance=2);

clear results
results(length(exp)) = struct(name=[], varsL=[], hasTimestampsL=[], varsR=[], hasTimestampsR=[], isValid=[]);
for iExp = 1:length(exp)
    results(iExp).name = exp(iExp).name;
    results(iExp).varsR = string(exp(iExp).vtdR.Properties.VariableNames)';
    results(iExp).hasTimestampsR = ismember('Timestamp', exp(iExp).vtdR.Properties.VariableNames);

    if ~isempty(exp(iExp).vtdL)
        results(iExp).varsL = string(exp(iExp).vtdL.Properties.VariableNames)';
        results(iExp).hasTimestampsL = ismember('Timestamp', exp(iExp).vtdL.Properties.VariableNames);
    else
        results(iExp).hasTimestampsL = false;
    end

    results(iExp).isValid = length(results(iExp).varsL) == 23 && length(results(iExp).varsR) == 23;
end
clear iExp

% Cull bad experiments (some are still pending DLC)
assert(isequal({exp.name}, {results.name}));
exp = exp([results.isValid]);
results = results([results.isValid]);
eu = [exp.eu];

% Rename variables
from = ["HandIpsiCam", "HandContraCam", "LegIpsiCam", "LegContraCam", "Jaw", "Tongue", "Spine"];
toL = ["HandL", "HandR", "FootL", "FootR", "Jaw", "Tongue", "Spine"];
toR = ["HandR", "HandL", "FootR", "FootL", "Jaw", "Tongue", "Spine"];
for iExp = 1:length(exp)
    if results(iExp).hasTimestampsL
        vtd = exp(iExp).vtdL;
        for i = 1:length(toL)
            if ismember(sprintf("%s_X", from(i)), results(iExp).varsL)
                vtd = renamevars(vtd, ...
                    [sprintf("%s_X", from(i)), sprintf("%s_Y", from(i)), sprintf("%s_Likelihood", from(i))], ...
                    [sprintf("%s_X", toL(i)), sprintf("%s_Y", toL(i)), sprintf("%s_Likelihood", toL(i))]);
            end
        end
        exp(iExp).vtdL = vtd;
    end
    if results(iExp).hasTimestampsR
        vtd = exp(iExp).vtdR;
        for i = 1:length(toR)
            if ismember(sprintf("%s_X", from(i)), results(iExp).varsR)
                vtd = renamevars(vtd, ...
                    [sprintf("%s_X", from(i)), sprintf("%s_Y", from(i)), sprintf("%s_Likelihood", from(i))], ...
                    [sprintf("%s_X", toR(i)), sprintf("%s_Y", toR(i)), sprintf("%s_Likelihood", toR(i))]);
            end
        end
        exp(iExp).vtdR = vtd;
    end
end

% Reset results
clearvars -except eu exp results xta p kinematics ROOTPATH

expIndices = cellfun(@(name) find(strcmpi(name, {exp.name}), 1, 'first'), {eu.ExpName});

%% Load data
% Load dip/rise triggered averages, the bootstraps contained within are kind of useless.
load(fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot\LickVsReach_DLC_dta_rta_25_100_200to800ms_units1to1225_0boots_20260613.mat"));

% Load bootstrapped per-cluster averages. Can skip next step unless you
% want to recluster/rebootstrap
% load(fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot\LickVsReach_DLC_miBoot_200to800ms_1225units_10000boots_20260614.mat")); % contains updated `p`

%% Make Press/Lick trials if they don't exist
clear hasTrials;
hasTrials.Press = arrayfun(@(eu) ~isempty(eu.Trials.Press), eu);
hasTrials.Lick = arrayfun(@(eu) ~isempty(eu.Trials.Lick), eu);

%%
euSub = eu(~hasTrials.Press);
    
tEu = euSub.alignTimestamps(["TIMEOUT_START", "WAITFORTOUCH", "LEVER_PRESSED", "LEVER_RELEASED", "LEVER_HELD", "LICK", "LICK_OFF", "REWARD_ON", "REWARD_OFF"], ...
    acRefEventName=["REWARD_ON"], euRefEventName=["RewardOn"], ...
    trialDurationTolerance=1);

%% Make trials
for iSession = 1:length(tEu)
    if isempty(tEu(iSession).TIMEOUT_START)
        continue
    end
    iEu = tEu(iSession).euIndices(1);
    pressTrials = Trial(euSub(iEu).EventTimes.TIMEOUT_START, euSub(iEu).EventTimes.Press, stopMode='first', exclude=euSub(iEu).EventTimes.Lick);
    lickTrials = Trial(euSub(iEu).EventTimes.TIMEOUT_START, euSub(iEu).EventTimes.Lick, stopMode='first', exclude=euSub(iEu).EventTimes.Press);
    for iEu = tEu(iSession).euIndices(:)'
        euSub(iEu).Trials.Press = pressTrials;
        euSub(iEu).Trials.Lick = lickTrials;
    end
end


%% Get trial aligned rasters
trialTypes = ["press", "lick"];
trialTypeDisplayNames = ["reach", "lick"];
trialTypeRefEventDisplayNames = ["bar-contact", "spout-contact"];
clear rd
for trialType = trialTypes
    rd.(trialType).dip(length(eu)) = struct(name='', trialType='', alignTo='', t=[], I=[], spikeIndex=[], duration=[], iti=[]);
    rd.(trialType).rise(length(eu)) = struct(name='', trialType='', alignTo='', t=[], I=[], spikeIndex=[], duration=[], iti=[]);
    for iEu = 1:length(eu)
        for dir = ["dip", "rise"]
            rd.(trialType).(dir)(iEu) = eu(iEu).getRasterData(char(trialType), window=[-4, 1], minTrialDuration=1, alignTo='stop', spikeTimes=xta.(dir)(iEu).t0);
        end
    end
    rd.(trialType).spike = eu.getRasterData(char(trialType), window=[-4, 1], minTrialDuration=2, alignTo='stop');
end

%% Save rasterdata (rd)
metaRasterData = rd;
save(fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot", sprintf("%s_metaRasterData.mat", datetime("now", Format="uuuuMMdd"))), 'metaRasterData');


%% Plot them
close all
exportPath = fullfile(ROOTPATH, "LickVsReach_DTA_RTA_boot\Figures", sprintf("%s_LickVsReach_DLC_dta_rta_%i_%i_%ito%ims_std", datetime("now", Format="uuuuMMdd"), 100*p.xta.dip.thresholdQuantile, 100*p.xta.dip.thresholdSubQuantile, p.spikeRes*1000*p.xta.dip.samples(1), p.spikeRes*1000*p.xta.rise.samples(2)), "Reach aligned dip raster");
if ~exist(exportPath, 'dir')
    mkdir(exportPath)
else
    rmdir(exportPath, 's')
    mkdir(exportPath)   
end

features = ["HandL", "HandR", "Jaw"];
xtaWindow = [-0.5, 0.5];
tLocal = xtaWindow(1):p.xta.res:xtaWindow(2);
scale = 1;

fig = figure(Units='inches', Position=[1, 1, 10, 6]);
tlp = tiledlayout(fig, 1, length(trialTypes), TileSpacing='compact', Padding='compact');
tl = gobjects(1, length(trialTypes));
ax = gobjects(3, length(trialTypes));
for iTrialType = 1:length(trialTypes)
    tl(iTrialType) = tiledlayout(tlp, 3, 1, TileSpacing='compact', Padding='compact');
    tl(iTrialType).Layout.Tile = iTrialType;
    for iAx = 1:3
        ax(iAx, iTrialType) = nexttile(tl(iTrialType));
    end
end

for iEu = 1:length(eu)
    % try
        colors = 'kbr';
        for iTrialType = 1:length(trialTypes)
            trialType = trialTypes(iTrialType);
            % title(tl(iTrialType), trialTypeDisplayNames(iTrialType), FontSize=10, FontWeight='bold');
            iAx = 0;
            for dir = ["spike", "dip", "rise"]
                iAx = iAx + 1;
                cla(ax(iAx, iTrialType));
    
                if isempty(rd.(trialType).(dir)(iEu).t)
                    continue
                end
    
                % Raster
                EphysUnit.plotRaster(ax(iAx, iTrialType), rd.(trialType).(dir)(iEu), xlim=[-4, 1], iti=true, onlyPlotSpikes=true, sz=1);
                legend(ax(iAx, iTrialType), 'off')
                hold(ax(iAx, iTrialType), 'on')
    
                % Dip-triggered movement trace
                if ismember(dir, ["dip", "rise"])
                    h = gobjects(3, 1);
                    for iDip = 1:length(rd.(trialType).(dir)(iEu).t)
                        iFeat = 0;
                        t0 = rd.(trialType).(dir)(iEu).t(iDip);
                        I = rd.(trialType).(dir)(iEu).I(iDip);
                        dur = p.spikeRes*single(xta.(dir)(iEu).duration(iDip));
                        plot(ax(iAx, iTrialType), [t0, t0], [I-0.5, I+0.5], 'k-')
                        c = colors(iAx);
                        plot(ax(1, iTrialType), t0, I, Color=c, Marker='o', MarkerSize=7.5)
                        plot(ax(1, iTrialType), [t0, t0+dur], [I, I], Color=c, Marker='o', MarkerSize=3.75, LineStyle='none')
                        % if dur > 0.201
                        %     text(ax(1), t0 + 0.5*dur, I-0.5, sprintf("%i", 1000*dur), Color=c, FontSize=6, HorizontalAlignment='center', VerticalAlignment='middle')
                        % end
                        for fn = features
                            iFeat = iFeat + 1;
                            c = getColor(iFeat, length(features), 0.7);
                            selT = isin(xta.(dir)(iEu).(fn).t, xtaWindow);
                            t = t0 + tLocal;
                            x = xta.(dir)(iEu).(fn).X(rd.(trialType).(dir)(iEu).spikeIndex(iDip), selT);
                            h(iFeat) = plot(ax(iAx, iTrialType), t, x.*scale + I,  LineWidth=1, Color=c, DisplayName=fn);
                        end
                    end
                    legend(h, Location='east', AutoUpdate=false);
                    yticks(ax(iAx, iTrialType), ax(1, iTrialType).YTick)
                end
    
                % Axes
                xline(ax(iAx, iTrialType), 0, '--', Color=[0.15, 0.15, 0.15, 0.5])
    
                % hold(ax(iAx), 'off')
                title(ax(iAx, iTrialType), sprintf("%ss aligned to %s", dir, trialTypeDisplayNames(iTrialType)))
            end
    
            xlabel(tl(iTrialType), sprintf("time from %s (s)", trialTypeRefEventDisplayNames(iTrialType)), FontSize=9)
            xlabel(ax, '')
        end

        title(tlp, sprintf("unit %i - %s", iEu, eu(iEu).getName()), Interpreter='none', FontWeight='bold');
        fontsize(fig, 9, 'points')
        print(fig, fullfile(exportPath, sprintf("reachAlignedDipRaster_unit_%03i", iEu)), '-dpng', '-r0')
    % catch
    %     warning("Could not plot unit %i", iEu)
    % end
end