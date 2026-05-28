
%% TODO:
% Individual units do not really always fit the grand cluster labels
% Try doing actual categorization, trial by trial, look for trials containing exactly ONE of the following:
% lick start, lick stop, lhand reach, lhand retract, rhand reach, rhand retract

% We'd probably need to do some sort of test to see if significant (do a
% bootstrap for each subcluster for each unit, only keep subcluster if it's a single movement)

if ~exist('expIndices', 'var')
    expIndices = [xta.dip.iExp];
end

p.mi.boot.features = ["Jaw", "HandL", "HandR", "Spine"];
p.mi.boot.clusters = 1:7;
p.mi.boot.alpha = 0.05;
p.mi.boot.nBoot = 1000;

features = p.mi.boot.features;
clusters = p.mi.boot.clusters;
alpha = p.mi.boot.alpha;
nBoot = p.mi.boot.nBoot;
windowPre = p.mi.windowPre;
windowPost = p.mi.windowPost;
assert(nBoot*alpha/length(features)/length(clusters) > 1, ...
    "nBoot=%i is not sufficient for alpha=%g + bon-ferroni correction for %i features x %i clusters.", nBoot, alpha, length(features), length(clusters))

rng(42)
if isempty(gcp('nocreate'))
    parpool('Processes');
end
tTic = tic();
ll = 0;
clear miBoot
nUnits = length(xta.dip);
miBoot(nUnits) = struct(dip=[], rise=[]);
for iUnit = 1:nUnits
    fprintf(repmat('\b', [1, ll]));
    ll = fprintf("Bootstrapping movement index by cluster... unit %i/%i... time elapsed (%.2fs)\n", iUnit, nUnits, toc(tTic));
    kine = kinematics(expIndices(iUnit));
    maxT = kine.HandR.t(end);
    for dir = ["dip", "rise"]
        miBootTemp = NaN(nBoot, length(clusters), length(features));
        for iClu = clusters
            nTrials = clusterSize(iUnit).(dir).n(iClu);
            if isnan(nTrials)
                continue
            end
            parfor iBoot = 1:nBoot
                xDiff = NaN(nTrials, length(features));
                for iFeat = 1:length(features)
                    fn = features(iFeat);
                    t = kine.(fn).t;
                    t0Boot = rand([nTrials, 1]) * maxT;
                    for iTrial = 1:nTrials
                        xPre = mean(kine.(fn).X(isin(t, t0Boot(iTrial) + windowPre)), 'all', 'omitnan');
                        xPost = mean(kine.(fn).X(isin(t, t0Boot(iTrial) + windowPost)), 'all', 'omitnan');
                        xDiff(iTrial, iFeat) = xPost - xPre;
                    end
                end
                % Avg across trials: Potentially problematic with NaNs, might need to average across trials then subtract;
                miBootTemp(iBoot, iClu, :) = mean(xDiff, 1, 'omitnan');
            end
        end
        miBoot(iUnit).(dir) = miBootTemp; % nBoot x clusters x features
    end
end

clear tTic ll iUnit kine maxT dir miBootTemp iClu nTrials iBoot xDiff iFeat fn t t0Boot iTrial xPre xPost xDiff

miObs(nUnits) = struct(dip=[], rise=[]);
for iUnit = 1:nUnits
    for dir = ["dip", "rise"]
        miObs(iUnit).(dir) = NaN(length(clusters), length(features));
        idx = mi.(dir)(iUnit).idx;
        for iClu = clusters
            nTrials = clusterSize(iUnit).(dir).n(iClu);
            if isnan(nTrials)
                continue
            end
            for iFeat = 1:length(features)
                fn = features(iFeat);
                miObs(iUnit).(dir)(iClu, iFeat) = mean(mi.(dir)(iUnit).(fn)(idx==iClu), 'all', 'omitnan');
            end
        end
    end
end
clear iUnit dir iFeat fn idx iClu

%% Save data
exportPath = fullfile("C:\SERVER\LickVsReach_DTA_RTA_boot", sprintf("LickVsReach_DLC_miBoot_%iunits_%iboots.mat", nUnits, p.mi.boot.nBoot));
save(exportPath, 'clusterSize', 'mi', 'miBoot', 'miObs', 'p', '-v7.3')
fprintf("Saved to %s\n", exportPath);

clear features clusters alpha nBoot windowPre windowPost

