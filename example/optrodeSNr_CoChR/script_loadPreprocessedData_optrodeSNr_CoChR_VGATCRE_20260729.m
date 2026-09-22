
%% Can skip all above, just do this to load
% Try running script_optrodeSNr_CoChR_VGATCre_20260729 again;
expDesc = "VGAT-Cre x SNr(AAV-flex-CoChR)";
eu = EphysUnit.load('C:\SERVER\Units\SNr_CoChR_VGATCre\SingleUnit_NonDuplicate_NonDrift_SNr_ValidVideos');

expNames = string({eu.ExpName}');
[uniqueExpNames, expToEuIndices, euToExpIndices] = unique(expNames);

load('C:\SERVER\Units\meta_SNr_CoChR_VGATCre_ValidVideos_withoutStimData.mat'); %'boot', 'c', 'eta', 'meta', 'kinematics', 'p'
load('C:\SERVER\Units\meta_SNr_CoChR_VGATCre_ValidVideos_justStimData.mat'); %'stimData'

% Load video tracking data
dlcResultsPath = 'C:\SERVER\DeepLabCut\Results\FourPawsJawTongueSpine_Emma';
exp = CompleteExperiment3(eu, cameras='lr', deeplabcutPath=dlcResultsPath);
exp.alignTimestamps(refEventNameArduino={'TIMEOUT_END'}, refEventNameEphys={'TimeoutOff'}, trialDurationTolerance=2);

clear results
results(length(exp)) = struct(name=[], varsL=[], hasTimestampsL=[], varsR=[], hasTimestampsR=[], isValid=[]);
for iExp = 1:length(exp)
    results(iExp).name = exp(iExp).name;
    results(iExp).varsR = string(exp(iExp).vtdR.Properties.VariableNames)';
    results(iExp).hasTimestampsR = ismember('Timestamp', exp(iExp).vtdR.Properties.VariableNames);

    results(iExp).varsL = string(exp(iExp).vtdL.Properties.VariableNames)';
    results(iExp).hasTimestampsL = ismember('Timestamp', exp(iExp).vtdL.Properties.VariableNames);

    results(iExp).isValid = length(results(iExp).varsL) == 23 && length(results(iExp).varsR) == 23;
end
clear iExp

% Rename variables
from = ["HandIpsiCam", "HandContraCam", "FootIpsiCam", "FootContraCam", "Jaw", "Tongue", "Spine"];
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
clear dlcResultsPath results from toL toR iExp vtd i

clear vs
vs = vertcat(who, 'vs');

% Load decoder data
load('C:\SERVER\Units\meta_SNr_CoChR_VGATCre_ValidVideos_decoderData.mat')
vs = who;