%% Load data
% load_ephysUnits;
read_reachDir_2tgt;
read_reachDir_4tgt;

%% Perm test SNr response AMPLITUDE for press vs lick
sel = c.hasLick & c.hasPress;
[boot.pressVsLick.h, boot.pressVsLick.p, boot.pressVsLick.ci, boot.pressVsLick.obs] = bootstrapAmplitude(eu(sel), ...
    'press', 'lick', responseWindowA=p.metaWindowPress, responseWindowB=p.metaWindowLick, baselineWindow=[-4, -2], ...
    allowedTrialDuration=[2, Inf], withReplacement=false, alpha=0.01);
c.isSelective.pressVsLick = false(1, length(eu));
c.isSelective.pressVsLick(sel) = boot.pressVsLick.h;

%% Perm test SNr response AMPLITUDE: 4tgt, contra vs. ipsi paw
N = arrayfun(@(eta) eta.N, trajCombined.eta, 'UniformOutput', false);
assert(minNumTrials == 4)
sel = N{1, 1} >= minNumTrials & N{2, 1} >= minNumTrials & N{3, 1} >= minNumTrials & N{4, 3} >= minNumTrials;


p.metaWindowReachDir = [-0.2, 0.1];
[boot.contraFrontVsIpsiFront4tgt.h, boot.contraFrontVsIpsiFront4tgt.p, boot.contraFrontVsIpsiFront4tgt.ci, boot.contraFrontVsIpsiFront4tgt.obs] = bootstrapAmplitude(euReachDir4Tgt(sel), ...
    trajCombined.trials{2, 1}(sel), trajCombined.trials{4, 3}(sel), correctionA=trajCombined.correction{2, 1}(sel), correctionB=trajCombined.correction{4, 3}(sel), ...
    responseWindowA=p.metaWindowReachDir, responseWindowB=p.metaWindowReachDir, baselineWindow=[-4, -2], ...
    allowedTrialDuration=[-Inf, Inf], withReplacement=false, alpha=0.01);
c.isSelective.contraFrontVsIpsiFront4tgt = false(1, length(euReachDir4Tgt));
c.isSelective.contraFrontVsIpsiFront4tgt(sel) = boot.contraFrontVsIpsiFront4tgt.h;

% Perm test SNr response AMPLITUDE: 4tgt, contra paw, out vs. in
N = arrayfun(@(eta) eta.N, trajCombined.eta, 'UniformOutput', false);
assert(minNumTrials == 4)
sel = N{1, 1} >= minNumTrials & N{2, 1} >= minNumTrials & N{3, 1} >= minNumTrials & N{4, 3} >= minNumTrials;

p.metaWindowReachDir = [-0.2, 0.1];
[boot.contraOutVsContraIn4tgt.h, boot.contraOutVsContraIn4tgt.p, boot.contraOutVsContraIn4tgt.ci, boot.contraOutVsContraIn4tgt.obs] = bootstrapAmplitude(euReachDir4Tgt(sel), ...
    trajCombined.trials{1, 1}(sel), trajCombined.trials{3, 1}(sel), correctionA=trajCombined.correction{1, 1}(sel), correctionB=trajCombined.correction{3, 1}(sel), ...
    responseWindowA=p.metaWindowReachDir, responseWindowB=p.metaWindowReachDir, baselineWindow=[-4, -2], ...
    allowedTrialDuration=[-Inf, Inf], withReplacement=false, alpha=0.01);
c.isSelective.contraOutVsContraIn4tgt = false(1, length(euReachDir4Tgt));
c.isSelective.contraOutVsContraIn4tgt(sel) = boot.contraOutVsContraIn4tgt.h;

% Perm test SNr response AMPLITUDE: 2tgt, contra paw, out vs. in
p.metaWindowReachDir = [-0.2, 0.1];
[boot.contraOutVsContraIn2tgt.h, boot.contraOutVsContraIn2tgt.p, boot.contraOutVsContraIn2tgt.ci, boot.contraOutVsContraIn2tgt.obs] = bootstrapAmplitude(euReachDir2Tgt, ...
    trajCombined2tgt.trials{1}, trajCombined2tgt.trials{2}, correctionA=trajCombined2tgt.correction{1}, correctionB=trajCombined2tgt.correction{2}, ...
    responseWindowA=p.metaWindowReachDir, responseWindowB=p.metaWindowReachDir, baselineWindow=[-4, -2], ...
    allowedTrialDuration=[-Inf, Inf], withReplacement=false, alpha=0.01);
c.isSelective.contraOutVsContraIn2tgt = boot.contraOutVsContraIn2tgt.h;
