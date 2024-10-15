%% Load data
load_ephysUnits;

%% Perm test SNr response amplitude differences between trial types
sel = c.hasLick & c.hasPress;
[boot.pressVsLick.h, boot.pressVsLick.p, boot.pressVsLick.ci, boot.pressVsLick.obs] = bootstrapAmplitude(eu(sel), ...
    'press', 'lick', responseWindowA=p.metaWindowPress, responseWindowB=p.metaWindowLick, ...
    allowedTrialDuration=[2, Inf], withReplacement=false, alpha=0.01);
c.isSelective.pressVsLick = false(1, length(eu));
c.isSelective.pressVsLick(sel) = boot.pressVsLick.h;
%%