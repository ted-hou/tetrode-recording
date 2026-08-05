eu = EphysUnit.load('C:\SERVER\Units\SNr_CoChR_VGATCre\SingleUnit_NonDuplicate_NonDrift_SNr', waveforms=false, spikecounts=false, spikerates=false);


nTrials = arrayfun(@(eu) length(eu.EventTimes.TimeoutOn), eu);
eu(nTrials==0) = [];

% [ac.ac, ac.euIndicesFirstInSession, ac.expIndices, ac.uniqueExpNames] = eu.loadArduinoConnection();
% eu.alignTimestamps(["WAITFORTOUCH", "TIMEOUT_START"], ac={ac.ac, ac.euIndicesFirstInSession, ac.expIndices, ac.uniqueExpNames}, acRefEventName="TIMEOUT_END", euRefEventName="TimeoutOff");

% Make press/lick trials

for iEu = 1:length(eu)
    selPress = eu(iEu).EventTimes.PressOff - eu(iEu).EventTimes.PressOn > 0.002;
    eu(iEu).Trials.Press = Trial(eu(iEu).EventTimes.TimeoutOn, eu(iEu).EventTimes.Press(selPress), stopMode='first', exclude=eu(iEu).EventTimes.Lick); 
    eu(iEu).Trials.Lick = Trial(eu(iEu).EventTimes.TimeoutOn, eu(iEu).EventTimes.Lick, stopMode='first', exclude=eu(iEu).EventTimes.Press(selPress));
end



%% Boot movement responses
c.hasPress = false(1, length(eu));
c.hasLick = false(1, length(eu));
for iEu = 1:length(eu)
    c.hasPress(iEu) = length(eu(iEu).Trials.Press) > 10;
    c.hasLick(iEu) = length(eu(iEu).Trials.Lick) > 10;
end
p.bootAlpha = 0.01;
p.nboot = 100000;
p.responseWindowPress = [-0.3, 0];
p.responseWindowLick = [-0.3, 0];
assert(isequal(p.responseWindowPress, [-0.3, 0]))
assert(isequal(p.responseWindowLick, [-0.3, 0]))
boot.press = struct('h', NaN(length(eu), 1), 'muDiffCI', NaN(length(eu), 2), 'muDiffObs', NaN(length(eu), 1));
boot.lick = struct('h', NaN(length(eu), 1), 'muDiffCI', NaN(length(eu), 2), 'muDiffObs', NaN(length(eu), 1));
[boot.press.h, boot.press.muDiffCI, boot.press.muDiffObs] = bootstrapMoveResponse( ...
    eu, 'press', nboot=p.nboot, alpha=p.bootAlpha, withReplacement=false, oneSided=false, ...
    responseWindow=p.responseWindowPress);
[boot.lick.h, boot.lick.muDiffCI, boot.lick.muDiffObs] = bootstrapMoveResponse( ...
    eu, 'lick', nboot=p.nboot, alpha=p.bootAlpha, withReplacement=false, oneSided=false, ...
    responseWindow=p.responseWindowLick);
fprintf(1, '\nAll done\n')

%% Report bootstraped movement response direction


% figure, histogram(boot.press.h)
c.isPressUp = boot.press.h' == 1 & c.hasPress;
c.isPressDown = boot.press.h' == -1 & c.hasPress;
c.isPressResponsive = (c.isPressUp | c.isPressDown);

% figure, histogram(boot.lick.h)
c.isLickUp = boot.lick.h' == 1 & c.hasLick;
c.isLickDown = boot.lick.h' == -1 & c.hasLick;
c.isLickResponsive = (c.isLickUp | c.isLickDown);

fprintf(1, ['%g good units, %g modulated for reach:\n' ...
    '\t%g (%.0f%%) are excited (p<%g);\n' ...
    '\t%g (%.0f%%) are inhibited (p<%g).\n'], ...
    nnz(c.hasPress), nnz(c.isPressResponsive), ...
    nnz(c.isPressUp), 100*nnz(c.isPressUp)/nnz(c.isPressResponsive), p.bootAlpha, ...
    nnz(c.isPressDown), 100*nnz(c.isPressDown)/nnz(c.isPressResponsive), p.bootAlpha);

fprintf(1, ['%g good units, %g modulated for lick:\n' ...
    '\t%g (%.0f%%) are excited (p<%g);\n' ...
    '\t%g (%.0f%%) are inhibited (p<%g).\n'], ...
    nnz(c.hasLick), nnz(c.isLickResponsive), ...
    nnz(c.isLickUp), 100*nnz(c.isLickUp)/nnz(c.isLickResponsive), p.bootAlpha, ...
    nnz(c.isLickDown), 100*nnz(c.isLickDown)/nnz(c.isLickResponsive), p.bootAlpha);

nTotal = nnz(c.isPressResponsive & c.isLickResponsive);
fprintf(1, ['%g good units, %g modulated for reach AND lick:\n' ...
    '\t%g (%.0f%%) are press-excited AND lick-excited;\n' ...
    '\t%g (%.0f%%) are press-inhibited AND lick-inhibited;\n' ...
    '\t%g (%.0f%%) are press-excited AND lick-inhibited;\n' ...
    '\t%g (%.0f%%) are press-inhibited AND lick-excited;\n'], ...
    nnz(c.hasPress & c.hasLick), nTotal, ...
    nnz(c.isPressUp & c.isLickUp), 100*nnz(c.isPressUp & c.isLickUp)/nTotal, ...
    nnz(c.isPressDown & c.isLickDown), 100*nnz(c.isPressDown & c.isLickDown)/nTotal, ...
    nnz(c.isPressUp & c.isLickDown), 100*nnz(c.isPressUp & c.isLickDown)/nTotal, ...
    nnz(c.isPressDown & c.isLickUp), 100*nnz(c.isPressDown & c.isLickUp)/nTotal)   

sel = c.hasPress & c.hasLick;
fprintf('05 Calculate: Of %i: %i (%i%%) showed modulation for BOTH, %i (%i%%) showed modulation for lick only, %i (%i%%) showed modulation for reach only, %i (%i%%) for neither.', ...
    nnz(sel), ...
    nnz(sel & c.isPressResponsive & c.isLickResponsive), round(nnz(sel & c.isPressResponsive & c.isLickResponsive)/nnz(sel)*100), ...
    nnz(sel & c.isLickResponsive & ~c.isPressResponsive), round(nnz(sel & c.isLickResponsive & ~c.isPressResponsive)/nnz(sel)*100), ...
    nnz(sel & c.isPressResponsive & ~c.isLickResponsive), round(nnz(sel & c.isPressResponsive & ~c.isLickResponsive)/nnz(sel)*100), ...
    nnz(sel & ~c.isPressResponsive & ~c.isLickResponsive), round(nnz(sel & ~c.isPressResponsive & ~c.isLickResponsive)/nnz(sel)*100) ...
    )

clear nTotal sel

save('C:\SERVER\Units\boot_SNr_CoChR_VGATCre.mat', 'boot', 'c', 'eta')
eu.save('C:\SERVER\Units\SNr_CoChR_VGATCre\SingleUnit_NonDuplicate_NonDrift_SNr')

%%
rd = eu.getRasterData('stim', window=[-1, 1], photoelectricBlankDuration=0.5e-3, minTrialDuration=0.75, maxTrialDuration=4);
if ~exist('E:\Figures\SNr_CoChR_VGATCre', 'dir')
    mkdir('E:\Figures\SNr_CoChR_VGATCre')
end

fig = figure(Units='inches', Position=[1, 1, 7, 7]);
ax = axes(fig);
for iEu = 1:length(eu)
    cla(ax)
    try
        EphysUnit.plotRaster(ax, rd(iEu), timeUnit='ms', xlim=[-1000, 4000]);
        legend(ax, 'off')
        print(fig, fullfile('E:\Figures\SNr_CoChR_VGATCre', sprintf('%s.png', eu(iEu).getName())), '-dpng', '-r0')
    end
end


%% Plot PSTH aligned to stim onset
clear artifacts
artifacts(1) = struct(event='StimOn', length=0.5, lengthUnit='ms', direction='right');
artifacts(2) = struct(event='StimOff', length=0.5, lengthUnit='ms', direction='right');
eta.stim = eu.getETA('count', 'stim', [-1, 1], resolution=0.020, alignTo='start', normalize=[-0.5, 0], artifacts=artifacts);
meta.stim = mean(eta.stim.X(:, isin(eta.stim.t, [0.05, 0.2])), 2, 'omitnan');
ax = axes(figure); hold(ax, 'on')
plot(ax, eta.stim.t, eta.stim.X(meta.stim>=0, :), Color=[0.8, 0.2, 0.2, 0.1])
plot(ax, eta.stim.t, eta.stim.X(meta.stim<0, :), Color=[0.2, 0.2, 0.8, 0.1])
plot(ax, eta.stim.t, mean(eta.stim.X, 1, 'omitnan'), Color='k', LineWidth=2)

%% Plot PETH aligned to reach
close all
eta.press = eu.getETA('count', 'press', [-4, 1], resolution=0.2, alignTo='start', normalize=[-4, -2], minTrialDuration=1);
eta.lick = eu.getETA('count', 'lick', [-4, 1], resolution=0.2, alignTo='start', normalize=[-4, -2], minTrialDuration=1);
meta.press = mean(eta.press.X(:, isin(eta.press.t, [-0.4, 0.2])), 2, 'omitnan');
meta.lick = mean(eta.lick.X(:, isin(eta.lick.t, [-0.4, 0.2])), 2, 'omitnan');

%%
ax = axes(figure); hold(ax, 'on')
plot(ax, eta.stim.t, eta.stim.X(c.isPressUp, :), Color=[0.8, 0.2, 0.2, 0.1])
plot(ax, eta.stim.t, eta.stim.X(c.isPressDown, :), Color=[0.2, 0.2, 0.8, 0.1])
plot(ax, eta.stim.t, mean(eta.stim.X, 1, 'omitnan'), Color='k', LineWidth=2)
plot(ax, eta.stim.t, mean(eta.stim.X(c.isPressUp, :), 1, 'omitnan'), Color=[0.8, 0.2, 0.2], LineWidth=2)
plot(ax, eta.stim.t, mean(eta.stim.X(c.isPressDown, :), 1, 'omitnan'), Color=[0.2, 0.2, 0.8], LineWidth=2)

%% Plot press PETH
ax = axes(figure); hold(ax, 'on')
plot(ax, eta.press.t, eta.press.X(c.isPressUp, :), Color=[0.8, 0.2, 0.2, 0.1])
plot(ax, eta.press.t, eta.press.X(c.isPressDown, :), Color=[0.2, 0.2, 0.8, 0.1])
plot(ax, eta.press.t, mean(eta.press.X, 1, 'omitnan'), Color='k', LineWidth=2)
xline(ax, 0, 'k-')
title(ax, 'press')
%%
ax = axes(figure); hold(ax, 'on')
plot(ax, eta.press.t, eta.press.X(meta.stim>+0.50, :), Color=[0.8, 0.2, 0.2, 0.1])
plot(ax, eta.press.t, eta.press.X(meta.stim<-0.5, :), Color=[0.2, 0.2, 0.8, 0.1])
plot(ax, eta.press.t, mean(eta.press.X(meta.stim>0, :), 1, 'omitnan'), Color=[0.8, 0.2, 0.2], LineWidth=1)
plot(ax, eta.press.t, mean(eta.press.X(meta.stim<0, :), 1, 'omitnan'), Color=[0.2, 0.2, 0.8], LineWidth=1)
xline(ax, 0, 'k-')
title(ax, 'press')
%%
ax = axes(figure); hold(ax, 'on')
plot(ax, eta.lick.t, eta.lick.X(c.isLickUp, :), Color=[0.8, 0.2, 0.2, 0.1])
plot(ax, eta.lick.t, eta.lick.X(c.isLickDown, :), Color=[0.2, 0.2, 0.8, 0.1])
plot(ax, eta.lick.t, mean(eta.lick.X, 1, 'omitnan'), Color='k', LineWidth=2)
xline(ax, 0, 'k-')
title(ax, 'lick')

ax = axes(figure); hold(ax, 'on')
sel = c.isPressResponsive;
scatter(ax, meta.press(sel), meta.stim(sel))
xline(ax, 0, 'k--')
yline(ax, 0, 'k--')
mdl = fitlm(meta.press(sel), meta.stim(sel));
plot(ax, mdl)

xlabel(ax, 'press')
ylabel(ax, 'stim')


ax = axes(figure); hold(ax, 'on')
sel = c.isLickResponsive;
scatter(ax, meta.lick(sel), meta.stim(sel))
xline(ax, 0, 'k--')
yline(ax, 0, 'k--')
mdl = fitlm(meta.lick(sel), meta.stim(sel));
plot(ax, mdl)

xlabel(ax, 'lick')
ylabel(ax, 'stim')



ax = axes(figure); hold(ax, 'on')
scatter(ax, meta.press, meta.lick)
xline(ax, 0, 'k--')
yline(ax, 0, 'k--')
mdl = fitlm(meta.press, meta.lick);
plot(ax, mdl)

xlabel(ax, 'press')
ylabel(ax, 'lick')