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
meta.stim = mean(eta.stim.X(:, isin(eta.stim.t, [0, 0.1])), 2, 'omitnan');
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

% Plot press PETH
ax = axes(figure); hold(ax, 'on')
plot(ax, eta.press.t, eta.press.X(meta.press>=0.5, :), Color=[0.8, 0.2, 0.2, 0.1])
plot(ax, eta.press.t, eta.press.X(meta.press<-0.25, :), Color=[0.2, 0.2, 0.8, 0.1])
plot(ax, eta.press.t, mean(eta.press.X, 1, 'omitnan'), Color='k', LineWidth=2)
title(ax, 'press')

ax = axes(figure); hold(ax, 'on')
plot(ax, eta.lick.t, eta.lick.X(meta.lick>=0.5, :), Color=[0.8, 0.2, 0.2, 0.1])
plot(ax, eta.lick.t, eta.lick.X(meta.lick<-0.25, :), Color=[0.2, 0.2, 0.8, 0.1])
plot(ax, eta.lick.t, mean(eta.lick.X, 1, 'omitnan'), Color='k', LineWidth=2)
title(ax, 'lick')

ax = axes(figure); hold(ax, 'on')
sel = meta.press >= 0.5 | meta.press <= -0.25;
scatter(ax, meta.press(sel), meta.stim(sel))
xline(ax, 0, 'k--')
yline(ax, 0, 'k--')
mdl = fitlm(meta.press(sel), meta.stim(sel));
plot(ax, mdl)

xlabel(ax, 'press')
ylabel(ax, 'stim')


ax = axes(figure); hold(ax, 'on')
sel = meta.lick >= 0.5 | meta.lick <= -0.25;
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