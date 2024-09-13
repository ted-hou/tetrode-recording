%% Load data
clear EMG
[EMG.exp, EMG.onset, EMG.isGoodTrial, EMG.eta, EMG.onsetThreshold, EMG.spuriousMovementThreshold] = read_emg();

%% Categorize trials based on EMG
spuriousMovementThreshold = 5;

normX = arrayfun(@(exp) exp.emg.touchAligned.normX, EMG.exp, UniformOutput=false);
normX = cat(1, normX{:});
t = EMG.exp(1).emg.touchAligned.t;

X = arrayfun(@(exp) exp.emg.touchAligned.X, EMG.exp, UniformOutput=false);
X = cat(1, X{:});

nTrials = size(normX, 1);
tThreshold = -0.5;
hMaxPre = max(normX(:, t <= tThreshold), [], 2, 'omitnan');
hMaxPost = max(normX(:, t > tThreshold), [], 2, 'omitnan');
isSimpleTrial = hMaxPre < spuriousMovementThreshold & hMaxPost >= spuriousMovementThreshold;
isComplexTrial = hMaxPre >= spuriousMovementThreshold & hMaxPost >= spuriousMovementThreshold;
isFlatTrial = hMaxPost < spuriousMovementThreshold;


%% Calculate neural ETA based on trial type.
clear eta
eta(length(EMG.exp)) = struct(simple=[], complex=[], flat=[]);
EMG.eta = eta;
EMG.eta(length(EMG.exp)) = struct(simple=[], complex=[], flat=[]);
csnTrials = cumsum([0, arrayfun(@(exp) length(exp.goodTrials), EMG.exp)]);
for iExp = 1:length(EMG.exp)
    selTrialInSession = csnTrials(iExp) + 1:csnTrials(iExp + 1);
    EMG.eta(iExp).simple = EMG.exp(iExp).eu.getETA('count', 'press_spontaneous', window=[-4, 0.5], resolution=0.1, normalize=[-4, -2], alignTo='stop', includeInvalid=true, ...
                trials=EMG.exp(iExp).goodTrials(isSimpleTrial(selTrialInSession)), correction=EMG.onset.emg(isSimpleTrial(selTrialInSession)));
    EMG.eta(iExp).complex = EMG.exp(iExp).eu.getETA('count', 'press_spontaneous', window=[-4, 0.5], resolution=0.1, normalize=[-4, -2], alignTo='stop', includeInvalid=true, ...
                trials=EMG.exp(iExp).goodTrials(isComplexTrial(selTrialInSession)), correction=EMG.onset.emg(isComplexTrial(selTrialInSession)));
    EMG.eta(iExp).flat = EMG.exp(iExp).eu.getETA('count', 'press_spontaneous', window=[-4, 0.5], resolution=0.1, normalize=[-4, -2], alignTo='stop', includeInvalid=true, ...
                trials=EMG.exp(iExp).goodTrials(isFlatTrial(selTrialInSession)), correction=EMG.onset.emg(isFlatTrial(selTrialInSession)));
end

etaSimple = [EMG.eta.simple];
etaComplex = [EMG.eta.complex];
etaFlat = [EMG.eta.flat];
etaMerged.simple = struct(X=cat(1, etaSimple.X), t=etaSimple(1).t, N=cat(1, etaSimple.N), D=cat(1, etaSimple.D));
etaMerged.complex = struct(X=cat(1, etaComplex.X), t=etaComplex(1).t, N=cat(1, etaComplex.N), D=cat(1, etaComplex.D));
etaMerged.flat = struct(X=cat(1, etaFlat.X), t=etaFlat(1).t, N=cat(1, etaFlat.N), D=cat(1, etaFlat.D));
clear etaSimple etaComplex etaFlat csnTrials iExp selTrialInSession


%% Figure S8 (EMG)
p.fontSize=9;
clear layout
layout.w = 8;
layout.h = 6;
layout.left.w = 3;
layout.mid.w = 3;
layout.right.w = 3;

p.lineWidth = 1.5;

close all
fig = figure(Units='inches', Position=[0, 0, layout.w, layout.h], DefaultAxesFontSize=p.fontSize);
layout.tl = tiledlayout(fig, 1, layout.left.w + layout.mid.w + layout.right.w, TileSpacing='loose');

layout.left.tl = tiledlayout(layout.tl, length(isSimpleTrial), 1, TileSpacing='compact', Padding='compact');
l = layout.left.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [1, layout.left.w];

layout.mid.tl = tiledlayout(layout.tl, 4, 1, TileSpacing='compact');
l = layout.mid.tl; l.Layout.Tile = 1 + layout.left.w; l.Layout.TileSpan = [1, layout.mid.w];

layout.right.tl = tiledlayout(layout.tl, 3, 1, TileSpacing='compact', Padding='compact');
l = layout.right.tl; l.Layout.Tile = 1 + layout.left.w + layout.mid.w; l.Layout.TileSpan = [1, layout.right.w];


% S8a. HeaEMG signal aligned to touch (good vs. bad)
ax = gobjects(3, 1);
ax(1) = nexttile(layout.left.tl, 1, [nnz(isSimpleTrial), 1]);
EphysUnit.plotETA(ax(1), struct(X=normX, t=t, N=ones(size(X, 1))), isSimpleTrial, XLim=[-4, 0], CLim=[0, 10], sortThreshold=0.25, signWindow=[-0.5, 0], sortWindow=[-2, 0]);
title(ax(1), sprintf('Simple trials (n=%i)', nnz(isSimpleTrial)))

ax(2) = nexttile(layout.left.tl, 1 + nnz(isSimpleTrial), [nnz(isComplexTrial), 1]);
EphysUnit.plotETA(ax(2), struct(X=normX, t=t, N=ones(size(X, 1))), isComplexTrial, XLim=[-4, 0], CLim=[0, 10], sortThreshold=0.25, signWindow=[-0.5, 0], sortWindow=[-2, 0]);
title(ax(2), sprintf('Complex trials (n=%i)', nnz(isComplexTrial)))

ax(3) = nexttile(layout.left.tl, 1 + nnz(isSimpleTrial) + nnz(isComplexTrial), [nnz(isFlatTrial), 1]);
EphysUnit.plotETA(ax(3), struct(X=normX, t=t, N=ones(size(X, 1))), isFlatTrial, XLim=[-4, 0], CLim=[0, 10], sortThreshold=0.25, signWindow=[-0.5, 0], sortWindow=[-2, 0]);
title(ax(3), sprintf('Flat trials (n=%i)', nnz(isFlatTrial)))

xlabel(ax, ''), ylabel(ax, '')
for i = 1:3
    colormap(ax(i), 'hot')
end

xlabel(layout.left.tl, 'Time to bar contact (s)', FontSize=p.fontSize)
ylabel(layout.left.tl, 'Trial', FontSize=p.fontSize)
colorbar(ax(1), 'off')
colorbar(ax(2), 'off')
hCb = colorbar(ax(3), Location='layout');
hCb.Layout.Tile = 'east';
hCb.Label.String = 'Normalized EMG (a.u.)';
% hCb.Label.FontSize = p.fontSize;

% S8b. Trial-averaged EMG signal alinged to touch
ax = nexttile(layout.mid.tl);
hold(ax, 'on')
plot(ax, t, mean(normX(isSimpleTrial, :), 1, 'omitnan'), 'r', Display=sprintf('Simple trials (n=%i)', nnz(isSimpleTrial)));
plot(ax, t, mean(normX(isComplexTrial, :), 1, 'omitnan'), 'b', Display=sprintf('Complex trials (n=%i)', nnz(isComplexTrial)));
plot(ax, t, mean(normX(isFlatTrial, :), 1, 'omitnan'), 'k', Display=sprintf('Flat trials (n=%i)', nnz(isFlatTrial)));
hold(ax, 'off')
ylim(ax, [-5, 15])
xlabel(ax, 'Time to bar contact (s)')
ylabel(ax, 'Normalized EMG (a.u.)')
legend(ax, Location='northoutside')

% S8c. Histogram of video detected onset times
edges = -4:1/30:0;
ax(1) = nexttile(layout.mid.tl);
histogram(ax(1), EMG.onset.video, edges, FaceColor='black');
title(ax(1), 'Video onset')
xlabel(ax(1), 'Time to bar contact (s)')
ylabel(ax(1), 'No. trials')

% S8d. Histogram of EMG detected onset times
ax(2) = nexttile(layout.mid.tl);
histogram(ax(2), EMG.onset.emg, edges, FaceColor='black');
title(ax(2), 'EMG onset')
xlabel(ax(2), 'Time to bar contact (s)')
ylabel(ax(2), 'No. trials')

xlim(ax(1:2), [-4, 0])

% S8f. Heatmap of neural activity, divided by trial type
ax = gobjects(3, 1);
ax(1) = nexttile(layout.right.tl);
[~, order, ~, latencySimple] = EphysUnit.plotETA(ax(1), etaMerged.simple, Clim=[-2, 2], SortWindow=[-2, 0], SignWindow=[-0.1, 0], sortThreshold=0.25, XLim=[-2, 0]);
title(ax(1), 'Simple trials')

ax(2) = nexttile(layout.right.tl);
EphysUnit.plotETA(ax(2), etaMerged.complex, Clim=[-2, 2], SortWindow=[-2, 0], sortThreshold=0.25, SignWindow=[-0.1, 0], order=order, XLim=[-2, 0]);
axTemp = axes(figure);
[~, order, ~, latencyComplex] = EphysUnit.plotETA(axTemp, etaMerged.complex, Clim=[-2, 2], SortWindow=[-2, 0], sortThreshold=0.25, SignWindow=[-0.1, 0], XLim=[-2, 0]);
close(axTemp.Parent)
title(ax(2), 'Complex trials')


ax(3) = nexttile(layout.right.tl);
EphysUnit.plotETA(ax(3), etaMerged.flat, Clim=[-2, 2], SortWindow=[-2, 0], sortThreshold=0.25, SignWindow=[-0.1, 0], order=order, XLim=[-2, 0]);
axTemp = axes(figure);
[~, order, ~, latencyFlat] = EphysUnit.plotETA(axTemp, etaMerged.flat, Clim=[-2, 2], SortWindow=[-2, 0], sortThreshold=0.25, SignWindow=[-0.1, 0], XLim=[-2, 0]);
close(axTemp.Parent)
title(ax(3), 'Flat trials')


xlabel(ax, ''), ylabel(ax, '')

xlabel(layout.right.tl, 'Time to EMG onset (s)', FontSize=p.fontSize)
ylabel(layout.right.tl, 'Unit', FontSize=p.fontSize)
colorbar(ax(1), 'off')
colorbar(ax(2), 'off')
hCb = colorbar(ax(3), Location='layout');
hCb.Layout.Tile = 'east';
hCb.Label.String = 'Normalized spike rate (a.u.)';
hCb.Label.FontSize = p.fontSize;

% S8e. Histogram of neural onset times
ax(3) = nexttile(layout.mid.tl);
hold(ax(3), 'on')
histogram(ax(3), latencySimple, -2:2/30:0, FaceColor='black');
% histogram(ax(3), latencyComplex, -2:2/30:0, FaceColor='red', FaceAlpha=0.1, EdgeColor='red', DisplayName='complex');
title(ax(3), 'Neural onset (simple trials)')
xlabel(ax(3), 'Time to EMG onset (s)')
xlim(ax(3), [-2, 0])
ylabel(ax(3), 'No. units')
% legend(ax(3), Location='north', Orientation='vertical')

copygraphics(fig, ContentType='vector')
