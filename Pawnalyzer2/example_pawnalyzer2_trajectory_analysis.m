
eu = EphysUnit.load('C:\SERVER\Units\acute_3cam_reach_direction_2tgts\SingleUnits_NonDuplicate');
for iEu = 1:length(eu)
    trials = eu(iEu).makeTrials('press_spontaneous2');
    goodTrials = trials(trials.duration() > 4);
    eu(iEu).Trials.PressSpontaneous = goodTrials;
end

% Group EphysUnits by session (euByExp)
pa = Pawnalyzer2(eu, refEvent='reward');
euByExp = arrayfun(@(x) x.eu, pa.exp, UniformOutput=false);
clear pa iEu trials goodTrials

pa = Pawnalyzer2(eu, refEvent='press');

%% Select good sessions
expNames = ["daisy23_20240626", "daisy23_20240627", "daisy24_20240626", "daisy25_20240701"];
selEu = ismember({eu.ExpName}, expNames);
eu = eu(selEu);
clearvars -except eu expNames

%% Make trajectories
pa = Pawnalyzer2(eu, refEvent='press');
pa.getClips(nFramesBefore=15, nFramesAfter=0, keepData=false, trials='PressSpontaneous');
pa.load('auto')

%% Analysis
%% Divide trajectories (nTargets x nSubclusters) and average into PETHs
close all
nDims = 8;
nSubClusters = 1;
nExp = pa.getLength('exp');
fig1 = figure();
ax1 = axes(fig1);
fig2 = figure(Units='normalized', OuterPosition=[0.1, 0.1, 0.7, 0.5]);
ax2 = gobjects(nExp, 1);
ax3 = gobjects(nExp, 1);
hold(ax1, 'on')
clear traj
traj(nExp) = struct(contra=[], ipsi=[], target=[], t=[], pca=[], kmeans=[]);
nt = 7;
for iExp = 2:nExp
    nFrames = pa.getLength('frame', exp=iExp, trial=1);
%     [traj(iExp).contra, traj(iExp).ipsi, traj(iExp).target, traj(iExp).t] = pa.getTrajectories(iExp, zero=nFrames - nt + 1);
    [traj(iExp).contra, traj(iExp).ipsi, traj(iExp).target, traj(iExp).t] = pa.getTrajectories(iExp);
    f = horzcat(traj(iExp).contra.x(:, end-nt+1:end), traj(iExp).contra.y(:, end-nt+1:end), traj(iExp).contra.z(:, end-nt+1:end));
    nTrials = size(f, 1);
    [traj(iExp).pca.coeff, traj(iExp).pca.score, ~, ~, traj(iExp).pca.explained, ~] = pca(f);
    traj(iExp).cluster.k = nSubClusters;
    targetNames = {'contra-in', 'contra-out'};
    [~, targetIndex] = ismember(traj(iExp).target, targetNames);
    targetIndex = targetIndex';
    traj(iExp).cluster.class = zeros(nTrials, 1);
    for iTarget = 1:2
        sel = targetIndex==iTarget;
%         traj(iExp).cluster.class(sel) = kmeans(traj(iExp).pca.score(sel, 1:nDims), nSubClusters, EmptyAction='drop');
        gm = fitgmdist(traj(iExp).pca.score(sel, 1:nDims), nSubClusters, RegularizationValue=0.01);
        traj(iExp).cluster.class(sel) = gm.cluster(traj(iExp).pca.score(sel, 1:nDims));
%         traj(iExp).cluster.class(sel) = kmeans(traj(iExp).pca.score(sel, 1:nDims), nSubClusters);
    end

%     plot(ax1, [0; cumsum(traj(iExp).pca.explained)], DisplayName=sprintf('iExp=%i', iExp));

%     ax2(iExp) = subplot(2, nExp, iExp);
%     hold(ax2(iExp), 'on')
%     colors = 'rgbm';
%     markers = ["o", "square", "diamond"];
%     h = gobjects(4, 1);
%     for iTarget = 1:2
%         h(iTarget) = scatter3(ax2(iExp), traj(iExp).pca.score(targetIndex==iTarget, 1), traj(iExp).pca.score(targetIndex==iTarget, 2), traj(iExp).pca.score(targetIndex==iTarget, 3), 15, colors(iTarget), 'filled', DisplayName=targetNames{iTarget});
%         for iSubCluster = 1:nSubClusters
%             sel = traj(iExp).cluster.class==iSubCluster & targetIndex==iTarget;
%             scatter3(ax2(iExp), traj(iExp).pca.score(sel, 1), traj(iExp).pca.score(sel, 2), traj(iExp).pca.score(sel, 3), 50, colors(iTarget), Marker=markers(iSubCluster), DisplayName=sprintf('Cluster %i',iSubCluster));
%         end
%     end
%     axis(ax2(iExp), 'image')
%     hold(ax2(iExp), 'off')
%     ax2(iExp).View = [15, 15];
%     if iExp == 1
%         legend(ax2(iExp), h, Location='northoutside', Orientation='horizontal')
%     end

    ax3(iExp) = subplot(1, nExp, iExp);
    xlabel(ax3(iExp), 'ML')
    ylabel(ax3(iExp), 'AP')
    zlabel(ax3(iExp), 'DV')
    hold(ax3(iExp), 'on')
    selFrames = nFrames - nt + 1:nFrames;
    h = gobjects(2, 1);
    for iTarget = 1:2
        for iSubCluster = 1:nSubClusters
            selTrials = targetIndex==iTarget & traj(iExp).cluster.class==iSubCluster;
            x = traj(iExp).contra.x(selTrials, selFrames);
            y = traj(iExp).contra.y(selTrials, selFrames);
            z = traj(iExp).contra.z(selTrials, selFrames);
            % assert(all(x(:, 1) == 0))
            % assert(all(y(:, 1) == 0))
            % assert(all(z(:, 1) == 0))
            h(iTarget) = plot3(ax3(iExp), median(x, 1, 'omitnan'), median(y, 1, 'omitnan'), median(z, 1, 'omitnan'), Color=colors(iTarget), DisplayName=targetNames{iTarget}, LineWidth=1.5);
            scatter3(ax3(iExp), median(x, 1, 'omitnan'), median(y, 1, 'omitnan'), median(z, 1, 'omitnan'), 2*(selFrames).^1.4, colors(iTarget), Marker=markers(iSubCluster), DisplayName=targetNames{iTarget})

            % x = traj(iExp).ipsi.x(selTrials, selFrames);
            % y = traj(iExp).ipsi.y(selTrials, selFrames);
            % z = traj(iExp).ipsi.z(selTrials, selFrames);
            % plot3(ax3(iExp), median(x, 1, 'omitnan'), median(y, 1, 'omitnan'), median(z, 1, 'omitnan'), Color=colors(iTarget), DisplayName=targetNames{iTarget}, LineWidth=1.5, LineStyle='--');
            % scatter3(ax3(iExp), median(x, 1, 'omitnan'), median(y, 1, 'omitnan'), median(z, 1, 'omitnan'), 2*(selFrames).^1.4, colors(iTarget), Marker=markers(iSubCluster), DisplayName=targetNames{iTarget})            
        end        
    end
    legend(h)
    hold(ax3(iExp), 'off')
    ax3(iExp).ZAxis.Direction = 'reverse';
    ax3(iExp).YAxis.Direction = 'reverse';
    axis(ax3(iExp), 'image')
    if ismember(pa.exp(iExp).animalName, {'desmond29'})
        ax3(iExp).XAxis.Direction = 'reverse';
    end
    ax3(iExp).View = [15, 15];
end
ylim(ax1, [0, 100])
xlabel(ax1, 'dims')
ylabel(ax1, 'var explained')
legend(ax1)
clear iExp f

% Plot PETH for each sub-clustered trajectories
% Calculate ETA
targetNames = {'contra-out', 'contra-front', 'contra-in', 'ipsi-front'};
for iExp = 1:nExp
    [~, targetIndex] = ismember(traj(iExp).target, targetNames);
    targetIndex = targetIndex';
    for iTarget = 1:2
        for iSubCluster = 1:nSubClusters
            sel = traj(iExp).cluster.class==iSubCluster & targetIndex==iTarget;
            traj(iExp).eta(iTarget, iSubCluster) = pa.exp(iExp).eu.getETA('count', 'press', window=[-3, 1], resolution=0.1, normalize=[-3, -1], alignTo='stop', includeInvalid=true, ...
                trials=pa.exp(iExp).eu(1).Trials.Press(sel));
        end
    end
end

% Plot ETA
% close all

for iExp = 1:nExp
    figure();
    i = 0;
    for iTarget = 1:2
        for iSubCluster = 1:nSubClusters
            i = i + 1;
            ax = subplot(nSubClusters, 4, i);
            if iSubCluster == 1 && iTarget == 1
                [~, order] = EphysUnit.plotETA(ax, traj(iExp).eta(iTarget, iSubCluster), clim=[-1.5, 1.5], sortWindow=[-1, 0]);
            else
                [~, order] = EphysUnit.plotETA(ax, traj(iExp).eta(iTarget, iSubCluster), clim=[-1.5, 1.5], order=order);
            end
            title(ax, sprintf('exp%i, tgt%i, clu%i (%i trials)', iExp, iTarget, iSubCluster, traj(iExp).eta(iTarget, iSubCluster).N(1)));
        end
    end
end

%% Divide trials by modulation sign/amplitude/latency and draw average trajectories
close all
nSubClusters = 3;
bin = cell(length(eu), 1);
for iEu = 1:length(eu)
    [x, t, d] = eu(iEu).getTrialAlignedData('count', [-3, 0], 'press', alignTo='stop', resolution=0.25, includeInvalid=true);
    % figure, plot(t, x'), hold on, plot(t, mean(x, 1), LineWidth=3)
    x0 = mean(x(:, t>=-0.5&t<=0), 2);
    edges = prctile(x0, 0:100/nSubClusters:100);
    [~, ~, bin{iEu}] = histcounts(x0, edges);
    if ~all(ismember(1:nSubClusters, bin{iEu}))
        warning('%i %s\n', iEu, num2str(ismember(1:nSubClusters, bin{iEu})))
    end
end

nEu = arrayfun(@(e) length(e.eu), pa.exp);
edges = [0, cumsum(nEu)];
expIndices = NaN(length(eu), 1);
for iExp = 1:length(pa.exp)
    expIndices(edges(iExp) + 1:edges(iExp + 1)) = iExp;
end

clear iEu nEu x0 edges x t d edges

colors = [1, 0, 0; 0, 1, 0; 0, 0, 1];
%% 2.2 This opens the user interface. If you close it by mistake run this to 
% restart. Your labelled data should not be affected unless you exit MATLAB

for iEu = 1:length(eu)
    iExp = expIndices(iEu);
    ax = axes(figure());
    h = gobjects(iBin, 1);
    hold(ax, 'on')
    for iBin = 1:nSubClusters
        selTrials = bin{iEu}==iBin;
        x = traj(iExp).contra.x(selTrials, selFrames);
        y = traj(iExp).contra.y(selTrials, selFrames);
        z = traj(iExp).contra.z(selTrials, selFrames);
        h(iBin) = plot3(ax, mean(x, 1, 'omitnan'), mean(y, 1, 'omitnan'), mean(z, 1, 'omitnan'), Color=colors(iBin, :), DisplayName=sprintf('div%i - %i trials', iBin, nnz(selTrials)), LineWidth=4);
        scatter3(ax, mean(x, 1, 'omitnan'), mean(y, 1, 'omitnan'), mean(z, 1, 'omitnan'), 2*(selFrames).^1.4, colors(iBin, :), Marker=markers(iBin))
        plot3(ax, x', y', z', Color=[colors(iBin, :), 0.1])
    end
    hold(ax, 'off')
    ax.ZAxis.Direction = 'reverse';
    ax.YAxis.Direction = 'reverse';
    axis(ax, 'image')
    if ismember(pa.exp(iExp).animalName, {'desmond29'})
        ax.XAxis.Direction = 'reverse';
    end
    ax.View = [15, 15];
    xlabel(ax, 'ML')
    ylabel(ax, 'AP')
    zlabel(ax, 'DV')
    legend(ax, h)
end
%% Finished manual labeling, now we combine front and side view and figure out the scaling
% close all
% fig = figure(Units='pixels', Position=[50, 200, 1500, 600]);
% slider = uicontrol(Parent=fig, Style='slider', value=16, min=1, max=16, Position=[287   142   986    20], SliderStep=[1/15, 1/15]);
% slider.Callback = @(src, evnt) updatePlot(pa, src, evnt);
% 
% function updatePlot(obj, src, evnt)
%     disp(src.Value)
%     for iExp = 1:5       
%         ax = subplot(1, 5, iExp);
%         cla(ax)
%         hold(ax, 'on')
%         xlabel(ax, 'X'), ylabel(ax, 'Y'), zlabel(ax, 'Z')
% 
%         [contra, ipsi, targetPos, t] = obj.getTrajectories(iExp);
% 
%         colors = ...
%         [...
%             1, 0, 0; ...
%             0, 1, 0; ...
%             0, 0, 1; ...
%             0.5, 0.5, 0.5 ...
%         ];
%         posNames = {'contra-out', 'contra-front', 'contra-in', 'ipsi-front'};
%         h = gobjects(1, 4);
%         selFrames = 1:src.Value;
%         for iPos = 1:4
%             relPos = posNames{iPos};
%             selTrials = targetPos==relPos;
%             plot3(ax, contra.x(selTrials, selFrames)', contra.y(selTrials, selFrames)', contra.z(selTrials, selFrames)', Color=[0.5*colors(iPos, :), 0.1], LineStyle='-', LineWidth=0.01);
% 
%             h(iPos) = plot3(ax, median(contra.x(selTrials, selFrames), 1, 'omitnan'), median(contra.y(selTrials, selFrames), 1, 'omitnan'), median(contra.z(selTrials, selFrames), 1, 'omitnan'), LineWidth=2.5, Color=colors(iPos, :), DisplayName=relPos);
%             scatter3(ax, median(contra.x(selTrials, selFrames), 1, 'omitnan'), median(contra.y(selTrials, selFrames), 1, 'omitnan'), median(contra.z(selTrials, selFrames), 1, 'omitnan'), 2*(selFrames).^1.4, colors(iPos, :), DisplayName=relPos)
%         end
%         ax.ZAxis.Direction = 'reverse';
%         if ismember(obj.exp(iExp).animalName, {'desmond29'})
%             ax.YAxis.Direction = 'reverse';
%         end
%         axis(ax, 'image')
%         ax.View = [15, 15];
%         if iExp == 1
%             legend(ax, h, Orientation='horizontal', Position=[0.330369750566468,0.765833336810274,0.382999994665384,0.059999998609225])
%         end
%     end
% 
% end