%% 
read_reachDir

%% Lever-4-pos
p.fontSize = 9;

clear layout
layout.w = 7;
layout.h = 5;
layout.top.h = 2;
layout.bottom.h = 3;
layout.top.left.w = 2;
layout.top.middle.w = 4;
layout.top.right.w = 4;

close all
fig = figure(Units='inches', Position=[0, 0, layout.w, layout.h]);
layout.tl = tiledlayout(fig, layout.top.h + layout.bottom.h, 1);
layout.top.tl = tiledlayout(layout.tl, 1, layout.top.left.w + layout.top.right.w + layout.top.right.w, TileSpacing='loose', Padding='loose');
l = layout.top.tl; l.Layout.Tile = 1; l.Layout.TileSpan = [layout.top.h, 1];

layout.bottom.tl = tiledlayout(layout.tl, 1, 4, TileSpacing='compact');
l = layout.bottom.tl; l.Layout.Tile = 1 + layout.top.h; l.Layout.TileSpan = [layout.bottom.h, 1];


% 7b. Grand average PETH
nTrials = zeros(4, 1);
for iExp = 1:length(expReachDir)
    trials = expReachDir(iExp).eu(1).getTrials('press');
    tTouch = [trials.Stop];
    motPos = expReachDir(iExp).eu(1).getMotorState(tTouch);
    for iPos = 1:4
        rawPos = posOrder(iExp, iPos);
        nTrials(iPos) = nTrials(iPos) + nnz(motPos == rawPos);
    end
end
ax = nexttile(layout.top.tl, 1 + layout.top.left.w + layout.top.middle.w, [1, layout.top.right.w]);
hold(ax, 'on')
for iPos = 1:4
    % plot(ax, etaMerged(iPos).t, mean(etaMerged(iPos).X, 1), Color=getColor(iPos, 4, 0.8), DisplayName=sprintf('%s (%i trials)', posNames{iPos}, nTrials(iPos)));
    plot(ax, etaReachDirMerged(iPos).t, mean(etaReachDirMerged(iPos).X, 1), Color=getColor(iPos, 4, 0.8), DisplayName=sprintf('%s', posNames{iPos}));
end
hold(ax, 'off')
legend(ax, Location='northwest')
ylabel(ax, 'Normalized spike rate (a.u.)')
xlabel(ax, 'Time to reach onset (s)')
% ylim(ax, [0, 150])
xlim(ax, [-4, 0])
title(ax, sprintf('Grand average (%i units)', length(euReachDir)))
fontsize(ax, p.fontSize, 'points')
fontname(ax, 'Arial')

clear iExp trials tTouch motPos iPos rawPos ax

% Plot population ETAs, 4 pos side by side
p.etaWindow = [-4, 0];
p.etaSortWindow = [-2, 0];
p.etaSignWindow = [-0.5, 0];
for iPos = 1:4
    ax = nexttile(layout.bottom.tl);
    if iPos == 1
        [~, order] = EphysUnit.plotETA(ax, etaReachDirMerged(iPos), event='reach onset', ...
            clim=[-2, 2], xlim=p.etaWindow, sortWindow=p.etaSortWindow, signWindow=p.etaSignWindow, ...
            sortThreshold=0.3, negativeSortThreshold=0.15);
        yticks(ax, unique([1, 50:50:size(etaReachDir, 1), size(etaReachDir, 1)]))
    else
        EphysUnit.plotETA(ax, etaReachDirMerged(iPos), event='reach onset', order=order, clim=[-1.5, 1.5], xlim=p.etaWindow);
        yticks(ax, []);
    end
    ylabel(ax, '');
    xlabel(ax, '');    
    title(ax, posNames{iPos})
    if iPos < 4
        colorbar(ax, 'off')
    else
        ax.Colorbar.Layout.Tile = 'east';
    end
    fontsize(ax, p.fontSize, 'points')
    fontname(ax, 'Arial')
end
xlabel(layout.bottom.tl, 'Time to reach onset (s)', FontSize=p.fontSize, FontName='Arial')
ylabel(layout.bottom.tl, 'Unit', FontSize=p.fontSize, FontName='Arial')


copygraphics(fig, ContentType='vector')

clear fig ax iPos order

%% 7b. Example trajectories from one session

close all

for iExp = 1
    ax = axes(figure(Units='inches', Position=[2, 2, 2, p.firstRowHeight]));
    hold(ax, 'on')
    h = gobjects(4, 1);
    for iPos = 1:4
        rawPos = posOrder(iExp, iPos);
        switch expReachDir(iExp).animalName
            case 'desmond29'        
                h(iPos) = plot(ax, trajectories(iExp).handContra.X(rawPos, :), trajectories(iExp).handContra.Y(rawPos, :), ...
                    DisplayName=posNames{iPos}, ...
                    Color=getColor(iPos, 4, 0.8), LineStyle='-');
                    plot(ax, trajectories(iExp).handContra.X(rawPos, end), trajectories(iExp).handContra.Y(rawPos, end), Color=getColor(iPos, 4, 0.8), Marker='.', MarkerSize=25)
            case {'desmond28', 'desmond30'}
                h(iPos) = plot(ax, -trajectories(iExp).handContra.X(rawPos, :), trajectories(iExp).handContra.Y(rawPos, :), ...
                    DisplayName=posNames{iPos}, ...
                    Color=getColor(iPos, 4, 0.8), LineStyle='-');
                    plot(ax, -trajectories(iExp).handContra.X(rawPos, end), trajectories(iExp).handContra.Y(rawPos, end), Color=getColor(iPos, 4, 0.8), Marker='.', MarkerSize=25)
        end
    end
    % for iPos = 1:4
    %     h(5) = plot(ax, mean(trajectories(iExp).jaw.X, 'all'), mean(trajectories(iExp).jaw.Y, 'all'), ...
    %         DisplayName='jaw', ...
    %         Color='k', Marker='o', MarkerSize=10, LineStyle='none');
    % end
    % axis(ax, 'image');
    ax.YDir = 'reverse';
%     xlim(ax, [0, 640])
%     ylim(ax, [0, 480])
    hold(ax, 'off')
    h = h(:);
    legend(ax, h(:), Interpreter='none', NumColumns=2, Position=[0.005156439690736,0.839620487870244,0.989583319673936,0.138020829918484])
    % xlabel(ax, 'X')
    % ylabel(ax, 'Y')
    xticks(ax, ax.XLim)
    xticklabels(ax, {'out', 'in'})
    yticks(ax, ax.YLim)
    yticklabels(ax, {'up', 'down'})
    xlim(ax, ax.XLim + [-10, 10])
    ylim(ax, ax.YLim + [-10, 10])
    % title(sprintf('%s', exp(iExp).name), Interpreter='none')
    fontsize(ax, p.fontSize, 'points')
    fontname(ax, 'Arial')

end

%% Plot difference matrix (upper diagonal)
% fig = figure();
% for i = 1:3
%     for j = i+1:4
%         ax = subplot(3, 4, (i-1)*4+j);
%         if i == 1 && j == 2
%             [~, order] = EphysUnit.plotETA(ax, etaDiff(i, j), event='touch', ...
%                 clim=[-2, 2], xlim=[-2, 2], sortWindow=[-1, 1], signWindow=[-0.5, 0.5], ...
%                 sortThreshold=1, negativeSortThreshold=0.5);
%         else
%             EphysUnit.plotETA(ax, etaDiff(i, j), event='touch', ...
%                 clim=[-2, 2], xlim=[-2, 2], sortWindow=[-1, 1], signWindow=[-0.5, 0.5], ...
%                 sortThreshold=1, negativeSortThreshold=0.5);
%         end
%     end
% end
% clear fig ax i j order

% %% 7b. Example trajectories from one session
% 
% close all
% 
% 
% for iExp = 1
%     ax = axes(figure(Units='inches', Position=[2, 2, 2, p.firstRowHeight]));
%     hold(ax, 'on')
%     h = gobjects(4, 1);
%     for iPos = 1:4
%         rawPos = posOrder(iExp, iPos);
%         switch expReachDir(iExp).animalName
%             case 'desmond29'        
%                 h(iPos) = plot(ax, trajectories(iExp).handContra.X(rawPos, :), trajectories(iExp).handContra.Y(rawPos, :), ...
%                     DisplayName=posNames{iPos}, ...
%                     Color=getColor(iPos, 4, 0.8), LineStyle='-');
%                     plot(ax, trajectories(iExp).handContra.X(rawPos, end), trajectories(iExp).handContra.Y(rawPos, end), Color=getColor(iPos, 4, 0.8), Marker='.', MarkerSize=25)
%             case {'desmond28', 'desmond30'}
%                 h(iPos) = plot(ax, -trajectories(iExp).handContra.X(rawPos, :), trajectories(iExp).handContra.Y(rawPos, :), ...
%                     DisplayName=posNames{iPos}, ...
%                     Color=getColor(iPos, 4, 0.8), LineStyle='-');
%                     plot(ax, -trajectories(iExp).handContra.X(rawPos, end), trajectories(iExp).handContra.Y(rawPos, end), Color=getColor(iPos, 4, 0.8), Marker='.', MarkerSize=25)
%         end
%     end
%     % for iPos = 1:4
%     %     h(5) = plot(ax, mean(trajectories(iExp).jaw.X, 'all'), mean(trajectories(iExp).jaw.Y, 'all'), ...
%     %         DisplayName='jaw', ...
%     %         Color='k', Marker='o', MarkerSize=10, LineStyle='none');
%     % end
%     % axis(ax, 'image');
%     ax.YDir = 'reverse';
% %     xlim(ax, [0, 640])
% %     ylim(ax, [0, 480])
%     hold(ax, 'off')
%     h = h(:);
%     legend(ax, h(:), Interpreter='none', NumColumns=2, Position=[0.005156439690736,0.839620487870244,0.989583319673936,0.138020829918484])
%     % xlabel(ax, 'X')
%     % ylabel(ax, 'Y')
%     xticks(ax, ax.XLim)
%     xticklabels(ax, {'out', 'in'})
%     yticks(ax, ax.YLim)
%     yticklabels(ax, {'up', 'down'})
%     xlim(ax, ax.XLim + [-10, 10])
%     ylim(ax, ax.YLim + [-10, 10])
%     % title(sprintf('%s', exp(iExp).name), Interpreter='none')
%     fontsize(ax, p.fontSize, 'points')
%     fontname(ax, 'Arial')
% 
% end
%  
% %% All
% close all
% 
% fig = figure(Units='inches', Position=[2, 2, 10, p.firstRowHeight]);
% 
% for iExp = 1:length(expReachDir)
%     ax = subplot(1, length(expReachDir), iExp);
%     hold(ax, 'on')
%     h = gobjects(3, 1);
%     for iPos = 1:3
%         rawPos = posOrder(iExp, iPos);
%         switch expReachDir(iExp).animalName
%             case 'desmond29'        
%                 h(iPos) = plot(ax, trajectories(iExp).handContra.X(rawPos, :), trajectories(iExp).handContra.Y(rawPos, :), ...
%                     DisplayName=sprintf('%s (n=%i)', posNames{iPos}, trajectories(iExp).handContra.n(rawPos)), ...
%                     Color=getColor(iPos, 4, 0.7), LineStyle='-');
%                     plot(ax, trajectories(iExp).handContra.X(rawPos, end), trajectories(iExp).handContra.Y(rawPos, end), Color=getColor(iPos, 4, 0.8), Marker='.', MarkerSize=25)
%             case {'desmond28', 'desmond30'}
%                 h(iPos) = plot(ax, -trajectories(iExp).handContra.X(rawPos, :), trajectories(iExp).handContra.Y(rawPos, :), ...
%                     DisplayName=sprintf('%s (n=%i)', posNames{iPos}, trajectories(iExp).handContra.n(rawPos)), ...
%                     Color=getColor(iPos, 4, 0.7), LineStyle='-');
%                     plot(ax, -trajectories(iExp).handContra.X(rawPos, end), trajectories(iExp).handContra.Y(rawPos, end), Color=getColor(iPos, 4, 0.8), Marker='.', MarkerSize=25)
%         end
%     end
%     % for iPos = 1:4
%     %     h(5) = plot(ax, mean(trajectories(iExp).jaw.X, 'all'), mean(trajectories(iExp).jaw.Y, 'all'), ...
%     %         DisplayName='jaw', ...
%     %         Color='k', Marker='o', MarkerSize=10, LineStyle='none');
%     % end
%     % axis(ax, 'image');
%     ax.YDir = 'reverse';
% %     xlim(ax, [0, 640])
% %     ylim(ax, [0, 480])
%     hold(ax, 'off')
%     h = h(:);
%     legend(ax, h(:), Interpreter='none', Location='northoutside')
%     % xlabel(ax, 'X')
%     % ylabel(ax, 'Y')
%     xticks(ax, ax.XLim)
%     xticklabels(ax, {'out', 'in'})
%     yticks(ax, ax.YLim)
%     yticklabels(ax, {'up', 'down'})
%     xlim(ax, ax.XLim + [-10, 10])
%     ylim(ax, ax.YLim + [-10, 10])
%     % title(sprintf('%s', exp(iExp).name), Interpreter='none')
%     fontsize(ax, p.fontSize, 'points')
%     fontname(ax, 'Arial')
% 
% end
%       

% %% 7c. Example units
% exampleUnitNames = { ...
%     'desmond29_20230608_Channel7_Unit1', ... % Up identical
%     'desmond29_20230614_Channel102_Unit1', ... % Up, big amplitube mod
%     'desmond29_20230608_Channel40_Unit1', ... % Up amplitude mod
%     'desmond29_20230608_Channel20_Unit1', ... % Down identical
%     'desmond29_20230608_Channel76_Unit1', ... % Sign change maybe?
%     'desmond30_20230615_Channel31_Unit1', ... % Sign change maybe?
%     };
% close all
% 
% fig = figure(Units='inches', Position=[0, 0, p.width - p.rightPanelWidth, p.secondRowHeight]);
% 
% for i = 1:length(exampleUnitNames)
%     iEu = find(strcmpi(euReachDir.getName, exampleUnitNames{i}));
%     trials = euReachDir(iEu).getTrials('press');
%     tTouch = [trials.Stop];
%     motPos = euReachDir(iEu).getMotorState(tTouch);
%     nTrials = zeros(4, 1);
%     iExp = expIndex(iEu);
%     etaEg(4) = struct('X', [], 't', [], 'N', [], 'D', []);
%     for iPos = 1:4
%         rawPos = posOrder(iExp, iPos);
%         nTrials(iPos) = nnz(motPos == rawPos);
%         etaEg(iPos) = euReachDir(iEu).getETA('count', 'press', window=[-4, 2], resolution=0.1, normalize='none', ...
%                 alignTo='stop', includeInvalid=true, trials=trials(motPos == rawPos));
%     end
% 
%     
%     ax = subplot(2, 3, i);
%     hold(ax, 'on')
%     for iPos = 1:4
%         % plot(ax, etaEg(iPos).t, etaEg(iPos).X / 0.1, Color=getColor(iPos, 4, 0.8), DisplayName=sprintf('%s (%i trials)', posNames{iPos}, nTrials(iPos)));
%         plot(ax, etaEg(iPos).t, etaEg(iPos).X / 0.1, Color=getColor(iPos, 4, 0.8), DisplayName=sprintf('%s', posNames{iPos}));
%     end
%     hold(ax, 'off')
%     % legend(ax, location='northwest')
% 
%     if i == 1
%         ylabel(ax, 'Spike rate (sp/s)', Position=[-5.344491194471983,-47.872268899958186,-1])
%         % legend(ax, Position=[0.12708334554592,0.940787484393089,0.777083322840433,0.05208333219505], Orientation='horizontal')
%     end
%     if i == 5
%         xlabel(ax, 'Time to bar contact (s)')
%     end
%     ylim(ax, [0, 150])
%     % title(ax, eu(iEu).getName, Interpreter='none')
%     fontsize(ax, p.fontSize, 'points')
%     fontname(ax, 'Arial')
% 
% end
% clear i iEu trials tTouch motPos iExp rawPos iPos

% %% 7c. META scatter for all units
% fig = figure(Units='inches', Position=[0, 0, (p.width - p.rightPanelWidth) / 2, p.firstRowHeight]);
% ax = axes(fig);
% hold(ax, 'on')
% [~, order] = sort(metaReachDir(2).pressDirectional);
% h = gobjects(4, 1);
% for iPos = 1:4
%     h(iPos) = plot(ax, 1:length(metaReachDir(iPos).pressDirectional), ...
%         metaReachDir(iPos).pressDirectional(order), '.', Color=getColor(iPos, 4), ...
%         DisplayName=posNames{iPos}, MarkerSize=4);
% end
% plot(ax, [1, length(metaReachDir(1).pressDirectional)], [0, 0], 'k--', DisplayName='y=0')
% legend(ax, h)
% xlabel(ax, 'Unit')
% ylabel(ax, 'Peri-touch response (a.u.)')
% fontsize(ax, p.fontSize, 'points')
% fontname(ax, 'Arial')
% xlim(ax, [0, length(metaReachDir(1).pressDirectional)])
% ax.Legend.Position = [0.143905121630677,0.71722468856705,0.35069443927043,0.257812493170301];
% 
% %% 7d. META "tuning curve" heatmap for all units
% clear etaTuning
% etaTuning = struct('X', [], 't', [], 'N', [], 'D', [], 'stats', []);
% etaTuning.X = horzcat(metaReachDir.pressDirectional);
% etaTuning.t = [1, 2, 3, 4];
% 
% fig = figure(Units='inches', Position=[0, 0, p.rightPanelWidth, p.secondRowHeight]);
% ax = axes(fig);
% [~, order] = EphysUnit.plotETA(ax, etaTuning, event='touch', ...
%             clim=[-1, 1], xlim=[0.5, 4.5], sortWindow=[1, 4], signWindow=[1, 4], ...
%             sortThreshold=0.5, negativeSortThreshold=0.25);
% xlabel(ax, '')
% xticks(ax, 1:4)
% xticklabels(ax, posNames)
% title(ax, '')
% copygraphics(fig)
