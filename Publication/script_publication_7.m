%% Lever-4-pos
p.fontSize = 8;
p.width = 7;
p.height = 5;
p.firstRowHeight = 2;
p.secondRowHeight = 3;
p.rightPanelWidth = 2;
p.leftPanelWidth = 5;
%% Load EU
eu = EphysUnit.load('C:\SERVER\Units\acute_3cam_reach_direction');

% Make CompleteExperiment3 from EU
exp = CompleteExperiment3(eu);
exp.alignTimestamps();
posOrder = zeros(length(exp), 4);
posNames = {'contra-out', 'contra-front', 'contra-in', 'ipsi-front'};
expIndex = zeros(length(eu), 1);

i = 1;
% Correct feature names
for iExp = 1:length(exp)
    expIndex(i:i+length(exp(iExp).eu)-1) = iExp;
    i = i + length(exp(iExp).eu);
    switch exp(iExp).animalName
        case {'desmond28', 'desmond30'}
            exp(iExp).vtdL = renamevars(exp(iExp).vtdL, ...
                {'handIpsiCam_X', 'handIpsiCam_Y', 'handIpsiCam_Likelihood', 'footIpsiCam_X', 'footIpsiCam_Y', 'footIpsiCam_Likelihood'}, ...
                {'handContra_X', 'handContra_Y', 'handContra_Likelihood', 'footContra_X', 'footContra_Y', 'footContra_Likelihood'});
            exp(iExp).vtdR = renamevars(exp(iExp).vtdR, ...
                {'handIpsiCam_X', 'handIpsiCam_Y', 'handIpsiCam_Likelihood', 'footIpsiCam_X', 'footIpsiCam_Y', 'footIpsiCam_Likelihood'}, ...
                {'handIpsi_X', 'handIpsi_Y', 'handIpsi_Likelihood', 'footIpsi_X', 'footIpsi_Y', 'footIpsi_Likelihood'});
            posOrder(iExp, :) = [4, 3, 2, 1];

        case 'desmond29'
            exp(iExp).vtdF = renamevars(exp(iExp).vtdF, ...
                {'handIpsi_X', 'handIpsi_Y', 'handIpsi_Likelihood', 'footIpsi_X', 'footIpsi_Y', 'footIpsi_Likelihood', ...
                'handContra_X', 'handContra_Y', 'handContra_Likelihood', 'footContra_X', 'footContra_Y', 'footContra_Likelihood'}, ...
                {'handContra_X', 'handContra_Y', 'handContra_Likelihood', 'footContra_X', 'footContra_Y', 'footContra_Likelihood', ...
                'handIpsi_X', 'handIpsi_Y', 'handIpsi_Likelihood', 'footIpsi_X', 'footIpsi_Y', 'footIpsi_Likelihood'});
            exp(iExp).vtdL = renamevars(exp(iExp).vtdL, ...
                {'handIpsiCam_X', 'handIpsiCam_Y', 'handIpsiCam_Likelihood', 'footIpsiCam_X', 'footIpsiCam_Y', 'footIpsiCam_Likelihood'}, ...
                {'handIpsi_X', 'handIpsi_Y', 'handIpsi_Likelihood', 'footIpsi_X', 'footIpsi_Y', 'footIpsi_Likelihood'});
            exp(iExp).vtdR = renamevars(exp(iExp).vtdR, ...
                {'handIpsiCam_X', 'handIpsiCam_Y', 'handIpsiCam_Likelihood', 'footIpsiCam_X', 'footIpsiCam_Y', 'footIpsiCam_Likelihood'}, ...
                {'handContra_X', 'handContra_Y', 'handContra_Likelihood', 'footContra_X', 'footContra_Y', 'footContra_Likelihood'});
            posOrder(iExp, :) = [1, 2, 3, 4];            
    end
end
clear iExp i

% Calculate ETA, grouped by lever pos
eu = [exp.eu];
clear eta

eta(length(eu), 4) = struct('X', [], 't', [], 'N', [], 'D', [], 'stats', []);
nTrials = zeros(length(eu), 4);

for iEu = 1:length(eu)
    trials = eu(iEu).getTrials('press');
    tTouch = [trials.Stop];
    motPos = eu(iEu).getMotorState(tTouch);
    nTrials(iEu, :) = histcounts(motPos, [0.5, 1.5, 2.5, 3.5, 4.5]);
    iExp = expIndex(iEu);
    for iPos = 1:4
        rawPos = posOrder(iExp, iPos);
        eta(iEu, iPos) = eu(iEu).getETA('count', 'press', window=[-4, 2], resolution=0.1, normalize=[-4,-1], ...
            alignTo='stop', includeInvalid=true, trials=trials(motPos == rawPos));
    end
    
    clear trials tTouch motPos iPos rawPos iExp
end
clear iEu

% Merge ETAs
clear etaMerged
etaMerged(1, 4) = struct('X', [], 't', [], 'N', [], 'D', [], 'stats', []);
for iPos = 1:4
    etaMerged(1, iPos).X = vertcat(eta(:, iPos).X);
    etaMerged(1, iPos).t = eta(1, iPos).t;
    etaMerged(1, iPos).N = NaN(length(eu), 1);
    etaMerged(1, iPos).D = NaN(length(eu), 1);
end

% Diff ETAs
clear etaDiff
etaDiff(3, 4) = struct('X', [], 't', [], 'N', [], 'D', [], 'stats', []);
for i = 1:3
    for j = i+1:4
        etaDiff(i, j).X = etaMerged(j).X - etaMerged(i).X;
        etaDiff(i, j).t = etaMerged(i).t;
        etaDiff(i, j).N = etaMerged(i).N;
        etaDiff(i, j).D = etaMerged(i).D;
    end
end

%% Calculate META
clear meta
meta(1, 4) = struct('pressDirectional', []);
for iPos = 1:4
    t = etaMerged(1, iPos).t;
    sel = t <= 0.5 & t >= -0.5;
    meta(iPos).pressDirectional = mean(etaMerged(iPos).X(:, sel), 2, 'omitnan');
end

%% 7b. Example trajectories from one session

close all


for iExp = 1
    ax = axes(figure(Units='inches', Position=[2, 2, 2, p.firstRowHeight]));
    hold(ax, 'on')
    h = gobjects(4, 1);
    for iPos = 1:4
        rawPos = posOrder(iExp, iPos);
        switch exp(iExp).animalName
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
 
%% All
close all

fig = figure(Units='inches', Position=[2, 2, 10, p.firstRowHeight]);

for iExp = 1:length(exp)
    ax = subplot(1, length(exp), iExp);
    hold(ax, 'on')
    h = gobjects(3, 1);
    for iPos = 1:3
        rawPos = posOrder(iExp, iPos);
        switch exp(iExp).animalName
            case 'desmond29'        
                h(iPos) = plot(ax, trajectories(iExp).handContra.X(rawPos, :), trajectories(iExp).handContra.Y(rawPos, :), ...
                    DisplayName=sprintf('%s (n=%i)', posNames{iPos}, trajectories(iExp).handContra.n(rawPos)), ...
                    Color=getColor(iPos, 4, 0.8), LineStyle='-');
                    plot(ax, trajectories(iExp).handContra.X(rawPos, end), trajectories(iExp).handContra.Y(rawPos, end), Color=getColor(iPos, 4, 0.8), Marker='.', MarkerSize=25)
            case {'desmond28', 'desmond30'}
                h(iPos) = plot(ax, -trajectories(iExp).handContra.X(rawPos, :), trajectories(iExp).handContra.Y(rawPos, :), ...
                    DisplayName=sprintf('%s (n=%i)', posNames{iPos}, trajectories(iExp).handContra.n(rawPos)), ...
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
    legend(ax, h(:), Interpreter='none', Location='northoutside')
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
        

%% 7b. Grand average
nTrials = zeros(4, 1);
for iExp = 1:length(exp)
    trials = exp(iExp).eu(1).getTrials('press');
    tTouch = [trials.Stop];
    motPos = exp(iExp).eu(1).getMotorState(tTouch);
    for iPos = 1:4
        rawPos = posOrder(iExp, iPos);
        nTrials(iPos) = nTrials(iPos) + nnz(motPos == rawPos);
    end
end
ax = axes(figure(Units='inches', Position=[0, 0, (p.width - p.rightPanelWidth) / 2, p.firstRowHeight]));
hold(ax, 'on')
for iPos = 1:4
    % plot(ax, etaMerged(iPos).t, mean(etaMerged(iPos).X, 1), Color=getColor(iPos, 4, 0.8), DisplayName=sprintf('%s (%i trials)', posNames{iPos}, nTrials(iPos)));
    plot(ax, etaMerged(iPos).t, mean(etaMerged(iPos).X, 1), Color=getColor(iPos, 4, 0.8), DisplayName=sprintf('%s', posNames{iPos}));
end
hold(ax, 'off')
legend(ax, location='northwest')
ylabel(ax, 'Normalized spike rate (a.u.)')
xlabel(ax, 'Time to bar contact (s)')
% ylim(ax, [0, 150])
title(ax, sprintf('Grand average (%i units)', length(eu)))
fontsize(ax, p.fontSize, 'points')
fontname(ax, 'Arial')


clear iExp trials tTouch motPos iPos rawPos ax


%% 7c. Example units
exampleUnitNames = { ...
    'desmond29_20230608_Channel7_Unit1', ... % Up identical
    'desmond29_20230614_Channel102_Unit1', ... % Up, big amplitube mod
    'desmond29_20230608_Channel40_Unit1', ... % Up amplitude mod
    'desmond29_20230608_Channel20_Unit1', ... % Down identical
    'desmond29_20230608_Channel76_Unit1', ... % Sign change maybe?
    'desmond30_20230615_Channel31_Unit1', ... % Sign change maybe?
    };
close all

fig = figure(Units='inches', Position=[0, 0, p.width - p.rightPanelWidth, p.secondRowHeight]);

for i = 1:length(exampleUnitNames)
    iEu = find(strcmpi(eu.getName, exampleUnitNames{i}));
    trials = eu(iEu).getTrials('press');
    tTouch = [trials.Stop];
    motPos = eu(iEu).getMotorState(tTouch);
    nTrials = zeros(4, 1);
    iExp = expIndex(iEu);
    etaEg(4) = struct('X', [], 't', [], 'N', [], 'D', []);
    for iPos = 1:4
        rawPos = posOrder(iExp, iPos);
        nTrials(iPos) = nnz(motPos == rawPos);
        etaEg(iPos) = eu(iEu).getETA('count', 'press', window=[-4, 2], resolution=0.1, normalize='none', ...
                alignTo='stop', includeInvalid=true, trials=trials(motPos == rawPos));
    end

    
    ax = subplot(2, 3, i);
    hold(ax, 'on')
    for iPos = 1:4
        % plot(ax, etaEg(iPos).t, etaEg(iPos).X / 0.1, Color=getColor(iPos, 4, 0.8), DisplayName=sprintf('%s (%i trials)', posNames{iPos}, nTrials(iPos)));
        plot(ax, etaEg(iPos).t, etaEg(iPos).X / 0.1, Color=getColor(iPos, 4, 0.8), DisplayName=sprintf('%s', posNames{iPos}));
    end
    hold(ax, 'off')
    % legend(ax, location='northwest')

    if i == 1
        ylabel(ax, 'Spike rate (sp/s)', Position=[-5.344491194471983,-47.872268899958186,-1])
        % legend(ax, Position=[0.12708334554592,0.940787484393089,0.777083322840433,0.05208333219505], Orientation='horizontal')
    end
    if i == 5
        xlabel(ax, 'Time to bar contact (s)')
    end
    ylim(ax, [0, 150])
    % title(ax, eu(iEu).getName, Interpreter='none')
    fontsize(ax, p.fontSize, 'points')
    fontname(ax, 'Arial')

end
clear i iEu trials tTouch motPos iExp rawPos iPos

%% 7c. META scatter for all units
fig = figure(Units='inches', Position=[0, 0, (p.width - p.rightPanelWidth) / 2, p.firstRowHeight]);
ax = axes(fig);
hold(ax, 'on')
[~, order] = sort(meta(2).pressDirectional);
h = gobjects(4, 1);
for iPos = 1:4
    h(iPos) = plot(ax, 1:length(meta(iPos).pressDirectional), ...
        meta(iPos).pressDirectional(order), '.', Color=getColor(iPos, 4), ...
        DisplayName=posNames{iPos}, MarkerSize=4);
end
plot(ax, [1, length(meta(1).pressDirectional)], [0, 0], 'k--', DisplayName='y=0')
legend(ax, h)
xlabel(ax, 'Unit')
ylabel(ax, 'Peri-touch response (a.u.)')
fontsize(ax, p.fontSize, 'points')
fontname(ax, 'Arial')
xlim(ax, [0, length(meta(1).pressDirectional)])
ax.Legend.Position = [0.143905121630677,0.71722468856705,0.35069443927043,0.257812493170301];

%% 7d. META "tuning curve" heatmap for all units
clear etaTuning
etaTuning = struct('X', [], 't', [], 'N', [], 'D', [], 'stats', []);
etaTuning.X = horzcat(meta.pressDirectional);
etaTuning.t = [1, 2, 3, 4];

fig = figure(Units='inches', Position=[0, 0, p.rightPanelWidth, p.secondRowHeight]);
ax = axes(fig);
[~, order] = EphysUnit.plotETA(ax, etaTuning, event='touch', ...
            clim=[-1, 1], xlim=[0.5, 4.5], sortWindow=[1, 4], signWindow=[1, 4], ...
            sortThreshold=0.5, negativeSortThreshold=0.25);
xlabel(ax, '')
xticks(ax, 1:4)
xticklabels(ax, posNames)
title(ax, '')
copygraphics(fig)

%% 7e. Decoding model for target pos

% 7a. Lever-4-pos diagram
DONE
% 7b. PETH (example units, 4 color)
DONE
% 7c. META scatter for all units
DONE
% 7d. Tuning curve heatmap for all units


%% Plot population ETAs, 4 pos side by side
% p.etaWindow = [-1, 1];
% p.etaSortWindow = [-0.5, 0.5];
% p.etaSignWindow = [-0.5, 0.5];
% close all
% fig = figure(Units='inches', Position=[0, 0, p.width, p.secondRowHeight]);
% for iPos = 1:4
%     ax = subplot(1, 4, iPos);
%     if iPos == 1
%         [~, order] = EphysUnit.plotETA(ax, etaMerged(iPos), event='touch', ...
%             clim=[-2, 2], xlim=p.etaWindow, sortWindow=p.etaSortWindow, signWindow=p.etaSignWindow, ...
%             sortThreshold=0.5, negativeSortThreshold=0.25);
%         ax.XLabel.Position = [3.830189632919599,247.4264953440326,1.000000000000014];
%     else
%         EphysUnit.plotETA(ax, etaMerged(iPos), event='touch', order=order, clim=[-2, 2], xlim=p.etaWindow);
%         ylabel(ax, '');
%         xlabel(ax, '');
%         yticks(ax, []);
%     end
%     title(ax, posNames{iPos})
%     if iPos < 4
%         colorbar(ax, 'off')
%     else
%         ax.Colorbar.Position = [0.916232638741429,0.10890052356021,0.011111111111111,0.814659685863874];
%     end
%     fontsize(ax, p.fontSize, 'points')
%     fontname(ax, 'Arial')
% end
% clear fig ax iPos order

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
