%% Lever-4-pos

%% Load EU
eu = EphysUnit.load('C:\SERVER\Units\acute_3cam_reach_direction');

%% Make CompleteExperiment3 from EU
exp = CompleteExperiment3(eu);
exp.alignTimestamps();
posOrder = zeros(length(exp), 4);
posNames = {'contraLat', 'contraInt', 'contraMed', 'IpsiInt'};
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

%% Calculate ETA, grouped by lever pos
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

%% Merge ETAs
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

% Calculate META
clear meta
meta(1, 4) = struct('pressDirectional', []);
for iPos = 1:4
    t = etaMerged(1, iPos).t;
    sel = t <= 0 & t >= -0.5;
    meta(iPos).pressDirectional = mean(etaMerged(iPos).X(:, sel), 2, 'omitnan');
end

%% Plot population ETAs, 4 pos side by side
close all
fig = figure();
for iPos = 1:4
    ax = subplot(1, 4, iPos);
    if iPos == 1
        [~, order] = EphysUnit.plotETA(ax, etaMerged(iPos), event='touch', ...
            clim=[-2, 2], xlim=[-2, 0], sortWindow=[-2, 0], signWindow=[-0.5, 0], ...
            sortThreshold=1, negativeSortThreshold=0.5);
    else
        EphysUnit.plotETA(ax, etaMerged(iPos), event='touch', order=order, ...
            clim=[-2, 2], xlim=[-2, 0], sortWindow=[-2, 0], signWindow=[-0.5, 0], ...
            sortThreshold=1, negativeSortThreshold=0.5);
    end
end
clear fig ax iPos order

fig = figure();
for i = 1:3
    for j = i+1:4
        ax = subplot(3, 4, (i-1)*4+j);
        if i == 1 && j == 2
            [~, order] = EphysUnit.plotETA(ax, etaDiff(i, j), event='touch', ...
                clim=[-2, 2], xlim=[-2, 0], sortWindow=[-2, 0], signWindow=[-0.5, 0], ...
                sortThreshold=1, negativeSortThreshold=0.5);
        else
            EphysUnit.plotETA(ax, etaDiff(i, j), event='touch', ...
                clim=[-2, 2], xlim=[-2, 0], sortWindow=[-2, 0], signWindow=[-0.5, 0], ...
                sortThreshold=1, negativeSortThreshold=0.5);
        end
    end
end
clear fig ax i j order

%% Plot META
fig = figure(Units='inches', Position=[0, 0, 3, 2]);
ax = axes(fig);
hold(ax, 'on')
[~, order] = sort(meta(2).pressDirectional);
h = gobjects(4, 1);
for iPos = 1:4
    h(iPos) = plot(ax, 1:length(meta(iPos).pressDirectional), ...
        meta(iPos).pressDirectional(order), '.', Color=getColor(iPos, 4), ...
        DisplayName=posNames{iPos});
end
plot(ax, [1, length(meta(1).pressDirectional)], [0, 0], 'k--', DisplayName='y=0')
legend(ax, h)
xlabel(ax, 'Unit')
ylabel(ax, 'Pre-reach response (a.u.)')
fontsize(ax, p.fontSize, 'points')
fontname(ax, 'Arial')
xlim(ax, [0, length(meta(1).pressDirectional)])

% 7a. Lever-4-pos diagram
% 7b. PETH (example units, 4 color)
% 7c. Tuning curve heatmap for all units