%% DONE: Make EphysUnit objects (ONLY RUN ONCE)
expNames = {...
    %"desmond28_20230530",...
	"desmond29_20230607",...
	"desmond29_20230608",...
	"desmond29_20230614",...
	"desmond30_20230615",...
	"desmond30_20230616"...
};

for iExp = 1:length(expNames)
    try
        clear tr ar bmr eu
        tr = TetrodeRecording.BatchLoadSimple(expNames{iExp}, true);

        % Make AR
        ar = AcuteRecording(tr, 'WT');
        bmr = ar.binMoveResponse(tr, 'Press', 'Window', [-1, 0], 'BaselineWindow', [-1, -0], 'Store', true);
        
        % Make EphysUnits from AR
        eu = EphysUnit(ar, tr=tr, savepath='C:\SERVER\Units\acute_3cam_reach_direction', cullITI=false, readWaveforms=false);
        
    catch ME
        warning('Error while processing file %g (%s)', iExp, expNames{iExp});
    end
end

%% Load EU
clear all
eu = EphysUnit.load('C:\SERVER\Units\acute_3cam_reach_direction');

%% Make CompleteExperiment3 from EU
exp = CompleteExperiment3(eu);
exp.alignTimestamps();

% Correct feature names
for iExp = 1:length(exp)
    switch exp(iExp).animalName
        case {'desmond28', 'desmond30'}
            exp(iExp).vtdL = renamevars(exp(iExp).vtdL, ...
                {'handIpsiCam_X', 'handIpsiCam_Y', 'handIpsiCam_Likelihood', 'footIpsiCam_X', 'footIpsiCam_Y', 'footIpsiCam_Likelihood'}, ...
                {'handContra_X', 'handContra_Y', 'handContra_Likelihood', 'footContra_X', 'footContra_Y', 'footContra_Likelihood'});
            exp(iExp).vtdR = renamevars(exp(iExp).vtdR, ...
                {'handIpsiCam_X', 'handIpsiCam_Y', 'handIpsiCam_Likelihood', 'footIpsiCam_X', 'footIpsiCam_Y', 'footIpsiCam_Likelihood'}, ...
                {'handIpsi_X', 'handIpsi_Y', 'handIpsi_Likelihood', 'footIpsi_X', 'footIpsi_Y', 'footIpsi_Likelihood'});
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
    end
end

%% Sanity check: get some video clips at lever touch
clear clips
nTrials = 5;
clips = cell(length(exp), nTrials);
for iExp = 1:length(exp)
    trialIndices = round(linspace(1, length(exp(iExp).eu(1).Trials.Press), nTrials+2));
    trialIndices = trialIndices(2:end-1);
    timestamps = [exp(iExp).eu(1).Trials.Press(trialIndices).Stop];
    clipsF = exp(iExp).getVideoClip(timestamps, side='f', bodyparts={'handIpsi', 'handContra', 'footIpsi', 'footContra', 'jaw', 'lever'}, numFramesBefore=30, numFramesAfter=30);

    switch exp(iExp).animalName
        case {'desmond28', 'desmond30'}
            clipsIpsi = exp(iExp).getVideoClip(timestamps, side='r', bodyparts={'handIpsi', 'footIpsi', 'jaw', 'lever'}, numFramesBefore=30, numFramesAfter=30);
            clipsContra = exp(iExp).getVideoClip(timestamps, side='l', bodyparts={'handContra', 'footContra', 'jaw', 'lever'}, numFramesBefore=30, numFramesAfter=30);
        case 'desmond29'
            clipsIpsi = exp(iExp).getVideoClip(timestamps, side='l', bodyparts={'handIpsi', 'footIpsi', 'jaw', 'lever'}, numFramesBefore=30, numFramesAfter=30);
            clipsContra = exp(iExp).getVideoClip(timestamps, side='r', bodyparts={'handContra', 'footContra', 'jaw', 'lever'}, numFramesBefore=30, numFramesAfter=30);
    end
    for iTrial = 1:length(timestamps)
        clips{iExp, iTrial} = horzcat(clipsF{iTrial}, clipsIpsi{iTrial}, clipsContra{iTrial});
    end
end
clear iExp timestamps clipsF clipsIpsi clipsContra trialIndices

% implay(clips{iSession, iTrial}), remember to maximize window
implay(clips{1, 1})

%% Plot ETA (normalized) for each EphysUnit, group by leverPos
eu = [exp.eu];
clear eta

for iEu = 1:length(eu)
    trials = eu(iEu).getTrials('press');
    tTouch = [trials.Stop];
    motPos = eu(iEu).getMotorState(tTouch);
    nTrials = histcounts(motPos, [0.5, 1.5, 2.5, 3.5, 4.5]);
    eta(4) = struct('X', [], 't', [], 'N', [], 'D', [], 'stats', []);
    for iPos = 1:4
        eta(iPos) = eu(iEu).getETA('count', 'press', window=[-4, 2], resolution=0.1, normalize=[-4,-1], ...
            alignTo='stop', includeInvalid=true, trials=trials(motPos == iPos));
    end
    
    ax = axes(figure(iEu));
    hold(ax, 'on')
    for iPos = 1:4
        plot(ax, eta(iPos).t, eta(iPos).X / 0.1, Color=getColor(iPos, 4, 0.8), DisplayName=sprintf('pos %i (%i trials)', iPos, nTrials(iPos)));
    end
    hold(ax, 'off')
    legend(ax, location='southwest')
    ylabel(ax, 'Z-scored \Deltaspike rate (a.u.)')
    xlabel(ax, 'Time to bar contact (s)')
    ylim(ax, [-30, 30])
    title(ax, eu(iEu).getName, Interpreter='none')
    
    print(ax.Parent, sprintf('C:\\SERVER\\Figures\\lever_4pos\\eta_norm\\%s', eu(iEu).getName), '-dpng')

    close(ax.Parent)

    clear trials tTouch motPos nTrials eta ax
end

%% Plot ETA for each EphysUnit, group by leverPos
eu = [exp.eu];
clear eta

for iEu = 1:length(eu)
    trials = eu(iEu).getTrials('press');
    tTouch = [trials.Stop];
    motPos = eu(iEu).getMotorState(tTouch);
    nTrials = histcounts(motPos, [0.5, 1.5, 2.5, 3.5, 4.5]);
    eta(4) = struct('X', [], 't', [], 'N', [], 'D', []);
    for iPos = 1:4
        eta(iPos) = eu(iEu).getETA('count', 'press', window=[-4, 2], resolution=0.1, normalize='none', ...
            alignTo='stop', includeInvalid=true, trials=trials(motPos == iPos));
    end
    
    ax = axes(figure(iEu));
    hold(ax, 'on')
    for iPos = 1:4
        plot(ax, eta(iPos).t, eta(iPos).X / 0.1, Color=getColor(iPos, 4, 0.8), DisplayName=sprintf('pos %i (%i trials)', iPos, nTrials(iPos)));
    end
    hold(ax, 'off')
    legend(ax, location='northwest')
    ylabel(ax, 'Spike rate (sp/s)')
    xlabel(ax, 'Time to bar contact (s)')
    ylim(ax, [0, 150])
    title(ax, eu(iEu).getName, Interpreter='none')
    
    print(ax.Parent, sprintf('C:\\SERVER\\Figures\\lever_4pos\\eta\\%s', eu(iEu).getName), '-dpng')

    close(ax.Parent)

    clear trials tTouch motPos nTrials eta ax
end

%% Gather 2d trajctories from front cam
p.resolution = 1/30;
p.window = [-1, 0];
p.features = {'handIpsi', 'handContra', 'jaw', 'lever'};
p.llThreshold = 0.90;

clear trajectories X Y L XX YY LL;

for iExp = 1:length(exp)
    for iFeat = 1:length(p.features)
        trials = exp(iExp).eu(1).getTrials('press');
        motPos = exp(iExp).eu(1).getMotorState([trials.Start]);
        nTrials = histcounts(motPos, [0.5, 1.5, 2.5, 3.5, 4.5]);
        t = p.window(1):p.resolution:p.window(2);
        X = zeros(4, length(t));
        Y = X;
        L = X;
        
        for iPos = 1:4
            [XX, YY, LL, ~] = exp(iExp).getTrajectoryByTrial('f', p.features{iFeat}, trials=trials(motPos==iPos), window=p.window, likelihoodThreshold=p.llThreshold);
            X(iPos, :) = mean(XX, 1, 'omitnan');
            Y(iPos, :) = mean(YY, 1, 'omitnan');
            L(iPos, :) = mean(LL, 1, 'omitnan');
        end

        [XAll, YAll, LLAll] = exp(iExp).getTrajectoryByTrial('f', p.features{iFeat}, trials=trials, window=p.window, likelihoodThreshold=p.llThreshold);

        trajectories(iExp).(p.features{iFeat}).X = X;
        trajectories(iExp).(p.features{iFeat}).Y = Y;
        trajectories(iExp).(p.features{iFeat}).L = L;
        trajectories(iExp).(p.features{iFeat}).XAll= XAll;
        trajectories(iExp).(p.features{iFeat}).YAll= YAll;
        trajectories(iExp).(p.features{iFeat}).LAll= LLAll;
        trajectories(iExp).(p.features{iFeat}).t = t;
        trajectories(iExp).(p.features{iFeat}).n = nTrials;
    end
end

%%
close all
for iExp = 1:length(exp)
    ax = axes(figure(Units='pixels', Position=[0 0 1024 768]));
    hold(ax, 'on')
    h = gobjects(4, 3);
    % for iPos = 1:4
    %     h(iPos, 1) = plot(ax, trajectories(iExp).handIpsi.X(iPos, :), trajectories(iExp).handIpsi.Y(iPos, :), ...
    %         DisplayName=sprintf('handIpsi %i (n=%i)', iPos, trajectories(iExp).handIpsi.n(iPos)), ...
    %         Color=getColor(iPos, 4, 0.7), LineStyle='-');
    %     plot(ax, trajectories(iExp).handIpsi.X(iPos, 1), trajectories(iExp).handIpsi.Y(iPos, 1), Color=getColor(iPos, 4, 0.7), Marker='o', MarkerSize=10)
    %     plot(ax, trajectories(iExp).handIpsi.X(iPos, end), trajectories(iExp).handIpsi.Y(iPos, end), Color=getColor(iPos, 4, 0.7), Marker='.', MarkerSize=25)
    % end
    for iPos = 1:4
        h(iPos, 2) = plot(ax, trajectories(iExp).handContra.X(iPos, :), trajectories(iExp).handContra.Y(iPos, :), ...
            DisplayName=sprintf('handContra %i (n=%i)', iPos, trajectories(iExp).handIpsi.n(iPos)), ...
            Color=getColor(iPos, 4, 0.7), LineStyle='-');
            plot(ax, trajectories(iExp).handContra.X(iPos, 1), trajectories(iExp).handContra.Y(iPos, 1), Color=getColor(iPos, 4, 0.8), Marker='o', MarkerSize=10)
            plot(ax, trajectories(iExp).handContra.X(iPos, end), trajectories(iExp).handContra.Y(iPos, end), Color=getColor(iPos, 4, 0.8), Marker='.', MarkerSize=25)
    end
    for iPos = 1:4
        h(iPos, 3) = plot(ax, trajectories(iExp).jaw.X(iPos, :), trajectories(iExp).jaw.Y(iPos, :), ...
            DisplayName=sprintf('jaw %i (n=%i)', iPos, trajectories(iExp).jaw.n(iPos)), ...
            Color=getColor(iPos, 4, 0.8), LineStyle=':');
    end
    axis(ax, 'image');
    ax.YDir = 'reverse';
%     xlim(ax, [0, 640])
%     ylim(ax, [0, 480])
    hold(ax, 'off')
    h = h(:, 2:3);
    legend(ax, h(:), Interpreter='none', Location='northwest')
    xlabel('X')
    ylabel('Y')
    title(sprintf('%s', exp(iExp).name), Interpreter='none')

    % print(ax.Parent, sprintf('C:\\SERVER\\Figures\\lever_4pos\\trajectories\\2d_%s', exp(iExp).name), '-dpng')
%     close(ax.Parent)
end
        

%% 2.1 Get training video clips for DLC
% Training set
trainClipsF = {};
trainClipsL = {};
trainClipsR = {};
trainTimesF = [];
trainTimesL = [];
trainTimesR = [];

for iExp = 1:length(exp)
    trials = exp(iExp).eu(1).getTrials('press');
    motPos = exp(iExp).eu(1).getMotorState([trials.Stop]);
    selTrials = [];
    for iPos = 1:4
        inPos = find(motPos == iPos);
        sel = randsample(length(inPos), 1, false);
        selTrials = vertcat(selTrials, trials(inPos(sel)));
    end
    [clipsF, tF] = exp(iExp).getVideoClip([selTrials.Stop], 'f', numFramesBefore=10, numFramesAfter=5, bodyParts={});
    [clipsL, tL] = exp(iExp).getVideoClip([selTrials.Stop], 'l', numFramesBefore=10, numFramesAfter=5, bodyParts={});
    [clipsR, tR] = exp(iExp).getVideoClip([selTrials.Stop], 'r', numFramesBefore=10, numFramesAfter=5, bodyParts={});
    trainClipsF = vertcat(trainClipsF, clipsF);
    trainClipsL = vertcat(trainClipsL, clipsL);
    trainClipsR = vertcat(trainClipsR, clipsR);
    tF = tF';
    tL = tL';
    tR = tR';
%     trainTimesF = vertcat(trainTimesF, tF(:));
%     trainTimesL = vertcat(trainTimesL, tL(:));
%     trainTimesR = vertcat(trainTimesR, tR(:));
end

trainMovieF = cat(4, trainClipsF{:});
trainMovieL = cat(4, trainClipsL{:});
trainMovieR = cat(4, trainClipsR{:});

% Write to disk
v = VideoWriter(sprintf('C:\\SERVER\\VideoTracking\\videos\\training_F_short.mp4'), 'MPEG-4');
open(v)
writeVideo(v, trainMovieF);
close(v);

v = VideoWriter(sprintf('C:\\SERVER\\VideoTracking\\videos\\training_L_short.mp4'), 'MPEG-4');
open(v)
writeVideo(v, trainMovieL);
close(v);

v = VideoWriter(sprintf('C:\\SERVER\\VideoTracking\\videos\\training_R_short.mp4'), 'MPEG-4');
open(v)
writeVideo(v, trainMovieR);
close(v);

%% 2.2 Get whole movies (only press trials)
% for iExp = 1:length(exp)
%     exp(iExp).eu(2:end) = [];
% end

clearvars -except eu exp sides iside
sides = 'FLR';
for iside = 1:3
    side = sides(iside);
    
    v = VideoWriter(sprintf('C:\\SERVER\\VideoTracking\\videos\\movie_%s.mp4', side), 'MPEG-4');
    v.FrameRate = 30;
    open(v)

    ephysTimes = [];
    motPos = [];
    trials = [];
    for iExp = 1:length(exp)
        theseTrials = exp(iExp).eu(1).getTrials('press');
        motPos = vertcat(motPos, exp(iExp).eu(1).getMotorState([theseTrials.Stop]));
        [clips, t] = exp(iExp).getVideoClip([theseTrials.Stop], lower(side), numFramesBefore=60, numFramesAfter=30, bodyParts={});
        clips = cat(4, clips{:});
        writeVideo(v, cat(4, clips));  
        t = t';
        ephysTimes = vertcat(ephysTimes, t(:));
        trials = vertcat(trials, theseTrials(:));
        clear clips t theseTrials
    end
    
    close(v);
    expNames = {exp.name};

    save(sprintf('movie_%s_meta.mat', side), 'ephysTimes', 'motPos', 'expNames', 'trials', '-mat', '-v7.3')


    clearvars -except eu exp sides iside
end

%% 2.3 Read data from DLC
[vtdF, meta, motPos, trials] = exp.readVideoTrackingDataShortConcatenated('C:\SERVER\VideoTracking\videos', 'movie_F', 2, 'f', assign=false);

% Correct feature names
for iExp = 1:length(exp)
    switch exp(iExp).animalName
        case {'desmond28', 'desmond30'}
            vtdF{iExp} = renamevars(vtdF{iExp}, ...
                {'handLeft_X', 'handLeft_Y', 'handLeft_Likelihood', 'handRight_X', 'handRight_Y', 'handRight_Likelihood', ...
                'handLeft_Smooth_X', 'handLeft_Smooth_Y', 'handLeft_Smooth_Likelihood', 'handRight_Smooth_X', 'handRight_Smooth_Y', 'handRight_Smooth_Likelihood'}, ...
                {'handContra_X', 'handContra_Y', 'handContra_Likelihood', 'handIpsi_X', 'handIpsi_Y', 'handIpsi_Likelihood', ...
                'handContra_Smooth_X', 'handContra_Smooth_Y', 'handContra_Smooth_Likelihood', 'handIpsi_Smooth_X', 'handIpsi_Smooth_Y', 'handIpsi_Smooth_Likelihood'});
        case 'desmond29'
            vtdF{iExp} = renamevars(vtdF{iExp}, ...
                {'handLeft_X', 'handLeft_Y', 'handLeft_Likelihood', 'handRight_X', 'handRight_Y', 'handRight_Likelihood', ...
                'handLeft_Smooth_X', 'handLeft_Smooth_Y', 'handLeft_Smooth_Likelihood', 'handRight_Smooth_X', 'handRight_Smooth_Y', 'handRight_Smooth_Likelihood'}, ...
                {'handIpsi_X', 'handIpsi_Y', 'handIpsi_Likelihood', 'handContra_X', 'handContra_Y', 'handContra_Likelihood', ...
                'handIpsi_Smooth_X', 'handIpsi_Smooth_Y', 'handIpsi_Smooth_Likelihood', 'handContra_Smooth_X', 'handContra_Smooth_Y', 'handContra_Smooth_Likelihood'});
    end
end


%% Gather 2d trajctories from front cam
p.resolution = 1/30;
p.window = [-1.5, 0.25];
p.features = {'handIpsi_Smooth', 'handContra_Smooth', 'jaw', 'lever'};
p.llThreshold = 0;

clear trajectories X Y L XX YY LL;

for iExp = 1:length(exp)
    for iFeat = 1:length(p.features)
        trials = exp(iExp).eu(1).getTrials('press');
        motPos = exp(iExp).eu(1).getMotorState([trials.Start]);
        nTrials = histcounts(motPos, [0.5, 1.5, 2.5, 3.5, 4.5]);
        t = p.window(1):p.resolution:p.window(2);
        X = zeros(4, length(t));
        Y = X;
        L = X;
        
        XX = cell(4, 1);
        YY = cell(4, 1);
        LL = cell(4, 1);
        for iPos = 1:4
            switch p.features{iFeat}
                case {'handIpsi', 'handContra', 'handIpsi_Smooth', 'handContra_Smooth'}
                    [XX{iPos}, YY{iPos}, LL{iPos}, ~] = exp(iExp).getTrajectoryByTrial('f', p.features{iFeat}, ...
                        trials=trials(motPos==iPos), window=p.window, likelihoodThreshold=p.llThreshold, ...
                        data=vtdF{iExp}, includeInvalid=true);
                otherwise
                    [XX{iPos}, YY{iPos}, LL{iPos}, ~] = exp(iExp).getTrajectoryByTrial('f', p.features{iFeat}, ...
                        trials=trials(motPos==iPos), window=p.window, likelihoodThreshold=p.llThreshold, ...
                        includeInvalid=true);
            end
            X(iPos, :) = median(XX{iPos}, 1, 'omitnan');
            Y(iPos, :) = median(YY{iPos}, 1, 'omitnan');
            L(iPos, :) = median(LL{iPos}, 1, 'omitnan');
        end
        
        switch p.features{iFeat}
            case {'handIpsi', 'handContra', 'handIpsi_Smooth', 'handContra_Smooth'}
                [XXAll, YYAll, LLAll, ~] = exp(iExp).getTrajectoryByTrial('f', p.features{iFeat}, ...
                    trials=trials, window=p.window, likelihoodThreshold=p.llThreshold, ...
                    data=vtdF{iExp}, includeInvalid=true);
            otherwise
                [XXAll, YYAll, LLAll, ~] = exp(iExp).getTrajectoryByTrial('f', p.features{iFeat}, ...
                    trials=trials, window=p.window, likelihoodThreshold=p.llThreshold, ...
                    includeInvalid=true);
        end

        trajectories(iExp).(p.features{iFeat}).X = X;
        trajectories(iExp).(p.features{iFeat}).Y = Y;
        trajectories(iExp).(p.features{iFeat}).L = L;
        trajectories(iExp).(p.features{iFeat}).XX = XX;
        trajectories(iExp).(p.features{iFeat}).YY = YY;
        trajectories(iExp).(p.features{iFeat}).LL = LL;
        trajectories(iExp).(p.features{iFeat}).XXAll = XXAll;
        trajectories(iExp).(p.features{iFeat}).YYAll = YYAll;
        trajectories(iExp).(p.features{iFeat}).LLAll = LLAll;
        trajectories(iExp).(p.features{iFeat}).t = t;
        trajectories(iExp).(p.features{iFeat}).n = nTrials;
    end
end
clear iExp iFeat trials motPos nTrials t X Y L iPos XX YY LL XXAll YYAll LLAll

%% Cluster 2d(front) trajectories
close all
K = 4;
for iExp = 1:5
%     XXc = cat(1, trajectories(iExp).handContra_Smooth.XX{:}); 
%     YYc = cat(1, trajectories(iExp).handContra_Smooth.YY{:});
%     XXi = cat(1, trajectories(iExp).handIpsi_Smooth.XX{:}); 
%     YYi = cat(1, trajectories(iExp).handIpsi_Smooth.YY{:});
    XXc = trajectories(iExp).handContra_Smooth.XXAll; 
    YYc = trajectories(iExp).handContra_Smooth.YYAll;
    XXi = trajectories(iExp).handIpsi_Smooth.XXAll; 
    YYi = trajectories(iExp).handIpsi_Smooth.YYAll;    
    F = [XXc, YYc, XXi, YYi];
    if (nnz(isnan(F)) > 0)
        FF = F';
        nanIndices = find(isnan(FF));
        FF(nanIndices) = FF(nanIndices-1);
        F = FF';
    end
    [~, score] = pca(F, NumComponents=10);
    idx = kmeans(score, K);
%     motPos = vertcat( ...
%         repmat(1, [trajectories(iExp).handContra_Smooth.n(1), 1]), ...
%         repmat(2, [trajectories(iExp).handContra_Smooth.n(2), 1]), ...
%         repmat(3, [trajectories(iExp).handContra_Smooth.n(3), 1]), ...
%         repmat(4, [trajectories(iExp).handContra_Smooth.n(4), 1]));
    trials = exp(iExp).eu(1).getTrials('press');
    motPos = exp(iExp).eu(1).getMotorState([trials.Start]);

    f = figure(Units='normalized', OuterPosition=[(iExp-1)*0.33, 0, 0.33, 1]);
    for iPos = 1:4
        ax = subplot(4, 1, iPos);
        hold(ax, 'on')
        h = [];
        for k = 1:K
            sel = idx==k & motPos==iPos;
            if nnz(sel) <= 2
                continue
            end
            h = [h, plot(mean(XXc(sel, :), 1), mean(YYc(sel, :), 1), DisplayName=sprintf('cluster %i (n=%i)', k, nnz(sel)), Color=getColor(k, K, 0.8))];
            plot(mean(XXc(sel, 1), 1), mean(YYc(sel, 1), 1), Marker='o', LineStyle='none', MarkerSize=3, Color=getColor(k, K, 0.8), DisplayName='Start')
            plot(mean(XXc(sel, end), 1), mean(YYc(sel, end), 1), Marker='.', LineStyle='none', MarkerSize=15, Color=getColor(k, K, 0.8), DisplayName='Stop')

            plot(mean(XXi(sel, :), 1), mean(YYi(sel, :), 1), LineStyle=':', DisplayName=sprintf('cluster %i (n=%i)', k, nnz(sel)), Color=getColor(k, K, 0.8));
            plot(mean(XXi(sel, 1), 1), mean(YYi(sel, 1), 1), Marker='o', LineStyle='none', MarkerSize=3, Color=getColor(k, K, 0.8), DisplayName='Start')
            plot(mean(XXi(sel, end), 1), mean(YYi(sel, end), 1), Marker='.', LineStyle='none', MarkerSize=15, Color=getColor(k, K, 0.8), DisplayName='Stop')
        end


        plot(ax, trajectories(iExp).lever.X(iPos, end), trajectories(iExp).lever.Y(iPos, end), ...
            DisplayName=sprintf('lever %i (n=%i)', iPos, trajectories(iExp).lever.n(iPos)), ...
            Color='black', LineStyle='none', Marker='o', MarkerSize=35);
        
        hold(ax, 'off')
        axis(ax, 'image')
        ax.YDir = 'reverse';
        legend(ax, h, Location='eastoutside')
        xlim(ax, [250, 500])
        ylim(ax, [200, 300])
        title(ax, sprintf('pos=%i', iPos))
        
        print(ax.Parent, sprintf('C:\\SERVER\\Figures\\lever_4pos\\eta_clustered\\trajectory_%s', exp(iExp).name), '-dpng')
    
    end

    % Calculate ETA for trajectory clusters
    for iEu = 1:length(exp(iExp).eu)
        thisEu = exp(iExp).eu(iEu);
        trials = thisEu.getTrials('press');
        motPos = thisEu.getMotorState([trials.Start]);
        eta(4, K) = struct('X', [], 't', [], 'N', [], 'D', []);
        for iPos = 1:4
            for k = 1:K
                selTrials = motPos==iPos & idx==k;
                if nnz(selTrials) > 0
                    eta(iPos, k) = thisEu.getETA('count', 'press', window=[-4, 2], resolution=0.1, normalize='none', ...
                        alignTo='stop', includeInvalid=true, trials=trials(selTrials));
                end
            end
        end
        f = figure(Units='normalized', OuterPosition=[(iExp-1)*0.33, 0, 0.33, 1]);
        for iPos = 1:4
            ax = subplot(4, 1, iPos);
            hold(ax, 'on')
            for k = 1:K
                selTrials = motPos==iPos & idx==k;
                if nnz(selTrials) > 3
                    plot(ax, eta(iPos, k).t, eta(iPos, k).X / 0.1, Color=getColor(k, K, 0.8), DisplayName=sprintf('cluster %i (%i trials)', k, nnz(selTrials)));
                end
            end
            hold(ax, 'off')
            ylim(ax, [0, 150])
            title(ax, sprintf('pos %i', iPos))
            legend(ax, location='northwest')
            if iPos==4
                xlabel(ax, 'Time to bar contact (s)')
                ylabel(ax, 'Spike rate (sp/s)')
            end
        end
        title(ax, eu(iEu).getName, Interpreter='none')
        
        print(f, sprintf('C:\\SERVER\\Figures\\lever_4pos\\eta_clustered\\%s', thisEu.getName), '-dpng')
    
        close(f)
    end
end

clear F FF XXc XXi YYc YYi K iExp nanIndices score idx trials motPos f ax h iPos iEu thisEu
%



%% Plot 2d (front) trajectories
close all
for iExp = 1:length(exp)
    ax = axes(figure(Units='pixels', Position=[0 0 1024 768]));
    hold(ax, 'on')
    h = gobjects(4, 3);
    for iPos = 1:4
        h(iPos, 1) = plot(ax, trajectories(iExp).handIpsi_Smooth.X(iPos, :), trajectories(iExp).handIpsi_Smooth.Y(iPos, :), ...
            DisplayName=sprintf('handIpsi %i (n=%i)', iPos, trajectories(iExp).handIpsi_Smooth.n(iPos)), ...
            Color=getColor(iPos, 4, 0.8), LineStyle='-');
        plot(ax, trajectories(iExp).handIpsi_Smooth.X(iPos, 1), trajectories(iExp).handIpsi_Smooth.Y(iPos, 1), Color=getColor(iPos, 4, 0.8), Marker='o', MarkerSize=20)
        plot(ax, trajectories(iExp).handIpsi_Smooth.X(iPos, end), trajectories(iExp).handIpsi_Smooth.Y(iPos, end), Color=getColor(iPos, 4, 0.8), Marker='.', MarkerSize=50)
    end
    for iPos = 1:4
        h(iPos, 2) = plot(ax, trajectories(iExp).handContra_Smooth.X(iPos, :), trajectories(iExp).handContra_Smooth.Y(iPos, :), ...
            DisplayName=sprintf('handContra %i (n=%i)', iPos, trajectories(iExp).handContra_Smooth.n(iPos)), ...
            Color=getColor(iPos, 4, 0.8), LineStyle='--');
            plot(ax, trajectories(iExp).handContra_Smooth.X(iPos, 1), trajectories(iExp).handContra_Smooth.Y(iPos, 1), Color=getColor(iPos, 4, 0.8), Marker='o', MarkerSize=20)
            plot(ax, trajectories(iExp).handContra_Smooth.X(iPos, end), trajectories(iExp).handContra_Smooth.Y(iPos, end), Color=getColor(iPos, 4, 0.8), Marker='.', MarkerSize=50)
    end
    for iPos = 1:4
        h(iPos, 3) = plot(ax, trajectories(iExp).lever.X(iPos, end), trajectories(iExp).lever.Y(iPos, end), ...
            DisplayName=sprintf('lever %i (n=%i)', iPos, trajectories(iExp).lever.n(iPos)), ...
            Color=getColor(iPos, 4, 0.8), LineStyle='none', Marker='x', MarkerSize=50);
    end
    for iPos = 1:4
        hJaw = plot(ax, trajectories(iExp).jaw.X(iPos, :), trajectories(iExp).jaw.Y(iPos, :), 'k.', ...
            DisplayName='jaw');
    end
    axis(ax, 'image');
    ax.YDir = 'reverse';
    hold(ax, 'off')
    legend(ax, [h(:); hJaw], Interpreter='none', Location='northwest')
    xlabel('X')
    ylabel('Y')
    title(sprintf('%s', exp(iExp).name), Interpreter='none')

    print(ax.Parent, sprintf('C:\\SERVER\\Figures\\lever_4pos\\trajectories\\2d_new_%s', exp(iExp).name), '-dpng')
end

clear iExp ax h iPos

%% Sanity check: get some video clips at lever touch
clear clips
nTrials = 10;

v = VideoWriter(sprintf('C:\\SERVER\\VideoTracking\\videos\\movie_%s.mp4', side), 'MPEG-4');
v.FrameRate = 30;
open(v)

clips = cell(length(exp), nTrials);
for iExp = 1:length(exp)
    trialIndices = round(linspace(1, length(exp(iExp).eu(1).Trials.Press), nTrials+2));
    trialIndices = trialIndices(2:end-1);
    timestamps = [exp(iExp).eu(1).Trials.Press(trialIndices).Stop];
    clipsF = exp(iExp).getVideoClip(timestamps, side='f', bodyparts={'handIpsi', 'handContra'}, numFramesBefore=30, numFramesAfter=30, vtd=vtdF{iExp});

    switch exp(iExp).animalName
        case {'desmond28', 'desmond30'}
            clipsIpsi = exp(iExp).getVideoClip(timestamps, side='r', bodyparts={'handIpsi', 'footIpsi', 'jaw', 'lever'}, numFramesBefore=30, numFramesAfter=30);
            clipsContra = exp(iExp).getVideoClip(timestamps, side='l', bodyparts={'handContra', 'footContra', 'jaw', 'lever'}, numFramesBefore=30, numFramesAfter=30);
        case 'desmond29'
            clipsIpsi = exp(iExp).getVideoClip(timestamps, side='l', bodyparts={'handIpsi', 'footIpsi', 'jaw', 'lever'}, numFramesBefore=30, numFramesAfter=30);
            clipsContra = exp(iExp).getVideoClip(timestamps, side='r', bodyparts={'handContra', 'footContra', 'jaw', 'lever'}, numFramesBefore=30, numFramesAfter=30);
    end
    for iTrial = 1:length(timestamps)
        clips{iExp, iTrial} = horzcat(clipsF{iTrial}, clipsIpsi{iTrial}, clipsContra{iTrial});
    end
end
clear iExp timestamps clipsF clipsIpsi clipsContra trialIndices

% implay(clips{iSession, iTrial}), remember to maximize window
implay(clips{1, 1})



%% 2.2 Get whole movies (only press trials)

v = VideoWriter(sprintf('C:\\SERVER\\VideoTracking\\videos\\movie_labelled_%s.mp4', 'F2'), 'MPEG-4');
v.FrameRate = 30;
open(v)

for iExp = 1:5
    theseTrials = exp(iExp).eu(1).getTrials('press');
    motPos = exp(iExp).eu(1).getMotorState([theseTrials.Stop]);
    for iPos = 1:4
        trialIndices = find(motPos==iPos);
        for i = floor(linspace(4, length(trialIndices)-4, 10))
            iTrial = trialIndices(i);
            [clip, ~] = exp(iExp).getVideoClip([theseTrials(iTrial).Stop], 'f', numFramesBefore=30, numFramesAfter=0, bodyparts={'handIpsi_Smooth', 'handContra_Smooth'}, ...
                vtd=vtdF{iExp}, file='C:\SERVER\VideoTracking\videos\movie_F.mp4');
            writeVideo(v, clip);  
        end
    end
end

close(v);


clear sides iside side v iExp theseTrials clip iTrial