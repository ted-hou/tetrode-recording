cropParams = [50, 75, 240, 320];
scale = 0.5;
nRows = 6;
nPowers = 3;
nDurations = 2;
clear videos clips moveTimes
videos(length(exp), nPowers, nDurations) = struct(press=[], lick=[]);
clips(length(exp), nPowers, nDurations) = struct(press=[], lick=[]);
moveTimes(length(exp), nPowers, nDurations) = struct(press=[], lick=[]);
for iExp = 1:length(exp)
    for iPower = 0:nPowers-1
        for iDuration = 1:nDurations
            for trialType = ["press", "lick"]
                if iPower == 0 % ctrl
                    sel = stimData(iExp).trialTypeCtrl==trialType;
                    t0 = stimData(iExp).stimCtrl(sel);
                else
                    sel = stimData(iExp).iPower==iPower & stimData(iExp).iDuration==iDuration & stimData(iExp).trialType==trialType;
                    t0 = stimData(iExp).stimOn(sel);
                end
                [clipsTemp, t] = exp(iExp).getVideoClip(t0, side='l', numFramesBefore=1*30, numFramesAfter=(3+5)*30, bodyParts={'HandL', 'Jaw'}, minLikelihood=0.5);
                for iClip = 1:length(clipsTemp)
                    clipsTemp{iClip} = clipsTemp{iClip}(cropParams(1):cropParams(1)+cropParams(3)-1, cropParams(2):cropParams(2)+cropParams(4)-1, :, :);
                end
            
                % Sort clips by movement latency
                psmh = stimData(iExp).psmh.stim.(trialType).N(sel, :);
                moveTime = Inf(size(psmh, 1), 1);
                for iTrial = 1:size(psmh, 1)
                    i0 = find(psmh(iTrial, :)>0, 1, 'first');
                    if ~isempty(i0)
                        moveTime(iTrial) = mean(stimData(iExp).psmh.stim.(trialType).edges(i0:i0+1));
                    end
                end
        
                % Save results
                clips(iExp).(trialType) = clipsTemp;
                moveTimes(iExp).(trialType) = moveTime;
            end
        end
    end
end

%% Combine sessions
for trialType = ["press", "lick"]
    CLIPS = arrayfun(@(clips) clips.(trialType), clips, UniformOutput=false);
    CLIPS = cat(1, CLIPS{:});
    MOVETIMES = arrayfun(@(moveTimes) moveTimes.(trialType), moveTimes, UniformOutput=false);
    MOVETIMES = cat(1, MOVETIMES{:});
    clips(length(exp) + 1).(trialType) = CLIPS; 
    moveTimes(length(exp) + 1).(trialType) = MOVETIMES; 
    clear CLIPS MOVETIMES
end
%% Make big video displaying all trials at once
for iExp = length(exp) + 1
    for trialType = ["press", "lick"]
        [~, trialOrder] = sort(moveTimes(iExp).(trialType), 'ascend');
        clipsTemp = clips(iExp).(trialType)(trialOrder);
        
        nCols = ceil(length(clipsTemp) / nRows);
        nFrames = size(clipsTemp{1}, 4);
        megaClip = zeros(cropParams(3)*nRows, cropParams(4)*nCols, 3, nFrames, 'uint8');
        i0 = 0;
        j0 = 0;
        for iClip = 1:length(clipsTemp)
            megaClip(i0+1:i0+cropParams(3), j0+1:j0+cropParams(4), :, :) = clipsTemp{iClip};
            j0 = j0 + cropParams(4);
            if j0 >= size(megaClip, 2)
                j0 = 0;
                i0 = i0 + cropParams(3);
            end
        end
        megaClip = imresize(megaClip, scale);
        videos(iExp).(trialType) = megaClip;

        % Writer videos to file
        exportPath = "E:\Figures\SNr_VGAT-Cre-CoChR\Videos";
        if ~exist(exportPath, 'dir')
            mkdir(exportPath)
        end
        if isempty(videos(iExp).(trialType))
            continue
        end
        fname = fullfile(exportPath, sprintf("%i_%s_%s.avi", iExp, stimData(iExp).name, trialType));
        fprintf('Writing video %s...\n', fname);
        v = VideoWriter(fname);
        v.FrameRate = 30;
        open(v)
        writeVideo(v, videos(iExp).(trialType))
        close(v)
    end
end