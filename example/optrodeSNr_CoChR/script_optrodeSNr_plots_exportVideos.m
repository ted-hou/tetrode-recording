cropParams = [50, 75, 240, 320];
scale = 0.5;
nRows = 6;
nPowers = 3;
nDurations = 2;
clear clips moveTimes
clips(length(exp), nPowers, nDurations) = struct(press=[], lick=[]);
moveTimes(length(exp), nPowers, nDurations) = struct(press=[], lick=[]);
for iExp = 1:length(exp)
    for iPower = 1:nPowers-1
        if iPower == 0
            nDurations = 1;
        else
            nDurations = 2;
        end
        for iDuration = 1:nDurations
            for trialType = ["press", "lick"]
                if iPower == 0 % ctrl
                    sel = stimData(iExp).trialTypeCtrl==trialType;
                    t0 = stimData(iExp).stimCtrl(sel);
                else
                    sel = stimData(iExp).iPower==iPower & stimData(iExp).iDuration==iDuration & stimData(iExp).trialType==trialType;
                    t0 = stimData(iExp).stimOn(sel);
                end
                fprintf("iExp=%i, iPower=%i, iDuration=%i, trialType=%s (n=%i)\n", iExp, iPower, iDuration, trialType, nnz(sel));

                [clipsTemp, t] = exp(iExp).getVideoClip(t0, side='l', numFramesBefore=1*30, numFramesAfter=(3+5)*30, bodyParts={'HandL', 'Jaw'}, minLikelihood=0.5);
                for iClip = 1:length(clipsTemp)
                    clipsTemp{iClip} = clipsTemp{iClip}(cropParams(1):cropParams(1)+cropParams(3)-1, cropParams(2):cropParams(2)+cropParams(4)-1, :, :);
                end
            
                % Sort clips by movement latency
                if iPower == 0
                    psmh = stimData(iExp).psmh.ctrl.(trialType).N(sel, :);
                    edges = stimData(iExp).psmh.ctrl.(trialType).edges;
                else
                    psmh = stimData(iExp).psmh.stim.(trialType).N(sel, :);
                    edges = stimData(iExp).psmh.stim.(trialType).edges;
                end
                moveTime = Inf(size(psmh, 1), 1);
                for iTrial = 1:size(psmh, 1)
                    i0 = find(psmh(iTrial, :)>0, 1, 'first');
                    if ~isempty(i0)
                        moveTime(iTrial) = mean(edges(i0:i0+1));
                    end
                end
        
                % Save results
                clips(iExp, iPower+1, iDuration).(trialType) = clipsTemp;
                moveTimes(iExp, iPower+1, iDuration).(trialType) = moveTime;
            end
        end
    end
end

%% Combine sessions
for iExp = 1:length(exp)
    for iPower = 1:nPowers-1
        if iPower == 0
            nDurations = 1;
        else
            nDurations = 2;
        end
        for iDuration = 1:nDurations
            for trialType = ["press", "lick"]
                CLIPS = arrayfun(@(clips) clips.(trialType), clips(1:length(exp), iPower+1, iDuration), UniformOutput=false);
                CLIPS = cat(1, CLIPS{:});
                MOVETIMES = arrayfun(@(moveTimes) moveTimes.(trialType), moveTimes(1:length(exp), iPower+1, iDuration), UniformOutput=false);
                MOVETIMES = cat(1, MOVETIMES{:});
                clips(length(exp) + 1, iPower+1, iDuration).(trialType) = CLIPS; 
                moveTimes(length(exp) + 1, iPower+1, iDuration).(trialType) = MOVETIMES; 
                clear CLIPS MOVETIMES
            end
        end
    end
end
%% Make big video displaying all trials at once
for iExp = length(exp) + 1
    for iPower = 1:nPowers-1
        if iPower == 0
            nDurations = 1;
        else
            nDurations = 2;
        end
        for iDuration = 1:nDurations
            for trialType = ["press", "lick"]
                try
                    [~, trialOrder] = sort(moveTimes(iExp, iPower+1, iDuration).(trialType), 'ascend');
                    clipsTemp = clips(iExp, iPower+1, iDuration).(trialType)(trialOrder);
                    
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
            
                    % Writer videos to file
                    exportPath = "C:\SERVER\Figures\SNr_VGAT-Cre-CoChR\Videos";
                    if ~exist(exportPath, 'dir')
                        mkdir(exportPath)
                    end
                    if isempty(videos(iExp).(trialType))
                        continue
                    end
                    fname = fullfile(exportPath, sprintf("%i_%s_iPower%i_iDuration%i_%s.avi", iExp, stimData(iExp).name, iPower, iDuration, trialType));
                    fprintf('Writing video %s...\n', fname);
                    v = VideoWriter(fname);
                    v.FrameRate = 30;
                    open(v)
                    writeVideo(v, megaClip)
                    close(v)
                catch
                    warning('Could not make video for "%i_%s_iPower%i_iDuration%i_%s.avi"', iExp, stimData(iExp).name, iPower, iDuration, trialType)
                end
            end
        end
    end
end