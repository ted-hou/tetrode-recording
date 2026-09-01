cropParams = [50, 75, 240, 320];
scale = 0.5;
nRows = 6;
nPowers = 3;
nDurations = 2;
maxCtrlClipsPerSession = ceil(6*12/length(exp));
clear clips moveTimes
clips(length(exp), nPowers, nDurations) = struct(press=[], lick=[]);
moveTimes(length(exp), nPowers, nDurations) = struct(press=[], lick=[]);
ll = 0;
tTicTotal = tic();
iter = 0;
for iExp = 1:length(exp)
    for iPower = 1:nPowers
        if iPower == 1
            nDurations = 1;
        else
            nDurations = 2;
        end
        for iDuration = 1:nDurations
            for trialType = ["press", "lick"]
                if iPower == 1 % ctrl
                    sel = find(stimData(iExp).trialTypeCtrl==trialType & stimData(iExp).stimCtrl >= 10 & stimData(iExp).stimCtrl < exp(iExp).vtdL.Timestamp(end)-10);
                    if length(sel) > maxCtrlClipsPerSession
                        sel = sel(randperm(length(sel), maxCtrlClipsPerSession));
                    end
                    t0 = stimData(iExp).stimCtrl(sel);
                else
                    sel = stimData(iExp).iPower==iPower-1 & stimData(iExp).iDuration==iDuration & stimData(iExp).trialType==trialType;
                    t0 = stimData(iExp).stimOn(sel);
                end
                iter = iter + 1;
                fprintf(repmat('\b', [1, ll]))
                pProgress = iter/(10*length(exp));
                ll = fprintf("iExp=%i, iPower=%i, iDuration=%i, trialType=%s (n=%i)... %.1f%%... %.1fs elapsed\n", iExp, iPower-1, iDuration, trialType, nnz(sel), 100*pProgress, toc(tTicTotal));

                % [clipsTemp, t] = exp(iExp).getVideoClip(t0, side='l', numFramesBefore=1*30, numFramesAfter=(3+5)*30, bodyParts={'HandL', 'Jaw'}, minLikelihood=0.5);
                [clipsTemp, t] = exp(iExp).getVideoClip(t0, side='l', numFramesBefore=1*30, numFramesAfter=(p.pulseDurations(iDuration)+3)*30, bodyParts={}, minLikelihood=0.5);
                for iClip = 1:length(clipsTemp)
                    % Crop
                    clipsTemp{iClip} = clipsTemp{iClip}(cropParams(1):cropParams(1)+cropParams(3)-1, cropParams(2):cropParams(2)+cropParams(4)-1, :, :);
                    % Downsample
                    clipsTemp{iClip} = imresize(clipsTemp{iClip}, scale);
                end
            
                % Get bar/spout contact latency
                if iPower == 1
                    psmh = stimData(iExp).psmh.ctrl.(trialType).N(sel, :);
                    edges = stimData(iExp).psmh.ctrl.(trialType).edges;
                    centers = (edges(1:end-1) + edges(2:end))/2;
                else
                    psmh = stimData(iExp).psmh.stim.(trialType).N(sel, :);
                    edges = stimData(iExp).psmh.stim.(trialType).edges;
                    centers = (edges(1:end-1) + edges(2:end))/2;
                end
                psmh = psmh(:, centers>=0);
                centers = centers(centers>=0);
                moveTime = Inf(size(psmh, 1), 1);
                for iTrial = 1:size(psmh, 1)
                    i0 = find(psmh(iTrial, :)>0, 1, 'first');
                    if ~isempty(i0)
                        moveTime(iTrial) = centers(i0);
                    end
                end

                % Annotate video with a red dot to indicate bar/spout contact has occured
                tLocal = t - t0;
                assert(length(clipsTemp) == length(moveTime))
                for iClip = 1:length(clipsTemp)
                    selFramesAfterMove = find(tLocal(iClip, :)>=moveTime(iClip));
                    for iFrame = selFramesAfterMove(:)'
                        clipsTemp{iClip}(:, :, :, iFrame) = insertShape(clipsTemp{iClip}(:, :, :, iFrame), "filled-rectangle", [1, 1, 10, 10], Color='red', Opacity=0.67);
                    end
                end

                % Save results
                clips(iExp, iPower, iDuration).(trialType) = clipsTemp;
                moveTimes(iExp, iPower, iDuration).(trialType) = moveTime;
            end
        end
    end
end

%% Combine sessions
for iExp = 1:length(exp)
    for iPower = 1:nPowers
        if iPower == 1
            nDurations = 1;
        else
            nDurations = 2;
        end
        for iDuration = 1:nDurations
            for trialType = ["press", "lick"]
                CLIPS = arrayfun(@(clips) clips.(trialType), clips(1:length(exp), iPower, iDuration), UniformOutput=false);
                CLIPS = cat(1, CLIPS{:});
                MOVETIMES = arrayfun(@(moveTimes) moveTimes.(trialType), moveTimes(1:length(exp), iPower, iDuration), UniformOutput=false);
                MOVETIMES = cat(1, MOVETIMES{:});
                clips(length(exp) + 1, iPower, iDuration).(trialType) = CLIPS; 
                moveTimes(length(exp) + 1, iPower, iDuration).(trialType) = MOVETIMES; 
                clear CLIPS MOVETIMES
            end
        end
    end
end
%% Make big video displaying all trials at once
for iExp = length(exp) + 1
    for iPower = 1:nPowers
        if iPower == 1
            nDurations = 1;
        else
            nDurations = 2;
        end
        for iDuration = 1:nDurations
            for trialType = ["press", "lick"]
                % try
                    fprintf('Sorting clips by movement time...')
                    [~, trialOrder] = sort(moveTimes(iExp, iPower, iDuration).(trialType), 'ascend');
                    clipsTemp = clips(iExp, iPower, iDuration).(trialType)(trialOrder);
                    fprintf('Done.\n')
                    
                    fprintf('Concatenating clips...')
                    nCols = ceil(length(clipsTemp) / nRows);
                    nFrames = size(clipsTemp{1}, 4);
                    h = size(clipsTemp{1}, 1);
                    w = size(clipsTemp{1}, 2);
                    megaClip = zeros(h*nRows, w*nCols, 3, nFrames, 'uint8');
                    i0 = 0;
                    j0 = 0;
                    for iClip = 1:length(clipsTemp)
                        megaClip(i0+1:i0+h, j0+1:j0+w, :, :) = clipsTemp{iClip};
                        j0 = j0 + w;
                        if j0 >= size(megaClip, 2)
                            j0 = 0;
                            i0 = i0 + h;
                        end
                    end
                    fprintf('Done.\n')
            
                    % Write videos to file
                    exportPath = "C:\SERVER\Figures\SNr_VGAT-Cre-CoChR\Videos";
                    if ~exist(exportPath, 'dir')
                        mkdir(exportPath)
                    end
                    if iPower == 1
                        powerDispName = 'ctrl';
                        dispName = sprintf("%s %s %s", trialType, powerDispName);
                        fname = fullfile(exportPath, sprintf("%i_%s_%s_%s.avi", iExp, stimData(iExp).name, trialType, powerDispName));
                    else
                        powerDispName = sprintf("%gmW", 1e3*p.laserPowers(iPower - 1));
                        dispName = sprintf("%s %s %is", trialType, powerDispName, p.pulseDurations(iDuration));
                        fname = fullfile(exportPath, sprintf("%i_%s_%s_%s_%is.avi", iExp, stimData(iExp).name, trialType, powerDispName, p.pulseDurations(iDuration)));
                    end
                    fprintf('Writing video %s...\n', fname);
                    v = VideoWriter(fname);
                    v.FrameRate = 30;
                    open(v)
                    writeVideo(v, megaClip)
                    close(v)

                    % Write a smaller cropped video with consitent aspect ratio for presentation
                    megaClip = megaClip(:, 1:min(w*12, size(megaClip, 2)), :, :);
                    % append some black before video starts
                    titleClip = zeros(size(megaClip, 1), size(megaClip, 2), size(megaClip, 3), 30, 'uint8');
                    for iFrame = 1:30
                        titleClip(:, :, :, iFrame) = insertText(titleClip(:, :, :, iFrame), size(titleClip, [2, 1])/2, dispName, FontSize=36, BoxOpacity=0, TextColor="white", AnchorPoint='center');
                    end
                    for iFrame = 1:size(megaClip, 4)
                        megaClip(:, :, :, iFrame) = insertText(megaClip(:, :, :, iFrame), size(megaClip, [2, 1])/2, dispName, FontSize=36, BoxOpacity=0, TextColor="white", AnchorPoint='center');
                        megaClip(:, :, :, iFrame) = insertText(megaClip(:, :, :, iFrame), [0, 0], sprintf("%.3fs", (iFrame-30-1)/30), FontSize=24, BoxOpacity=0, TextColor="white", AnchorPoint='LeftTop');
                    end
                    megaClip = cat(4, titleClip, megaClip);
                    v = VideoWriter(sprintf("%s_cropped", fname));
                    v.FrameRate = 30;
                    open(v)
                    writeVideo(v, megaClip)
                    close(v)
                % catch
                %     warning('Could not make video for "%i_%s_iPower%i_iDuration%i_%s.avi"', iExp, stimData(iExp).name, iPower, iDuration, trialType)
                % end
            end
        end
    end
end