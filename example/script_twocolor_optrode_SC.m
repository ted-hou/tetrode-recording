

%% Remove duplicates (Slow)
eu = EphysUnit.load('\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\Units\TwoColor_SC\Batch2', waveforms=false, spikecounts=false, spikerates=false);
%
euAll = eu;

% Remove drift, low spike rate units
clear c

% Remove drift
c.isDrifting = detectDriftingUnits(eu, smoothWindow=300, tolerance=0.05, spikeRateThreshold=0.5, includeITI=true);

% Filter by spike rate
msr = arrayfun(@(eu) eu.SpikeRateStats.median, eu);
eu = eu(~c.isDrifting);

eu = eu.removeMultiUnits(cullZeros=true);
[eu, isDuplicate] = eu.removeDuplicates(0.7);


% Make spontaneous reach/lick trials
% p.minTrialLength = 4;
% for iEu = 1:length(eu)
%     eu(iEu).Trials.Press = eu(iEu).makeTrials('press_spontaneous_clean', minSpontaneousTrialDuration=p.minTrialLength);
%     eu(iEu).Trials.Lick = eu(iEu).makeTrials('lick_spontaneous_clean', minSpontaneousTrialDuration=p.minTrialLength);
% end

% eu.save('\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\Units\TwoColor_SC\SingleUnit_NonDuplicate_NonDrift_SNr')
% 
% %% Fix extra red stim event for bad session
% for iEu = find(string({eu.ExpName}) == "desmond38_20250403")
%     assert(isscalar(iEu))
%     assert(eu(iEu).EventTimes.LaserModRedOn(1) - eu(iEu).EventTimes.LaserModRedOff(1) < 1e-6)
%     iBadStimOn = find(eu(iEu).EventTimes.StimOn == eu(iEu).EventTimes.LaserModRedOn(1));
%     iBadStimOff = find(eu(iEu).EventTimes.StimOff == eu(iEu).EventTimes.LaserModRedOff(1));
%     eu(iEu).EventTimes.StimOn(iBadStimOn) = [];
%     eu(iEu).EventTimes.StimOff(iBadStimOff) = [];
%     eu(iEu).EventTimes.LaserModRedOn(1) = [];
%     eu(iEu).EventTimes.LaserModRedOff(1) = [];
%     eu(iEu).Trials.Stim = eu(iEu).makeTrials('stim');
% end
% 
% 
% eu.save('\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\Units\TwoColor_SNr_SCRetro\SingleUnit_NonDuplicate_NonDrift_SNr')

%% Get medial vs. lateral trials, save to EU (do only once)
p.minTrialLength = 4;
[~, ia, ~] = unique({eu.ExpName});
for iEu = 1:length(eu)
    switch eu(iEu).getAnimalName()
        case {'desmond38', 'daisy26'}
            implantSide = 'R';
        case 'desmond39'
            implantSide = 'L';
        otherwise 
            error()
    end


    pressTrials = Trial(eu(iEu).EventTimes.PressOff, eu(iEu).EventTimes.PressOn, 'first');
    pressTrials = pressTrials(pressTrials.duration() >= p.minTrialLength);
    pressTimes = [pressTrials.Stop];
    lickTrials = Trial(eu(iEu).EventTimes.LickOff, eu(iEu).EventTimes.LickOn, 'first');
    lickTrials = lickTrials(lickTrials.duration() >= p.minTrialLength);
    lickTimes = [lickTrials.Stop];

    leftOn = eu(iEu).EventTimes.CueLeftOn;
    leftOff = eu(iEu).EventTimes.CueLeftOff;
    rightOn = eu(iEu).EventTimes.CueRightOn;
    rightOff = eu(iEu).EventTimes.CueRightOff;

    if length(leftOn) == length(leftOff) + 1
        leftOff(end + 1) = Inf;
    end
    if length(rightOn) == length(rightOff) + 1
        rightOff(end + 1) = Inf;
    end
    assert(length(leftOn) == length(leftOff))
    assert(length(rightOn) == length(rightOff))
    assert(all(leftOn - leftOff <= 0))
    assert(all(rightOn - rightOff <= 0))

    cueTrialsLeft = Trial(leftOn, leftOff, advancedValidation=false);
    cueTrialsRight = Trial(rightOn, rightOff, advancedValidation=false);
    isLeft = cueTrialsLeft.inTrial(pressTimes);
    isRight = cueTrialsRight.inTrial(pressTimes);

    eu(iEu).Trials.CueLeft = cueTrialsLeft;
    eu(iEu).Trials.CueRight = cueTrialsRight;
    eu(iEu).Trials.PressLeft = pressTrials(isLeft);
    eu(iEu).Trials.PressRight = pressTrials(isRight);

    % Press trials must occur when one LED was on
    sel = cueTrialsLeft.inTrial(pressTimes) | cueTrialsRight.inTrial(pressTimes);
    if nnz(sel) < length(sel)
        eu(iEu).Trials.Press = pressTrials(sel);
        if ismember(iEu, ia)
            fprintf('Removed %i (of %i) press trials because they occured when both Cue LEDs were off.\n', length(sel) - nnz(sel), length(sel));
        end
    else
        eu(iEu).Trials.Press = pressTrials;
    end

    % Lick trials must occur when both LEDs were off
    sel = ~cueTrialsLeft.inTrial(lickTimes) & ~cueTrialsRight.inTrial(lickTimes);
    if nnz(sel) < length(sel)
        eu(iEu).Trials.Lick = lickTrials(sel);
        if ismember(iEu, ia)
            fprintf('Removed %i (of %i) lick trials because they occured when Cue LED(s) was on.\n', length(sel) - nnz(sel), length(sel));
        end
    else
        eu(iEu).Trials.Lick = lickTrials;
    end

    switch implantSide
        case 'L'
            eu(iEu).Trials.PressSpontaneousMedial = eu(iEu).Trials.PressLeft;
            eu(iEu).Trials.PressSpontaneousLateral = eu(iEu).Trials.PressRight;
        case 'R'
            eu(iEu).Trials.PressSpontaneousMedial = eu(iEu).Trials.PressRight;
            eu(iEu).Trials.PressSpontaneousLateral = eu(iEu).Trials.PressLeft;
    end
end

% Count number of trials by session
arrayfun(@(eu) fprintf('%s: %i med, %i lat, %i lick;\n', eu.ExpName, length(eu.Trials.PressSpontaneousMedial), length(eu.Trials.PressSpontaneousLateral), length(eu.Trials.Lick)), eu(ia));
%%
clear ia sel implantSide cueTrialsLeft cueTrialsRight isLeft isRight leftOn leftOff pressTrials pressTimes lickTrials lickTimes iEu
eu.save('\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\Units\TwoColor_SC\SingleUnit_NonDuplicate_NonDrift_SC')

%%
clear, clc
eu = EphysUnit.load('\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\Units\TwoColor_SC\SingleUnit_NonDuplicate_NonDrift_SC', waveforms=false, spikecounts=false, spikerates=false);

%% Make Raster
clear rd
rd.stim = eu.getRasterData('stimtwocolor', window=[-0.1, 0.15], durErr=1e-3, shutterDelay=0, photoelectricBlankDuration=0.5e-3);
rd.press = eu.getRasterData('press', window=[-4, 0], alignTo='stop');
rd.pressMed = eu.getRasterData('press_spontaneous_medial', window=[-4, 0], alignTo='stop');
rd.pressLat = eu.getRasterData('press_spontaneous_lateral', window=[-4, 0], alignTo='stop');
rd.lick = eu.getRasterData('lick', window=[-4, 0], alignTo='stop');

% Make ETA
clear eta
eta.pressRaw = eu.getETA('count', 'press', [-4, 0], resolution=0.1, alignTo='stop', includeInvalid=false, normalize='none');
eta.pressMedRaw = eu.getETA('count', 'press_spontaneous_medial', [-4, 0], resolution=0.1, alignTo='stop', includeInvalid=false, normalize='none');
eta.pressLatRaw = eu.getETA('count', 'press_spontaneous_lateral', [-4, 0], resolution=0.1, alignTo='stop', includeInvalid=false, normalize='none');
eta.lickRaw = eu.getETA('count', 'lick', [-4, 0], resolution=0.1, alignTo='stop', includeInvalid=false, normalize='none');
eta.press = eu.getETA('count', 'press', [-4, 0], resolution=0.1, alignTo='stop', includeInvalid=false, normalize=[-4, -2]);
eta.pressMed = eu.getETA('count', 'press_spontaneous_medial', [-4, 0], resolution=0.1, alignTo='stop', includeInvalid=false, normalize=[-4, -2]);
eta.pressLat = eu.getETA('count', 'press_spontaneous_lateral', [-4, 0], resolution=0.1, alignTo='stop', includeInvalid=false, normalize=[-4, -2]);
eta.lick = eu.getETA('count', 'lick', [-4, 0], resolution=0.1, alignTo='stop', includeInvalid=false, normalize=[-4, -2]);
eta.pressRaw.X = eta.pressRaw.X ./ 0.1;
eta.pressMedRaw.X = eta.pressMedRaw.X ./ 0.1;
eta.pressLatRaw.X = eta.pressLatRaw.X ./ 0.1;
eta.lickRaw.X = eta.lickRaw.X ./ 0.1;

%% ETA Stim
p.isiBaselineWindow = [-0.1, 0];
p.stimBluePowers = [25 50]*1e-6; 
p.stimRedPowers = [25 50 100 500, 2000, 16000]*1e-6;
p.stimBlueDurations = [20]*1e-3;
p.stimRedDurations = [20]*1e-3;

p.isiWindow = [-0.4, 0.4];
p.isiRes = 1e-3;
p.xlim.stim = [-0.1, 0.3];
p.xlim.move = [-4, 0];
p.rasterSzStim = 1;
p.rasterSzMove = 1;

close all
XBlue = cell(length(eu), 1);
XRed = cell(length(eu), 1);
for iEu = 1:length(eu)
    groupsBlue = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'duration'}, selectBy=struct(power=p.stimBluePowers, duration=p.stimBlueDurations, location=[], wavelength=[470, 473]));
    groupsRed = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'duration'}, selectBy=struct(power=p.stimRedPowers, duration=p.stimRedDurations, location=[], wavelength=635));

    % rd = eu(iEu).getRasterData('stimtwocolor', p.isiWindow, trials=[groupsRed.trials], alignTo='start', shutterDelay=0, sort=false, photoelectricBlankDuration=0.5e-3);
    % EphysUnit.plotRaster(rd)

    [isi, t] = eu(iEu).getMeanPEISI('stimtwocolor', [groupsBlue.trials], window=p.isiWindow, resolution=p.isiRes, photoelectricBlankDuration=0.5e-3);
    selBaseline = t>p.isiBaselineWindow(1) & t<p.isiBaselineWindow(2);
    normSR = (1./isi - mean(1./isi(:, selBaseline), 'omitnan')) ./ std(1./isi(:, selBaseline), 0, 2, 'omitnan');
    XBlue{iEu} = normSR;

    [isi, t] = eu(iEu).getMeanPEISI('stimtwocolor', [groupsRed.trials], window=p.isiWindow, resolution=p.isiRes, ...
        photoelectricBlankDuration=0.5e-3);
    selBaseline = t>p.isiBaselineWindow(1) & t<p.isiBaselineWindow(2);
    normSR = (1./isi - mean(1./isi(:, selBaseline), 'omitnan')) ./ std(1./isi(:, selBaseline), 0, 2, 'omitnan');
    XRed{iEu} = normSR;

    % fprintf('nan=%i, nan=%i\n', nnz(isnan(XBlue{iEu})), nnz(isnan(XRed{iEu})))
end

eta.stimBlue = struct(X=cat(1, XBlue{:}), t=t, N=[], D=[], stats=[]);
eta.stimRed = struct(X=cat(1, XRed{:}), t=t, N=[], D=[], stats=[]);

clear XBlue XRed iEu groupsBlue groupsRed isi t selBaseline normSR

%% Calculate META
clear meta
p.metaWindow = [-0.1, 0];
p.posRespThreshold = 2;
p.negRespThreshold = -1;
t = eta.press.t;
meta.press = mean(eta.press.X(:, t>=p.metaWindow(1) & t<=p.metaWindow(2)), 2, 'omitnan');
meta.pressMed = mean(eta.pressMed.X(:, t>=p.metaWindow(1) & t<=p.metaWindow(2)), 2, 'omitnan');
meta.pressLat = mean(eta.pressLat.X(:, t>=p.metaWindow(1) & t<=p.metaWindow(2)), 2, 'omitnan');
meta.lick = mean(eta.lick.X(:, t>=p.metaWindow(1) & t<=p.metaWindow(2)), 2, 'omitnan');
t = eta.pressRaw.t;
meta.pressRaw = mean(eta.pressRaw.X(:, t>=p.metaWindow(1) & t<=p.metaWindow(2)), 2, 'omitnan');
meta.pressMedRaw = mean(eta.pressMedRaw.X(:, t>=p.metaWindow(1) & t<=p.metaWindow(2)), 2, 'omitnan');
meta.pressLatRaw = mean(eta.pressLatRaw.X(:, t>=p.metaWindow(1) & t<=p.metaWindow(2)), 2, 'omitnan');
meta.lickRaw = mean(eta.lickRaw.X(:, t>=p.metaWindow(1) & t<=p.metaWindow(2)), 2, 'omitnan');
clear t

p.metaWindowStim = [0.005, 0.015];
p.posRespThresholdStim = 2;
p.negRespThresholdStim = -2;

t = eta.stimBlue.t;
meta.stimBlue = mean(eta.stimBlue.X(:, t>=p.metaWindowStim(1) & t<=p.metaWindowStim(2)), 2, 'omitnan');
t = eta.stimRed.t;
meta.stimRed = mean(eta.stimRed.X(:, t>=p.metaWindowStim(1) & t<=p.metaWindowStim(2)), 2, 'omitnan');
clear t
% 
% c.isPressUp =         meta.press >= p.posRespThreshold;
% c.isPressDown =       meta.press <= p.negRespThreshold;
% c.isPressResponsive = c.isPressUp | c.isPressDown;
% c.isLickUp =          meta.lick >= p.posRespThreshold;
% c.isLickDown =        meta.lick <= p.negRespThreshold;
% c.isLickResponsive =  c.isLickUp | c.isLickDown;
% c.isLickUnresponsiveButUp = ~c.isLickResponsive & meta.lick > 0;
% c.isLickUnresponsiveButDown = ~c.isLickResponsive & meta.lick < 0;
% c.isPressUnresponsiveButUp = ~c.isPressResponsive & meta.press > 0;
% c.isPressUnresponsiveButDown = ~c.isPressResponsive & meta.press < 0;

c.isStimBlueUp = meta.stimBlue >= p.posRespThresholdStim;
c.isStimBlueDown = meta.stimBlue <= p.negRespThresholdStim;
c.isStimRedUp = meta.stimRed >= p.posRespThresholdStim;
c.isStimRedDown = meta.stimRed <= p.negRespThresholdStim;

c.isStimBlueUpRedUpThereforeChrimsonMaybe = c.isStimBlueUp & c.isStimRedUp;
c.isStimBlueUpRedNotUpThereforeCoChrMaybe = c.isStimBlueUp & ~c.isStimRedUp;
c.isStimBlueNotUpRedUpThereforeChrimsonMaybe = ~c.isStimBlueUp & c.isStimRedUp;
c.isStimBlueNotUpRedNotUp = ~c.isStimBlueUp & ~c.isStimRedUp;

%% Boot response dir
p.bootAlpha = 0.05;
p.nboot = 100000;
p.metaWindow = [-0.5, 0];
boot.press = struct('h', NaN(length(eu), 1), 'muDiffCI', NaN(length(eu), 2), 'muDiffObs', NaN(length(eu), 1));
boot.pressMed = struct('h', NaN(length(eu), 1), 'muDiffCI', NaN(length(eu), 2), 'muDiffObs', NaN(length(eu), 1));
boot.pressLat = struct('h', NaN(length(eu), 1), 'muDiffCI', NaN(length(eu), 2), 'muDiffObs', NaN(length(eu), 1));
boot.lick = struct('h', NaN(length(eu), 1), 'muDiffCI', NaN(length(eu), 2), 'muDiffObs', NaN(length(eu), 1));
[boot.press.h, boot.press.muDiffCI, boot.press.muDiffObs] = bootstrapMoveResponse( ...
    eu, 'press', nboot=p.nboot, alpha=p.bootAlpha, withReplacement=false, oneSided=false, ...
    responseWindow=p.metaWindow, allowedTrialDuration=[0, Inf]);
[boot.pressMed.h, boot.pressMed.muDiffCI, boot.pressMed.muDiffObs] = bootstrapMoveResponse( ...
    eu, 'press_spontaneous_medial', nboot=p.nboot, alpha=p.bootAlpha, withReplacement=false, oneSided=false, ...
    responseWindow=p.metaWindow, allowedTrialDuration=[0, Inf]);
[boot.pressLat.h, boot.pressLat.muDiffCI, boot.pressLat.muDiffObs] = bootstrapMoveResponse( ...
    eu, 'press_spontaneous_lateral', nboot=p.nboot, alpha=p.bootAlpha, withReplacement=false, oneSided=false, ...
    responseWindow=p.metaWindow, allowedTrialDuration=[0, Inf]);
[boot.lick.h, boot.lick.muDiffCI, boot.lick.muDiffObs] = bootstrapMoveResponse( ...
    eu, 'lick', nboot=p.nboot, alpha=p.bootAlpha, withReplacement=false, oneSided=false, ...
    responseWindow=p.metaWindow, allowedTrialDuration=[0, Inf]);
fprintf(1, '\nAll done\n')

% Report bootstraped movement response direction
assert(nnz(isnan(boot.lick.h)) == 0)
assert(nnz(isnan(boot.press.h)) == 0)

figure, histogram(boot.press.h)
c.isPressUp = boot.press.h == 1;
c.isPressDown = boot.press.h == -1;
c.isPressResponsive = c.isPressUp | c.isPressDown;
c.isPressUnresponsiveButUp = ~c.isPressResponsive & meta.press > 0;
c.isPressUnresponsiveButDown = ~c.isPressResponsive & meta.press < 0;

c.isPressMedUp = boot.pressMed.h == 1;
c.isPressMedDown = boot.pressMed.h == -1;
c.isPressMedResponsive = c.isPressMedUp | c.isPressMedDown;
c.isPressMedUnresponsiveButUp = ~c.isPressMedResponsive & meta.pressMed > 0;
c.isPressMedUnresponsiveButDown = ~c.isPressMedResponsive & meta.pressMed < 0;

c.isPressLatUp = boot.pressLat.h == 1;
c.isPressLatDown = boot.pressLat.h == -1;
c.isPressLatResponsive = c.isPressLatUp | c.isPressLatDown;
c.isPressLatUnresponsiveButUp = ~c.isPressLatResponsive & meta.pressLat > 0;
c.isPressLatUnresponsiveButDown = ~c.isPressLatResponsive & meta.pressLat < 0;

figure, histogram(boot.lick.h)
c.isLickUp = boot.lick.h == 1;
c.isLickDown = boot.lick.h == -1;
c.isLickResponsive = c.isLickUp | c.isLickDown;
c.isLickUnresponsiveButUp = ~c.isLickResponsive & meta.lick > 0;
c.isLickUnresponsiveButDown = ~c.isLickResponsive & meta.lick < 0;

fprintf(1, ['%g/%g (%.0f%%) units (press responsive):\n' ...
    '\t%g (%.0f%%) are excited (p<%g);\n' ...
    '\t%g (%.0f%%) are inhibited (p<%g).\n'], ...
    nnz(c.isPressResponsive), length(eu), 100*nnz(c.isPressResponsive)/length(eu), ...
    nnz(c.isPressUp), 100*nnz(c.isPressUp)/nnz(c.isPressResponsive), p.bootAlpha, ...
    nnz(c.isPressDown), 100*nnz(c.isPressDown)/nnz(c.isPressResponsive), p.bootAlpha);

fprintf(1, ['%g/%g (%.0f%%) units (press medial responsive):\n' ...
    '\t%g (%.0f%%) are excited (p<%g);\n' ...
    '\t%g (%.0f%%) are inhibited (p<%g).\n'], ...
    nnz(c.isPressMedResponsive), length(eu), 100*nnz(c.isPressMedResponsive)/length(eu), ...
    nnz(c.isPressMedUp), 100*nnz(c.isPressMedUp)/nnz(c.isPressMedResponsive), p.bootAlpha, ...
    nnz(c.isPressMedDown), 100*nnz(c.isPressMedDown)/nnz(c.isPressMedResponsive), p.bootAlpha);

fprintf(1, ['%g/%g (%.0f%%) units (press lateral responsive):\n' ...
    '\t%g (%.0f%%) are excited (p<%g);\n' ...
    '\t%g (%.0f%%) are inhibited (p<%g).\n'], ...
    nnz(c.isPressLatResponsive), length(eu), 100*nnz(c.isPressLatResponsive)/length(eu), ...
    nnz(c.isPressLatUp), 100*nnz(c.isPressLatUp)/nnz(c.isPressLatResponsive), p.bootAlpha, ...
    nnz(c.isPressLatDown), 100*nnz(c.isPressLatDown)/nnz(c.isPressLatResponsive), p.bootAlpha);

fprintf(1, ['%g/%g (%.0f%%) units (lick responsive):\n' ...
    '\t%g (%.0f%%) are excited (p<%g);\n' ...
    '\t%g (%.0f%%) are inhibited (p<%g).\n'], ...
    nnz(c.isLickResponsive), length(eu), 100*nnz(c.isLickResponsive)/length(eu), ...
    nnz(c.isLickUp), 100*nnz(c.isLickUp)/nnz(c.isLickResponsive), p.bootAlpha, ...
    nnz(c.isLickDown), 100*nnz(c.isLickDown)/nnz(c.isLickResponsive), p.bootAlpha);

nTotal = nnz(c.isPressResponsive & c.isLickResponsive);
fprintf(1, ['%g/%g (%.0f%%) units (press and lick responsive):\n' ...
    '\t%g (%.0f%%) are press-excited AND lick-excited;\n' ...
    '\t%g (%.0f%%) are press-inhibited AND lick-inhibited;\n' ...
    '\t%g (%.0f%%) are press-excited AND lick-inhibited;\n' ...
    '\t%g (%.0f%%) are press-inhibited AND lick-excited;\n'], ...
    nTotal, length(eu), 100*nTotal/length(eu), ...
    nnz(c.isPressUp & c.isLickUp), 100*nnz(c.isPressUp & c.isLickUp)/nTotal, ...
    nnz(c.isPressDown & c.isLickDown), 100*nnz(c.isPressDown & c.isLickDown)/nTotal, ...
    nnz(c.isPressUp & c.isLickDown), 100*nnz(c.isPressUp & c.isLickDown)/nTotal, ...
    nnz(c.isPressDown & c.isLickUp), 100*nnz(c.isPressDown & c.isLickUp)/nTotal)   

nTotal = nnz(c.isPressMedResponsive & c.isPressLatResponsive);
fprintf(1, ['%g/%g (%.0f%%) units (press-med and press-lat responsive):\n' ...
    '\t%g (%.0f%%) are press-med-excited AND press-lat-excited;\n' ...
    '\t%g (%.0f%%) are press-med-inhibited AND press-lat-inhibited;\n' ...
    '\t%g (%.0f%%) are press-med-excited AND press-lat-inhibited;\n' ...
    '\t%g (%.0f%%) are press-med-inhibited AND press-lat-excited;\n'], ...
    nTotal, length(eu), 100*nTotal/length(eu), ...
    nnz(c.isPressMedUp & c.isPressLatUp), 100*nnz(c.isPressMedUp & c.isPressLatUp)/nTotal, ...
    nnz(c.isPressMedDown & c.isPressLatDown), 100*nnz(c.isPressMedDown & c.isPressLatDown)/nTotal, ...
    nnz(c.isPressMedUp & c.isPressLatDown), 100*nnz(c.isPressMedUp & c.isPressLatDown)/nTotal, ...
    nnz(c.isPressMedDown & c.isPressLatUp), 100*nnz(c.isPressMedDown & c.isPressLatUp)/nTotal)   

fprintf('Calculate: Of %i: %i (%i%%) showed modulation for BOTH, %i (%i%%) showed modulation for lick only, %i (%i%%) showed modulation for reach only, %i (%i%%) for neither.', ...
    length(eu), ...
    nnz(c.isPressResponsive & c.isLickResponsive), round(nnz(c.isPressResponsive & c.isLickResponsive)/length(eu)*100), ...
    nnz(c.isLickResponsive & ~c.isPressResponsive), round(nnz(c.isLickResponsive & ~c.isPressResponsive)/length(eu)*100), ...
    nnz(c.isPressResponsive & ~c.isLickResponsive), round(nnz(c.isPressResponsive & ~c.isLickResponsive)/length(eu)*100), ...
    nnz(~c.isPressResponsive & ~c.isLickResponsive), round(nnz(~c.isPressResponsive & ~c.isLickResponsive)/length(eu)*100) ...
    )

clear nTotal

save('C:\SERVER\Units\meta_TwoColor_SC.mat', 'c', 'p', 'boot', 'rd', 'eta', 'meta')

%% LDA to decode reach vs. lick using opto-tagged neural populations
clear pLDA;
pLDA.window = [-4, 0];
pLDA.baselineWindow = [-4, -1];
pLDA.responseWindow = [-0.2, -0];
pLDA.responseWindow2tgt = [-0.2, -0];
% pLDA.minTrialDuration = 2;
pLDA.minNumUnits = 10;
pLDA.res = 0.1;
pLDA.nBoot = 10000;
pLDA.kFold = 5;

GROUPNAME = {'CoChR', 'ChrimsonR', 'ChrimsonRClean', 'Neither'};
SEL = {c.isStimBlueUpRedNotUpThereforeCoChrMaybe, c.isStimBlueNotUpRedUpThereforeChrimsonMaybe | c.isStimBlueUpRedUpThereforeChrimsonMaybe, c.isStimBlueNotUpRedUpThereforeChrimsonMaybe, c.isStimBlueNotUpRedNotUp};
clear likelihood nUnits goodExpNames isGoodExp
clear SR RESP

for iGroup = 1:length(GROUPNAME)
    groupName = GROUPNAME{iGroup};
    % Select sessions with enough press/lick trials, enough units
    allUnitIndices = find(SEL{iGroup});
    [goodExpNames.(groupName), ~, ic] = unique(string({eu(allUnitIndices).ExpName}));
    nUnits.(groupName) = histcounts(ic, 1:max(ic)+1);
    goodExpNames.(groupName) = goodExpNames.(groupName)(nUnits.(groupName) >= pLDA.minNumUnits);
    isGoodExp.(groupName) = nUnits.(groupName) >= pLDA.minNumUnits;
    nUnits.(groupName) = nUnits.(groupName)(nUnits.(groupName) >= pLDA.minNumUnits);
    
    fprintf('Selected %i sessions with >%i units.\n', length(goodExpNames.(groupName)), pLDA.minNumUnits)
    fprintf('\tUnits per session: %s\n', num2str(nUnits.(groupName)));
    
    % Extract data
    t = pLDA.window(1):pLDA.res:pLDA.window(2);
    t = (t(1:end-1) + t(2:end))./2;
    clear sr resp
    sr(length(goodExpNames.(groupName))) = struct(press=[], lick=[], t=[]);
    resp(length(goodExpNames.(groupName))) = struct(press=[], lick=[], baseline=[]);
    for iExp = 1:length(goodExpNames.(groupName))
        unitIndicesInExp = allUnitIndices(strcmpi(string({eu(allUnitIndices).ExpName}), goodExpNames.(groupName){iExp}));
        unitIndicesInExp = unitIndicesInExp(:)';
        pressTrials = eu(unitIndicesInExp(1)).Trials.Press;
        lickTrials = eu(unitIndicesInExp(1)).Trials.Lick;
        sr(iExp).press = NaN(length(pressTrials), length(t), length(unitIndicesInExp));
        sr(iExp).lick = NaN(length(lickTrials), length(t), length(unitIndicesInExp));
        for i = 1:length(unitIndicesInExp)
            iEu = unitIndicesInExp(i);
            [sr(iExp).press(:, :, i), ~] = eu(iEu).getTrialAlignedData('count', pLDA.window, 'press', trials=pressTrials, alignTo='stop', resolution=pLDA.res, includeInvalid=false);
            [sr(iExp).lick(:, :, i), ~] = eu(iEu).getTrialAlignedData('count', pLDA.window, 'lick', trials=lickTrials, alignTo='stop', resolution=pLDA.res, includeInvalid=false);
        end
        sr(iExp).press = sr(iExp).press ./ pLDA.res;    
        sr(iExp).lick = sr(iExp).lick ./ pLDA.res;
        sr(iExp).t = t;
    
        selT = t>=pLDA.baselineWindow(1) & t<=pLDA.baselineWindow(2);
        
        mu = mean(sr(iExp).press(:, selT, :), [1, 2], 'omitnan');
        sd = std(sr(iExp).press(:, selT, :), 0, [1, 2], 'omitnan');
        sr(iExp).press = (sr(iExp).press - mu) ./ sd;
    
        mu = mean(sr(iExp).lick(:, selT, :), [1, 2], 'omitnan');
        sd = std(sr(iExp).lick(:, selT, :), 0, [1, 2], 'omitnan');
        sr(iExp).lick = (sr(iExp).lick - mu) ./ sd;
    
    
        resp(iExp).press = squeeze(mean(sr(iExp).press(:, t>=pLDA.responseWindow(1) & t<=pLDA.responseWindow(2), :), 2, 'omitnan'));
        resp(iExp).lick = squeeze(mean(sr(iExp).lick(:, t>=pLDA.responseWindow(1) & t<=pLDA.responseWindow(2), :), 2, 'omitnan'));
        resp(iExp).baseline = cat(1, ...
            squeeze(mean(sr(iExp).press(:, t>=pLDA.baselineWindow(1) & t<=pLDA.baselineWindow(2), :), 2, 'omitnan')), ...
            squeeze(mean(sr(iExp).lick(:, t>=pLDA.baselineWindow(1) & t<=pLDA.baselineWindow(2), :), 2, 'omitnan')) ...
            );
    
        % Remove NaNs
        selPress = all(~isnan(resp(iExp).press), 2);
        selLick = all(~isnan(resp(iExp).lick), 2);
        selBaseline = all(~isnan(resp(iExp).baseline), 2);
    
        resp(iExp).press = resp(iExp).press(selPress, :);
        resp(iExp).lick = resp(iExp).lick(selLick, :);
        resp(iExp).baseline = resp(iExp).baseline(selBaseline, :);
    
        sr(iExp).press = sr(iExp).press(selPress, :, :);
        sr(iExp).lick = sr(iExp).lick(selLick, :, :);

        SR.(groupName)(iExp) = sr(iExp);
        RESP.(groupName)(iExp) = resp(iExp);
    end
    
    clear t iExp unitIndicesInExp pressTrials lickTrials i iEu selT mu sd selPress selLick selBaseline
    
    % Fit LDA
    likelihood.(groupName)(length(sr)) = struct(press=[], lick=[], baseline=[]);
    t = sr(1).t;
    
    for iExp = 1:length(sr)
        nPress = size(resp(iExp).press, 1);
        nLick = size(resp(iExp).lick, 1);
        nBaseline = size(resp(iExp).baseline, 1);
        nTrials = nPress + nLick;
    
        % Fit model using response window
        X = vertcat(resp(iExp).press, resp(iExp).lick, resp(iExp).baseline);
        Y = vertcat(repmat("press", [nPress, 1]), repmat("lick", [nLick, 1]), repmat("baseline", [nBaseline, 1]));
        mdl = fitcdiscr(X, Y, Prior='empirical', CrossVal='on', KFold=pLDA.kFold);
    
        % Predict full timecourse using fitted model
        likelihood.(groupName)(iExp).press = NaN(nTrials, length(t));
        likelihood.(groupName)(iExp).lick = NaN(nTrials, length(t));
        likelihood.(groupName)(iExp).baseline = NaN(nTrials, length(t));
        likelihood.(groupName)(iExp).trueLabel = Y;
        likelihood.(groupName)(iExp).t = t;
    
        assert(length(mdl.ClassNames) == 3)
        [~, iClass] = ismember(["press", "lick", "baseline"], mdl.ClassNames);
        for i = 1:length(t)
            XPress = squeeze(sr(iExp).press(:, i, :));
            XLick = squeeze(sr(iExp).lick(:, i, :));
            XAll = vertcat(XPress, XLick);
    
            score = NaN(nPress + nLick, 3);
    
            for iFold = 1:pLDA.kFold
                testIndices = mdl.Partition.test(iFold);
                testIndices = testIndices(1:nTrials)';
                [~, score(testIndices, 1:3)] = mdl.Trained{iFold}.predict(XAll(testIndices, :));
            end
    
            likelihood.(groupName)(iExp).press(:, i) = score(:, iClass(1));
            likelihood.(groupName)(iExp).lick(:, i) = score(:, iClass(2));
            likelihood.(groupName)(iExp).baseline(:, i) = score(:, iClass(3));
        end
        likelihood.(groupName)(iExp).df = likelihood.(groupName)(iExp).press - likelihood.(groupName)(iExp).lick; % df = press - lick
    end
    clear iExp nPress nLick nTrials X Y mdl i t XPress XLick score iFold testIndices trainIndices
end

%% LDA with bootstrap to get confidence intervals
% Fit LDA
t = sr(1).t;
pLDA.nBoot = 10000;
    
rng(42);

pool = parpool();
for iGroup = 1:length(GROUPNAME)
    groupName = GROUPNAME{iGroup};
    dfBoot = arrayfun(@(llh) llh.df, likelihood.(groupName), UniformOutput=false);
    dfBoot = cat(1, dfBoot{:});
    dfBoot = NaN([size(dfBoot), pLDA.nBoot]);
    pressBoot = dfBoot;
    lickBoot = dfBoot;
    resp = RESP.(groupName);
    sr = SR.(groupName);

    parfor iBoot = 1:pLDA.nBoot
        % fprintf('%i\n', iBoot);
        df = cell(length(sr), 1);
        press = df;
        lick = df;
        for iExp = 1:length(sr)
            nPress = size(resp(iExp).press, 1);
            nLick = size(resp(iExp).lick, 1);
            nBaseline = size(resp(iExp).baseline, 1);
            nTrials = nPress + nLick;
        
            % Fit model using response window
            X = vertcat(resp(iExp).press, resp(iExp).lick);
            X = X(randperm(size(X, 1)), :);
            X = vertcat(X, resp(iExp).baseline);
            Y = vertcat(repmat("press", [nPress, 1]), repmat("lick", [nLick, 1]), repmat("baseline", [nBaseline, 1]));
            mdl = fitcdiscr(X, Y, Prior='empirical', CrossVal='on', KFold=pLDA.kFold);
        
            % Predict full timecourse using fitted model
            press{iExp} = NaN(nTrials, length(t));
            lick{iExp} = NaN(nTrials, length(t));
            assert(length(mdl.ClassNames) == 3)
            [~, iClass] = ismember(["press", "lick", "baseline"], mdl.ClassNames);
            for i = 1:length(t)
                XPress = squeeze(sr(iExp).press(:, i, :));
                XLick = squeeze(sr(iExp).lick(:, i, :));
                XAll = vertcat(XPress, XLick);
    
                score = NaN(nPress + nLick, 3);
        
                for iFold = 1:pLDA.kFold
                    testIndices = mdl.Partition.test(iFold);
                    testIndices = testIndices(1:nTrials)';
                    [~, score(testIndices, 1:3)] = mdl.Trained{iFold}.predict(XAll(testIndices, :));
                end
    
                press{iExp}(:, i) = score(:, iClass(1));
                lick{iExp}(:, i) = score(:, iClass(2));
            end
            df{iExp} = press{iExp} - lick{iExp}; % df = press - lick
        end
        dfBoot(:, :, iBoot) = cat(1, df{:});
        pressBoot(:, :, iBoot) = cat(1, press{:});
        lickBoot(:, :, iBoot) = cat(1, lick{:});
    end
    

    DFBOOT.(groupName) = dfBoot;
    PRESSBOOT.(groupName) = pressBoot;
    LICKBOOT.(groupName) = lickBoot;
    clear iExp nPress nLick nTrials X Y mdl i press lick XPress XLick score likelihoodBoot df iBoot
end

delete(pool)
clear pool iGroup groupName dfBoot pressBoot lickBoot

save('C:\SERVER\Units\lda_twocolor_optrode_SC_fullBootData_20250602.mat', 'DFBOOT', 'PRESSBOOT', 'LICKBOOT', 'pLDA', 'RESP', 'SR', '-v7.3')

% Quick summary (99% CI, mean) of bootstrap for lick vs reach
clear dfBootStats pressBootStats lickBootStats;
clear DFBOOTSTATS PRESSBOOTSTATS LICKBOOTSTATS;

for iGroup = 1:length(GROUPNAME)
    groupName = GROUPNAME{iGroup};
    dfBoot = DFBOOT.(groupName);
    pressBoot = PRESSBOOT.(groupName);
    lickBoot = LICKBOOT.(groupName);

    sr = SR.(groupName);
    resp = RESP.(groupName);

    Y = cell(length(sr), 1);
    for iExp = 1:length(sr)
        nPress = size(resp(iExp).press, 1);
        nLick = size(resp(iExp).lick, 1);
        Y{iExp} = vertcat(repmat("press", [nPress, 1]), repmat("lick", [nLick, 1]));
    end
    Y = cat(1, Y{:});
    
    isPress = Y == "press";
    isLick = Y == "lick";
    
    dfBootStats.press.X = transpose(squeeze(mean(dfBoot(isPress, :, :), 1, 'omitnan')));
    dfBootStats.press.mu = mean(dfBootStats.press.X, 1, 'omitnan');
    dfBootStats.press.ci = quantile(dfBootStats.press.X, [0.005, 0.995], 1);
    
    dfBootStats.lick.X = transpose(squeeze(mean(dfBoot(isLick, :, :), 1, 'omitnan')));
    dfBootStats.lick.mu = mean(dfBootStats.lick.X, 1, 'omitnan');
    dfBootStats.lick.ci = quantile(dfBootStats.lick.X, [0.005, 0.995], 1);
    
    dfBootStats.all.X = transpose(squeeze(mean(dfBoot, 1, 'omitnan')));
    dfBootStats.all.mu = mean(dfBootStats.all.X, 1, 'omitnan');
    dfBootStats.all.ci = quantile(dfBootStats.all.X, [0.005, 0.995], 1);

    pressBootStats.press.X = transpose(squeeze(mean(pressBoot(isPress, :, :), 1, 'omitnan')));
    pressBootStats.press.mu = mean(pressBootStats.press.X, 1, 'omitnan');
    pressBootStats.press.ci = quantile(pressBootStats.press.X, [0.005, 0.995], 1);

    pressBootStats.lick.X = transpose(squeeze(mean(pressBoot(isLick, :, :), 1, 'omitnan')));
    pressBootStats.lick.mu = mean(pressBootStats.lick.X, 1, 'omitnan');
    pressBootStats.lick.ci = quantile(pressBootStats.lick.X, [0.005, 0.995], 1);

    pressBootStats.all.X = transpose(squeeze(mean(pressBoot, 1, 'omitnan')));
    pressBootStats.all.mu = mean(pressBootStats.all.X, 1, 'omitnan');
    pressBootStats.all.ci = quantile(pressBootStats.all.X, [0.005, 0.995], 1);

    lickBootStats.press.X = transpose(squeeze(mean(lickBoot(isPress, :, :), 1, 'omitnan')));
    lickBootStats.press.mu = mean(lickBootStats.press.X, 1, 'omitnan');
    lickBootStats.press.ci = quantile(lickBootStats.press.X, [0.005, 0.995], 1);

    lickBootStats.lick.X = transpose(squeeze(mean(lickBoot(isLick, :, :), 1, 'omitnan')));
    lickBootStats.lick.mu = mean(lickBootStats.lick.X, 1, 'omitnan');
    lickBootStats.lick.ci = quantile(lickBootStats.lick.X, [0.005, 0.995], 1);

    lickBootStats.all.X = transpose(squeeze(mean(lickBoot, 1, 'omitnan')));
    lickBootStats.all.mu = mean(lickBootStats.all.X, 1, 'omitnan');
    lickBootStats.all.ci = quantile(lickBootStats.all.X, [0.005, 0.995], 1);

    DFBOOTSTATS.(groupName) = dfBootStats;
    PRESSBOOTSTATS.(groupName) = pressBootStats;
    LICKBOOTSTATS.(groupName) = lickBootStats;
end

clear iGroup Y iExp nPress nLick isPress isLick dfBoot dfBootStats

t = SR.(GROUPNAME{1}).t;

save('C:\SERVER\Units\lda_twocolor_optrode_SC_20250602.mat', 'pLDA', 'likelihood', 'SR', 'RESP', 't', 'goodExpNames', 'nUnits', 'DFBOOTSTATS', 'PRESSBOOTSTATS', 'LICKBOOTSTATS')
%% Plot results (individual sessions)
close all
xl = [-0.5, 0];
fig = figure(Units='inches', Position=[0.2 0.2 10 10]);
tlp = tiledlayout(fig, length(GROUPNAME), 1, TileSpacing='tight', Padding='tight');

for iGroup = 1:length(GROUPNAME)
    groupName = GROUPNAME{iGroup};
    t = likelihood.(groupName)(1).t;
    tl = tiledlayout(tlp, 1, 4, TileSpacing='tight', Padding='tight');
    tl.Layout.Tile = iGroup;
    iExp = 0;
    for i = 1:length(isGoodExp.(groupName))
        ax = nexttile(tl);

        if ~isGoodExp.(groupName)(i)
            ax.Visible = false;
            continue;
        end

        iExp = iExp + 1;

        isPress = likelihood.(groupName)(iExp).trueLabel == "press";
        isLick = likelihood.(groupName)(iExp).trueLabel == "lick";

        hold(ax, 'on')
        h = gobjects(2, 1);
        X.(groupName)(iExp).truePress = mean(likelihood.(groupName)(iExp).press(isPress, :), 1, 'omitnan');
        X.(groupName)(iExp).trueLick = mean(likelihood.(groupName)(iExp).lick(isLick, :), 1, 'omitnan');
        X.(groupName)(iExp).falsePress = mean(likelihood.(groupName)(iExp).press(isLick, :), 1, 'omitnan');
        X.(groupName)(iExp).falseLick = mean(likelihood.(groupName)(iExp).lick(isPress, :), 1, 'omitnan');
        h(1) = plot(ax, t, X.(groupName)(iExp).truePress, 'red', LineWidth=1.5, DisplayName=sprintf('true reach (n=%i)', nnz(isPress)));
        h(2) = plot(ax, t, X.(groupName)(iExp).trueLick, 'blue', LineWidth=1.5, DisplayName=sprintf('true lick (n=%i)', nnz(isLick)));
        h(3) = plot(ax, t, X.(groupName)(iExp).falsePress, 'red', LineStyle='--', LineWidth=1.5, DisplayName=sprintf('false reach (n=%i)', nnz(isPress)));
        h(4) = plot(ax, t, X.(groupName)(iExp).falseLick, 'blue', LineStyle='--', LineWidth=1.5, DisplayName=sprintf('false lick (n=%i)', nnz(isLick)));
        hold(ax, 'off')
        title(ax, sprintf('%s (%i units)', goodExpNames.(groupName)(iExp), nUnits.(groupName)(iExp)), Interpreter='none')

        ylim(ax, [0, 1])
        xticks(ax, [-0.5 -0.2 0])
        xlim(ax, xl)
        % xline(ax, 0, 'k:')
        % yline(ax, 0, 'k:')
    
        legend(ax, h, Location='southwest', Orientation='vertical')
    
        fontsize(ax, 9, 'points')
    end
    ylabel(tl, groupName, FontSize=9)
end
xlabel(tlp, 'Time to contact (s)', FontSize=9)
ylabel(tlp, 'p(reach) - p(lick)', FontSize=9)


% Plot LDA results, average across sessions
for groupName = GROUPNAME
    groupName = groupName{:};
    XMean.(groupName).truePress = mean(vertcat(X.(groupName).truePress), 1, 'omitnan');
    XMean.(groupName).trueLick = mean(vertcat(X.(groupName).trueLick), 1, 'omitnan');
    XMean.(groupName).falsePress = mean(vertcat(X.(groupName).falsePress), 1, 'omitnan');
    XMean.(groupName).falseLick = mean(vertcat(X.(groupName).falseLick), 1, 'omitnan');
end

fig = figure(Units='inches', Position=[0.2 0.2 8 3]);
tl = tiledlayout(fig, 1, length(GROUPNAME), TileSpacing='tight', Padding='loose');

for iGroup = 1:length(GROUPNAME)
    groupName = GROUPNAME{iGroup};
    t = likelihood.(groupName)(1).t;
    dfBootStats = DFBOOTSTATS.(groupName);
    pressBootStats = PRESSBOOTSTATS.(groupName);
    lickBootStats = LICKBOOTSTATS.(groupName);

    ax = nexttile(tl);
    hold(ax, 'on')
    h = gobjects(4, 1);
    h(1) = plot(ax, t, XMean.(groupName).truePress, 'red', LineWidth=1.5, DisplayName='true reach');
    h(2) = plot(ax, t, XMean.(groupName).trueLick, 'blue', LineWidth=1.5, DisplayName='true lick');
    h(3) = plot(ax, t, XMean.(groupName).falsePress, 'red', LineStyle=':', LineWidth=1.5, DisplayName='false reach');
    h(4) = plot(ax, t, XMean.(groupName).falseLick, 'blue', LineStyle=':', LineWidth=1.5, DisplayName='false lick');
    % patch(ax, [t, flip(t)], [dfBootStats.all.ci(1, :), flip(dfBootStats.all.ci(2, :))], 'black', EdgeColor='black', FaceAlpha=0.15, EdgeAlpha=0.5, DisplayName='99% CI');
    patch(ax, [t, flip(t)], [pressBootStats.press.ci(1, :), flip(pressBootStats.press.ci(2, :))], 'red', EdgeColor='red', FaceAlpha=0.2, EdgeAlpha=0.5, DisplayName='99% CI');
    % patch(ax, [t, flip(t)], [pressBootStats.lick.ci(1, :), flip(pressBootStats.lick.ci(2, :))], 'blue', EdgeColor='blue', FaceAlpha=0.2, EdgeAlpha=0.5, DisplayName='99% CI');
    patch(ax, [t, flip(t)], [lickBootStats.lick.ci(1, :), flip(lickBootStats.lick.ci(2, :))], 'blue', EdgeColor='blue', FaceAlpha=0.2, EdgeAlpha=0.5, DisplayName='99% CI');
    patch(ax, [pLDA.responseWindow, flip(pLDA.responseWindow)], [-1, -1, 1, 1], 'yellow', FaceAlpha=0.1, EdgeAlpha=0.5, DisplayName='training');
    hold(ax, 'off')
    ylim(ax, [0, 1])
    xticks(ax, [-0.5 -0.2 0])
    xlim(ax, xl)
    % xline(ax, 0, 'k:')
    % yline(ax, 0, 'k:')
    title(ax, groupName, FontSize=9)

    fontsize(ax, 9, 'points')
end
xlabel(tl, 'Time to contact (s)', FontSize=9)
ylabel(tl, 'p(reach) - p(lick)', FontSize=9)
lgd = legend(ax, h, Location='southwest', Orientation='horizontal');
lgd.Layout.Tile = 'north';



clear fig tl t iExp isPress isLick ax h dfBootStats
% clear allUnitIndices goodExpNames ic nUnits t SEL GROUPNAME iGroup groupName

%% Figure 1. Sort by optotagging (blue/both/red/neither), then sub-divide by move response (decrese/increase/unmodulated)
sel = strcmpi({eu.ExpName}, 'daisy26_20250519');

p.xlim.stimETA = [-50, 50];
p.xlim.moveETA = [-2, 0];

groupVarStim = NaN(length(eu), 1);
groupVarStim(c.isStimBlueUpRedNotUpThereforeCoChrMaybe) = 0; % CoChR
groupVarStim(c.isStimBlueUpRedUpThereforeChrimsonMaybe) = 10; % Lots of Chrimson or Chrimson+CoChR double label?
groupVarStim(c.isStimBlueNotUpRedUpThereforeChrimsonMaybe) = 20; % Low titer of Chrimson or distant Chrimson-cell?
groupVarStim(c.isStimBlueNotUpRedNotUp) = 30; % No opto response


% Subgroups: 
% press-dec & lick-dec, 
% press-dec & ~lick-dec, 
% ~press-dec & lick-dec, 
% ~press-dec & ~lick-dec 
selLastGroup = ~c.isPressResponsive & ~c.isLickResponsive;
groupVarStim(c.isPressDown  & c.isLickDown  & ~selLastGroup) = groupVarStim(c.isPressDown  & c.isLickDown  & ~selLastGroup) + 0;
groupVarStim(c.isPressDown  & ~c.isLickDown & ~selLastGroup) = groupVarStim(c.isPressDown  & ~c.isLickDown & ~selLastGroup) + 1;
groupVarStim(~c.isPressDown & c.isLickDown  & ~selLastGroup) = groupVarStim(~c.isPressDown & c.isLickDown  & ~selLastGroup) + 2;
groupVarStim(~c.isPressDown & ~c.isLickDown & ~selLastGroup) = groupVarStim(~c.isPressDown & ~c.isLickDown & ~selLastGroup) + 3;
groupVarStim(selLastGroup) = groupVarStim(selLastGroup) + 4;

fig = figure(Units="inches", Position=[0 0 10.5, 7.5]);
ax = arrayfun(@(i) subplot(1, 4, i), 1:4);

[~, order] = EphysUnit.plotETA(ax(3), eta.stimBlue, sel, ...
    clim=[-10, 10], xlim=p.xlim.stimETA, hidecolorbar=true, timeUnit='ms', ...
    sortWindow=[0, 0.05], signWindow=[0.01, 0.02], ...
    sortThreshold=6, negativeSortThreshold=3, onsetDirection='forward', onsetPattern=[0 1 1], ...
    sortGroup=groupVarStim(sel));
EphysUnit.plotETA(ax(4), eta.stimRed, sel, order=order, ...
    clim=[-10, 10], xlim=p.xlim.stimETA, hidecolorbar=true, timeUnit='ms');

EphysUnit.plotETA(ax(1), eta.press, sel, order=order, ...
    clim=[-1.5, 1.5], xlim=p.xlim.moveETA, sortWindow=[-3, 0], signWindow=[-0.3, 0], ...
    sortThreshold=0.25, negativeSortThreshold=0.25, hidecolorbar=true);
EphysUnit.plotETA(ax(2), eta.lick, sel, order=order, ...
    clim=[-1.5, 1.5], xlim=p.xlim.moveETA, hidecolorbar=true);

title(ax(1), 'Reach')
title(ax(2), 'Lick')
title(ax(3), sprintf('470-473 nm\n%g-%g uW\n%g-%g ms', 1e6*min(p.stimBluePowers), 1e6*max(p.stimBluePowers), 1e3*min(p.stimBlueDurations), 1e3*max(p.stimBlueDurations)))
title(ax(4), sprintf('635 nm\n%g-%g uW\n%g-%g ms', 1e6*min(p.stimRedPowers), 1e6*max(p.stimRedPowers), 1e3*min(p.stimRedDurations), 1e3*max(p.stimRedDurations)))
for iAx = 1:2
    applyCustomColormap(ax(iAx), [-2, 2], hlim=[0.375, 0, 0, -0.375], llim=[0.125, 0.5, 0.5, 0.25], hpwr=.5, lpwr=1, h0=0.33);
end
for iAx = 3:4
    applyCustomColormap(ax(iAx), [-10, 10], hlim=[0.375, 0, 0, -0.375], llim=[0.125, 0.5, 0.5, 0.25], hpwr=.5, lpwr=1, h0=0.33);
end

xline(ax(3), 0, 'k--')
xline(ax(4), 0, 'k--')

xlabel(ax(1), 'Time to bar contact (s)')
xlabel(ax(2), 'Time to spout contact (s)')
xlabel(ax(3:4), 'Time from opto on (ms)')


N = histcounts(groupVarStim(sel), -5:10:35);
for iAx = 1:4
    yline(ax(iAx), cumsum(N(1:end-1)) + 1, 'k', LineWidth=3);
end

NTen = cumsum([0, N]);
for iTen = 0:3
    N = histcounts(groupVarStim(sel), iTen*10 + [-0.5, 2.5, 3.5, 4.5]);
    for iAx = 1:4
        yline(ax(iAx), NTen(iTen + 1) + cumsum(N) + 1, 'k--', LineWidth=1);
    end
end

copygraphics(fig, BackgroundColor='none', ContentType='vector')

% clear fig ax order iAx N

%%
fig = figure;
ax = subplot(1, 4, 1); hold(ax, 'on')
sel = c.isStimBlueUpRedNotUpThereforeCoChrMaybe;
plot(ax, eta.pressRaw.t, mean(eta.pressRaw.X(sel, :), 1, 'omitnan'), Color='red', LineWidth=1.5, DisplayName='Reach');
plot(ax, eta.lickRaw.t, mean(eta.lickRaw.X(sel, :), 1, 'omitnan'), Color='blue', LineWidth=1.5, DisplayName='Lick');
title(ax, 'Blue')
ylim(ax, [0, 60])


ax = subplot(1, 4, 2); hold(ax, 'on')
sel = c.isStimBlueUpRedUpThereforeChrimsonMaybe;
plot(ax, eta.pressRaw.t, mean(eta.pressRaw.X(sel, :), 1, 'omitnan'), Color='red', LineWidth=1.5, DisplayName='Reach');
plot(ax, eta.lickRaw.t, mean(eta.lickRaw.X(sel, :), 1, 'omitnan'), Color='blue', LineWidth=1.5, DisplayName='Lick');
title(ax, 'Blue & red')
ylim(ax, [0, 60])


ax = subplot(1, 4, 3); hold(ax, 'on')
sel = c.isStimBlueNotUpRedUpThereforeChrimsonMaybe;
plot(ax, eta.pressRaw.t, mean(eta.pressRaw.X(sel, :), 1, 'omitnan'), Color='red', LineWidth=1.5, DisplayName='Reach');
plot(ax, eta.lickRaw.t, mean(eta.lickRaw.X(sel, :), 1, 'omitnan'), Color='blue', LineWidth=1.5, DisplayName='Lick');
title(ax, 'Red only')
ylim(ax, [0, 60])


ax = subplot(1, 4, 4); hold(ax, 'on')
sel = c.isStimBlueNotUpRedNotUp;
plot(ax, eta.pressRaw.t, mean(eta.pressRaw.X(sel, :), 1, 'omitnan'), Color='red', LineWidth=1.5, DisplayName='Reach');
plot(ax, eta.lickRaw.t, mean(eta.lickRaw.X(sel, :), 1, 'omitnan'), Color='blue', LineWidth=1.5, DisplayName='Lick');
title(ax, 'Neither')

ylim(ax, [0, 60])

%%
msr = arrayfun(@(eu) eu.SpikeRateStats.median, eu);
xyl = [-25, 75];
sz = 10;
fig = figure;
ax = axes(fig); hold(ax, 'on')
sel = c.isStimBlueUpRedNotUpThereforeCoChrMaybe;
scatter(ax, meta.pressRaw(sel) - msr(sel)', meta.lickRaw(sel) - msr(sel)', sz, 'blue');
sel = c.isStimBlueNotUpRedUpThereforeChrimsonMaybe;
scatter(ax, meta.pressRaw(sel) - msr(sel)', meta.lickRaw(sel) - msr(sel)', sz, 'red');
sel = c.isStimBlueUpRedUpThereforeChrimsonMaybe;
scatter(ax, meta.pressRaw(sel) - msr(sel)', meta.lickRaw(sel) - msr(sel)', sz, 'black');
xlabel('\DeltaSR_{reach} (sp/s)')
ylabel('\DeltaSR_{lick} (sp/s)')
xline(ax, 0, 'k--')
yline(ax, 0, 'k--')
axis(ax, 'equal')
plot(ax, xyl, xyl, 'k--')
xlim(ax, xyl)
ylim(ax, xyl)


%% Figure 2. Sort by optotagging (blue/both/red/neither), then sub-divide by move response (decrese/increase/unmodulated)
sel = true(size(eu));

p.xlim.stimETA = [-50, 50];
p.xlim.moveETA = [-2, 0];

groupVarStim = NaN(length(eu), 1);
groupVarStim(c.isStimBlueUpRedNotUpThereforeCoChrMaybe) = 0; % CoChR
groupVarStim(c.isStimBlueUpRedUpThereforeChrimsonMaybe) = 10; % Lots of Chrimson or Chrimson+CoChR double label?
groupVarStim(c.isStimBlueNotUpRedUpThereforeChrimsonMaybe) = 20; % Low titer of Chrimson or distant Chrimson-cell?
groupVarStim(c.isStimBlueNotUpRedNotUp) = 30; % No opto response


% Subgroups: 
% press-dec & lick-dec, 
% press-dec & ~lick-dec, 
% ~press-dec & lick-dec, 
% ~press-dec & ~lick-dec 
selLastGroup = ~c.isPressResponsive & ~c.isLickResponsive;
groupVarStim(c.isPressDown  & c.isLickDown  & ~selLastGroup) = groupVarStim(c.isPressDown  & c.isLickDown  & ~selLastGroup) + 0;
groupVarStim(c.isPressDown  & ~c.isLickDown & ~selLastGroup) = groupVarStim(c.isPressDown  & ~c.isLickDown & ~selLastGroup) + 1;
groupVarStim(~c.isPressDown & c.isLickDown  & ~selLastGroup) = groupVarStim(~c.isPressDown & c.isLickDown  & ~selLastGroup) + 2;
groupVarStim(~c.isPressDown & ~c.isLickDown & ~selLastGroup) = groupVarStim(~c.isPressDown & ~c.isLickDown & ~selLastGroup) + 3;
groupVarStim(selLastGroup) = groupVarStim(selLastGroup) + 4;

fig = figure(Units="inches", Position=[0 0 10.5, 7.5]);
ax = arrayfun(@(i) subplot(1, 4, i), 1:4);

[~, order] = EphysUnit.plotETA(ax(2), eta.lick, sel, ...
    clim=[-1.5, 1.5], xlim=p.xlim.moveETA, hidecolorbar=true, ...
    sortWindow=[-3, 0], signWindow=[-0.3, 0], ...
    sortThreshold=0.25, negativeSortThreshold=0.25, ...
    sortGroup=groupVarStim(sel));
EphysUnit.plotETA(ax(1), eta.press, sel, order=order, ...
    clim=[-1.5, 1.5], xlim=p.xlim.moveETA, hidecolorbar=true, ...
    sortWindow=[-3, 0], signWindow=[-0.3, 0], ...
    sortThreshold=0.25, negativeSortThreshold=0.25, ...
    sortGroup=groupVarStim(sel));
EphysUnit.plotETA(ax(3), eta.stimBlue, sel, order=order, ...
    clim=[-10, 10], xlim=p.xlim.stimETA, hidecolorbar=true, timeUnit='ms');
EphysUnit.plotETA(ax(4), eta.stimRed, sel, order=order, ...
    clim=[-10, 10], xlim=p.xlim.stimETA, hidecolorbar=true, timeUnit='ms');

title(ax(1), 'Reach')
title(ax(2), 'Lick')
title(ax(3), sprintf('470-473 nm\n%g-%g uW\n%g-%g ms', 1e6*min(p.stimBluePowers), 1e6*max(p.stimBluePowers), 1e3*min(p.stimBlueDurations), 1e3*max(p.stimBlueDurations)))
title(ax(4), sprintf('635 nm\n%g-%g uW\n%g-%g ms', 1e6*min(p.stimRedPowers), 1e6*max(p.stimRedPowers), 1e3*min(p.stimRedDurations), 1e3*max(p.stimRedDurations)))
for iAx = 1:2
    applyCustomColormap(ax(iAx), [-2, 2], hlim=[0.375, 0, 0, -0.375], llim=[0.125, 0.5, 0.5, 0.25], hpwr=.5, lpwr=1, h0=0.33);
end
for iAx = 3:4
    applyCustomColormap(ax(iAx), [-10, 10], hlim=[0.375, 0, 0, -0.375], llim=[0.125, 0.5, 0.5, 0.25], hpwr=.5, lpwr=1, h0=0.33);
end

xline(ax(3), 0, 'k--')
xline(ax(4), 0, 'k--')

xlabel(ax(1), 'Time to bar contact (s)')
xlabel(ax(2), 'Time to spout contact (s)')
xlabel(ax(3:4), 'Time from opto on (ms)')


N = histcounts(groupVarStim(sel), -5:10:35);
for iAx = 1:4
    yline(ax(iAx), cumsum(N(1:end-1)) + 1, 'k', LineWidth=3);
end

NTen = cumsum([0, N]);
for iTen = 0:3
    N = histcounts(groupVarStim(sel), iTen*10 + [-0.5, 2.5, 3.5, 4.5]);
    for iAx = 1:4
        yline(ax(iAx), NTen(iTen + 1) + cumsum(N) + 1, 'k--', LineWidth=1);
    end
end

copygraphics(fig, BackgroundColor='none', ContentType='vector')

% clear fig ax order iAx N



%% Figure 4. Sort by peri-movement response (lick-inc/both-inc/reach-inc/neither-inc) and then subdivide by optotagging (blue/both/red/neither)
p.xlim.stimETA = [-50, 50];
p.xlim.moveETA = [-2, 0];

% Major groups: 
% press-dec & lick-dec, 
% press-dec & ~lick-dec, 
% ~press-dec & lick-dec, 
% ~press-dec & ~lick-dec 
groupVar = NaN(length(eu), 1);
groupVar(c.isPressResponsive & ~c.isLickResponsive) = 0;
groupVar(c.isPressResponsive & c.isLickResponsive) = 10;
groupVar(~c.isPressResponsive & c.isLickResponsive) = 20;
groupVar(~c.isPressResponsive & ~c.isLickResponsive) = 30;
assert(nnz(isnan(groupVar)) == 0)

% Subgroups: % BlueOn, % Red On, Blue On, % Red On, % Neither On
groupVar(c.isStimBlueUpRedNotUpThereforeCoChrMaybe) = groupVar(c.isStimBlueUpRedNotUpThereforeCoChrMaybe) + 0; % CoChR
groupVar(c.isStimBlueUpRedUpThereforeChrimsonMaybe) = groupVar(c.isStimBlueUpRedUpThereforeChrimsonMaybe) + 1; % Lots of Chrimson or Chrimson+CoChR double label?
groupVar(c.isStimBlueNotUpRedUpThereforeChrimsonMaybe) = groupVar(c.isStimBlueNotUpRedUpThereforeChrimsonMaybe) + 2; % Low titer of Chrimson or distant Chrimson-cell?
groupVar(c.isStimBlueNotUpRedNotUp) = groupVar(c.isStimBlueNotUpRedNotUp) + 3; % No opto response


fig = figure(Units="inches", Position=[0 0 11, 7.5]);
ax = arrayfun(@(i) subplot(1, 6, i), 1:6);

EphysUnit.plotETA(ax(1), eta.press, ...
    clim=[-1.5, 1.5], xlim=p.xlim.moveETA, sortWindow=[-3, 0], signWindow=[-0.3, 0], ...
    sortThreshold=0.25, negativeSortThreshold=0.25, hidecolorbar=true);
EphysUnit.plotETA(ax(2), eta.lick, ...
    clim=[-1.5, 1.5], xlim=p.xlim.moveETA, sortWindow=[-3, 0], signWindow=[-0.3, 0], ...
    sortThreshold=0.25, negativeSortThreshold=0.25, hidecolorbar=true);

[~, order] = EphysUnit.plotETA(ax(3), eta.press, sortGroup=groupVar, ...
    clim=[-1.5, 1.5], xlim=p.xlim.moveETA, sortWindow=[-3, 0], signWindow=[-0.3, 0], ...
    sortThreshold=0.25, negativeSortThreshold=0.25, hidecolorbar=true);
EphysUnit.plotETA(ax(4), eta.lick, order=order, ...
    clim=[-1.5, 1.5], xlim=p.xlim.moveETA, hidecolorbar=true);

EphysUnit.plotETA(ax(5), eta.stimBlue, order=order, ...
    clim=[-10, 10], xlim=p.xlim.stimETA, hidecolorbar=true, timeUnit='ms');
EphysUnit.plotETA(ax(6), eta.stimRed, order=order, ...
    clim=[-10, 10], xlim=p.xlim.stimETA, hidecolorbar=true, timeUnit='ms');

title(ax([1, 3]), 'Reach')
title(ax([2, 4]), 'Lick')
title(ax(5), sprintf('470-473 nm\n%g-%g uW\n%g-%g ms', 1e6*min(p.stimBluePowers), 1e6*max(p.stimBluePowers), 1e3*min(p.stimBlueDurations), 1e3*max(p.stimBlueDurations)))
title(ax(6), sprintf('635 nm\n%g-%g uW\n%g-%g ms', 1e6*min(p.stimRedPowers), 1e6*max(p.stimRedPowers), 1e3*min(p.stimRedDurations), 1e3*max(p.stimRedDurations)))
for iAx = 1:4
    applyCustomColormap(ax(iAx), [-2, 2], hlim=[0.375, 0, 0, -0.375], llim=[0.125, 0.5, 0.5, 0.25], hpwr=.5, lpwr=1, h0=0.33);
end
for iAx = 5:6
    applyCustomColormap(ax(iAx), [-10, 10], hlim=[0.375, 0, 0, -0.375], llim=[0.125, 0.5, 0.5, 0.25], hpwr=.5, lpwr=1, h0=0.33);
    xline(ax(iAx), 0, 'k--')
end


xlabel(ax([1, 3]), 'Time to bar contact (s)')
xlabel(ax([2, 4]), 'Time to spout contact (s)')
xlabel(ax([5, 6]), 'Time from opto on (ms)')


N = histcounts(groupVar, [-5, 5, 15, 25, 35]);
for iAx = 3:6
    yline(ax(iAx), cumsum(N(1:end-1)) + 1, 'k', LineWidth=3);
end

NTen = cumsum([0, N]);
for iTen = 0 : length(NTen)-1
    N = histcounts(groupVar, iTen*10 + [-0.5, 0.5, 1.5, 2.5, 3.5]);
    for iAx = 3:6
        yline(ax(iAx), NTen(iTen + 1) + cumsum(N) + 1, 'k--', LineWidth=1);
    end
end


copygraphics(fig, BackgroundColor='none', ContentType='vector')
