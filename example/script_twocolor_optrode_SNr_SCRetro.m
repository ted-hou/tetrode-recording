

%% Remove duplicates (Slow)
eu = EphysUnit.load('\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\Units\TwoColor_SNr_SCRetro\Batch2', waveforms=false, spikecounts=false, spikerates=false, animalNames={'desmond38'});
%
euAll = eu;

eu = eu.removeMultiUnits(cullZeros=true);
[eu, isDuplicate] = eu.removeDuplicates(0.7);
% eu.save('\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\Units\TwoColor_SNr_SCRetro\SingleUnit_NonDuplicate')

% Remove drift, low spike rate units
clear c
% eu = EphysUnit.load('\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\Units\TwoColor_SNr_SCRetro\SingleUnit_NonDuplicate', waveforms=false, spikecounts=false, spikerates=false);

% Remove drift
c.isDrifting = detectDriftingUnits(eu, smoothWindow=300, tolerance=0.05, spikeRateThreshold=15, includeITI=true);


% Filter by spike rate
msr = arrayfun(@(eu) eu.SpikeRateStats.median, eu);
p.minSpikeRate = 15;
c.isSNr = msr >= p.minSpikeRate;

eu = eu(c.isSNr & ~c.isDrifting);

% eu.save('\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\Units\TwoColor_SNr_SCRetro\SingleUnit_NonDuplicate_NonDrift_SNr')

% Make spontaneous reach/lick trials
% p.minTrialLength = 4;
% for iEu = 1:length(eu)
%     eu(iEu).Trials.Press = eu(iEu).makeTrials('press_spontaneous_clean', minSpontaneousTrialDuration=p.minTrialLength);
%     eu(iEu).Trials.Lick = eu(iEu).makeTrials('lick_spontaneous_clean', minSpontaneousTrialDuration=p.minTrialLength);
% end

eu.save('\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\Units\TwoColor_SNr_SCRetro\SingleUnit_NonDuplicate_NonDrift_SNr')

%% Fix extra red stim event for bad session
for iEu = find(string({eu.ExpName}) == "desmond38_20250403")
    assert(isscalar(iEu))
    assert(eu(iEu).EventTimes.LaserModRedOn(1) - eu(iEu).EventTimes.LaserModRedOff(1) < 1e-6)
    iBadStimOn = find(eu(iEu).EventTimes.StimOn == eu(iEu).EventTimes.LaserModRedOn(1));
    iBadStimOff = find(eu(iEu).EventTimes.StimOff == eu(iEu).EventTimes.LaserModRedOff(1));
    eu(iEu).EventTimes.StimOn(iBadStimOn) = [];
    eu(iEu).EventTimes.StimOff(iBadStimOff) = [];
    eu(iEu).EventTimes.LaserModRedOn(1) = [];
    eu(iEu).EventTimes.LaserModRedOff(1) = [];
    eu(iEu).Trials.Stim = eu(iEu).makeTrials('stim');
end


eu.save('\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\Units\TwoColor_SNr_SCRetro\SingleUnit_NonDuplicate_NonDrift_SNr')

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
    end

    % Lick trials must occur when both LEDs were off
    sel = ~cueTrialsLeft.inTrial(lickTimes) & ~cueTrialsRight.inTrial(lickTimes);
    if nnz(sel) < length(sel)
        eu(iEu).Trials.Lick = lickTrials(sel);
        if ismember(iEu, ia)
            fprintf('Removed %i (of %i) lick trials because they occured when Cue LED(s) was on.\n', length(sel) - nnz(sel), length(sel));
        end
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

clear ia sel implantSide cueTrialsLeft cueTrialsRight isLeft isRight leftOn leftOff pressTrials pressTimes lickTrials lickTimes iEu
eu.save('\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\Units\TwoColor_SNr_SCRetro\SingleUnit_NonDuplicate_NonDrift_SNr')

%%
clear, clc
eu = EphysUnit.load('\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\Units\TwoColor_SNr_SCRetro\SingleUnit_NonDuplicate_NonDrift_SNr', waveforms=false, spikecounts=false, spikerates=false);

% Make Raster
clear rd
rd.stim = eu.getRasterData('stimtwocolor', window=[-0.1, 0.4], durErr=1e-3, shutterDelay=0, photoelectricBlankDuration=0.5e-3);
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

% ETA Stim
p.isiBaselineWindow = [-0.1, 0];
p.stimBluePowers = [25, 50, 100, 500]*1e-6; 
p.stimRedPowers = [25, 50, 100, 500, 2000, 8000, 16000]*1e-6;
p.stimBlueDurations = [10, 20]*1e-3;
p.stimRedDurations = [10, 20]*1e-3;

p.isiWindow = [-0.4, 0.4];
p.isiRes = 1e-3;
p.xlim.stim = [-0.1, 0.3];
p.xlim.move = [-4, 0];
p.path = 'C:\SERVER\Figures\TwoColor_SNr_SCRetro';
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

% Calculate META
clear meta
p.metaWindow = [-0.3, 0];
p.posRespThreshold = 0.5;
p.negRespThreshold = -0.25;
t = eta.press.t;
meta.press = mean(eta.press.X(:, t>=p.metaWindow(1) & t<=p.metaWindow(2)), 2, 'omitnan');
meta.pressMed = mean(eta.pressMed.X(:, t>=p.metaWindow(1) & t<=p.metaWindow(2)), 2, 'omitnan');
meta.pressLat = mean(eta.pressLat.X(:, t>=p.metaWindow(1) & t<=p.metaWindow(2)), 2, 'omitnan');
meta.lick = mean(eta.lick.X(:, t>=p.metaWindow(1) & t<=p.metaWindow(2)), 2, 'omitnan');
t = eta.pressRaw.t;
meta.pressRaw = mean(eta.press.X(:, t>=p.metaWindow(1) & t<=p.metaWindow(2)), 2, 'omitnan');
meta.pressMedRaw = mean(eta.pressMedRaw.X(:, t>=p.metaWindow(1) & t<=p.metaWindow(2)), 2, 'omitnan');
meta.pressLatRaw = mean(eta.pressLatRaw.X(:, t>=p.metaWindow(1) & t<=p.metaWindow(2)), 2, 'omitnan');
meta.lickRaw = mean(eta.lick.X(:, t>=p.metaWindow(1) & t<=p.metaWindow(2)), 2, 'omitnan');
clear t

p.metaWindowStim = [0.005, 0.050];
p.posRespThresholdStim = 2;
p.negRespThresholdStim = -1;

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

% Boot response dir
p.bootAlpha = 0.05;
p.nboot = 100000;
p.metaWindow = [-0.3, 0];
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

%% Report bootstraped movement response direction
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

save('C:\SERVER\Units\meta_TwoColor_SNr_SCRetro.mat', 'c', 'p', 'boot', 'rd', 'eta', 'meta')


%% Combined PEISI and Stim Raster
% Stim Rasters

fig = figure(Units='inches', Position=[0, 0, 12, 8]);
clear layout

layout.w = [3, 4]; % Stim, reach/lick
layout.h = [3, 3, 4]; % Raster, peisi/peth

layout.tl = tiledlayout(fig, sum(layout.h), sum(layout.w), TileSpacing='tight', Padding='tight', TileIndexing='columnmajor');
layout.ax = gobjects(length(layout.h), length(layout.w));
layout.ax(1, 1) = nexttile(layout.tl, [sum(layout.h(1:2)), layout.w(1)]);
layout.ax(2, 1) = nexttile(layout.tl, [layout.h(3), layout.w(1)]);
layout.ax(1, 2) = nexttile(layout.tl, [layout.h(1), layout.w(2)]);
layout.ax(2, 2) = nexttile(layout.tl, [layout.h(2), layout.w(2)]);
layout.ax(3, 2) = nexttile(layout.tl, [layout.h(3), layout.w(2)]);

if ~exist(p.path, 'dir')
    mkdir(p.path)
end

for iEu = 1:length(eu)
    % try
        cla(layout.ax(1:2, 1))
        cla(layout.ax(:, 2))
        % Make Raster
        ax = layout.ax(1, 1);
        EphysUnit.plotRaster(ax, rd.stim(iEu), xlim=p.xlim.stim*1e3, sz=p.rasterSzStim, timeUnit='ms');
        ax.Legend.Location = 'northeast';
        ax.Legend.FontSize = 6;
        title(ax, eu(iEu).getName(), Interpreter="none")

        % Make PE-ISI
        ax = layout.ax(2, 1);
        groups = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power'});
        isi = NaN(length(groups), length(p.isiWindow(1):p.isiRes:p.isiWindow(2)));
        deltaSR = isi;
        for iGrp = 1:length(groups)
            [isi(iGrp, :), t] = eu(iEu).getMeanPEISI('stimtwocolor', groups(iGrp).trials, window=p.isiWindow, resolution=p.isiRes, ...
                photoelectricBlankDuration=0.5e-3);
            deltaSR(iGrp, :) = 1./isi(iGrp, :) - mean(1./isi(iGrp, t<0), 'omitnan');
        end
        imagesc(ax, 1e3*t, [], deltaSR)
        xline(ax, 0)
        xlim(ax, p.xlim.stim*1e3)
        clim(ax, [-50, 50])
        ax.YAxisLocation = 'right';
        colormap(ax, 'turbo')
        h = colorbar(ax, 'westoutside');
        h.Label.String = '\Deltasp/s';
        yticks(ax, 1:length(groups));
        yticklabels(ax, {groups.label})
        xlabel(ax, 'Time from opto onset (ms)')
        title(ax, 'Stim PETH (from ISI)')

        % Make Raster Press
        trialTypes = ["press", "lick"];
        colors = ["red", "blue"];
        for iTrialType = 1:2
            trialType = trialTypes(iTrialType);
            ax = layout.ax(iTrialType, 2);
            EphysUnit.plotRaster(ax, rd.(trialType)(iEu), xlim=p.xlim.move, sz=p.rasterSzMove);
            title(ax, trialType, Interpreter="none")

            ax = layout.ax(3, 2);
            hold(ax, 'on')
            fieldname = sprintf('%sRaw', trialType);
            plot(ax, eta.(fieldname).t, eta.(fieldname).X(iEu, :), LineWidth=1.5, Color=colors(iTrialType), DisplayName=trialTypes(iTrialType))
            xlim(ax, p.xlim.move)
            title(ax, 'Move PETH')
            xlabel(ax, 'Time to bar/spout contact (s)')
            ylabel(ax, 'sp/s')
            legend(ax, Location='northwest')
        end

        print(fig, sprintf('%s\\%s.png', p.path, eu(iEu).getName()), '-dpng')
    % end
end

clear fig layout iEu ax groups isi deltaSR iGrp h trialTypes colors iTrialType trialType

%% Figure 1. Sort by optotagging (blue/both/red/neither), then sub-divide by move response (decrese/increase/unmodulated)
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

[~, order] = EphysUnit.plotETA(ax(3), eta.stimBlue, ...
    clim=[-10, 10], xlim=p.xlim.stimETA, hidecolorbar=true, timeUnit='ms', ...
    sortWindow=[0, 0.05], signWindow=[0.01, 0.02], ...
    sortThreshold=6, negativeSortThreshold=3, onsetDirection='forward', onsetPattern=[0 1 1], ...
    sortGroup=groupVarStim);
EphysUnit.plotETA(ax(4), eta.stimRed, order=order, ...
    clim=[-10, 10], xlim=p.xlim.stimETA, hidecolorbar=true, timeUnit='ms');

EphysUnit.plotETA(ax(1), eta.press, order=order, ...
    clim=[-1.5, 1.5], xlim=p.xlim.moveETA, sortWindow=[-3, 0], signWindow=[-0.3, 0], ...
    sortThreshold=0.25, negativeSortThreshold=0.25, hidecolorbar=true);
EphysUnit.plotETA(ax(2), eta.lick, order=order, ...
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


N = histcounts(groupVarStim, -5:10:35);
for iAx = 1:4
    yline(ax(iAx), cumsum(N(1:end-1)) + 1, 'k', LineWidth=3);
end

NTen = cumsum([0, N]);
for iTen = 0:3
    N = histcounts(groupVarStim, iTen*10 + [-0.5, 2.5, 3.5, 4.5])
    for iAx = 1:4
        yline(ax(iAx), NTen(iTen + 1) + cumsum(N) + 1, 'k--', LineWidth=1);
    end
end

copygraphics(fig, BackgroundColor='none', ContentType='vector')

% clear fig ax order iAx N


%% Figure2. Sort by peri-movement response (incongruent/congruent)
p.xlim.stimETA = [-50, 50];
p.xlim.moveETA = [-2, 0];

groupVar = NaN(length(eu), 1);
groupVar(c.isPressUnresponsiveButDown & c.isLickUp) = 0;
groupVar(c.isPressUnresponsiveButDown & c.isLickUnresponsiveButUp) = 0;
groupVar(c.isPressDown & c.isLickUp) = 0;
groupVar(c.isPressDown & c.isLickUnresponsiveButUp) = 0;
groupVar(c.isPressUp & c.isLickUnresponsiveButDown) = 1;
groupVar(c.isPressUp & c.isLickDown) = 1;
groupVar(c.isPressUnresponsiveButUp & c.isLickUnresponsiveButDown) = 1;
groupVar(c.isPressUnresponsiveButUp & c.isLickDown) = 1;
groupVar(c.isPressUnresponsiveButDown & c.isLickUnresponsiveButDown) = 2;
groupVar(c.isPressUnresponsiveButDown & c.isLickDown) = 2;
groupVar(c.isPressDown & c.isLickUnresponsiveButDown) = 2;
groupVar(c.isPressDown & c.isLickDown) = 2;
groupVar(c.isPressUp & c.isLickUp) = 3;
groupVar(c.isPressUp & c.isLickUnresponsiveButUp) = 3;
groupVar(c.isPressUnresponsiveButUp & c.isLickUp) = 3;
groupVar(c.isPressUnresponsiveButUp & c.isLickUnresponsiveButUp) = 3;


fig = figure(Units="inches", Position=[0 0 12.5, 7.5]);
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


N = histcounts(groupVar, -0.5:2:3.5);
for iAx = 3:6
    yline(ax(iAx), cumsum(N(1:end-1)) + 1, 'k', LineWidth=3);
end

copygraphics(fig, BackgroundColor='none', ContentType='vector')


%% Figure 3. Sort by peri-movement response (incongruent/congruent) and then subdivide by optotagging (blue/both/red/neither)
p.xlim.stimETA = [-50, 50];
p.xlim.moveETA = [-2, 0];

groupVar = NaN(length(eu), 1);
groupVar(c.isPressUnresponsiveButDown & c.isLickUp) = 0;
groupVar(c.isPressUnresponsiveButDown & c.isLickUnresponsiveButUp) = 0;
groupVar(c.isPressDown & c.isLickUp) = 0;
groupVar(c.isPressDown & c.isLickUnresponsiveButUp) = 0;
groupVar(c.isPressUp & c.isLickUnresponsiveButDown) = 0;
groupVar(c.isPressUp & c.isLickDown) = 0;
groupVar(c.isPressUnresponsiveButUp & c.isLickUnresponsiveButDown) = 0;
groupVar(c.isPressUnresponsiveButUp & c.isLickDown) = 0;
groupVar(c.isPressUnresponsiveButDown & c.isLickUnresponsiveButDown) = 1;
groupVar(c.isPressUnresponsiveButDown & c.isLickDown) = 1;
groupVar(c.isPressDown & c.isLickUnresponsiveButDown) = 1;
groupVar(c.isPressDown & c.isLickDown) = 1;
groupVar(c.isPressUp & c.isLickUp) = 1;
groupVar(c.isPressUp & c.isLickUnresponsiveButUp) = 1;
groupVar(c.isPressUnresponsiveButUp & c.isLickUp) = 1;
groupVar(c.isPressUnresponsiveButUp & c.isLickUnresponsiveButUp) = 1;

groupVar = groupVar * 10;

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


N = histcounts(groupVar, [-5, 5, 15]);
for iAx = 3:6
    yline(ax(iAx), cumsum(N(1:end-1)) + 1, 'k', LineWidth=3);
end

NTen = cumsum([0, N]);
for iTen = 0 : 1
    N = histcounts(groupVar, iTen*10 + (-0.5:1:3.5));
    for iAx = 3:6
        yline(ax(iAx), NTen(iTen + 1) + cumsum(N) + 1, 'k--', LineWidth=1);
    end
end

copygraphics(fig, BackgroundColor='none', ContentType='vector')


%% Figure 4. Sort by peri-movement response (1+dec/1+dec/neither-dec) and then subdivide by optotagging (blue/both/red/neither)
p.xlim.stimETA = [-50, 50];
p.xlim.moveETA = [-2, 0];

% Major groups: 
% press-dec & lick-dec, 
% press-dec & ~lick-dec, 
% ~press-dec & lick-dec, 
% ~press-dec & ~lick-dec 
groupVar = NaN(length(eu), 1);
groupVar(c.isPressDown | c.isLickDown) = 0;
groupVar(~(c.isPressDown | c.isLickDown) & (c.isPressUp | c.isLickUp)) = 10;
groupVar(~c.isPressResponsive & ~c.isLickResponsive) = 20;
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


N = histcounts(groupVar, [-5, 5, 15, 25]);
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
