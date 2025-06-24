function [eu, rd, eta, meta, p, c, boot] = read_SNr_SCRetro(varargin)
parser = inputParser();
parser.addParameter('eu', [])
parser.addParameter('metaPath', '')
parser.addParameter('readBootButRecalculateMeta', true, @islogical)
parser.addParameter('metaSavePath', 'C:\\SERVER\\Units\\meta_TwoColor_SNr_SCRetro_%s.mat')
parser.addParameter('stimBluePowers', [100, 500, 2000]*1e-6, @isnumeric)
parser.addParameter('stimRedPowers', [2000, 8000, 16000]*1e-6, @isnumeric)
parser.addParameter('stimBluePowersMustContain', ["*1e6==100", "*1e6==500", "*1e6==2000"], @isstring) % Blue stim powers must contain each of these elements
parser.addParameter('stimRedPowersMustContain', ["*1e6==2000", "*1e6>=8000"], @isstring) % Red stim powers must contain each of these elements, e.g. ["*1e-6==2000", "*1e-6>=8000"] means must have 2mW train and 8(and above)mW train
parser.addParameter('special_groupRedPowersOver8000', false, @islogical)
parser.addParameter('stimBlueDurations', [10, 20]*1e-3, @isnumeric)
parser.addParameter('stimRedDurations', [10, 20]*1e-3, @isnumeric)
parser.parse(varargin{:})
metaPath = parser.Results.metaPath;
readBootButRecalculateMeta = parser.Results.readBootButRecalculateMeta;
metaSavePath = sprintf(parser.Results.metaSavePath, char(datetime('today'), 'yyyyMMdd'));
stimBluePowersMustContain = parser.Results.stimBluePowersMustContain;
stimRedPowersMustContain = parser.Results.stimRedPowersMustContain;

if ~isempty(parser.Results.eu)
    eu = parser.Results.eu;
else
    eu = EphysUnit.load('\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\Units\TwoColor_SNr_SCRetro\SingleUnit_NonDuplicate_NonDrift_SNr', waveforms=false, spikecounts=false, spikerates=false);
end

if ~isempty(metaPath) % 'C:\SERVER\Units\meta_TwoColor_SNr_SCRetro.mat'
    load(metaPath)
    if ~readBootButRecalculateMeta
        return
    end
end

clear rd
clear eta

% ETA Stim
p.isiBaselineWindow = [-0.1, 0];
p.stimBluePowers = parser.Results.stimBluePowers; 
p.stimRedPowers = parser.Results.stimRedPowers;
p.stimBlueDurations = parser.Results.stimBlueDurations;
p.stimRedDurations = parser.Results.stimRedDurations;
p.stimBluePowersMustContain = stimBluePowersMustContain;
p.stimRedPowersMustContain = stimRedPowersMustContain;

p.isiWindow = [-0.4, 0.4];
p.isiRes = 1e-3;
p.xlim.stim = [-0.1, 0.3];
p.xlim.move = [-4, 0];
p.rasterSzStim = 1;
p.rasterSzMove = 1;

% First filter out units without requested stim conditions
[~, expToEuIndex, euToExpIndex] = unique(string({eu.ExpName}));
selExp = true(length(expToEuIndex), 1);
if ~isempty(stimBluePowersMustContain)
    for iExp = 1:length(expToEuIndex)
        thisEu = eu(expToEuIndex(iExp));
        groups = thisEu.groupTwoColorStimTrials({'wavelength', 'power', 'duration'}, selectBy=struct(power=p.stimBluePowers, duration=p.stimBlueDurations, location=[], wavelength=[470, 473]));
        powers = [groups.power];
        for iCond = 1:length(stimBluePowersMustContain)
            found = any(eval(sprintf('powers%s', stimBluePowersMustContain(iCond))));
            if ~found
                fprintf("Exp %s has been exluded, has blue powers [%s], needs [%s].\n", thisEu.ExpName, string(powers).join(", "), stimBluePowersMustContain.join(", "));
                break
            end
        end
        selExp(iExp) = selExp(iExp) && found;
    end
end

if ~isempty(stimRedPowersMustContain)
    for iExp = 1:length(expToEuIndex)
        thisEu = eu(expToEuIndex(iExp));
        groups = thisEu.groupTwoColorStimTrials({'wavelength', 'power', 'duration'}, selectBy=struct(power=p.stimRedPowers, duration=p.stimRedDurations, location=[], wavelength=[590, 635]));
        powers = [groups.power];
        for iCond = 1:length(stimRedPowersMustContain)
            found = any(eval(sprintf('powers%s', stimRedPowersMustContain(iCond))));
            if ~found
                fprintf("Exp %s has been exluded, has red powers [%s], needs [%s].\n", thisEu.ExpName, string(powers).join(", "), stimRedPowersMustContain.join(", "));
                break
            end
        end
        selExp(iExp) = selExp(iExp) && found;
    end
end
selEu = ismember(euToExpIndex, find(selExp));
eu = eu(selEu);


%%
X = cell(length(p.stimBluePowers) + 1, 1);
for i = 1:length(X)
    X{i} = cell(length(eu), 1);
end
XBlue = dictionary([string(p.stimBluePowers*1e6), "all"]', X);


X = cell(length(p.stimRedPowers) + 1, 1);
for i = 1:length(X)
    X{i} = cell(length(eu), 1);
end
% assert(parser.Results.special_groupRedPowersOver8000, "I do not have time to write proper code. Please fix this for expandability")
redPowers = p.stimRedPowers;
redPowers(redPowers >= 8000e-6) = [];
X = X(1:end-1);
XRed = dictionary([string(redPowers*1e6), "8000+", "all"]', X);


%%

for iEu = 1:length(eu)
    groupsBlue = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'duration'}, selectBy=struct(power=[100, 500, 2000]*1e-6, duration=p.stimBlueDurations, location=[], wavelength=[470, 473]));
    groupsRed = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'duration'}, selectBy=struct(power=[2000, 8000, 16000]*1e-6, duration=p.stimRedDurations, location=[], wavelength=635));


    for iGrp = 1:length(groupsBlue) + 1
        if iGrp <= length(groupsBlue)
            trials = groupsBlue(iGrp).trials;
            key = string(groupsBlue(iGrp).power*1e6);
        else
            trials = [groupsBlue.trials];
            key = "all";
        end
        
        [isi, t] = eu(iEu).getMeanPEISI('stimtwocolor', trials, window=p.isiWindow, resolution=p.isiRes, photoelectricBlankDuration=0.5e-3);
        selBaseline = t>p.isiBaselineWindow(1) & t<p.isiBaselineWindow(2);
        normSR = (1./isi - mean(1./isi(:, selBaseline), 'omitnan')) ./ std(1./isi(:, selBaseline), 0, 2, 'omitnan');
        X = XBlue(key);
        if length(X) == 1
            X = X{1};
        end
        X{iEu} = normSR;
        XBlue(key) = {X};
    end

    assert(nnz(ismember([8000, 16000]*1e-6, [groupsRed.power])) == 1, "I do not have time to write proper code. Please fix this for expandability")
    for iGrp = 1:length(groupsRed) + 1
        if iGrp <= length(groupsRed)
            trials = groupsRed(iGrp).trials;
            if groupsRed(iGrp).power*1e6 < 8000
                key = string(groupsRed(iGrp).power*1e6);
            else
                key = "8000+";
            end
        else
            trials = [groupsRed.trials];
            key = "all";
        end
        
        [isi, t] = eu(iEu).getMeanPEISI('stimtwocolor', trials, window=p.isiWindow, resolution=p.isiRes, photoelectricBlankDuration=0.5e-3);
        selBaseline = t>p.isiBaselineWindow(1) & t<p.isiBaselineWindow(2);
        normSR = (1./isi - mean(1./isi(:, selBaseline), 'omitnan')) ./ std(1./isi(:, selBaseline), 0, 2, 'omitnan');
        X = XRed(key);
        if length(X) == 1
            X = X{1};
        end
        X{iEu} = normSR;
        XRed(key) = {X};
    end

    % fprintf('nan=%i, nan=%i\n', nnz(isnan(XBlue{iEu})), nnz(isnan(XRed{iEu})))
end
%%
eta.stimBlue = configureDictionary("string", "struct");
eta.stimRed = configureDictionary("string", "struct");

for key = XBlue.keys
    X = XBlue(key);
    X = X{1};
    eta.stimBlue(key) = struct(X=cat(1, X{:}), t=t, N=[], D=[], stats=[]);
end

for key = XRed.keys
    X = XRed(key);
    X = X{1};
    eta.stimRed(key) = struct(X=cat(1, X{:}), t=t, N=[], D=[], stats=[]);
end
%%
clear XBlue XRed iEu groupsBlue groupsRed isi t selBaseline normSR

% Make Raster
rd.stim = eu.getRasterData('stimtwocolor', window=[-0.1, 0.4], durErr=1e-3, shutterDelay=0, photoelectricBlankDuration=0.5e-3);
rd.press = eu.getRasterData('press', window=[-4, 0], alignTo='stop');
rd.pressMed = eu.getRasterData('press_spontaneous_medial', window=[-4, 0], alignTo='stop');
rd.pressLat = eu.getRasterData('press_spontaneous_lateral', window=[-4, 0], alignTo='stop');
rd.lick = eu.getRasterData('lick', window=[-4, 0], alignTo='stop');

% Make ETA
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

meta.stimBlue = configureDictionary("string", "cell");
meta.stimRed = configureDictionary("string", "cell");
for key = eta.stimBlue.keys'
    t = eta.stimBlue(key).t;
    X = eta.stimBlue(key).X;
    meta.stimBlue(key) = {mean(X(:, t>=p.metaWindowStim(1) & t<=p.metaWindowStim(2)), 2, 'omitnan')};
end
for key = eta.stimRed.keys'
    t = eta.stimRed(key).t;
    X = eta.stimRed(key).X;
    meta.stimRed(key) = {mean(X(:, t>=p.metaWindowStim(1) & t<=p.metaWindowStim(2)), 2, 'omitnan')};
end
clear t X
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

% c.isStimBlueUp = meta.stimBlue >= p.posRespThresholdStim;
% c.isStimBlueDown = meta.stimBlue <= p.negRespThresholdStim;
% c.isStimRedUp = meta.stimRed >= p.posRespThresholdStim;
% c.isStimRedDown = meta.stimRed <= p.negRespThresholdStim;
% 
% c.isStimBlueUpRedUpThereforeChrimsonMaybe = c.isStimBlueUp & c.isStimRedUp;
% c.isStimBlueUpRedNotUpThereforeCoChrMaybe = c.isStimBlueUp & ~c.isStimRedUp;
% c.isStimBlueNotUpRedUpThereforeChrimsonMaybe = ~c.isStimBlueUp & c.isStimRedUp;
% c.isStimBlueNotUpRedNotUp = ~c.isStimBlueUp & ~c.isStimRedUp;

% Boot response dir
if isempty(metaPath) && ~readBootButRecalculateMeta
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
end

%% Report bootstraped movement response direction
assert(nnz(isnan(boot.lick.h)) == 0)
assert(nnz(isnan(boot.press.h)) == 0)

% Sometimes we booted all the EphysUnits and stored it. Now we slice it
% with selEu
if ~isempty(metaPath) && readBootButRecalculateMeta
    for field = ["press", "lick", "pressMed", "pressLat"]
        if isfield(boot, field) && size(boot.(field).h, 1) == length(selEu)
            boot.(field).h = boot.(field).h(selEu, :);
            boot.(field).muDiffCI = boot.(field).muDiffCI(selEu, :);
            boot.(field).muDiffObs = boot.(field).muDiffObs(selEu, :);
        end
    end
end

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

save(metaSavePath, 'c', 'p', 'boot', 'rd', 'eta', 'meta')
