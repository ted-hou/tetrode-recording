%% Find osci lick units through circular statistics

%%
Fs = 100;            % Sampling frequency
% eta.anyLickRaw = eu.getETA('count', 'anylick', [-0.25, 0.25], resolution=1/Fs, normalize='none', alignTo='stop', includeInvalid=true);
% t = eta.anyLickRaw.t; 
% x = eta.anyLickRaw.X'*Fs;
% x = normalize(x, 1, 'zscore', 'robust');
% eta.anyLickNorm = eta.anyLickRaw;
% eta.anyLickNorm.X = x';
% 
% eta.firstLickRaw = eu.getETA('count', 'firstlick', [0.01, 0.51], resolution=1/Fs, normalize='none', alignTo='stop', includeInvalid=true);
% t = eta.firstLickRaw.t; 
% x = eta.firstLickRaw.X'*Fs;
% x = normalize(x, 1, 'zscore', 'robust');
% eta.firstLickNorm = eta.firstLickRaw;
% eta.firstLickNorm.X = x';
% 
% meta.firstLick = mean(eta.firstLickRaw.X(:, t < 0.02), 2, 'omitnan');

%% Find lick pairs defining the start and stop of trials
maxLickInterval = 0.2;
minLickInterval = 0.05;
minBoutCycles = 4;
maxBoutCycles = 4;

[~, expEuIndices, euExpIndices] = unique({eu.ExpName});

trialsCircLick = cell(length(expEuIndices), 1);
trialsCircLickBaseline = cell(length(expEuIndices), 1);
bouts = trialsCircLick;
for iExp = 1:length(expEuIndices)
    trialsCircLick{iExp} = eu(expEuIndices(iExp)).getTrials('circlick');
    trialsCircLick{iExp} = trialsCircLick{iExp}(trialsCircLick{iExp}.duration <= maxLickInterval & trialsCircLick{iExp}.duration >= minLickInterval);
    trialsCircLickBaseline{iExp} = eu(expEuIndices(iExp)).makeTrials('circlickbaseline');
    bouts{iExp} = eu(expEuIndices(iExp)).getTrials('lickbout', sorted=false, minInterval=minLickInterval, maxInterval=maxLickInterval, ...
        minBoutCycles=minBoutCycles, maxBoutCycles=maxBoutCycles);
end
trialsCircLick = trialsCircLick(euExpIndices);
trialsCircLickBaseline = trialsCircLickBaseline(euExpIndices);
bouts = bouts(euExpIndices);

durations = cellfun(@(x) x.duration, trialsCircLick(expEuIndices), UniformOutput=false);
boutCycles = cellfun(@(t) sum(~t.isEmpty, 2), bouts, UniformOutput=false);

tl = tiledlayout(figure(), 1, 2);
ax = nexttile(tl);
interLickInterval = cat(2, durations{:});
histogram(ax, cat(2, interLickInterval))
xlabel(ax, 'inter-lick interval (s)')
ylabel(ax, 'count')

ax = nexttile(tl);
histogram(ax, cat(1, boutCycles{:}))

clear ax iExp tl

%% Calculate circular PETH, where each lick occurs at 0 (or 2pi)
eta.circLick = eu.getETA('count', 'circlick', window=[0, 2*pi], resolution=2*pi/30, normalize='none', trials=trialsCircLick);
% eta.circLick = eu.getETA('count', 'circlick', window=[0, 2*pi], resolution=2*pi/30, normalize='none', trials=cellfun(@(trials) trials(:), bouts, UniformOutput=false));

% Control circlick using randomly selected
eta.circLickBaseline = eu.getETA('count', 'circlick', window=[0, 2*pi], resolution=2*pi/30, normalize='none', trials=trialsCircLickBaseline);

eta.lickBout = eu.getETA('count', 'lickbout', window=[0, 2*pi*maxBoutCycles], resolution=2*pi/30, normalize='none', trials=bouts);

% Units are spike counts per bin averaged across trials. Hard to convert to sp/s so we won't do
% it. Values are non-negative.

%% Bootstrap to find the significance of average Z vector magnitudes (shuffle bins, not trials)
clear bootCLick
[bootCLick.magH, bootCLick.magCI] = bootCircLick(eta.circLick, alpha=0.01, nBoot=100000, replace=false, seed=42);
[bootCLick.magHBaseline, bootCLick.magCIBaseline] = bootCircLick(eta.circLickBaseline, alpha=0.01, nBoot=100000, replace=false, seed=42);

c.isLickOsci = bootCLick.magH(:)';
c.isBaselineOsci = bootCLick.magHBaseline(:)';
c.isLick = c.isLickOsci & ~c.isBaselineOsci;

%% Calculate peri-lick lick frequency histograms for any first lick in trial
[~, expEuIndices] = unique({eu.ExpName});
FsLick = 500;
lickHistEdges = 0.01:1/FsLick:2;
lickHistCenters = 0.5*(lickHistEdges(2:end) + lickHistEdges(1:end-1));

lickHistCounts = zeros(size(lickHistCenters));
lickHistNLicks = 0;
for iExp = 1:length(expEuIndices)
    iEu = expEuIndices(iExp);
    firstLickTimes = [eu(iEu).makeTrials('firstlick').Stop];
    allLickTimes = eu(iEu).EventTimes.Lick;
    for iLick = 1:length(firstLickTimes)
        edgesGlobal = firstLickTimes(iLick) + lickHistEdges;
        n = histcounts(allLickTimes, edgesGlobal);
        lickHistCounts = lickHistCounts + n;
        lickHistNLicks = lickHistNLicks + 1;
    end
end
lickHist.firstLick.count = lickHistCounts;
lickHist.firstLick.pdf = lickHistCounts ./ lickHistNLicks;
lickHist.firstLick.t = lickHistCenters;
lickHist.firstLick.edges = lickHistEdges;

clear expEuIndices FsLick lickHistEdges lickHistCenters lickHistNLicks lickHistCounts iExp iEu firstLickTimes allLickTimes iLick edgesGlobal n

% Calculate peri-correct-lick lick frequency histograms for osci cells
sel = c.hasLick & c.hasPress & c.isLick;
[~, expEuIndices] = unique({eu(sel).ExpName});
FsLick = 500;
lickHistEdges = 0.01:1/FsLick:2;
lickHistCenters = 0.5*(lickHistEdges(2:end) + lickHistEdges(1:end-1));

lickHistCounts = zeros(size(lickHistCenters));
lickHistNLicks = 0;
for iExp = 1:length(expEuIndices)
    iEu = expEuIndices(iExp);
    trials = eu(iEu).getTrials('lick');
    trials = trials(trials.duration() >= 4);
    firstLickTimes = [trials.Stop];
    allLickTimes = eu(iEu).EventTimes.Lick;
    for iLick = 1:length(firstLickTimes)
        edgesGlobal = firstLickTimes(iLick) + lickHistEdges;
        n = histcounts(allLickTimes, edgesGlobal);
        lickHistCounts = lickHistCounts + n;
        lickHistNLicks = lickHistNLicks + 1;
    end
end

lickHist.correctLickOsci.count = lickHistCounts;
lickHist.correctLickOsci.pdf = lickHistCounts ./ lickHistNLicks;
lickHist.correctLickOsci.t = lickHistCenters;
lickHist.correctLickOsci.edges = lickHistEdges;

clear sel trials expEuIndices FsLick lickHistEdges lickHistCenters lickHistNLicks lickHistCounts iExp iEu firstLickTimes allLickTimes iLick edgesGlobal n

function [magH, magCI] = bootCircLick(eta, varargin)
    parser = inputParser();
    parser.addRequired('eta', @isstruct);
    parser.addParameter('alpha', 0.01, @isnumeric);
    parser.addParameter('nBoot', 10000, @isnumeric);
    parser.addParameter('replace', true, @islogical) % true for bootstrap, false for wda
    parser.addParameter('seed', 42, @isnumeric)
    parser.parse(eta, varargin{:});
    r = parser.Results;
    X = r.eta.X;
    t = r.eta.t;
    alpha = r.alpha;
    nBoot = r.nBoot;
    replace = r.replace;
    seed = r.seed;
    rng(seed);

    magH = false(size(X, 1), 1);
    magCI = zeros(size(X, 1), 2);

    nUnits = size(X, 1);
    nBins = length(t);
    for iUnit = 1:nUnits
        if replace
            I = randi(nBins, [nBoot, nBins]);
        else
            I = zeros(nBoot, nBins);
            for iBoot = 1:nBoot
                I(iBoot, :) = randperm(nBins);
            end
        end
        x = X(iUnit, :);
        x(1) = mean(x([2, 30])); % First bin has lick artifact usually
        zObs = mean(x.*exp(t*1i));
        zRand = mean(x(I).*exp(t*1i), 2);
        
        magObs = abs(zObs);
        magRand = abs(zRand);
        magCI(iUnit, :) = quantile(magRand, [0, 1 - alpha]);
        magH(iUnit) = magObs > magCI(iUnit, 2);
    end
end