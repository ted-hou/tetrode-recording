%% Find osci lick units through circular statistics

%%
% eta.anyLickLeftRaw = eu.getETA('count', 'anylick', [-2, 2], resolution=1/Fs, normalize='none', alignTo='stop', includeInvalid=true);
% t = eta.anyLickLeftRaw.t; 
% x = eta.anyLickLeftRaw.X'*Fs;
% x = normalize(x, 1, 'zscore', 'robust');
% eta.anyLickLeftNorm = eta.anyLickLeftRaw;
% eta.anyLickLeftNorm.X = x';

%% Find lick pairs defining the start and stop of trials
maxLickInterval = 0.2;
minLickInterval = 0.075;

[~, expEuIndices, euExpIndices] = unique({eu.ExpName});

trials = cell(length(expEuIndices), 1);
for iExp = 1:length(expEuIndices)
    trials{iExp} = eu(expEuIndices(iExp)).getTrials('circlick');
    trials{iExp} = trials{iExp}(trials{iExp}.duration <= maxLickInterval & trials{iExp}.duration >= minLickInterval);
end
trials = trials(euExpIndices);

durations = cellfun(@(x) x.duration, trials(expEuIndices), UniformOutput=false);

ax = axes(figure());
histogram(ax, cat(2, durations{:}))
xlabel(ax, 'inter-lick interval (s)')
ylabel(ax, 'count')

clear ax iExp

%% Calculate circular PETH, where each lick occurs at 0 (or 2pi)
eta.circLick = eu.getETA('count', 'circlick', window=[0, 2*pi], resolution=2*pi/30, normalize='none', trials=trials);
% Units are spike counts per bin averaged across trials. Hard to convert to sp/s so we won't do
% it. Values are non-negative.

% Bootstrap to find the significance of average Z vector magnitudes (shuffle bins, not trials)

[magH, magCI] = bootCircLick(eta.circLick, alpha=0.01, nBoot=100000, replace=true);

%% Plot everything
close all

meanLickInterval = mean(cat(2, durations{:}));
meanBinWidth = meanLickInterval / length(eta.circLick.t);

eta.circLick.Z = eta.circLick.X.*exp(eta.circLick.t*1i)./meanBinWidth;
eta.circLick.Z(:, 1) = mean(eta.circLick.X(:, [2, 30]), 2).*exp(eta.circLick.t(1)*1i)./meanBinWidth;
meanZ = mean(eta.circLick.Z, 2);
% sel = 1:length(eu);
% sel = abs(meanZ) > 2;
% sel = c.isLick;
sel = magH;

tl = tiledlayout(figure(), 1, 2);
ax = nexttile(tl);
histogram(ax, angle(meanZ(sel)) ./ pi, -1:1/30:1)
xlabel(ax, 'phase (\pi)')
ylabel(ax, 'count')
ax = nexttile(tl);
histogram(ax, abs(meanZ(sel)), 50)
xlabel(ax, 'magnitude (a.u.)')
ylabel(ax, 'count')
title(tl, sprintf('Distribution of phase and magnitude for %i units', nnz(sel)))

phaseCirc = angle(meanZ(sel));
[~, I] = sort(phaseCirc);


eta.circLickNorm = eta.circLick;
eta.circLickNorm.X = normalize(eta.circLick.X, 2, 'zscore', 'robust');

tl = tiledlayout(figure(Units='inches', Position=[0, 0, 8, 5]), 1, 4);

ax = nexttile(tl);
EphysUnit.plotETA(ax, eta.circLickNorm, sel, order=I, clim=[-2, 2], xlim=[0, 2]*pi, hidecolorbar=true);
xticks (ax, [0, 1, 2] .* pi);
xticklabels(ax, ["0", "\pi", "2\pi"]);
xlabel(ax, 'Lick phase')
title(ax, '')
ylabel(ax, '')


ax = nexttile(tl);
EphysUnit.plotETA(ax, eta.anyLickNorm, sel, order=I, clim=[-2, 2], xlim=[-0.25, 0.25], hidecolorbar=true);
xlabel(ax, 'Time to lick')
title(ax, '')
ylabel(ax, '')

ax = nexttile(tl);
EphysUnit.plotETA(ax, eta.anyLickRightNorm, sel, order=I, clim=[-2, 2], xlim=[0, 0.5], hidecolorbar=true);
xlabel(ax, 'Time to lick')
title(ax, '')
ylabel(ax, '')

ax = nexttile(tl);
EphysUnit.plotETA(ax, eta.firstLickNorm, sel, order=I, clim=[-2, 2], xlim=[0.01, 0.5]);
xlabel(ax, 'Time to first lick')
ax.Colorbar.Layout.Tile = 'east';
title(ax, '')
ylabel(ax, '')

ylabel(tl, 'Unit', FontSize=p.fontSize, FontName='Arial')
title(tl, sprintf('%i lick entrained units (p<0.01)', nnz(sel)), FontSize=p.fontSize, FontName='Arial')

% Plot eta
% The grad student uses Euler's formular to convert angle (t, between 0 and 2pi) and magnitude (spike counts) to cartesian coords:
% (mag*exp(theta*i) = mag*cos(theta) + mag*sin(theta)*i;
% Then real and imaginary parts are x and y coords, respectively.
% close all
clear i
edges = 0:2*pi/30:2*pi;
fig = figure(Units='inches', Position=[0, 0, 8, 8]);
tl = tiledlayout(fig, 8, 8);
indices = find(magH(:)');
indices = indices(I);

for iEu = indices(1:64)
    ax = nexttile(tl);
    z = eta.circLick.Z(iEu, :);
    hold(ax, 'on')
    patch(ax, real(z), imag(z), 'black', FaceAlpha=0.25, EdgeColor='black')
    plot(ax, 10*[0, real(mean(z))], 10*[0, imag(mean(z))], 'k', LineWidth=2)
    lims = [-1, 1] * max(abs(z));
    plot(ax, [0, 0], lims, 'k:')
    plot(ax, lims, [0, 0], 'k:')
    xlim(ax, lims)
    ylim(ax, lims)
    axis(ax, 'equal')
    % xlim(ax, [-0.5, 0.5])
    % ylim(ax, [-0.5, 0.5])
end


function [magH, magCI] = bootCircLick(eta, varargin)
    parser = inputParser();
    parser.addRequired('eta', @isstruct);
    parser.addParameter('alpha', 0.01, @isnumeric);
    parser.addParameter('nBoot', 10000, @isnumeric);
    parser.addParameter('replace', true, @islogical) % true for bootstrap, false for wda
    parser.parse(eta, varargin{:});
    r = parser.Results;
    X = r.eta.X;
    t = r.eta.t;
    alpha = r.alpha;
    nBoot = r.nBoot;
    replace = r.replace;

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
