Fs = 100;            % Sampling frequency    
eta.firstLickRaw = eu.getETA('count', 'firstlick', [0.01, 0.51], resolution=1/Fs, normalize='none', alignTo='stop', includeInvalid=true);
t = eta.firstLickRaw.t; 
x = eta.firstLickRaw.X'*Fs;
x = normalize(x, 1, 'zscore', 'robust');
eta.firstLickNorm = eta.firstLickRaw;
eta.firstLickNorm.X = x';

meta.firstLick = mean(eta.firstLickRaw.X(:, t < 0.02), 2, 'omitnan');

%% FFT to find cells with 8Hz modulation
Fs = 100;
t = eta.firstLickRaw.t;  
x = eta.firstLickRaw.X'*Fs;
x = normalize(x, 1, 'zscore', 'robust');
nUnits = size(x, 2);

Y = fft(x);
L = length(t); % length of signal;
freq = Fs/L*(0:L-1);
magnitude = abs(Y/L);
Y(magnitude<1e-6) = 0;
phase = angle(Y);

theta = 0.37;%prctile(magnitude(freq==8, :), 80);
isLick = magnitude(freq==8, :) > theta;
c.isLick = isLick;
nnz(isLick)

fig = figure(Units='inches', Position=[0, 0, 8, 2]);
tl = tiledlayout(fig, 1, 4);
ax = nexttile(tl);
histogram(ax, magnitude(freq==8, :), 100)

ax = nexttile(tl);
histogram(ax, phase(freq==8, c.isLick), 50)
xticks((-1:1)*pi)
xticklabels(["-\pi", "0", "\pi"])


ax = nexttile(tl);
plot(ax, eta.firstLickRaw.t, mean(x(:, isLick), 2))
title(ax, sprintf('theta=%g, N = %g (%.1f%%)', theta, nnz(isLick), 100*nnz(isLick)/length(eu)))
xlabel(ax, 'Time to any lick (s)')
ylabel(ax, 'Normalized spike rate (pop avg, a.u.)')

ax = nexttile(tl);
scatter(ax, meta.lick(sel), meta.firstLick(sel), 5, 'black')
axis(ax, 'equal')
xlabel(ax, 'Pre-lick responses (a.u.)')
ylabel(ax, 'First-lick response(a.u.)')

%% Plot some ETA heatmaps for funsies
tl = tiledlayout(figure(Units='inches', Position=[0, 0, 6, 5]), 1, 2);

[~, I] = sort(phase(freq==8, c.isLick));

ax = nexttile(tl);
[~, ~] = EphysUnit.plotETA(ax, eta.anyLickNorm, c.isLick, order=I, ...
    clim=[-2, 2], xlim=[-0.25, 0.25], hidecolorbar=true);
xlabel(ax, 'Time to lick')
title(ax, '')
ylabel(ax, '')

ax = nexttile(tl);
[~, ~] = EphysUnit.plotETA(ax, eta.firstLickNorm, c.isLick, order=I, ...
    clim=[-2, 2], xlim=[0.01, 0.5]);
xlabel(ax, 'Time to first lick')
ax.Colorbar.Layout.Tile = 'east';
title(ax, '')
ylabel(ax, '')

ylabel(tl, 'Unit', FontSize=p.fontSize, FontName='Arial')
title(tl, sprintf('%i lick entrained units (FFT)', nnz(c.isLick)), FontSize=p.fontSize, FontName='Arial')

%% Calculate peri-lick lick frequency histograms
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
lickHistCountsNorm = lickHistCounts ./ lickHistNLicks;

Y = fft(lickHistCounts);
L = length(lickHistCounts);             % Length of signal

ax = axes(figure);
plot(ax, FsLick/L*(0:L-1), abs(Y/L),"LineWidth",3)
xlim(ax, [0, 32])
title(ax, "Complex Magnitude of fft Spectrum")
xlabel(ax, "f (Hz)")
ylabel(ax, "|fft(X)|")

% clear iEu iLick n expEuIndices edgesGlobal firstLickTimes allLickTimes
% 
% clear t x i Y Fs FsLick T L P2 P1 f P8 P6 P10 P16 P14 P18 theta relTheta isLick ax P P2 P1 signWindow sortWindow plotWindow latency
