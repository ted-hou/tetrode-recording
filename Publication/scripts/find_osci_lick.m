Fs = 100;            % Sampling frequency    
eta.firstLickRaw = eu.getETA('count', 'firstlick', [0, 0.5], resolution=1/Fs, normalize='none', alignTo='stop', includeInvalid=true);
t = eta.firstLickRaw.t; 
x = eta.firstLickRaw.X'*Fs;
x = normalize(x, 1, 'zscore', 'robust');
eta.firstLickNorm = eta.firstLickRaw;
eta.firstLickNorm.X = x';

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

fig = figure(Units='inches', Position=[0, 0, 6, 2]);
tl = tiledlayout(fig, 1, 3);
ax = nexttile(tl);
histogram(ax, magnitude(freq==8, :), 100)

ax = nexttile(tl);
histogram(ax, phase(freq==8, :), 50)
xticks((-1:1)*pi)
xticklabels(["-\pi", "0", "\pi"])


theta = 0.37;%prctile(magnitude(freq==8, :), 80);
isLick = magnitude(freq==8, :) > theta;
c.isLick = isLick;
nnz(isLick)
ax = nexttile(tl);
plot(ax, eta.firstLickRaw.t, mean(x(:, isLick), 2))
title(ax, sprintf('theta=%g, N = %g (%.1f%%)', theta, nnz(isLick), 100*nnz(isLick)/length(eu)))
xlabel(ax, 'Time to any lick (s)')
ylabel(ax, 'Normalized spike rate (pop avg, a.u.)')
% 
% ax = subplot(1, 3, 2);
% hold(ax, 'on')
% histogram(ax, P8, 50, FaceColor='none');
% ylim(ax, 'manual')
% plot(ax, [theta, theta], ax.YLim, 'r--');
% xlabel(ax, 'Power at 8 Hz')
% ylabel(ax, '# units')
% 
% ax = subplot(1, 3, 3); hold(ax, 'on')               
% L = length(t);             % Length of signal
% f = Fs*(0:(L/2))/L;
% 
% P = zeros(26, nnz(isLick));
% for i = find(isLick)
%     Y = fft(x(:, i));
%     P2 = abs(Y/L);
%     P1 = P2(1:L/2+1);
%     P1(2:end-1) = 2*P1(2:end-1);
%     P(:, i) = P1;
% end
% plot(ax, f, P1)
% title(ax, 'Single-Sided Amplitude Spectrum of X(t)')
% xlabel(ax, 'f (Hz)')
% ylabel(ax, '|P1(f)|')
% xticks(ax, [0, 4, 8, 12])
% 
% t = eta.firstLickNorm.t;
% meta.firstLickNorm = transpose(mean(eta.firstLickNorm.X(:, t >= -0.05 & t < 0), 2, 'omitnan'));

%
% signWindow = [-0.13, -0.01];
% sortWindow = [-0.13, -0.01];
% plotWindow = [-0.25, 0.25];
% EphysUnit.plotETA(eta.anyLickNorm, [], xlim=plotWindow, clim=[-2 2], sortWindow=sortWindow, signWindow=signWindow, sortThreshold=0.5, negativeSortThreshold=Inf); title('Lick ETA')
% ax = EphysUnit.plotETA(eta.anyLickNorm, isLick, xlim=plotWindow, clim=[-2 2], sortWindow=sortWindow, signWindow=signWindow, sortThreshold=0.5, negativeSortThreshold=Inf); title('Lick ETA')
% ax.Parent.Position(3) = 0.25;
% ax.Parent.Position(4) = 0.4;%1*nnz(isLick)/nnz(c.hasLick);
% xlabel(ax, 'Time to lick (s)')
% title(ax, 'Lick ETA (ITI)')
% 
% 
% % ax = EphysUnit.plotETA(eta.anyLickNorm, c.isLick & c.isLickResponsive, xlim=plotWindow, clim=[-4 4], sortWindow=sortWindow, signWindow=signWindow, sortThreshold=0.5, negativeSortThreshold=Inf); title('Lick ETA')
% % ax.Parent.Position(3) = 0.25;
% % ax.Parent.Position(4) = 1*nnz(c.isLick & c.isLickResponsive)/nnz(c.hasLick);
% % xlabel(ax, 'Time relative to lick (s)')
% % title(ax, 'Lick ETA (ITI)')
% 
% [ax, ~, ~, latency] = EphysUnit.plotDoubleETA(eta.anyLickNorm, eta.lick, c.isLick & c.isLickResponsive, 'Lick (ITI)', 'First Lick', xlim=[-4,0], clim=[-2, 2], sortWindow=sortWindow, signWindow=signWindow, sortThreshold=0.5, negativeSortThreshold=Inf); 
% xlabel(ax(1), 'Time to lick (s)'); xlabel(ax(2), 'Time to first lick (s)')
% xlim(ax(1), [-0.25, 0.25])
% ax(1).FontSize=12;
% ax(2).FontSize=12;
% ax(1).Parent.Position(4) = 0.4;%1*nnz(c.isLick & c.isLickResponsive)/nnz(c.isLickResponsive & c.hasPress);
% 
% 
% %
% fprintf(1, 'Out of %g units, %g (%.1f%%) has oscilatory lick-related activity.\n', length(eu), nnz(c.isLick), nnz(c.isLick) / length(eu) * 100);
% fprintf(1, 'Out of %g lick-responsive units, %g (%.1f%%) has oscilatory lick-related activity.\n', nnz(c.isLickResponsive), nnz(c.isLick & c.isLickResponsive), nnz(c.isLick & c.isLickResponsive) / nnz(c.isLick) * 100)
% fprintf(1, '\t%g (%.1f%%) are lick activated, %g(%.1f%%) are suppressed.\n', nnz(c.isLick & c.isLickUp), nnz(c.isLick & c.isLickUp) / nnz(c.isLick & c.isLickResponsive) * 100, ...
%         nnz(c.isLick & c.isLickDown), nnz(c.isLick & c.isLickDown) / nnz(c.isLick & c.isLickResponsive) * 100)
% 
% 

%% Calculate peri-lick lick frequency histograms
[~, expEuIndices] = unique({eu.ExpName});
FsLick = 500;
lickHistEdges = 0:1/FsLick:0.5;
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
