close all
Fs = 100;            % Sampling frequency    


eta.anyLickRightRaw = eu.getETA('count', 'anylick', [-0, 0.5], resolution=1/Fs, normalize='none');
t = eta.anyLickRightRaw.t; 
x = eta.anyLickRightRaw.X'*100;
x = normalize(x, 1, 'zscore', 'robust');
eta.anyLickRightNorm = eta.anyLickRightRaw;
eta.anyLickRightNorm.X = x';

%%
ax = axes();
hold(ax, 'on')
clear P8 P6 P10 P16 P14 P18
for i = 1:size(x, 2)
    Y = fft(x(:, i));                
    T = 1/Fs;             % Sampling period       
    L = length(t);             % Length of signal
    t = (0:L-1)*T;        % Time vector
    
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    plot(ax, f,P1) 
    title('Single-Sided Amplitude Spectrum of X(t)')
    xlabel('f (Hz)')
    ylabel('|P1(f)|')
    P8(i) = P1(f==8);
    P6(i) = P1(f==6);
    P10(i) = P1(f==10);
    P16(i) = P1(f==16);
    P14(i) = P1(f==14);
    P18(i) = P1(f==18);
end

theta = 0.5;
relTheta = 0;
isLick = P8 > P6 + relTheta & P8 > P10 + relTheta & P16 > P14 + relTheta & P16 > P18 + relTheta & P8 > theta;
c.isLick = isLick;
nnz(isLick)
figure()
ax = subplot(1, 2, 1);
plot(t, x(:, isLick))
title(sprintf('N = %g (%.1f%%)', nnz(isLick), 100*nnz(isLick)/length(eu)))

ax = subplot(1, 2, 2); hold(ax, 'on')               
T = 1/Fs;             % Sampling period       
L = length(t);             % Length of signal
t = (0:L-1)*T;        % Time vector
f = Fs*(0:(L/2))/L;

P = zeros(26, nnz(isLick));
for i = find(isLick)
    Y = fft(x(:, i));
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    P(:, i) = P1;
end
plot(f, P1)
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

t = eta.anyLickNorm.t;
meta.anyLickNorm = transpose(mean(eta.anyLickNorm.X(:, t >= -0.05 & t < 0), 2, 'omitnan'));

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
clear t x i Y Fs T L P2 P1 f P8 P6 P10 P16 P14 P18 theta relTheta isLick ax P P2 P1 signWindow sortWindow plotWindow latency
