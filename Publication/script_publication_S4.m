
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
sel = magH(:)' & c.hasLick(:)' & c.hasPress(:)';

tl = tiledlayout(figure(Units='inches', Position=[0, 0, 2, 3]), 2, 1);
ax = nexttile(tl);
histogram(ax, angle(meanZ(sel)), (-1:1/10:1).*pi)
xlabel(ax, 'Phase')
% ylabel(ax, 'count')
xticks (ax, (-1:1:1) .* pi);
xticklabels(ax, arrayfun(@(x) sprintf('%g\\pi', x), -1:1:1, UniformOutput=false));
set(ax, FontSize=p.fontSize, fontName='Arial');


ax = nexttile(tl);
histogram(ax, abs(meanZ(sel)), 0:20)
xlabel(ax, 'Magnitude (sp/s)')
% ylabel(ax, 'count')
set(ax, FontSize=p.fontSize, fontName='Arial');

ylabel(tl, 'Count', FontSize=p.fontSize, fontName='Arial')

phaseCirc = angle(meanZ(sel));
[~, I] = sort(phaseCirc);

%%
eta.circLickNorm = eta.circLick;
eta.circLickNorm.X = normalize(eta.circLick.X, 2, 'zscore', 'robust');

eta.lickBoutNorm = eta.lickBout;
eta.lickBoutNorm.X = normalize(eta.lickBout.X, 2, 'zscore', 'robust');

tl = tiledlayout(figure(Units='inches', Position=[0, 0, 6, 3]), 1, 3);

ax = nexttile(tl);
EphysUnit.plotETA(ax, eta.circLickNorm, sel, order=I, clim=[-2, 2], xlim=[0, 2]*pi, hidecolorbar=true);
xticks (ax, [0, 1, 2] .* pi);
xticklabels(ax, arrayfun(@(x) sprintf('%i\\pi', x), 0:2, UniformOutput=false));
xlabel(ax, 'Lick phase')
title(ax, '')
ylabel(ax, '')


ax = nexttile(tl);
EphysUnit.plotETA(ax, eta.lickBoutNorm, sel, order=I, clim=[-2, 2], xlim=[0, 2*pi*maxBoutCycles], hidecolorbar=true);
xlabel(ax, 'Lick phase')
xticks (ax, (0:2:8).*pi);
xticklabels(ax, arrayfun(@(x) sprintf('%i\\pi', x), 0:2:8, UniformOutput=false));
xlabel(ax, 'Lick phase')
title(ax, '')
ylabel(ax, '')

ax = nexttile(tl);
EphysUnit.plotETA(ax, eta.anyLickNorm, sel, order=I, clim=[-2, 2], xlim=[-0.25, 0.25], hidecolorbar=false);
xlabel(ax, 'Time to any lick (s)')
ax.Colorbar.Layout.Tile = 'east';
title(ax, '')
ylabel(ax, '')

% ax = nexttile(tl);
% EphysUnit.plotETA(ax, eta.firstLickNorm, sel, order=I, clim=[-2, 2], xlim=[0.01, 0.5]);
% xlabel(ax, 'Time to first lick (s)')
% ax.Colorbar.Layout.Tile = 'east';
% title(ax, '')
% ylabel(ax, '')



ylabel(tl, 'Unit', FontSize=p.fontSize, FontName='Arial')
title(tl, sprintf('%i lick entrained units (p<0.01)', nnz(sel)), FontSize=p.fontSize, FontName='Arial')

%% Plot eta
% The grad student uses Euler's formular to convert angle (t, between 0 and 2pi) and magnitude (spike counts) to cartesian coords:
% (mag*exp(theta*i) = mag*cos(theta) + mag*sin(theta)*i;
% Then real and imaginary parts are x and y coords, respectively.
% close all
clear i
edges = 0:2*pi/30:2*pi;
fig = figure(Units='inches', Position=[0, 0, 8, 8]);
w = ceil(sqrt(nnz(sel))); h = w;
tl = tiledlayout(fig, h, w, TileSpacing='tight', Padding='tight');
indices = find(sel);
indices = indices(I);

for iEu = indices
    ax = nexttile(tl);
    z = eta.circLick.Z(iEu, :);
    hold(ax, 'on')
    patch(ax, real(z), imag(z), 'black', FaceAlpha=0, EdgeColor='black')
    plot(ax, 10*[0, real(mean(z))], 10*[0, imag(mean(z))], 'k', LineWidth=2)
    lims = [-1, 1] * max(abs(z));
    plot(ax, [0, 0], lims, 'k:')
    plot(ax, lims, [0, 0], 'k:')
    xlim(ax, lims)
    ylim(ax, lims)
    axis(ax, 'equal')
    axis(ax, 'off')
    % xlim(ax, [-0.5, 0.5])
    % ylim(ax, [-0.5, 0.5])
end

