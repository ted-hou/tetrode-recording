
%5d. Peri-lick lick prob histogram
ax = nexttile(layout.bottom.right.tl, [layout.bottom.right.middle.h, 1]);
histogram(ax, BinEdges=lickHist.firstLick.edges, BinCounts=lickHist.firstLick.count, Normalization='probability', EdgeColor='none', FaceColor='black', FaceAlpha=1)
xlim(ax, [0, 1/9*4]); % Assuming 9 Hz, 5 licks (4 cycles)
xticks(ax, 0:0.2:0.5)
yticks(ax, [])
xlabel(ax, 'Time from first lick (s)')
ylabel(ax, ' lick\newlineprob', HorizontalAlignment='center')
fontsize(ax, p.fontSize, 'points')
hold(ax, 'off')
hLetter = text(ax, 0, 0, 'e', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.25, ax.Position(4) + 0.25, 0];


% 5e. Reach vs. Lick (scatter), Reach vs. Lick (osci subset, scatter), Pie-chart
sz = 7;
AX = gobjects(1, 2);
% 1. Reach vs. Lick scatter
ax = nexttile(layout.bottom.right.tl, [layout.bottom.right.bottom.h, 1]);
hold(ax, 'on')
sel = c.hasLick & c.hasPress;
x = meta.lick;
y = meta.press;
subselResp = sel & c.isLick;
subselNone = sel & ~c.isLick;
h = gobjects(2, 1);
h(2) = scatter(ax, x(subselNone), y(subselNone), sz, 'black', 'filled', Marker='o', MarkerFaceAlpha=0.5, MarkerEdgeAlpha=0.5, DisplayName=sprintf('others (%i)', nnz(subselNone)));   
h(1) = scatter(ax, x(subselResp), y(subselResp), sz, [0, 0.5, 0.5], 'filled', Marker='o', MarkerFaceAlpha=0.5, MarkerEdgeAlpha=0.5, DisplayName=sprintf('lick-entrained (%i)', nnz(subselResp)));

plot(ax, [-10, 10], [0, 0], 'k:');
plot(ax, [0, 0], [-10, 10], 'k:');
plot(ax, [-10, 10], [-10, 10], 'k:')

xlabel(ax, 'Peri-lick activity (a.u.)')
ylabel(ax, 'Peri-reach activity (a.u.)')
% axis(ax, 'equal')
xlim(ax, [-2, 5])
ylim(ax, [-2, 5])
% legend(h, Location='north')
fontsize(ax, p.fontSize, 'points')

hLetter = text(ax, 0, 0, 'f', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
ax.Units = 'inches';
hLetter.HorizontalAlignment = 'right';
hLetter.VerticalAlignment = 'top';
hLetter.Position = [-0.25, ax.Position(4) + 0.25, 0];


% 5f

eta.circLick.Z = eta.circLick.X.*exp(eta.circLick.t*1i);
eta.circLick.Z(:, 1) = mean(eta.circLick.X(:, [2, 30]), 2).*exp(eta.circLick.t(1)*1i);
meanZ = mean(eta.circLick.Z, 2);
eta.lickBoutNorm = eta.lickBout;
eta.lickBoutNorm.X = normalize(eta.lickBout.X, 2, 'zscore', 'robust');

maxBoutCycles = 4;
sel = c.hasPress & c.hasLick & c.isLick; 
phase = angle(meanZ(sel));
phase(phase < 0) = phase(phase < 0) + 2*pi;
[sortedPhase, I] = sort(phase);
[~, ~] = EphysUnit.plotETA(ax(4), eta.lickBoutNorm, sel, order=I, ...
    clim=[-2, 2], xlim=[0, 2*pi*maxBoutCycles], hidecolorbar=true);
xticks (ax(4), (0:2:8).*pi);
xticklabels(ax(4), [{'0'}, arrayfun(@(x) sprintf('%i\\pi', x), 2:2:8, UniformOutput=false)]);
title(ax(4), 'Lick-entrained')
ylabel(ax(4), 'Unit')
xlabel(ax(4), 'Lick phase')
xlim(ax(4), [0, 2*pi*maxBoutCycles])
ylim(ax(4), [0, nnz(sel)+1])
yt = 0:100:nnz(sel);
yt(1) = 1;
if round(yt(end)./100) == round(nnz(sel)./100)
    yt(end) = nnz(sel);
else
    yt(end + 1) = nnz(sel);
end
yt = unique(yt);
yticks(ax(4), yt)
fontsize(ax, p.fontSize, 'points')
fontname(ax, 'Arial')
axc = ax(4);