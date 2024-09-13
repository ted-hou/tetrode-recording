function [x, y] = ax2figPosition(ax, fig, x, y)
assert(isgraphics(ax, 'axes'))
assert(isgraphics(fig, 'figure'))

if ax.Parent ~= fig
    warning('Fig is not parent of Ax. Unexpected behavior is afoot.')
end

figUnits = fig.Units;
axUnits = ax.Units;

ax.Units = 'normalized';
fig.Units = 'normalized';

xl = ax.XLim;
yl = ax.YLim;

normX = (x - xl(1))./diff(xl);
normY = (y - yl(1))./diff(yl);

x0 = ax.Position(1);
y0 = ax.Position(2);
w = ax.Position(3);
h = ax.Position(4);

x = x0 + w.*normX;
y = y0 + h.*normY;

ax.Units = axUnits;
fig.Units = figUnits;