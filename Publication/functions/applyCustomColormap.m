function cm=applyCustomColormap(ax, varargin)

p = inputParser();
p.addRequired('ax', @(x) isgraphics(x, 'Axes'))
p.addOptional('clim', [])
p.addParameter('hlim', [0.35, 0, 0, -0.35])
p.addParameter('h0', 0.35)
p.addParameter('llim', [0.125, 0.5, 0.5, 0.25])
p.addParameter('hpwr', 0.5)
p.addParameter('lpwr', 1)
p.parse(ax, varargin{:})
r = p.Results;
ax = r.ax;

if isempty(r.clim)
    r.clim = ax.CLim;
end

nPos = abs(r.clim(2))*1000;
nNeg = abs(r.clim(1))*1000;
nColors = nPos + nNeg;
assert(mod(nPos, 1) == 0)
assert(mod(nNeg, 1) == 0)
cm = ones(nColors, 3);
cm(:, 1) = [linspace(r.hlim(1).^(1/r.hpwr), r.hlim(2).^(1/r.hpwr), nNeg).^r.hpwr, -linspace(abs(r.hlim(3)).^(1/r.hpwr), abs(r.hlim(4)).^(1/r.hpwr), nPos).^r.hpwr] + r.h0;
cm(:, 3) = [linspace(r.llim(1).^(1/r.lpwr), r.llim(2).^(1/r.lpwr), nNeg), linspace(r.llim(3).^(1/r.lpwr), r.llim(4).^(1/r.lpwr), nPos)].^r.lpwr;

cm(cm>1)=1;
cm(cm<0)=0;
cm = hsl2rgb(cm);
colormap(ax, cm)
clim(ax, r.clim)