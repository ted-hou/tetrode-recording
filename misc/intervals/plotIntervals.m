function plotIntervals(ax, intervals, color, h)
if nargin < 3
    color = 'red';
end
if nargin < 4
    h = [-1, 1];
end
hold(ax, 'on')
for i = 1:size(intervals, 2)
    a1 = intervals(1, i);
    a2 = intervals(2, i);
    patch(ax, [a1, a2, a2, a1], h([1, 1, 2, 2]), color, EdgeColor=color, FaceColor=color, FaceAlpha=0.25, EdgeAlpha=1);
end
hold(ax, 'off')