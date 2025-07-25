function c = getColor(i, n, maxHue, varargin)
    p = inputParser();
    p.addOptional('s', 1);
    p.addOptional('l', 0.4);
    p.parse(varargin{:});
    s = p.Results.s;
    l = p.Results.l;
    if nargin < 3
        if n <= 4
            c = 'rgbm';
            c = c(i);
            return
        end
        maxHue = 0.8;
    end
    if length(i) > 1
        c = zeros(length(i), 3);
        for j = 1:length(i)
            c(j, :) = getColor(i(j), n, maxHue);
        end
        return
    else
        c = hsl2rgb([maxHue*(i-1)./(n-1), s, l]);
    end
end