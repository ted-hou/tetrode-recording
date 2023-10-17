function c = getColor(i, n, maxHue)
    if nargin < 3
        if n <= 4
            c = 'rgbm';
            c = c(i);
            return
        end
        maxHue = 0.8;
    end
    c = hsl2rgb([maxHue*(i-1)./(n-1), 1, 0.4]);
end