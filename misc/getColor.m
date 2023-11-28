function c = getColor(i, n, maxHue)
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
        c = hsl2rgb([maxHue*(i-1)./(n-1), 1, 0.4]);
    end
end