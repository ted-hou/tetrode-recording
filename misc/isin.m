function b = isin(x, window, inclusive)
    if nargin < 3
        inclusive = true;
    end
    if inclusive
        b = x>=window(1) & x<=window(2);
    else
        b = x>window(1) & x<window(2);
    end
end
