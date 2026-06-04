function varargout = isin(x, window, inclusive, assumeAscending)
% b = isin(x, window, inclusive, false)
% faster and returns only start/stop indices
    if nargin < 3
        inclusive = true;
    end
    if nargin < 4
        assumeAscending = false;
    end
    if ~assumeAscending
        if inclusive
            b = x>=window(1) & x<=window(2);
        else
            b = x>window(1) & x<window(2);
        end
        varargout = {b};
    % Try to optimize this, assuming x is sorted (ascending), as this is
    % usually a list of timestamps.
    else
        iStart = binarySearchStartIndex(x, window(1), inclusive);
        iStop = binarySearchStopIndex(x, window(2), inclusive);
        if isempty(iStart) || isempty(iStop)
            varargout = {[], []};
        else
            varargout = {iStart, iStop};
        end
    end
end

% Binary search for first element greater than (or equal to) k
% assume array X is sorted ascending
function i = binarySearchStartIndex(X, k, inclusive)
    % skip search if last element <(<=) threshold
    if (inclusive && X(end) < k) || (~inclusive && X(end) <= k)
        i = [];
        return
    end

    a = 1; 
    b = length(X); 
    while a < b 
        mid = floor((a + b) / 2); 
        if X(mid) <= k
            a = mid + 1; 
        else 
            b = mid; 
        end 
    end 
    if (a < length(X) && X(a) <= k) 
          a = a + 1; 
    end
    if inclusive && a > 1 && X(a - 1) == k
        i = a - 1;
    else
        i = a;
    end
end 

% Binary search for last element smaller than (or equal to) k
% assume array X is sorted ascending
function i = binarySearchStopIndex(X, k, inclusive) 
    % skip search if last element >(>=) threshold
    if (inclusive && X(1) > k) || (~inclusive && X(1) >= k)
        i = [];
        return
    end

    a = 1; 
    b = length(X); 
    while a < b 
        mid = floor((a + b) / 2); 
        if X(mid) >= k
            b = mid; 
        else 
            a = mid + 1; 
        end 
    end 
    if (a < length(X) && X(a) < k) 
          a = a + 1; 
    end
    if inclusive && X(a) <= k
        i = a;
    else
        i = a - 1;
    end
end