function C = intervalUnion(A, B)
% C = intervalUnion(A): A are a 2xn (or 2xm) arrays where each 
% column contains an interval [a1, a2] where a1<=a2. a1s do not need to be
% sorted
% Returns C: a 2xn array of intervals that represent the union of all intervals in A.
if nargin < 2
    B = [];
else
    assert(isnumeric(B) && size(B, 1) == 2)
end

assert((isnumeric(A) && size(A, 1) == 2) || isempty(A))
A = [A, B];

if isempty(A) && isempty(B)
    C = [];
    return
end

assert(all(A(1, :) <= A(2, :)), 'Invalid intervals: input A does not satisfy a1 <= a2.');
A(:, A(1, :) == A(2, :)) = []; % Remove intervals with zero length;

if isempty(A) && isempty(B)
    C = [];
    return
end

tokens = [-ones(1, size(A, 2)); ones(1, size(A, 2))]; % negative one for opening, positive one for closing, an equal number of opening and closing brackets represent a disjoint interval
tokens = tokens(:)';
A = A(:)';

% Sort all the values (isOpening labels goes with the sorting)
[A, order] = sort(A, 'ascend');
tokens = tokens(order);

assert(tokens(1)==-1 && tokens(end)==1, 'Sorted edges should still start with an opening and with a closing edge.')

stack = cumsum(tokens);

assert(stack(1)==-1 && stack(end)==0)

isOpening = strfind([0, stack], [0, -1]);
isClosing = strfind(stack, [-1, 0]) + 1;

C = [A(isOpening); A(isClosing)];