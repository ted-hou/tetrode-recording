iTest = 0;
iTest = iTest + 1;
[iStart, iStop] = isin(1:10, [2, 4], false, true);
assert(iStart == 3 && iStop == 3, "%i iStart=%i iStop=%i", iTest, iStart, iStop)

iTest = iTest + 1;
[iStart, iStop] = isin(1:10, [2.1, 3.9], false, true);
assert(iStart == 3 && iStop == 3, "%i iStart=%i iStop=%i", iTest, iStart, iStop)

iTest = iTest + 1;
[iStart, iStop] = isin(1:10, [10, 21], false, true);
assert(isempty(iStart) && isempty(iStop), "%i iStart=%i iStop=%i", iTest, iStart, iStop)

iTest = iTest + 1;
[iStart, iStop] = isin(1:10, [0, 1], false, true);
assert(isempty(iStart) && isempty(iStop), "%i iStart=%i iStop=%i", iTest, iStart, iStop)

iTest = iTest + 1;
[iStart, iStop] = isin(1:10, [20, 1], false, true);
assert(isempty(iStart) && isempty(iStop), "%i iStart=%i iStop=%i", iTest, iStart, iStop)


[iStart, iStop] = isin(1:10, [0.5, 1.5], false, true);
assert(iStart==1 && iStop==1, "%i iStart=%i iStop=%i", iTest, iStart, iStop)


[iStart, iStop] = isin(1:10, [0.5, 1], false, true);
assert(isempty(iStart) && isempty(iStop), "%i iStart=%i iStop=%i", iTest, iStart, iStop)








iTest = iTest + 1;
[iStart, iStop] = isin(1:10, [2, 4], true, true);
assert(iStart == 2 && iStop == 4, "%i iStart=%i iStop=%i", iTest, iStart, iStop)

iTest = iTest + 1;
[iStart, iStop] = isin(1:10, [2.1, 3.9], true, true);
assert(iStart == 3 && iStop == 3, "%i iStart=%i iStop=%i", iTest, iStart, iStop)

iTest = iTest + 1;
[iStart, iStop] = isin(1:10, [10, 21], true, true);
assert(iStart == 10 && iStop == 10, "%i iStart=%i iStop=%i", iTest, iStart, iStop)

iTest = iTest + 1;
[iStart, iStop] = isin(1:10, [9.8, 10.1], true, true);
assert(iStart == 10 && iStop == 10, "%i iStart=%i iStop=%i", iTest, iStart, iStop)

iTest = iTest + 1;
[iStart, iStop] = isin(1:10, [9.8, 25], true, true);
assert(iStart == 10 && iStop == 10, "%i iStart=%i iStop=%i", iTest, iStart, iStop)

iTest = iTest + 1;
[iStart, iStop] = isin(1:10, [9.8, -9], true, true);
assert(isempty(iStart) && isempty(iStop), "%i iStart=%i iStop=%i", iTest, iStart, iStop)

[iStart, iStop] = isin(1:10, [0, 1], true, true);
assert(iStart == 1 && iStop == 1, "%i iStart=%i iStop=%i", iTest, iStart, iStop)

[iStart, iStop] = isin(1:10, [20, 1], true, true);
assert(isempty(iStart) && isempty(iStop), "%i iStart=%i iStop=%i", iTest, iStart, iStop)

[iStart, iStop] = isin(1:10, [0.5, 1.5], true, true);
assert(iStart==1 && iStop==1, "%i iStart=%i iStop=%i", iTest, iStart, iStop)

[iStart, iStop] = isin(1:10, [0.5, 1], true, true);
assert(iStart==1 && iStop==1, "%i iStart=%i iStop=%i", iTest, iStart, iStop)






