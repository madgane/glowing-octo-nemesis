
function [entryLocs] = findMaxEntryLocs(gMatrix)

[iRow,~] = size(gMatrix);
lValue = max(max(gMatrix));
[lIndices] = find(lValue == gMatrix);

foundLength = length(lIndices);
entryLocs = zeros(foundLength,1);

for iIndex = 1:foundLength
    entryLocs(iIndex,1) = mod(lIndices(iIndex,1) - 1,iRow) + 1;
    entryLocs(iIndex,1) = entryLocs(iIndex,1) + sqrt(-1) * (floor((lIndices(iIndex,1) - 1) / iRow) + 1);
end


