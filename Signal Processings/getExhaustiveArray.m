
function [userArray] = getExhaustiveArray(nUsers,tUsers)

cIndex = 1;
combFact = factorial(nUsers) / (factorial(tUsers) * factorial(nUsers - tUsers));
userArray = zeros(combFact,tUsers);
for iGroup = 1:2^nUsers
    
    binVec = dec2bin(iGroup,nUsers);
    binSum = find(sum(binVec,1) == 49);
    if length(binSum) == tUsers
        userArray(cIndex,:) = binSum;
        cIndex = cIndex + 1;
    end
end
