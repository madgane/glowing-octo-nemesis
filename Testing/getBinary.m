
function [matrixA] = getBinary(lengthA,combA)

rIndex = 1;
for iIndex = 1:(2^(lengthA) - 1)
    
    arrayA = zeros(lengthA,1);
    A = dec2bin(iIndex,lengthA);
    arrayA(:) = A(:);arrayA = (arrayA == 49) * 1;
    
    if (sum(arrayA) == combA)
        matrixA(:,rIndex) = arrayA;
        rIndex = rIndex + 1;
    end    
end

end
