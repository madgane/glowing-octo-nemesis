function [Ainv] = recursiveInverse(A)

[mSize,~] = size(A);

if mod(mSize,2) == 0
    
else
end


for iIndex = 1 : 4 : mSize
    
    S1 = iIndex;
    E1 = S1 + 1;
    S2 = E1 + 1;
    E2 = S2 + 1;
    
    U1 = triangleInv_2(A(S1:E1,S1:E1));
    U2 = triangleInv_2(A(S2:E2,S2:E2));
    U12 = -U1 * A(S1:E1,S2:E2) * U2;
    
    Ainv(S1:E1,S1:E1) = U1;
    Ainv(S2:E2,S2:E2) = U2;
    Ainv(S1:E1,S2:E2) = U12;
    
end

end

function [Xinv] = triangleInv_2(X)

    Xinv = zeros(size(X));
    Xinv(1,1) = 1 / X(1,1);
    Xinv(2,2) = 1 / X(2,2);
    Xinv(1,2) = -Xinv(1,1) * X(1,2) * Xinv(2,2);  
    
end
