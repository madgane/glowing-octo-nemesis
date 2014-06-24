
function [U1, D1, V1] = iterativeSVD_testFile(X,nIter)

iIter = 1;
X = X';

[nR,nC] = size(X);
Q = cell(nIter,1);R = cell(nIter,1);

while iIter <= nIter
    
    if iIter == 1
        [Q{iIter,1}, R{iIter,1}] = performQRDecomposition(X);
    else
        [Q{iIter,1}, R{iIter,1}] = performQRDecomposition(R{iIter - 1,1}');
    end
    iIter = iIter + 1;
end

U = eye(nC);V = eye(nR);

for iIter = 1:nIter
    if mod(iIter - 1,2) == 0
        V = V * Q{iIter,1};
    else
        U = Q{iIter,1}' * U;
    end
end

if mod(nIter,2)
    D = R{nIter,1}';
else
    D = R{nIter,1};
end

V1 = V;U1 = U';D1 = D;
    
end

function [Q,R] = performQRDecomposition(A)

[nRows, nCols] = size(A);

R = A;
Q = eye(nRows);

for iCol = 1:min(nCols,nRows)
    
    colVec = [1 ; zeros(nRows - iCol,1)];
    
    vecA = R(iCol:end,iCol);
    
    vecNormA = norm(vecA);
    vecB = vecNormA * colVec;
    
    I = eye(nRows - iCol + 1);
    vecW = vecA - vecB;
    
    H = I - vecW * vecW' * inv(vecW' * vecA);
        
    H = blkdiag(eye(iCol - 1),H);    
    Q = Q * H';R = H * R;

end

end
