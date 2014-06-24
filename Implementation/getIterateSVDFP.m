
function [U1, D1, V1, fU1, fD1, fV1] = getIterateSVDFP(X,nIter,FP,doubleOutput)

global FSD;
FSD = (2^(FP)) - 1;

iIter = 1;
fX = packed_fp(X);
X = X';fX = cplx_herm(fX);

[nR,nC] = size(X);
Q = cell(nIter,1);R = cell(nIter,1);
fQ = cell(nIter,1);fR = cell(nIter,1);

while iIter <= nIter
    
    if iIter == 1
        [Q{iIter,1}, R{iIter,1}, fQ{iIter,1}, fR{iIter,1}] = performQRDecomposition(X,fX,'false');
    else
        [Q{iIter,1}, R{iIter,1}, fQ{iIter,1}, fR{iIter,1}] = performQRDecomposition(R{iIter - 1,1}',cplx_herm(fR{iIter - 1,1}),'false');
    end
    iIter = iIter + 1;
end

U = eye(nC);V = eye(nR);
fU = packed_fp(U);fV = packed_fp(V);

for iIter = 1:nIter
    if mod(iIter - 1,2) == 0
        V = V * Q{iIter,1};
        fV = cplx_mult(fV,fQ{iIter,1});
    else
        U = Q{iIter,1}' * U;
        fU = cplx_mult(cplx_herm(fQ{iIter,1}),fU);
    end
end

if mod(nIter,2)
    D = R{nIter,1}';
    fD = cplx_herm(fR{nIter,1});
else
    D = R{nIter,1};
    fD = fR{nIter,1};
end

V1 = V;U1 = U';D1 = D;
fV1 = fV;fU1 = cplx_herm(fU);fD1 = fD;

if strcmp(doubleOutput,'true')
    fV1 = double(fV1) / FSD;fU1 = double(fU1) / FSD;fD1 = double(fD1) / FSD;
end    
    
end

function [Q,R,fQ,fR] = performQRDecomposition(A,fA,doubleFormat)

global FSD;

if nargin == 1
    fA = packed_fp(A);
    doubleFormat = 'false';
end

[nRows, nCols] = size(A);

R = A;fR = fA;
Q = eye(nRows);fQ = packed_fp(Q);

for iCol = 1:min(nCols,nRows)
    
    colVec = [1 ; zeros(nRows - iCol,1)];fcolVec = packed_fp(colVec);
    
    vecA = R(iCol:end,iCol);fvecA = fR(iCol:end,iCol);
    
    vecNormA = norm(vecA);fvecNormA = fp_fake_norm(fvecA);
    vecB = vecNormA * colVec;fvecB = cplx_mult(fcolVec,fvecNormA);
    
    I = eye(nRows - iCol + 1);fI = packed_fp(I);
    vecW = vecA - vecB;fvecW = cplx_diff(fvecA,fvecB);
    
    H = I - vecW * vecW' * inv(vecW' * vecA);
    
    WWH = cplx_mult(fvecW,cplx_herm(fvecW));
    WHA = cplx_mult(cplx_herm(fvecW),fvecA);
    WHAI = cplx_inv(WHA);
    WWH_WHAI = cplx_mult(WWH,WHAI);
    fH = cplx_diff(fI,WWH_WHAI);
    
    H = blkdiag(eye(iCol - 1),H);
    fH = int16(blkdiag(double(packed_fp(eye(iCol - 1))),double(fH)));
    
    Q = Q * H';R = H * R;
    fQ = cplx_mult(fQ,cplx_herm(fH));fR = cplx_mult(fH,fR);
    
end

if strcmp(doubleFormat,'true')
    fQ = double(fQ) / FSD;
    fR = double(fR) / FSD;
end

end

% FP implementation of basic operations

function [outMatrix] = cplx_mult(x,y)

global FSD;

[nRowsX,nColsX] = size(x);
[nRowsY,nColsY] = size(y);
outMatrix = int16(zeros(nRowsX,nColsY));

if nColsX ~= nRowsY
    if (nRowsY - nColsY) * nRowsY * nColsY == 0
        x32 = int32(x);y32 = int32(y);
        tempFrac = int32(FSD);
        
        for iRow = 1:nRowsX
            for iCol = 1:nColsX
                temp32_r = (real(x32(iRow,iCol)) * real(y32) - imag(x32(iRow,iCol)) * imag(y32));
                temp32_i = (real(x32(iRow,iCol)) * imag(y32) + imag(x32(iRow,iCol)) * real(y32));                
                temp32_r = int16(temp32_r / tempFrac);
                temp32_i = int16(temp32_i / tempFrac);
                outMatrix(iRow,iCol) = int16(complex(temp32_r,temp32_i));
            end
        end
        return;
    else
        display('Mismatched Matrices !');
    end
end

x32 = int32(x);y32 = int32(y);
overFlowProt = int32(2^(ceil(log2(nColsX))));
tempFrac = int32(int32(FSD) / overFlowProt);

for iRow = 1:nRowsX
    for iCol = 1:nColsY
        temp32_r = int32(0);
        temp32_i = int32(0);
        for rIndex = 1:nColsX
            temp32_r = temp32_r + (real(x32(iRow,rIndex)) * real(y32(rIndex,iCol)) - imag(x32(iRow,rIndex)) * imag(y32(rIndex,iCol))) / overFlowProt;
            temp32_i = temp32_i + (real(x32(iRow,rIndex)) * imag(y32(rIndex,iCol)) + imag(x32(iRow,rIndex)) * real(y32(rIndex,iCol))) / overFlowProt;
        end
        
        temp32_r = int16(temp32_r / tempFrac);
        temp32_i = int16(temp32_i / tempFrac);
        outMatrix(iRow,iCol) = int16(complex(temp32_r,temp32_i));
    end
end

end

function [out_val] = packed_fp(x)

global FSD;

out_val = x;
out_val = int16(out_val * FSD);

end

function [out_val] = cplx_norm(x)

global FSD;

x32 = int32(x);
nRows = length(x);
norm_r = int32(0);
norm_i = int32(0);
maxShift = 2^ceil(log2(nRows));

for iRow = 1:nRows
    norm_r = norm_r + int32((real(x32(iRow,1)) * real(x32(iRow,1))) / maxShift);
    norm_i = norm_i + int32((imag(x32(iRow,1)) * imag(x32(iRow,1))) / maxShift);
end

norm_r = int32((int64(norm_r) + int64(norm_i)) / 2);

norm_r = fp_sqrt(norm_r);
out_val = complex(int16(norm_r * maxShift * 2 / (FSD)),int16(0));

end

function [out_val] = fp_sqrt(x)

out_val = int32(sqrt(double(x)));

end

function [out_val] = fp_fake_norm(x)

xD = double(x);
out_val = int16(norm(xD));

end

function [out_val] = cplx_herm(x)

xD = double(x);
out_val = int16(xD');

end

function [out_val] = cplx_inv(x)

global FSD;

x32_r = int32(real(x));
x32_i = int32(imag(x));

num32_r = x32_r * FSD / 2;
num32_i = -x32_i * FSD / 2;

den32 = (x32_r * x32_r / 2) + (x32_i * x32_i / 2);

den32 = den32 / FSD;
out_val = complex(int16(num32_r / den32),int16(num32_i / den32));

end

function [out_val] = cplx_sum(x,y)

out_val = complex(real(x) + real(y),imag(x) + imag(y));

end

function [out_val] = cplx_diff(x,y)

out_val = complex(real(x) - real(y),imag(x) - imag(y));

end
