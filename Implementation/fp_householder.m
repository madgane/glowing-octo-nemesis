function [Q,R] = fp_householder(A,FP)

global FSD;

FSD = (2^(FP)) - 1;

fA = packed_fp(A);

% Compute the QR decomposition of an m-by-n matrix A using
% Householder transformations.

[m,n] = size(A);

R = A; fR = fA;

for j = 1:min(m,n)-1
    
    I = eye((m+1-j));
    Q(:,:,j) = eye(m);
    
    fI = packed_fp(I);
    fQ(:,:,j) = packed_fp(eye(m));
    
    normx = norm(R(j:end,j));
    fnormx = fp_fake_norm(fR(j:end,j));
    
    u = R(j:end,j);
    u(1) = u(1)- exp(i*angle(R(j,j))) * normx;
    
    fu = fR(j:end,j);
    fu(1) = int16(double(fu(1)) - double(cplx_mpy(packed_fp(exp(1i*angle(double(fR(j,j))))),fnormx)));
    
    %v = u/||u||
    
    v = u/norm(u);
    fv = cplx_div_int(fu,fp_fake_norm(fu));
    
    %Q = I -(1 + (x'v)/(v'x))*vv'
    
    fvm = cplx_mpy(cplx_conj(fu),fv);
    fvu  = cplx_mpy(cplx_conj(fv),fu);
    fvv = cplx_mpy(fv,cplx_conj(fv));
    
    fQ(j:end,j:end,j) = fI - cplx_mpy(( packed_fp(1) + (fvm / fvu)),fvv);
    fR(j:end,j:end) = cplx_mpy(fQ(j:end,j:end,j) * fR(j:end,j:end));
    
    Q(j:end,j:end,j) = I -(1+((u'*v)/(v'*u)))*v*v';
    R(j:end,j:end) = Q(j:end,j:end,j)*R(j:end,j:end);
    
end

Q_tmp = Q(:,:,1);
fQ_tmp = fQ(:,:,1);

for j = 2:min(m,n)-1
    Q_tmp = Q_tmp * Q(:,:,j);
    fQ_tmp = cplx_mpy(fQ_tmp,fQ(:,:,j));
end

Q = Q_tmp;
R = (Q' * A);

end









% FP implementation of basic operations

function [outMatrix] = cplx_mpy(x,y)

global FSD;

[nRowsX,nColsX] = size(x);
[nRowsY,nColsY] = size(y);
outMatrix = zeros(nRowsX,nColsY);

if nColsX ~= nRowsY
    warning(1,'Dimension mismatch');
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
        outMatrix(iRow,iCol) = complex(temp32_r,temp32_i);
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
    norm_r = norm_r + int32((real(x32(iRow,1)) .* real(x32(iRow,1))) ./ maxShift);
    norm_i = norm_i + int32((imag(x32(iRow,1)) .* imag(x32(iRow,1))) ./ maxShift);
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

function [out_val] = cplx_div_int(x,v)

global FSD;

v = int32(v);
x32 = int32(x);
x32_r = real(x32);x32_i = imag(x32);

x32_r = int16((x32_r * FSD) / v);
x32_i = int16((x32_i * FSD) / v);

out_val = complex(x32_r,x32_i);

end

function [out_val] = cplx_conj(x)

xD = double(x);
out_val = int16(xD');

end

