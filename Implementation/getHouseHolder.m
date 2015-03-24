
function [Q,R] = getHouseHolder(varargin)

if nargin == 1
    A = varargin{1};
    algoType = 0;
else
    A = varargin{1};
    algoType = varargin{2};
end

if (algoType == 1)
    
    
    R = A;
    [M,N] = size(A);
    
    for iCol = 1:min(M,N)-1
        
        q = eye(min(M,N));
        X = R(iCol:end,iCol);
        U = X;
        
        normX = norm(X);
        U(1,1) = U(1,1) - normX;
        V = U;
        
        normU = norm(V)^2;
        fproductA = (1 + ((X' * V)/(V' * X))) / normU;
        
        q(iCol:end,iCol:end) = eye(length(V)) - fproductA * (V * V');
        R = q * R;
        
        if iCol == 1
            Q = q;
        else
            Q = q * Q;
        end
        
    end
    
else
    
    F = 2^12  - 1;
    QDisp = @(x)(display(x));
    QDispY = @(x,y)(display(x));
    
    [M,N] = size(A);
    
    R = A;
    Q = eye(M,M);
    pMatrix = zeros(M,M);
    wVec = zeros(M,1);
    
    for iCol = 1:N
        
        bProduct = 0;
        aProduct = abs(R(iCol,iCol))^2;
        QDisp(aProduct);
        
        for iRow = (iCol + 1):M
            
            wVec(iRow,1) = R(iRow,iCol);
            xProduct = abs(R(iRow,iCol))^2;
            bProduct = bProduct + xProduct;
            
        end
        
        xProduct = aProduct + bProduct;
        xProduct = sqrt(xProduct);
        
        wVec(iCol,1) = R(iCol,iCol) - xProduct;
        aProduct = abs(wVec(iCol,1))^2 + bProduct
        temp = wVec(iCol,1) * conj(R(iCol,iCol)) + bProduct
        xProduct = real(temp * conj(temp))
        
        temp1 = temp * temp;
        aProduct = xProduct * aProduct
        
        temp1 = (temp1 + xProduct);
        temp = -temp1 / aProduct
        
        if iCol == 1
            xRow = N;
        else
            xRow = iCol - 1;
            pMatrix(xRow,xRow) = 1;
        end
        
        for mRow = iCol:M
            
            pMatrix(xRow,mRow) = 0;
            pMatrix(mRow,xRow) = 0;
            temp2 = temp * wVec(mRow,1)
            
            for mCol = iCol : N
                
                temp1 = temp2 * conj(wVec(mCol,1))
                pMatrix(mRow,mCol) = temp1;
                
            end
            
            pMatrix(mRow,mRow) = 1 + pMatrix(mRow,mRow);
            
        end
        
        QDispY(pMatrix,15);
        
        if iCol == N
            R = (pMatrix * R)';
            Q = pMatrix * Q;
        else
            R = pMatrix * R;
            Q = pMatrix * Q;
        end
        
    end
    
end

