function [Qmatrix, Rmatrix] = getGivensRotation(X)

[nRows,nCols] = size(X);

Xr = X;
Qmatrix = eye(nRows);

for iCol = 1:nCols
    for iRow = iCol:nRows
        if iCol == iRow
            RotMatrix = eye(nRows);
            thetaOne = -angle(Xr(iRow,iCol));
            RotMatrix(iRow,iCol) = exp(sqrt(-1) * thetaOne);
        else
            RotMatrix = eye(nRows);
            thetaOne = atan(abs(Xr(iRow,iCol))/abs(Xr(iCol,iCol)));
            thetaTwo = -angle(Xr(iCol,iCol));thetaThree = -angle(Xr(iRow,iCol));
            RotMatrix(iCol,iCol) = cos(thetaOne) * exp(sqrt(-1) * thetaTwo);
            RotMatrix(iRow,iRow) = cos(thetaOne) * exp(sqrt(-1) * thetaThree);
            RotMatrix(iRow,iCol) = -sin(thetaOne) * exp(sqrt(-1) * thetaTwo);
            RotMatrix(iCol,iRow) = sin(thetaOne) * exp(sqrt(-1) * thetaThree);
        end
        
        Xr = RotMatrix * Xr;
        Qmatrix = RotMatrix * Qmatrix;
        
    end
end

Rmatrix = Xr;
Qmatrix = Qmatrix';

end