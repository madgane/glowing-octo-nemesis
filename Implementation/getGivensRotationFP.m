% function [Qmatrix, Rmatrix, varargout] = getGivensRotationFP(X,wordlength,fraclength)
% 
% [nRows,nCols] = size(X);
% Xr = sfi(X,wordlength,fraclength);
% Qmatrix = sfi(eye(nRows),wordlength,fraclength);
% 
% for iCol = 1:nCols
%     for iRow = iCol:nRows
%         if iCol == iRow
%             RotMatrix = sfi(eye(nRows),wordlength,fraclength);
%             thetaOne = -angle(Xr.data(iRow,iCol));
%             RotMatrix(iRow,iCol) = sfi(exp(sqrt(-1) * thetaOne));
%         else
%             RotMatrix = eye(nRows);
%             thetaOne = atan(abs(Xr.data(iRow,iCol))/abs(Xr.data(iCol,iCol)));
%             thetaTwo = -mod(angle(Xr.data(iCol,iCol)),pi);thetaThree = -mod(angle(Xr.data(iRow,iCol)),pi);
%             RotMatrix(iCol,iCol) = quantize(cos(thetaOne) * sfi(exp(sqrt(-1) * thetaTwo)),1,wordlength,fraclength);
%             RotMatrix(iRow,iRow) = quantize(cos(thetaOne) * sfi(exp(sqrt(-1) * thetaThree)),1,wordlength,fraclength);
%             RotMatrix(iRow,iCol) = quantize(-sin(thetaOne) * sfi(exp(sqrt(-1) * thetaTwo)),1,wordlength,fraclength);
%             RotMatrix(iCol,iRow) = quantize(sin(thetaOne) * sfi(exp(sqrt(-1) * thetaThree)),1,wordlength,fraclength);
%         end
%         
%         Xr = quantize(RotMatrix * Xr,1,wordlength,fraclength);
%         Qmatrix = quantize(RotMatrix * Qmatrix,1,wordlength,fraclength);
%         
%     end
% end
% 
% if nargout > 2
%     varargout = cell(1,2);
%     varargout{1,1} = Qmatrix';
%     varargout{1,2} = Xr;
% end
% 
% Rmatrix = Xr.data;
% Qmatrix = Qmatrix.data';

[nRows,nCols] = size(X);
Xr = sfi(X,wordlength,fraclength);
Qmatrix = sfi(eye(nRows),wordlength,fraclength);

for iCol = 1:nCols
    for iRow = iCol:nRows
        if iCol == iRow
            RotMatrix = sfi(eye(nRows),wordlength,fraclength);
            thetaOne = -angle(Xr.data(iRow,iCol));
            RotMatrix(iRow,iCol) = sfi(exp(sqrt(-1) * thetaOne));
        else
            RotMatrix = eye(nRows);
            thetaOne = atan(abs(Xr.data(iRow,iCol))/abs(Xr.data(iCol,iCol)));
            thetaTwo = -mod(angle(Xr.data(iCol,iCol)),pi);thetaThree = -mod(angle(Xr.data(iRow,iCol)),pi);
            RotMatrix(iCol,iCol) = quantize(cos(thetaOne) * sfi(exp(sqrt(-1) * thetaTwo)),1,wordlength,fraclength);
            RotMatrix(iRow,iRow) = quantize(cos(thetaOne) * sfi(exp(sqrt(-1) * thetaThree)),1,wordlength,fraclength);
            RotMatrix(iRow,iCol) = quantize(-sin(thetaOne) * sfi(exp(sqrt(-1) * thetaTwo)),1,wordlength,fraclength);
            RotMatrix(iCol,iRow) = quantize(sin(thetaOne) * sfi(exp(sqrt(-1) * thetaThree)),1,wordlength,fraclength);
        end
        
        Xr = quantize(RotMatrix * Xr,1,wordlength,fraclength);
        Qmatrix = quantize(RotMatrix * Qmatrix,1,wordlength,fraclength);
        
    end
end

if nargout > 2
    varargout = cell(1,2);
    varargout{1,1} = Qmatrix';
    varargout{1,2} = Xr;
end

Rmatrix = Xr.data;
Qmatrix = Qmatrix.data';
