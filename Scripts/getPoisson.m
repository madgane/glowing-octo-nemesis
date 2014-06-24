
function [poisRndVar] = getPoisson(lda,nRows,nCols)

poisRndVar = zeros(nRows,nCols);

for iRow = 1:nRows    
    cLambda = lda(1,iRow);
    for iCol = 1:nCols
        U = rand;j = 0;p = exp(-cLambda);F = p;
        while U > F
            p = cLambda * p / (j + 1);
            F = F + p;
            j = j + 1;
        end
        poisRndVar(iRow,iCol) = j;
    end
end
        
