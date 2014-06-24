
function [Pmatrix, Pavg] = performWFAlgorithm(PinvMatrix,Pt,Gains)

iIteration = 0;
maxIteration = 1e5;

if nargin == 2
    Ppow = diag(PinvMatrix' * PinvMatrix);
else
    Ppow = pinv(Gains)';
end

if isempty(PinvMatrix)
    Pavg = 0;
    Pmatrix = PinvMatrix;
    return;
end

Pbound(1,1) = 1e20;Pbound(2,1) = 0;Psum = 0;
Ppow = Ppow .* (Ppow > 0) + 1e25 .* (Ppow == 0);

while (1)    
   
    Pavg = mean(Pbound) - Ppow;
    Pavg = Pavg .* (Pavg > 0);
    Psum = sum(Pavg);
        
    if Psum >= Pt
        Pbound(1,1) = mean(Pbound);
    else
        Pbound(2,1) = mean(Pbound);
    end
    
    if abs(Psum - Pt) < 1e-5
        break;
    end  
    
    if iIteration > maxIteration
        break;
    else
        iIteration = iIteration + 1;
    end    
    
end

% Pavg = ones(size(Pavg)) * Pt / length(Pavg);

Pmatrix = PinvMatrix;
avgPWR = (1./sqrt(diag(Pmatrix' * Pmatrix)));
for iCol = 1:length(avgPWR)
    if (avgPWR(iCol,1) ~= Inf)
        Pmatrix(:,iCol) = Pmatrix(:,iCol) * avgPWR(iCol,1);
    else
        Pmatrix(:,iCol) = Pmatrix(:,iCol) * 0;
    end
end

Pmatrix = Pmatrix * diag(sqrt(Pavg));

end
