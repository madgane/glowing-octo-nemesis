function [mu_star Pwr] = bisectionEstimateMU(X,C,pTotal)

epsilon = 1e-5;

iIter = 0;
maxIter = 1e5;

[D L] = eig(X);
phiX = D' * C * D;

muK = 0;
Pwr = real(diag(phiX)) ./ (real(diag(L)) + muK).^2;
Pwr = (Pwr >= 0) .* Pwr;

if sum(Pwr) <= (pTotal + epsilon)
    continueAgain = 0;
else
    continueAgain = 1;
    muMax = 1e10; muMin = 0;muK = muMax;
end

while continueAgain

    muK = (muMax +  muMin) * 0.5;
    Pwr = real(diag(phiX)) ./ (real(diag(L)) + muK).^2;
    Pwr = (Pwr >= 0) .* Pwr;
    
    if sum(Pwr) < pTotal
        muMax = muK;
    else
        muMin = muK;
    end

    iIter = iIter + 1;
    if abs(sum(Pwr) - pTotal) < epsilon
        continueAgain = 0;
    end
    
    if iIter > maxIter
        continueAgain = 0;
    end
    
end

if sum(Pwr < 0)
    display('Warning ! Negative Powers !');
end

mu_star = muK;





