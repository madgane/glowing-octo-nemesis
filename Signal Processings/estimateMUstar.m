function [mu_star Pwr] = estimateMUstar(X,H,U,W,alpha,pTotal)

epsilon = 1e-4;
nK = length(H);
[D L] = eig(X);

phiMat = 0;
for iK = 1:nK
    phiMat = phiMat + alpha(iK,1) * H(:,:,iK)' * U(:,:,iK) * W(:,:,iK)^2 * U(:,:,iK)' * H(:,:,iK);
end

phiX = D' * phiMat * D;
muMax = 1e2;muMin = 0;

while (1)

    muK = (muMax +  muMin) * 0.5;
    Pwr = abs(diag(phiX)) ./ (abs(diag(L)) + muK).^2;
    Pwr = (abs(Pwr) >= 0) .* abs(Pwr);
    
    if sum(Pwr) < pTotal
        muMax = muK;
    else
        muMin = muK;
    end

    sum(Pwr)
    if abs(sum(Pwr) - pTotal) < epsilon
        break;
    end
    
end

mu_star = muK;





