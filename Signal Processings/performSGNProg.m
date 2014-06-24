function [WSRMaximization_t,PwrVector_n,SNRconstraints_n] = performSGNProg(WSRMaximization_t,SNRconstraints)

sgPhi = 0.2;
oldNorm = 5e10;
P = WSRMaximization_t.P;
H = WSRMaximization_t.H;
W = WSRMaximization_t.W;
nLayers = WSRMaximization_t.nLayers;
pFactor = WSRMaximization_t.pFactor;

while (1)
    
    sgApprox = pFactor .* (SNRconstraints ./ (SNRconstraints + 1));
    
    cvx_begin gp
    
        expressions SgPWR IfPWR ovPWR Xpwr;
        variables Gamma(nLayers,1) pwrFactor(nLayers,1);
        maximize (prod(Gamma.^sgApprox))

        subject to

            (1 - sgPhi) * SNRconstraints <= Gamma;
            -(1 + sgPhi) * SNRconstraints <= Gamma;

            for iLayer = 1:nLayers

                Xpwr = W{iLayer,1}' * H{iLayer,1} * P(:,iLayer);
                SgPWR = pwrFactor(iLayer,1) * norm(Xpwr)^2;
                IfPWR = 0;
                for kLayer = 1:nLayers
                    if kLayer ~= iLayer
                        Xpwr = W{iLayer,1}' * H{iLayer,1} * P(:,kLayer);
                        IfPWR = IfPWR + pwrFactor(kLayer,1) * norm(Xpwr)^2;
                    end
                end

                Gamma(iLayer,1) * SgPWR^(-1) * (1 + IfPWR) <= 1;

            end

            ovPWR = 0;
            for iLayer = 1:nLayers
                ovPWR = pwrFactor(iLayer,1) * norm(P(:,iLayer))^2 + ovPWR;
            end
            
            ovPWR <= WSRMaximization_t.Pt;
            
    cvx_end
    
    newNorm = norm(((Gamma - SNRconstraints) ./ SNRconstraints),1);
    
    if (oldNorm - newNorm) < 1e-2
        break;
    else
        oldNorm = newNorm;
        SNRconstraints = Gamma;
    end
    
end

PwrVector_n = pwrFactor;
SNRconstraints_n = Gamma;
WSRMaximization_t.cvx_status.sgn = cvx_status;

end