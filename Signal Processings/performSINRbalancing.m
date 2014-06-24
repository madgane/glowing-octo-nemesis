
function [SINRbalancing_t] = performSINRbalancing(SINRbalancing_t)

kUsers = SINRbalancing_t.nUsers;
pFactor = SINRbalancing_t.pFactor;
H = SINRbalancing_t.H;
W = SINRbalancing_t.W;
Pt = SINRbalancing_t.Pt;

[nReceive,nTransmit] = size(H{1,1});
GAMMA_MIN_G = 0;GAMMA_MAX_G = Pt * 100;
Gamma_min = GAMMA_MIN_G;Gamma_max = GAMMA_MAX_G;

while (1)

Gamma = (Gamma_min + Gamma_max) / 2;
qGamma = sqrt(1 + (1 ./ (Gamma.^pFactor)));

cvx_begin
  
    variable X(nTransmit,kUsers) complex;
    minimize 0
    subject to

        for iUser = 1:kUsers
            {[X' * H{iUser,1}' * W{iUser,1} ;1], ...
                qGamma(iUser,1) * W{iUser,1}' * H{iUser,1} * X(:,iUser)} <In> complex_lorentz(kUsers + 1);
        end
        
        norm(X(:)) <= sqrt(Pt);
       
cvx_end

if strcmp(cvx_status,'Solved')
    Gamma_min = Gamma;
else
    Gamma_max = Gamma;
end

if strcmp(cvx_status,'Solved')
    for iUser = 1:kUsers
        W{iUser,1} = calculateMMSEWvector(H,X,iUser);
    end
end

if strcmp(cvx_status,'Solved')
    if abs(Gamma_max - Gamma_min) < 1e-3
        break;
    end
end

end

SINRbalancing_t.W = W;
SINRbalancing_t.P = sqrt(Pt / trace(X' * X)) * X;

end




