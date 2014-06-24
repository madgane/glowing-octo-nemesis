
function [Pmatrix Pavg] = performQueuedWF(PinvMatrix,Pt,Q)

pAllocationType = 'WSR_ZF';

Q = Q / log2(exp(1));
[~,nUsers] = size(PinvMatrix);
pwrUsers = diag(PinvMatrix' * PinvMatrix);

if sum(Q) == 0
    Pmatrix = zeros(size(PinvMatrix));
    return;
end

switch pAllocationType
    
    case 'WSR_ZF'
        
        Pbound(1,1) = 1e20;Pbound(2,1) = 0;Psum = 0;
        pwrUsers = pwrUsers .* (pwrUsers > 0) + 1e25 .* (pwrUsers == 0);
        
        while (1)
            
            Pavg = mean(Pbound) * Q - pwrUsers;
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
            
        end
        
        cvxP = Pavg;
        
    case 'DifferentialQueues'
        
        cvx_begin
        
        variable cvxP(nUsers,1)
        variable t(nUsers,1)
        minimize(norm(t,2))
        
        subject to
        norm(cvxP,1) <= Pt;
        max(Q - log(1 + cvxP ./ pwrUsers),0) <= t;
        
        cvx_end
        
end

Pavg = cvxP;
Pmatrix = PinvMatrix;
avgPWR = 1./sqrt(pwrUsers);

for iCol = 1:length(avgPWR)
    if (avgPWR(iCol,1) ~= Inf)
        Pmatrix(:,iCol) = Pmatrix(:,iCol) * avgPWR(iCol,1);
    else
        Pmatrix(:,iCol) = Pmatrix(:,iCol) * 0;
    end
end

Pmatrix = Pmatrix * diag(sqrt(Pavg));
