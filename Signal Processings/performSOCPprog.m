function [WSRMaximization_t] = performSOCPprog(WSRMaximization_t,SNRconstraints)

P = WSRMaximization_t.P;
H = WSRMaximization_t.H;
W = WSRMaximization_t.W;
nLayers = WSRMaximization_t.nLayers;

[nTransmit,~] = size(P);
qGamma = sqrt(1 + (1./SNRconstraints));

cvx_begin

    variable X(nTransmit,nLayers) complex;
    variable rho
    
    minimize (rho)
    
    subject to
    
        for iLayer = 1:nLayers
            {[X' * H{iLayer,1}' * W{iLayer,1} ;1], ...
                qGamma(iLayer,1) * W{iLayer,1}' * H{iLayer,1} * X(:,iLayer)} <In> complex_lorentz(nLayers + 1);
        end
        
        norm(X(:)) <= sqrt(WSRMaximization_t.Pt) * rho;

cvx_end

WSRMaximization_t.P = X / rho;
WSRMaximization_t.cvx_status.socp = cvx_status;

end