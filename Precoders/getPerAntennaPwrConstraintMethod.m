function [SimParams,SimStructs] = getPerAntennaPwrConstraintMethod(SimParams,SimStructs)

maxIters = 5000;
iterCapacity = zeros(maxIters,1);

for iBand = 1:SimParams.nBands
    for iBase = 1:SimParams.nBases
        
        rankH = min(SimParams.nTxAntenna,SimParams.nRxAntenna);
        cUser = SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1};
        
        re_iterate = 1;
        cLambda = inv(diag(ones(SimParams.nTxAntenna,1)));
        H = SimStructs.linkChan{iBase,iBand}(:,:,cUser(1,1));
        perAntPwr = SimStructs.baseStruct{iBase,1}.sPower(1,iBand) / SimParams.nTxAntenna;
        powerProfile = repmat(perAntPwr,SimParams.nTxAntenna,1);
        
        %         powerProfile = rand(SimParams.nTxAntenna,1);
        %         powerProfile = SimStructs.baseStruct{iBase,1}.sPower(1,iBand) * powerProfile / sum(powerProfile);
        
        switch SimParams.weightedSumRateMethod
            
            case 'IterApproach'
                
                iIter = 0;
                alpha = 1;
                H = H * sqrt(alpha);
                
                while re_iterate
                    
                    sqCLambda = sqrt(cLambda);
                    Htilde = H * sqCLambda;
                    [~,Lt,Vt] = svd(Htilde);
                    
                    if SimParams.nRxAntenna == 1
                        LM = 1./Lt(1,1).^2;
                    else
                        LM = 1./diag(Lt).^2;
                    end
                    
                    LtM = ones(SimParams.nTxAntenna,1); LtM(1:length(LM),1) = LM;
                    LtM = LtM .* (LtM <= 1) + ones(SimParams.nTxAntenna,1) .* (LtM > 1);
                    
                    Psi = real(diag(Vt * diag(LtM) * Vt'));
                    nLambda = (1 / alpha) * powerProfile + Psi .* diag(cLambda);
                    
                    EM = zeros(SimParams.nTxAntenna,1);
                    EM(1:length(LM)) = max(1 - LM,0);
                    
                    Stilde = Vt * diag(EM) * Vt';
                    S = diag(nLambda)^0.5 * Stilde * diag(nLambda)^0.5;
                    
                    
                    [V,E] = eig(S);V = V * sqrt(E);
                    deviation = real(trace(diag(nLambda)^-1 * diag((powerProfile - diag(S)))));
                    if abs(deviation) <= 1e-4
                        re_iterate = 0;
                    else
                        cLambda = diag(nLambda);
                    end
                    
                    iIter = iIter + 1;
                    iterCapacity(iIter,1) = real(log(det(eye(rankH) + H * S * H')));
                    
                end
                
            case 'MSApproach'
                
                iIter = 0;                
                while re_iterate
                    
                    sqCLambda = sqrt(cLambda);
                    Htilde = H * sqCLambda;
                    [~,Lt,Vt] = svd(Htilde);
                    
                    if SimParams.nRxAntenna == 1
                        LM = 1./Lt(1,1).^2;
                    else
                        LM = 1./diag(Lt).^2;
                    end
                    
                    LtM = ones(SimParams.nTxAntenna,1); LtM(1:length(LM),1) = LM;
                    LtM = LtM .* (LtM <= 1) + ones(SimParams.nTxAntenna,1) .* (LtM > 1);
                    
                    Psi = real(diag(Vt * diag(LtM) * Vt'));                   
                    F = ones(SimParams.nTxAntenna,1) - Psi;
                    randK = repmat(min(Psi),SimParams.nTxAntenna,1);
                    
                    Psi_Plus = randK + F;Psi_Minus = randK;
                    nLambda = (diag(Psi_Plus)^-1) * (powerProfile + diag(Psi_Minus) * diag(cLambda));
                    
                    EM = zeros(SimParams.nTxAntenna,1);
                    EM(1:length(LM)) = max(1 - LM,0);
                    
                    Stilde = Vt * diag(EM) * Vt';
                    S = diag(nLambda)^0.5 * Stilde * diag(nLambda)^0.5;                    
                    
                    [V,E] = eig(S);V = V * sqrt(E);
                    deviation = real(trace(diag(nLambda)^-1 * diag((powerProfile - diag(S)))));
                    if abs(deviation) <= 1e-4
                        re_iterate = 0;
                    else
                        cLambda = diag(nLambda);
                    end
                    
                    iIter = iIter + 1;
                    iterCapacity(iIter,1) = real(log(det(eye(rankH) + H * S * H')));
                    
                end
                
                
            case 'SVDApproach'
                
                [~,D,V] = svd(H);
                if SimParams.nRxAntenna == 1
                    D = [D(1,1) ; zeros(SimParams.nTxAntenna - rankH,1)];
                else
                    D = [diag(D) ; zeros(SimParams.nTxAntenna - rankH,1)];
                end
                
                [V,~] = performWFAlgorithm(V,sum(powerProfile),D);
                
            case 'OptApproach'
                
                cvx_begin sdp
                
                expression A
                variable S(SimParams.nTxAntenna,SimParams.nTxAntenna) complex hermitian
                
                maximize(log_det(eye(rankH) + H * S * H'))
                
                subject to
                
                S >= 0;
                diag(S) <= powerProfile;
                
                cvx_end
                
                [V, D] = eig(S);
                V = V * sqrt(D);
                
            case 'AltOptApproach'
                                
                Szz = diag(powerProfile)^(-1);
                Ht = H * Szz^(-0.5);
                
                cvx_begin sdp quiet
                
                variable t
                variable Sxx(SimParams.nTxAntenna,SimParams.nTxAntenna) complex hermitian
                
                maximize(t);
                
                log_det(Ht.' * Sxx * Ht + eye(SimParams.nTxAntenna)) >= t;
                
                subject to
                
                real(trace(Sxx)) <= 1;
                                
                cvx_end
                
                cvx_begin sdp quiet

                variable t
                variable Szz(SimParams.nTxAntenna,SimParams.nTxAntenna) complex hermitian

                minimize(t);
                
                log_det(H.' * Sxx * H + Szz) - log_det(Szz) >= t;
                
                subject to
                
                trace(Szz * diag(powerProfile)) <= 1;
                
                
                
                cvx_begin sdp
                
                cvx_end
                
                
                
        end
        
        SimStructs.baseStruct{iBase,1}.P{iBand,1} = V;
        SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1} = repmat(cUser,1,SimParams.nTxAntenna);
        SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1} = (1:SimParams.nTxAntenna);
        
        if SimParams.iSNR == 1
            if SimParams.iDrop == 1
                SimParams.Debug.tempResource{1,1} = zeros(maxIters,length(SimParams.snrIndex));
            end
        end
        
        if SimParams.iDrop == 1
            SimParams.Debug.tempResource{1,1}(:,SimParams.iSNR) = iterCapacity;
        else
            SimParams.Debug.tempResource{1,1}(:,SimParams.iSNR) = iterCapacity + SimParams.Debug.tempOne(:,SimParams.iSNR);
        end
        
%         hold all;
%         plot(iterCapacity(iterCapacity ~= 0),'o-.');keyboard;
        
    end
end


end
