
function [SimParams,SimStructs] = getWSRMMatrix(SimParams,SimStructs)

nSB = SimParams.nBands;
nBS = SimParams.nBases;
cH = SimStructs.linkChan;

for iSB = 1:nSB
    
    currentUsers = cell(nBS,1);
    for iBase = 1:SimParams.nBases
        currentUsers{iBase,1} = SimStructs.baseStruct{iBase,1}.assignedUsers{iSB,1};
    end
    
    nUsers = 0;
    Queues = zeros(nUsers,nBS);
    for iBase = 1:nBS
        for iUser = 1:length(currentUsers{iBase,1})
            cUser = currentUsers{iBase,1}(iUser,1);
            Queues(iUser,iBase) = SimStructs.userStruct{cUser,1}.trafficStats.backLogPkt;
        end
        
        nUsers = max(nUsers,length(currentUsers{iBase,1}));
    end
    
    switch SimParams.weightedSumRateMethod
            
        case 'WSRMApproach'
            
            qWeights = Queues;
            re_iterate = 1;cvx_hist = 100;
            m_k_n_o = -rand(nUsers,nBS) * 10;            
            
            while re_iterate
                
                cvx_begin
                
                expression m_k_n(nUsers,nBS)
                variable M_k_n(SimParams.nTxAntenna,nUsers,nBS) complex
                variables t_k_n(nUsers,nBS) b_k_n(nUsers,nBS) g_k_n(nUsers,nBS)
                variables obj_var(nBS,1) obj_func_var
                
                maximize(obj_func_var);
                
                subject to
                
                for iBase = 1:nBS
                    sum(qWeights(:,iBase) .* t_k_n(:,iBase)) >= obj_var(iBase,1);
                end
                
                obj_func_var <= sum(obj_var);
                
                for iBase = 1:nBS
                    for iUser = 1:length(currentUsers{iBase,1})
                        
                        if_vector = sqrt(SimParams.N);
                        cUser = currentUsers{iBase,1}(iUser,1);
                        
                        for jBase = 1:nBS
                            cCH = cH{jBase,iSB}(:,:,cUser);
                            if iBase ~= jBase                                
                                for jUser = 1:length(currentUsers{jBase,1})
                                    if_vector = [if_vector ; cCH * M_k_n(:,jUser,jBase)];
                                end
                            else
                                for jUser = 1:length(currentUsers{jBase,1})
                                    if jUser ~= iUser
                                        if_vector = [if_vector ; cCH * M_k_n(:,jUser,jBase)];
                                    end
                                end
                            end
                        end
                        
                        norm(if_vector,2) <= sqrt(b_k_n(iUser,iBase));
                        log(1 + g_k_n(iUser,iBase)) >= t_k_n(iUser,iBase);
                        
                        cCH = cH{iBase,iSB}(:,:,cUser);
                        imag(cCH * M_k_n(:,iUser,iBase)) == 0;
                        m_k_n(iUser,iBase) = g_k_n(iUser,iBase) - b_k_n(iUser,iBase);
                        
                        4 * real(cCH * M_k_n(:,iUser,iBase)) + m_k_n_o(iUser,iBase)^2 + 2 * m_k_n_o(iUser,iBase) * (m_k_n(iUser,iBase) - m_k_n_o(iUser,iBase)) ...
                            >= (g_k_n(iUser,iBase) + b_k_n(iUser,iBase))^2;
                                                                        
                    end
                    
                    {vec(M_k_n(:,:,iBase)),sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iSB))} == complex_lorentz(length(vec(M_k_n(:,:,iBase))));
                end
                
                cvx_end
                
                if strcmp(cvx_status,'Solved')
                    m_k_n_o = g_k_n - b_k_n;
                    if abs(cvx_optval - cvx_hist) <= 1e-3
                        re_iterate = 0;
                    else
                        cvx_hist = cvx_optval;
                    end
                else
                    m_k_n_o = m_k_n_o * 2;
                end
                
            end
 
        otherwise
            
            display('Undefined Optimization Approach !');
            
    end
    
    for iBase = 1:nBS
        SimStructs.baseStruct{iBase,1}.P{iSB,1} = M_k_n(:,:,iBase);
    end

    
end

