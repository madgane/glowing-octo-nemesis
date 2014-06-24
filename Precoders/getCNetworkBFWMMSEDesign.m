
function [SimParams SimStructs] = getCNetworkBFWMMSEDesign(SimParams,SimStructs)

iIter = 0;
maxIter = 1e4;
epsilonCheck = min(1e-4,max(SimParams.sPower)^(-2));
nStreams = min(SimParams.maxRank,SimParams.nRxAntenna);

SumCapacity = cell(SimParams.nBands,1);

for iBand = 1:SimParams.nBands
   
    continueAgain = 1;
    W = cell(SimParams.nUsers,SimParams.nBases);
    U = cell(SimParams.nUsers,SimParams.nBases);
    V = cell(SimParams.nUsers,SimParams.nBases);
    linkChannel = SimStructs.linkChan;
    
    for iBase = 1:SimParams.nBases
        for iUser = 1:SimParams.nUsers
            V{iUser,iBase} = complex(ones(SimParams.nTxAntenna,nStreams),ones(SimParams.nTxAntenna,nStreams));
            V{iUser,iBase} = sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand) / (SimParams.nUsers / SimParams.nBases)) * V{iUser,iBase} / trace(V{iUser,iBase}' * V{iUser,iBase});
        end
    end
    
    while continueAgain
        
        % U Matrix calculation
        
        for iBase = 1:SimParams.nBases
            for iUser = 1:SimParams.nUsers  
                
                J = eye(SimParams.nRxAntenna) * SimParams.N;
                
                for jBase = 1:SimParams.nBases
                    for jUser = 1:SimParams.nUsers
                        HV = linkChannel{jBase,iBand}(:,:,iUser) * V{jUser,jBase};
                        J = J + HV * HV';
                    end
                end
                
                H = linkChannel{iBase,iBand}(:,:,iUser);
                
                HdVd = H * V{iUser,iBase};
                U{iUser,iBase} = J \ HdVd;
                W{iUser,iBase} = inv(eye(nStreams) - U{iUser,iBase}' * H * V{iUser,iBase});
                                
            end
        end
        
        for iBase = 1:SimParams.nBases            
            linkedUsers = (1:SimParams.nUsers)';
            
            Isum = 0;Dsum = 0;
            for jBase = 1:SimParams.nBases
                for iUser = 1:SimParams.nUsers
                    cUser = SimStructs.userStruct{iUser,1};
                    H_HU = linkChannel{iBase,iBand}(:,:,iUser)' * U{iUser,jBase};
                    Isum = Isum + cUser.weighingFactor * H_HU * W{iUser,jBase} * H_HU';
                    if jBase == iBase
                        W_2 = W{iUser,iBase} * W{iUser,iBase};
                        Dsum = Dsum + cUser.weighingFactor^2 * H_HU * W_2 * H_HU';
                    end
                end
            end
            
            mu_star = bisectionEstimateMU(Isum,Dsum,SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
            Isum = Isum + mu_star * eye(SimParams.nTxAntenna);
            
            Iinv = pinv(Isum);
            for iUser = 1:length(linkedUsers)
                cIndex = linkedUsers(iUser,1);
                cUser = SimStructs.userStruct{cIndex,1};
                V{cIndex,iBase} = cUser.weighingFactor * Iinv * linkChannel{iBase,iBand}(:,:,cIndex)' * U{cIndex,iBase} * W{cIndex,iBase};
            end
            
        end
        
        if ~iIter
            continueAgain = 1;
        else
            currDeviation = 0;
            for iBase = 1:SimParams.nBases
                for iUser = 1:SimParams.nUsers
                    currDeviation = currDeviation + abs(log(det(W{iUser,iBase})) - log(det(W_prev{iUser,iBase})));
                end
            end
            
            if currDeviation < epsilonCheck
                continueAgain = 0;
            end
            
            if iIter > maxIter
                continueAgain = 0;
                display('Lack of Convergence !');
            end
        end
        
        W_prev = W;
        iIter = iIter + 1;
        SumCapacity{iBand,1} = [SumCapacity{iBand,1} ; performMockReception(SimParams,SimStructs,V,iBand)];
        
    end
    
    % Assigning the V and U to the corresponding users
    
    Pmatrix = cell(SimParams.nBases,1);
    linkUsers = cell(SimParams.nBases,1);
    
    for iUser = 1:SimParams.nUsers
        
        txPower = zeros(SimParams.nBases,1);
                
        for iBase = 1:SimParams.nBases
            txPower(iBase,1) = trace(V{iUser,iBase} * V{iUser,iBase}');
        end
        
        [~,linkedBS] = max(txPower);        
        linkUsers{linkedBS,1} = [linkUsers{linkedBS,1} ; iUser];   
        Pmatrix{linkedBS,1} = [Pmatrix{linkedBS,1} V{iUser,linkedBS}];
        SimStructs.userStruct{iUser,1}.W{iBand,1} = U{iUser,linkedBS};
        SimStructs.userStruct{iUser,1}.baseNode = linkedBS;
        SimStructs.userStruct{iUser,1}.neighNode = find(linkedBS ~= 1:SimParams.nBases);
             
    end
    
    for iBase = 1:SimParams.nBases
        SimStructs.baseStruct{iBase}.P{iBand,1} = Pmatrix{iBase,1};
        SimStructs.baseStruct{iBase}.assignedUsers{iBand,1} = linkUsers{iBase,1};
    end    
    
end

