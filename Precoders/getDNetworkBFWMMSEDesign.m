
function [SimParams SimStructs] = getDNetworkBFWMMSEDesign(SimParams,SimStructs)

iIter = 0;
maxIter = 1e4;
epsilonCheck = min(1e-4,max(SimParams.sPower)^(-2));
nStreams = min(SimParams.maxRank,SimParams.nRxAntenna);

for iBand = 1:SimParams.nBands
    
    continueAgain = 1;
    W = cell(SimParams.nUsers,1);
    U = cell(SimParams.nUsers,1);
    V = cell(SimParams.nUsers,SimParams.nBases);
    Haug = cell(SimParams.nUsers,SimParams.nBases);
    
    linkChannel = SimStructs.linkChan;
    Oprecoder = ones(SimParams.nTxAntenna,SimParams.maxRank);
    
    for iUser = 1:SimParams.nUsers
        for iBase = 1:SimParams.nBases
            V{iUser,iBase} = complex(Oprecoder,Oprecoder);
            V{iUser,iBase} = sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand) / (SimParams.nUsers / SimParams.nBases)) * V{iUser,iBase} / trace(V{iUser,iBase}' * V{iUser,iBase});
            
            Haug{iUser,iBase} = linkChannel{iBase,iBand}(:,:,iUser);
        end
    end
    
    while continueAgain
        
        % U Matrix calculation
        
        for iUser = 1:SimParams.nUsers

            J = eye(SimParams.nRxAntenna) * SimParams.N;
            
            for jUser = 1:SimParams.nUsers
                for iBase = 1:SimParams.nBases
                    HV = Haug{iUser,iBase} * V{jUser,iBase};
                    J = J + HV * HV';
                end
            end
            
            Utemp = 0;Wtemp = 0;
            for iBase = 1:SimParams.nBases
                Wtemp = Wtemp + Haug{iUser,iBase} * V{iUser,iBase};
                Utemp = Utemp + inv(J) * Haug{iUser,iBase} * V{iUser,iBase};                
            end
            
            U{iUser,1} = inv(J) * Utemp;
            W{iUser,1} = inv(eye(nStreams) - U{iUser,1}' * Wtemp);
            
        end
        
        for iBase = 1:SimParams.nBases
            
            Isum = 0;Dsum = 0;
            for iUser = 1:SimParams.nUsers
                cUser = SimStructs.userStruct{iUser,1};
                H_HU = Haug{iUser,iBase}' * U{iUser,1};
                Isum = Isum + cUser.weighingFactor * H_HU * W{iUser,1} * H_HU';
                W_2 = W{iUser,1} * W{iUser,1};
                Dsum = Dsum + cUser.weighingFactor^2 * H_HU * W_2 * H_HU';
            end
            
            mu_star = bisectionEstimateMU(Isum,Dsum,SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
            Isum = Isum + mu_star * eye(SimParams.nTxAntenna);

            Iinv = pinv(Isum);
            for iUser = 1:SimParams.nUsers
                cUser = SimStructs.userStruct{iUser,1};
                V{iUser,iBase} = cUser.weighingFactor * Iinv * Haug{iUser,iBase}' * U{iUser,1} * W{iUser,1};
            end
            
        end
                
        clc;
        [V{1,1} V{1,2} V{2,1} V{2,2}]
        
        if ~iIter
            continueAgain = 1;
        else
            currDeviation = 0;
            for iUser = 1:SimParams.nUsers
                currDeviation = currDeviation + abs(log(det(W{iUser,1})) - log(det(W_prev{iUser,1})));
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
    end
    
    % Assigning the V and U to the corresponding users
    
    linkUsers = 1:SimParams.nUsers;
    for iBase = 1:SimParams.nBases
        
        P = [];aUsers = [];
        sI = (iBase - 1) * SimParams.nTxAntenna + 1;
        eI = sI + SimParams.nTxAntenna - 1;
        SimStructs.baseStruct{iBase,1}.linkedUsers = linkUsers';
        
        for iUser = 1:SimParams.nUsers
            P = [P V{iUser,iBase}];
            aUsers = [aUsers ; repmat(iUser,nStreams,1)];
        end
        
        SimStructs.baseStruct{iBase,1}.P{iBand,1} = P;
        SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1} = aUsers;
        
    end
    
    for iUser = 1:SimParams.nUsers
        SimStructs.userStruct{iUser,1}.W{iBand,1} = U{iUser,1};
        SimStructs.userStruct{iUser,1}.baseNode = 1:SimParams.nBases;
        SimStructs.userStruct{iUser,1}.neighNode = [];
    end
    
end

