
function [SimParams,SimStructs] = getDWeightedMMSEDesign(SimParams,SimStructs)

iIter = 0;
maxIter = 1e4;
epsilonCheck = min(1e-4,max(SimParams.sPower)^(-2));
nStreams = min(SimParams.maxRank,SimParams.nRxAntenna);

iterationPlot = 0;
SumCapacity = cell(SimParams.nBands,1);

for iBand = 1:SimParams.nBands

    continueAgain = 1;
    W = cell(SimParams.nUsers,1);
    U = cell(SimParams.nUsers,1);
    V = cell(SimParams.nUsers,1);
    linkChannel = SimStructs.linkChan;
    
    Ui = cell(SimParams.nUsers,1);
    
    for iUser = 1:SimParams.nUsers
        V{iUser,1} = complex(ones(SimParams.nTxAntenna,nStreams),ones(SimParams.nTxAntenna,nStreams));
        V{iUser,1} = sqrt(max(SimParams.sPower) / (SimParams.nUsers / SimParams.nBases)) * V{iUser,1} / trace(V{iUser,1}' * V{iUser,1});
    end
    
    while continueAgain
        
        for iBase = 1:SimParams.nBases
            
            cBase = SimStructs.baseStruct{iBase,1};
            linkedUsers = cBase.linkedUsers;kUsers = length(linkedUsers);
            
            % U Matrix calculation
            
            for iUser = 1:kUsers
                
                cUser = linkedUsers(iUser,1);
                J = eye(SimParams.nRxAntenna) * SimParams.N;
                
                for jUser = 1:kUsers
                    cjUser = linkedUsers(jUser,1);
                    HV = linkChannel{iBase,iBand}(:,:,cUser) * V{cjUser,1};
                    J = J + HV * HV';
                end
                
                H = linkChannel{iBase,iBand}(:,:,cUser);
                
                HdVd = H * V{cUser,1};
                U{cUser,1} = inv(J) * HdVd;
                W{cUser,1} = inv(eye(nStreams) - U{cUser,1}' * H * V{cUser,1});
                Ui{cUser,1} = J - HdVd * HdVd' - eye(SimParams.nRxAntenna) * SimParams.N;
                
            end            
            
            Isum = 0;Dsum = 0;
            for iUser = 1:kUsers
                userIndex = linkedUsers(iUser,1);
                cUser = SimStructs.userStruct{userIndex,1};
                H_HU = linkChannel{iBase,iBand}(:,:,userIndex)' * U{userIndex,1};
                Isum = Isum + cUser.weighingFactor * H_HU * W{userIndex,1} * H_HU';
                if cUser.baseNode == iBase
                    W_2 = W{userIndex,1} * W{userIndex,1};
                    Dsum = Dsum + cUser.weighingFactor^2 * H_HU * W_2 * H_HU';
                end
            end
            
            mu_star = bisectionEstimateMU(Isum,Dsum,SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
            Isum = Isum + mu_star * eye(SimParams.nTxAntenna);
            
            Iinv = pinv(Isum);
            for iUser = 1:length(linkedUsers)
                cIndex = linkedUsers(iUser,1);
                cUser = SimStructs.userStruct{cIndex,1};
                V{cIndex,1} = cUser.weighingFactor * Iinv * linkChannel{iBase,iBand}(:,:,cIndex)' * U{cIndex,1} * W{cIndex,1};
            end
            
        end
        
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
        SumCapacity{iBand,1} = [SumCapacity{iBand,1} ; performMockReception(SimParams,SimStructs,V,iBand)];
        
    end
    
    % Assigning the V and U to the corresponding users
    
    for iBase = 1:SimParams.nBases
        cBase = SimStructs.baseStruct{iBase,1};
        
        cBase.P{iBand,1} = [];
        cBase.assignedUsers{iBand,1} = [];
        cBase.assignedStreams{iBand,1} = [];
        
        for iUser = 1:length(cBase.linkedUsers)
            cUserIndex = cBase.linkedUsers(iUser,1);
            cUser = SimStructs.userStruct{cUserIndex,1};
            cBase.P{iBand,1} = [cBase.P{iBand,1} , V{cUserIndex,1}];
            cUser.W{iBand,1} = U{cUserIndex,1};
            
            xStreams = (1:nStreams)';
            cBase.assignedUsers{iBand,1} = [cBase.assignedUsers{iBand,1} ; repmat(cUserIndex,length(xStreams),1)];
            cBase.assignedStreams{iBand,1} = [cBase.assignedStreams{iBand,1} ; xStreams];
            SimStructs.userStruct{cUserIndex,1} = cUser;
        end
        
        SimStructs.baseStruct{iBase,1} = cBase;
    end
    
end

if iterationPlot
    figure(10);hold all;
    plot(SumCapacity{1,1});
end
