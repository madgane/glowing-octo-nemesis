
function [SimParams SimStructs] = getStrWeightedMMSEDesign(SimParams,SimStructs)

iIter = 0;
maxIter = 1e4;
epsilonCheck = min(1e-4,max(SimParams.sPower)^(-2));
nStreams = min(SimParams.maxRank,SimParams.nRxAntenna);

SumCapacity = cell(SimParams.nBands,1);
if isfield(SimParams.Debug.privateExchanges,'takeOverBand')
    bandRange = SimParams.Debug.privateExchanges.takeOverBand;
else
    bandRange = 1:SimParams.nBands;
    SimParams.Debug.tempResource{2,1} = cell(SimParams.nUsers,1);
end            
        
for iBand = bandRange

    continueAgain = 1;
    W = cell(SimParams.nUsers,1);
    U = cell(SimParams.nUsers,1);
    V = cell(SimParams.nUsers,1);
    linkChannel = SimStructs.linkChan;
    
    combUsers = [];
    gCombUsers = [];
    gCombStreams = [];
    for iBase = 1:SimParams.nBases
        combUsers = [combUsers ; unique(SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1})];
        gCombUsers = [gCombUsers ; SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1}];
        gCombStreams = [gCombStreams ; SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1}];
    end
    
    for iUser = 1:length(combUsers)
        
        combUIndex = combUsers(iUser,1);
        cUser = SimStructs.userStruct{combUIndex,1};
        [~,~,Z] = svd(linkChannel{cUser.baseNode,iBand}(:,:,combUIndex));
        
        uLocs = find(gCombUsers == combUIndex);
        uStreams = gCombStreams(uLocs,1);
        
        V{combUIndex,1} = zeros(SimParams.nTxAntenna,nStreams);
        for iStream = 1:length(uStreams)
            V{combUIndex,1}(:,uStreams(iStream,1)) = complex(ones(SimParams.nTxAntenna,1),ones(SimParams.nTxAntenna,1));
        end
        
        V{combUIndex,1} = sqrt(max(SimParams.sPower) / (SimParams.nUsers / SimParams.nBases)) * V{combUIndex,1} / (SimParams.nTxAntenna * SimParams.nRxAntenna * 2);
        
    end
    
    while continueAgain
        
        % U Matrix calculation
        
        for iUser = 1:length(combUsers)
            
            combUIndex = combUsers(iUser,1);
            cUser = SimStructs.userStruct{combUIndex,1};
            J = eye(SimParams.nRxAntenna) * SimParams.N;
            
            for jUser = 1:length(combUsers)
                combIFUIndex = combUsers(jUser,1);
                ifUser = SimStructs.userStruct{combIFUIndex,1};
                HV = linkChannel{ifUser.baseNode,iBand}(:,:,combUIndex) * V{combIFUIndex,1};
                J = J + HV * HV';
            end
            
            H = linkChannel{cUser.baseNode,iBand}(:,:,combUIndex);
            
            U{combUIndex,1} = J \ (H * V{combUIndex,1});
            W{combUIndex,1} = inv(eye(nStreams) - U{combUIndex,1}' * H * V{combUIndex,1});
            
        end
        
        for iBase = 1:SimParams.nBases
            
            cBase = SimStructs.baseStruct{iBase,1};
            linkedUsers = unique(cBase.assignedUsers{iBand,1});
            
            Isum = 0;Dsum = 0;
            for iUser = 1:length(combUsers)
                combUIndex = combUsers(iUser,1);
                cUser = SimStructs.userStruct{combUIndex,1};
                H_HU = linkChannel{iBase,iBand}(:,:,combUIndex)' * U{combUIndex,1};
                Isum = Isum + cUser.weighingFactor * H_HU * W{combUIndex,1} * H_HU';
                if cUser.baseNode == iBase
                    W_2 = W{combUIndex,1} * W{combUIndex,1};
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
            for iUser = 1:length(combUsers)
                combUIndex = combUsers(iUser,1);
                currDeviation = currDeviation + abs(log(det(W{combUIndex,1})) - log(det(W_prev{combUIndex,1})));
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
        [sumCap,qDev,rVec] = performMockReception(SimParams,SimStructs,V,U,W,iBand);
        
        SumCapacity{iBand,1} = [SumCapacity{iBand,1} ; sumCap];
        for iUser = 1:SimParams.nUsers
            if ~isfield(SimParams.Debug.privateExchanges,'takeOverBand')
                SimParams.Debug.tempResource{2,SimParams.iDrop}{iUser,1} = [SimParams.Debug.tempResource{2,SimParams.iDrop}{iUser,1} rVec(iUser,1)];
            end
            SimParams.Debug.tempResource{3,SimParams.iDrop}{iUser,1} = [SimParams.Debug.tempResource{3,SimParams.iDrop}{iUser,1} qDev(iUser,1)];
            SimParams.Debug.tempResource{4,SimParams.iDrop}{iUser,iBand} = [SimParams.Debug.tempResource{4,SimParams.iDrop}{iUser,iBand} (rVec(iUser,1))];
        end
        
    end
    
    % Assigning the V and U to the corresponding users
    
    for iBase = 1:SimParams.nBases
        
        cBase = SimStructs.baseStruct{iBase,1};
        assignedUsers = cBase.assignedUsers{iBand,1};
        
        cBase.P{iBand,1} = [];
        uAssignedUsers = unique(assignedUsers);
        agUser = zeros(nStreams * length(uAssignedUsers),1);
        
        for iUser = 1:length(uAssignedUsers)
            
            cUserIndex = uAssignedUsers(iUser,1);
            cUser = SimStructs.userStruct{cUserIndex,1};
            cBase.P{iBand,1} = [cBase.P{iBand,1} , V{cUserIndex,1}];
            
            sI = (iUser - 1) * nStreams + 1;
            eI = sI + nStreams - 1;
            agUser(sI:eI,1) = cUserIndex;
            
            cUser.W{iBand,1} = U{cUserIndex,1};
            SimStructs.userStruct{cUserIndex,1} = cUser;
        end
        
        cBase.assignedUsers{iBand,1} = agUser;
        SimStructs.baseStruct{iBase,1} = cBase;
        
    end
    
end
