
function [SimParams SimStructs] = getPreWeightedMMSEDesign(SimParams,SimStructs)

iIter = 0;
maxIter = 1e4;
epsilonCheck = min(1e-4,max(SimParams.sPower)^(-2));
nStreams = min(SimParams.maxRank,SimParams.nRxAntenna);

for iBand = 1:SimParams.nBands

    continueAgain = 1;
    W = cell(SimParams.nUsers,1);
    U = cell(SimParams.nUsers,1);
    V = cell(SimParams.nUsers,1);
    linkChannel = SimStructs.linkChan;
    
    combUsers = [];
    for iBase = 1:SimParams.nBases
        combUsers = [combUsers ; unique(SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1})];
    end
    
    userN = cell(length(combUsers),1);
    for iUser = 1:length(combUsers)
        combUIndex = combUsers(iUser,1);
        cUser = SimStructs.userStruct{combUIndex,1};
        [~,~,Z] = svd(linkChannel{cUser.baseNode,iBand}(:,:,combUIndex));
        
        Chn = linkChannel{cUser.baseNode,iBand}(:,:,combUIndex);
        chnInstant = SimParams.elapsedFBDuration(combUIndex) * SimParams.sampTime;
        userN{iUser,1} = abs(1 - besselj(0,(2 * pi * SimParams.userDoppler(combUIndex,1) * chnInstant))) * trace(Chn*Chn') * eye(SimParams.nRxAntenna);

        V{combUIndex,1} = Z(:,1:SimParams.maxRank);
        V{combUIndex,1} = sqrt(max(SimParams.sPower) / (SimParams.nUsers / SimParams.nBases)) * V{combUIndex,1};
    end
    
    while continueAgain
        
        % U Matrix calculation
        
        for iUser = 1:length(combUsers)
            
            combUIndex = combUsers(iUser,1);
            cUser = SimStructs.userStruct{combUIndex,1};
            J = eye(SimParams.nRxAntenna) * SimParams.N + userN{iUser,1} * SimParams.robustNoise;
            
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
        
    end
    
    % Assigning the V and U to the corresponding users
    
    for iBase = 1:SimParams.nBases
        
        cBase = SimStructs.baseStruct{iBase,1};
        assignedUsers = cBase.assignedUsers{iBand,1};
        agUser = zeros(nStreams * length(assignedUsers),1);
        cBase.P{iBand,1} = [];
        
        for iUser = 1:length(assignedUsers)
            
            cUserIndex = assignedUsers(iUser,1);
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
