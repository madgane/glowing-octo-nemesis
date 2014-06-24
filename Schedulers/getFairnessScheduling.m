function [SimParams,SimStructs] = getFairnessScheduling(SimParams,SimStructs)

charScheduling = char(SimParams.SchedType);
uscore_index = find(charScheduling == '_');

if ~isempty(uscore_index)
    fairnessType = charScheduling(uscore_index(1,1) + 1:end);
    switch fairnessType
        case {'BF','AF'}
            schedulingTypePF = 0;
        case 'PF'
            schedulingTypePF = 1;
        case 'SP'
            schedulingTypePF = 2;
    end
else
    display('Missing PF metric evaluation !');
end

for iBand = 1:SimParams.nBands
    
    for iBase = 1:SimParams.nBases
        
        cBase = SimStructs.baseStruct{iBase,1};
        
        linkedUsers = cBase.linkedUsers;
        kUsers = length(linkedUsers);
        
        iIndex = 0;
        fairnessMetric = zeros(kUsers,1);
        linkChannel  = SimStructs.linkChan{iBase,iBand}(:,:,linkedUsers);
        
        X = zeros(SimParams.nTxAntenna,kUsers * SimParams.maxRank);
        
        for iUser = 1:kUsers
            
            [U,D,~] = svd(linkChannel(:,:,iUser));
            avgServiceRate = SimStructs.userStruct{linkedUsers(iUser,1)}.PFmetric;
            currentWeight = SimStructs.userStruct{linkedUsers(iUser,1)}.weighingFactor;
            avgArrivalRate = SimStructs.userStruct{linkedUsers(iUser,1)}.trafficConfig.avgArrRate;
            
            for iRank = 1:SimParams.maxRank
                iIndex = iIndex + 1;
                currentRate = log2(1 + D(iRank,iRank).^2 * SimStructs.baseStruct{iBase,1}.sPower(1,iBand) / (SimParams.N * SimParams.maxRank));
                
                switch fairnessType
                    case 'BF'
                        fairnessMetric(iIndex,1) = (currentRate * currentWeight) / avgServiceRate;
                    case 'AF'
                        fairnessMetric(iIndex,1) = (currentRate * avgArrivalRate * currentWeight) / avgServiceRate;
                end
                
                if SimParams.queueWt
                    X(:,iIndex) = (U(:,iRank)' * linkChannel(:,:,iUser)).' * (currentWeight);
                else
                    X(:,iIndex) = (U(:,iRank)' * linkChannel(:,:,iUser)).' * sign(currentWeight);
                end
                
            end
            
        end
        
        switch schedulingTypePF
            
            case 0
                
                [~,sortI] = sort(fairnessMetric,'descend');
                currentUsers = floor((sortI(1:SimParams.muxRank,1) - 1) / SimParams.maxRank) + 1;
                currentStreams = mod((sortI(1:SimParams.muxRank,1) - 1),SimParams.maxRank) + 1;
                SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1} = linkedUsers(currentUsers,1);
                SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1} = currentStreams;
                
            case 1
                
                [~,sortI] = sort(fairnessMetric,'descend');
                percentileLimit = floor(0.5 * length(fairnessMetric));
                
                G = [];
                sortI = sortI(1:percentileLimit,1);
                
                ppVolume = zeros(percentileLimit,1);
                schedUsers = zeros(SimParams.muxRank,1);
                schedStreams = zeros(SimParams.muxRank,1);
                
                for iRank = 1:SimParams.muxRank
                    
                    for iUser = 1:percentileLimit
                        U = [G X(:,sortI(iUser,1))];
                        ppVolume(iUser,1) = real(det(U' * U));
                    end
                    
                    [~,sortJ] = sort(ppVolume,'descend');
                    currentIndex = sortI(sortJ(1,1),1);
                    
                    G = [G X(:,currentIndex)];
                    schedUsers(iRank,1) = floor((currentIndex - 1)/SimParams.maxRank) + 1;
                    schedStreams(iRank,1) = mod((currentIndex - 1),SimParams.maxRank) + 1;
                    
                end
                
                SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1} = linkedUsers(schedUsers,1);
                SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1} = schedStreams;
                
            case 2
                
                N = eye(SimParams.nTxAntenna);
                schedUsers = zeros(SimParams.muxRank,1);                
                nullProj = zeros(SimParams.muxRank * kUsers,1);
                for iRank = 1:SimParams.muxRank
                    
                    for iUser = 1:kUsers * SimParams.maxRank
                        currentUser = floor((iUser - 1)/SimParams.maxRank) + 1;
                        weightFactor = SimStructs.userStruct{linkedUsers(currentUser,1)}.weighingFactor;
                        currentRate = log2(1 + norm(N * X(:,iUser)).^2 * SimStructs.baseStruct{iBase,1}.sPower(1,iBand) / (SimParams.N * SimParams.maxRank));
                        nullProj(iUser,1) = currentRate * weightFactor;
                    end
                    
                    [~,sortI] = sort(nullProj,'descend');
                    schedUsers(iRank,1) = sortI(1,1);
                    
                    schedIndices = schedUsers(schedUsers ~= 0,1);
                    N = eye(SimParams.nTxAntenna) - X(:,schedIndices) * pinv(X(:,schedIndices)' * X(:,schedIndices)) * X(:,schedIndices)';
                end
                
                schedStreams = mod((schedUsers - 1),SimParams.maxRank) + 1;
                schedUsers = floor((schedUsers - 1)./SimParams.maxRank) + 1;    
                
                SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1} = linkedUsers(schedUsers,1);
                SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1} = schedStreams;
                
        end
                
    end
    
end

end
