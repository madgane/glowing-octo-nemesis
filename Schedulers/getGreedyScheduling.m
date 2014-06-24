function [SimParams,SimStructs] = getGreedyScheduling(SimParams,SimStructs)

for iBase = 1:SimParams.nBases
    
    uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
    kUsers = length(uIndices);
    
    schedUsers = zeros(SimParams.nTxAntenna,1);
    schedStreams = zeros(SimParams.muxRank,1);
    
    for iBand = 1:SimParams.nBands
        
        iIndex = 0;augE = [];
        xLocs = zeros(kUsers * SimParams.maxRank,2);
        eH = SimStructs.linkChan{iBase,iBand}(:,:,uIndices);
        
        for iUser = 1:kUsers
            cUser = uIndices(iUser,1);
            
            [U,~,~] = svd(eH(:,:,iUser));
            if SimParams.queueWt
                M = U' * eH(:,:,iUser) * (SimStructs.userStruct{cUser,1}.weighingFactor);
            else
                M = U' * eH(:,:,iUser) * sign(SimStructs.userStruct{cUser,1}.weighingFactor);
            end
            for iRank = 1:SimParams.maxRank
                iIndex = iIndex + 1;
                augE = [augE norm(M(iRank,:))];
                xLocs(iIndex,:) = [iUser iRank];
            end
        end
        
        [~,sortA] = sort(augE,'descend');
        for iRank = 1:SimParams.muxRank
            schedUsers(iRank,1) = xLocs(sortA(1,iRank),1);
            schedStreams(iRank,1) = xLocs(sortA(1,iRank),2);
        end
        
        SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1} = uIndices(schedUsers);
        SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1} = schedStreams;
        
    end
    
end