function [SimParams,SimStructs] = getLeakageMatrix(SimParams,SimStructs)

for iBase = 1:SimParams.nBases

    for iBand = 1:SimParams.nBands
    
        augH = [];
        Q = zeros(SimParams.muxRank,1);
        pickUsers = SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1};
        pickStreams = SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1};
        kUsers = length(pickUsers);
        
        for iUser = 1:kUsers
            cUser = pickUsers(iUser,1);cStream = pickStreams(iUser,1);
            [W,~,~] = svd(SimStructs.linkChan{iBase,iBand}(:,:,cUser));
            augH = [augH ; W(:,cStream)' * SimStructs.linkChan{iBase,iBand}(:,:,cUser)];
            SimStructs.userStruct{cUser,1}.W{iBand,1}(:,cStream) = W(:,cStream);
            Q(iUser,1) = SimStructs.userStruct{cUser,1}.trafficStats.backLogPkt;
        end

        eigLoading = zeros(kUsers,1);
        eP = zeros(SimParams.nTxAntenna,kUsers);
        
        for iUser = 1:kUsers
            uChan = augH(iUser,:);
            ifChan = augH((iUser ~= 1:kUsers),:);
            
            matP = pinv(SimParams.N * eye(SimParams.nTxAntenna) + (ifChan' * ifChan)) * (uChan' * uChan);
            [V,D] = eig(matP);[m,i] = max(abs(diag(D)));
            eP(:,iUser) = V(:,i);
            eigLoading(iUser,1) = m;
        end
        
        if SimParams.queueWt == 2
            [SimStructs.baseStruct{iBase}.P{iBand,1}] = performQueuedWF(eP,SimStructs.baseStruct{iBase,1}.sPower(1,iBand),Q);
        else
            [SimStructs.baseStruct{iBase}.P{iBand,1}] = performQueuedWF(eP,SimStructs.baseStruct{iBase,1}.sPower(1,iBand),eigLoading);
        end
        
    end
    
end
   
end
