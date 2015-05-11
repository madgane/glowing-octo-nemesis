function [SimParams,SimStructs] = getZFMatrix(SimParams,SimStructs)

debugCode = 1;

for iBase = 1:SimParams.nBases

    for iBand = 1:SimParams.nBands
    
        augH = [];
        Q = zeros(SimParams.muxRank,1);
        pickUsers = SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1};
        pickStreams = SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1};
        for iUser = 1:length(pickUsers)
            cUser = pickUsers(iUser,1);cStream = pickStreams(iUser,1);
            [W,~,~] = svd(SimStructs.linkChan{iBase,iBand}(:,:,cUser));
            augH = [augH ; W(:,cStream)' * SimStructs.linkChan{iBase,iBand}(:,:,cUser)];
            SimStructs.userStruct{cUser,1}.W{iBand,1}(:,cStream) = W(:,cStream);
            Q(iUser,1) = SimStructs.userStruct{cUser,1}.trafficStats.backLogPkt;
        end

        eP = (pinv(augH' * augH)) * augH';
        switch SimParams.queueWt
            case 1
                X = performWFAlgorithm(eP,SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
            case 2
                X = performQueuedWF(eP,SimStructs.baseStruct{iBase,1}.sPower(1,iBand),Q);
            otherwise
                X = performWFAlgorithm(eP,SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
        end
        
        if (debugCode == 1)
            P = zeros(SimParams.nTxAntenna,SimParams.maxRank,length(SimStructs.baseStruct{iBase,1}.linkedUsers));
            for iUser = 1:length(pickUsers)
                P(:,pickStreams(iUser,1),pickUsers(iUser,1)) = X(:,iUser);
            end
            SimStructs.baseStruct{iBase,1}.P{iBand,1} = P;
        else
            SimStructs.baseStruct{iBase,1}.P{iBand,1} = X;
        end
        
        
    end
    
end
 
if (debugCode == 1)
    [SimParams, SimStructs] = getSkipScheduling(SimParams,SimStructs);
    [SimParams, SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE');
end

end
