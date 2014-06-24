function [SimParams,SimStructs] = getBFMatrix(SimParams,SimStructs)

for iBase = 1:SimParams.nBases

    for iBand = 1:SimParams.nBands
    
        augH = [];
        pickUsers = SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1};
        pickStreams = SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1};
    
        for iUser = 1:length(pickUsers)
            cUser = pickUsers(iUser,1);cStream = pickStreams(iUser,1);
            [W,D,V] = svd(SimStructs.linkChan{iBase,iBand}(:,:,cUser));
            
            augH = [augH , V(:,cStream) * D(cStream,cStream)];
            SimStructs.userStruct{cUser,1}.W{iBand,1} = W(:,cStream);
        end

        eP = augH;
        [SimStructs.baseStruct{iBase}.P{iBand,1}] = performWFAlgorithm(eP,SimStructs.baseStruct{iBase,1}.sPower(1,iBand),1./diag(eP' * eP));
%         [SimParams SimStructs] = updateMMSEreceiver(SimParams,SimStructs,1);
        
    end
    
end
   
end
