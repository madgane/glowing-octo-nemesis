function [SimParams,SimStructs] = getTDMZFMatrix(SimParams,SimStructs)

for iBase = 1:SimParams.nBases

    for iBand = 1:SimParams.nBands
    
        augH = [];
        pickUsers = SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1};
        pickStreams = SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1};
    
        for iUser = 1:length(pickUsers)
            cUser = pickUsers(iUser,1);cStream = pickStreams(iUser,1);
            [W,~,~] = svd(SimStructs.linkChan{iBase,iBand}(:,:,cUser));
            augH = [augH ; W(:,cStream)' * SimStructs.linkChan{iBase,iBand}(:,:,cUser)];
            SimStructs.userStruct{cUser,1}.W{iBand,1} = W(:,cStream);
        end

        eP = pinv(augH);
        if mod((SimParams.iDrop - 1),SimParams.nBases) == (iBase - 1)
            effPower = SimStructs.baseStruct{iBase,1}.sPower(1,iBand) * SimParams.nBases;
            [SimStructs.baseStruct{iBase}.P{iBand,1}] = performWFAlgorithm(eP,effPower);
        else
            SimStructs.baseStruct{iBase}.P{iBand,1} = eP * 0;
        end       
        
    end
    
end
   
end
