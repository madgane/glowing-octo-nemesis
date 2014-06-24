function [SimParams,SimStructs] = updateZFprecoders(SimParams,SimStructs,performWF)

for iBand = 1:SimParams.nBands
    for iBase = 1:SimParams.nBases
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
        if performWF
            [SimStructs.baseStruct{iBase}.P{iBand,1}] = performWFAlgorithm(eP,SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
        else
            for iUser = 1:length(pickUsers)
                eP(:,1) = eP(:,1) / norm(eP(:,1));
            end
            eP = eP * sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand) / length(pickUsers));
            [SimStructs.baseStruct{iBase}.P{iBand,1}] = eP;
        end
    end
end
