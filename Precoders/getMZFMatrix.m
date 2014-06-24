function [SimParams,SimStructs] = getMZFMatrix(SimParams,SimStructs)

epsilon = 0.20;

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
        [X,Pwr] = performWFAlgorithm(eP,SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
        [sIndexP] = find((Pwr / max(Pwr)) >= epsilon);
        
        augX = zeros(size(augH));
        for ivUser = 1:length(sIndexP)
            ivcUser = sIndexP(ivUser);
            augX(ivUser,:) = augH(ivcUser,:);
        end
        
        eP = pinv(augX);
        [SimStructs.baseStruct{iBase}.P{iBand,1}] = performWFAlgorithm(eP,SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
        SimStructs.baseStruct{iBase}.assignedUsers{iBand,1} = pickUsers(sIndexP);
        SimStructs.baseStruct{iBase}.assignedStreams{iBand,1} = pickStreams(sIndexP);
        
    end
    
end
   
end


