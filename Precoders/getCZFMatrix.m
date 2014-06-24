function [SimParams,SimStructs] = getCZFMatrix(SimParams,SimStructs)

for iBand = 1:SimParams.nBands

    assignedUsers = cell(SimParams.nBases,1);
    assignedStreams = cell(SimParams.nBases,1);
    
    for iBase = 1:SimParams.nBases
        assignedUsers{iBase,1} = SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1};
        assignedStreams{iBase,1} = SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1};
    end
    
    Haug = cell(SimParams.nBases,1);
    Umatrix = cell(SimParams.nBases,SimParams.maxRank);
    
    for iBase = 1:SimParams.nBases
        for iUser = 1:length(assignedUsers{iBase,1})
            cUser = SimStructs.userStruct{assignedUsers{iBase,1}(iUser,1),1};
            H = SimStructs.linkChan{iBase,iBand}(:,:,assignedUsers{iBase,1}(iUser,1));
            
            [U,~,~] = svd(H);
            Umatrix{iBase,iUser} = U;
            M = U' * H * sign(cUser.weighingFactor);
            Haug{iBase,1} = [Haug{iBase,1} ; M(assignedStreams{iBase,1}(iUser,1),:)];
            SimStructs.userStruct{assignedUsers{iBase,1}(iUser,1),1}.W{iBand,1} = Umatrix{iBase,iUser}(:,assignedStreams{iBase,1}(iUser,1));
        end
    end
    
    for iBase = 1:SimParams.nBases
        for jBase = 1:SimParams.nBases
            if jBase ~= iBase
                for iUser = 1:length(assignedUsers{jBase,1})
                    cUser = SimStructs.userStruct{assignedUsers{jBase,1}(iUser,1),1};
                    H = SimStructs.linkChan{iBase,iBand}(:,:,assignedUsers{jBase,1}(iUser,1));
                    M = Umatrix{jBase,iUser}' * H * sign(cUser.weighingFactor);
                    Haug{iBase,1} = [Haug{iBase,1} ; M(assignedStreams{jBase,1}(iUser,1),:)];
                end
            end
        end
        
        eP = pinv(Haug{iBase,1});eP = eP(:,1:length(assignedUsers{iBase,1}));
        [SimStructs.baseStruct{iBase,1}.P{iBand,1}] = performWFAlgorithm(eP,SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
        
    end       
        
end

