function [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,rxType,bsIndices)

preLogue;
if nargin ~= 4
    bsIndices = 1:nBases;
end

pertNoise = 1e-20;
W = cell(nUsers,nBands);

for iBand = 1:nBands
    
    switch rxType
        
        case 'Ones'
            
            for iBase = bsIndices
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    W{cUser,iBand} = ones(SimParams.nRxAntenna,SimParams.maxRank);
                end
            end
            
        case 'BF'
            
            for iBase = bsIndices
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    [U,~,~] = svd(cH{iBase,iBand}(:,:,cUser));
                    W{cUser,iBand} = U;
                end
            end
            
        case 'MMSE'
            
            for iBase = bsIndices
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    H = cH{iBase,iBand}(:,:,cUser);
                    for iLayer = 1:SimParams.maxRank
                        R = eye(SimParams.nRxAntenna);
                        for jBase = bsIndices
                            for jUser = 1:usersPerCell(jBase,1)
                                R = R + H * SimStructs.baseStruct{iBase,1}.P{iBand,1}(:,:,jUser) * SimStructs.baseStruct{iBase,1}.P{iBand,1}(:,:,jUser)' * H';
                            end
                        end
                        W{cUser,iBand} = R \ (H * SimStructs.baseStruct{iBase,1}.P{iBand,1}(:,iLayer,iUser)) + pertNoise;
                    end
                end
            end
    end
    
end

for iBand = 1:nBands
    for iBase = bsIndices
        for iUser = 1:usersPerCell(iBase,1);
            cUser = cellUserIndices{iBase,1}(iUser,1);
            SimStructs.userStruct{cUser,1}.W{iBand,1} = W{cUser,iBand};
            SimStructs.userStruct{cUser,1}.pW{iBand,1} = W{cUser,iBand};
        end
    end
end
