function [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,rxType,bsIndices)

cH = SimStructs.linkChan;
nBases = SimParams.nBases;
nBands = SimParams.nBands;
nUsers = SimParams.nUsers;
maxRank = SimParams.maxRank;
cellUserIndices = cell(nBases,1);
usersPerCell = zeros(nBases,1);

if nargin ~= 4
    bsIndices = 1:nBases;
end

for iBase = 1:nBases
    for iBand = 1:nBands
        cellUserIndices{iBase,1} = [cellUserIndices{iBase,1} ; SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1}];
    end
    cellUserIndices{iBase,1} = unique(cellUserIndices{iBase,1});
    usersPerCell(iBase,1) = length(cellUserIndices{iBase,1});
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
                    for iLayer = 1:maxRank
                        R = eye(SimParams.nRxAntenna);
                        for jBase = bsIndices
                            for jUser = 1:usersPerCell(jBase,1)
                                R = R + H * SimStructs.baseStruct{iBase,1}.P{iBand,1}(:,:,jUser) * SimStructs.baseStruct{iBase,1}.P{iBand,1}(:,:,jUser)' * H';
                            end
                        end
                        W{cUser,iBand}(:,iLayer) = R \ (H * SimStructs.baseStruct{iBase,1}.P{iBand,1}(:,iLayer,iUser)) + pertNoise;
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
