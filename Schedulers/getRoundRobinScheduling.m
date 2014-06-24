function [SimParams,SimStructs] = getRoundRobinScheduling(SimParams,SimStructs)

for iBase = 1:SimParams.nBases

    uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
    kUsers = length(uIndices);
        
    for iBand = 1:SimParams.nBands
        
        sortIndex = ((SimParams.iDrop - 1) * SimParams.muxRank) + (0:SimParams.muxRank - 1);
        sortIndex = unique(mod(sortIndex,kUsers) + 1);
        SimStructs.baseStruct{iBase}.assignedUsers{iBand,1} = uIndices(sortIndex);
        SimStructs.baseStruct{iBase}.assignedStreams{iBand,1} = ones(SimParams.muxRank,1);

    end
end
end

