function [SimParams,SimStructs] = getRouletteWheelScheduling(SimParams,SimStructs)

for iBase = 1:SimParams.nBases

    uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
    kUsers = length(uIndices);
        
    for iBand = 1:SimParams.nBands
        
        sortIndex = unique(randi(kUsers,SimParams.muxRank,1));
        SimStructs.baseStruct{iBase}.assignedUsers{iBand,1} = uIndices(sortIndex);
        SimStructs.baseStruct{iBase}.assignedStreams{iBand,1} = ones(SimParams.muxRank,1);

    end
end
end
