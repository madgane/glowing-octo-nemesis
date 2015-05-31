function [SimParams,SimStructs] = getSkipScheduling(SimParams,SimStructs)

for iBase = 1:SimParams.nBases

    maxRank = SimParams.maxRank;
    uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
    uIndices = repmat(uIndices',maxRank,1);uIndices = uIndices(:);
    
    uStreams = 1:maxRank;
    uStreams = repmat(uStreams',(length(uIndices) / maxRank),1);
    
    for iBand = 1:SimParams.nBands
        SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1} = uIndices;
        SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1} = uStreams;
    end
end
end

