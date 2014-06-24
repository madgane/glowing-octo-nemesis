function [SimParams,SimStructs] = getSkipScheduling(SimParams,SimStructs)

for iBase = 1:SimParams.nBases

    maxRank = SimParams.maxRank;
    uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
    kUsers = length(uIndices);
    
    uIndices = conv(ones(maxRank,1),upsample(uIndices,maxRank));
    uIndices = uIndices(1:maxRank * kUsers);
            
    for iBand = 1:SimParams.nBands
        SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1} = uIndices;
        SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1} = ones(kUsers * maxRank,1);
    end
end
end

