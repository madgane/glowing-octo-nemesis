function [SimParams,SimStructs] = getNetworkScheduling(SimParams,SimStructs)

linkedUsers = cell(SimParams.nBases,1);
assUsers = cell(SimParams.nBases,SimParams.nBands);
assStreams = cell(SimParams.nBases,SimParams.nBands);

for iUser = 1:SimParams.nUsers
    SimStructs.userStruct{iUser,1}.baseNode = [];
    SimStructs.userStruct{iUser,1}.neighNode = [];
end

for iBand = 1:SimParams.nBands
    
    Haug = [];iIndex = 0;
    uLocs = zeros(SimParams.nUsers * SimParams.maxRank,2);
    
    for iUser = 1:SimParams.nUsers
        H = [];
        for iBase = 1:SimParams.nBases
            H = [H SimStructs.linkChan{iBase,iBand}(:,:,iUser)];
        end
        
        [U,~,~] = svd(H);M = U' * H;
        M = M.' * sign(SimStructs.userStruct{iUser,1}.weighingFactor);
        Haug = [Haug M];
        
        for iStream = 1:SimParams.maxRank
            iIndex = iIndex + 1;
            uLocs(iIndex,:) = [iUser iStream];
        end        
    end
    
    [~,~,sortI] = qr(Haug,'vector');
    
    modLocs = uLocs(sortI,:);
    maxMuxRank = min((SimParams.nTxAntenna * SimParams.nBases),(SimParams.nUsers * SimParams.nRxAntenna));
    
    for iRank = 1:maxMuxRank
        
        cUser = modLocs(iRank,1);cRank = modLocs(iRank,2);
        SimStructs.userStruct{cUser,1}.baseNode = [1:SimParams.nBases];
        
        for iBase = 1:SimParams.nBases
            linkedUsers{iBase,1} = [linkedUsers{iBase,1} ; cUser];
            assUsers{iBase,iBand} = [assUsers{iBase,iBand} ; cUser];
            assStreams{iBase,iBand} = [assStreams{iBase,iBand} ; cRank];            
        end
        
    end 
    
end

for iBase = 1:SimParams.nBases
    for iBand = 1:SimParams.nBands
        SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1} = assUsers{iBase,iBand};
        SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1} = assStreams{iBase,iBand};
    end
    
    SimStructs.baseStruct{iBase,1}.linkedUsers = unique(linkedUsers{iBase,1});
end

