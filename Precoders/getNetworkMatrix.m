function [SimParams,SimStructs] = getNetworkMatrix(SimParams,SimStructs)

for iBand = 1:SimParams.nBands

    activeUsers = SimStructs.baseStruct{1,1}.assignedUsers{iBand,1};
    activeRanks = SimStructs.baseStruct{1,1}.assignedStreams{iBand,1};
    
    Haug = [];
    for iUser = 1:length(activeUsers)
        
        H = [];
        cUser = activeUsers(iUser,1);cRank = activeRanks(iUser,1);
        for iBase = 1:SimParams.nBases
           H = [H SimStructs.linkChan{iBase,iBand}(:,:,cUser)];
        end
        
        [U,~,~] = svd(H);M = U' * H;
        M = M * sign(SimStructs.userStruct{cUser,1}.weighingFactor);
        Haug = [Haug ; M(cRank,:)];
        
    end
    
    eP = pinv(Haug);
    eP = performWFAlgorithm(eP,SimParams.sPower(1,iBand) * SimParams.nBases);
    for iBase = 1:SimParams.nBases
        sI = (iBase - 1) * SimParams.nTxAntenna + 1;
        eI = sI + SimParams.nTxAntenna - 1;eK = eP(sI:eI,:);
        [SimStructs.baseStruct{iBase,1}.P{iBand,1}] = eK;
    end
    
end
    
    
    