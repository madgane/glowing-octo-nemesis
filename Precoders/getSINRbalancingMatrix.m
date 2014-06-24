function [SimParams,SimStructs] = getSINRbalancingMatrix(SimParams,SimStructs)

SINRbalancing_t.nUsers = SimParams.muxRank;
SINRbalancing_t.H = cell(SimParams.muxRank,1);
SINRbalancing_t.W = cell(SimParams.muxRank,1);
SINRbalancing_t.pFactor = ones(SimParams.muxRank,1);

for iBase = 1:SimParams.nBases
    
    uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
    kUsers = length(uIndices);
    
    pF = zeros(kUsers,1);
    for iUser = 1:kUsers
        pF(iUser,1) = SimStructs.userStruct{uIndices(iUser,1)}.PFmetric;
    end
    
    for iBand = 1:SimParams.nBands
        
        SINRbalancing_t.Pt = SimStructs.baseStruct{iBase,1}.sPower(1,iBand);
        
        W = cell(kUsers,1);
        eG = zeros(SimParams.maxRank,kUsers);
        eH = SimStructs.linkChan{iBase,iBand}(:,:,uIndices);
        for iUser = 1:kUsers
            d = svd(eH(:,:,iUser));
            eG(1:length(d),iUser) = log2(1 + d.^2 * SimStructs.baseStruct{iBase,1}.sPower(1,iBand)) / pF(iUser,1);
            [W{iUser,1}, ~, ~] = svd(eH(:,:,iUser));
        end
        
        [~, sortIndex] = sort(eG(:),'descend');sortIndex = sortIndex - 1;
        
        for iRank = 1:SimParams.maxRank
            if (iRank > SimParams.muxRank)
                break;
            end
            cUser = floor(sortIndex(iRank,1) / SimParams.maxRank) + 1;
            cStream = mod(sortIndex(iRank,1),SimParams.maxRank) + 1;
            
            SINRbalancing_t.H{iRank,1} = eH(:,:,cUser);
            SINRbalancing_t.W{iRank,1} = W{cUser,1}(:,cStream);
        end
        
        [SINRbalancing_t] = performSINRbalancing(SINRbalancing_t);
        SimStructs.baseStruct{iBase}.P(:,:,iBand) = SINRbalancing_t.P;
        SimStructs.baseStruct{iBase}.allocPattern(:,iBand) = sortIndex;
        SimStructs.baseStruct{iBase}.allocGains(:,:,iBand) = eG;
        
        for iRank = 1:SimParams.muxRank
            cUser = floor(sortIndex(iRank,1) / SimParams.maxRank) + 1;
            cStream = mod(sortIndex(iRank,1),SimParams.maxRank) + 1;            
            SimStructs.userStruct{uIndices(cUser,1)}.W(:,cStream,iBand) = SINRbalancing_t.W{iRank,1};
        end
        
    end
    
end

end
