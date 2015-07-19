function [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,rxType,bsIndices)

stringCells = strsplit(rxType,'_');
if length(stringCells) == 1
    cH = SimStructs.linkChan;
else
    cH = SimStructs.prevChan;
end

proLogue;
rxType = stringCells{1};
updateHistory = 'true';

if nargin ~= 4
    bsIndices = 1:nBases;
end

pertNoise = 1e-40;
W = cell(nUsers,nBands);
M0 = cell(SimParams.nBases,1);

switch rxType
    
    case 'Ones'
        for iBand = 1:nBands
            for iBase = bsIndices
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    W{cUser,iBand} = ones(SimParams.nRxAntenna,SimParams.maxRank);
                end
            end
        end
        
    case 'Random'
        for iBand = 1:nBands
            for iBase = bsIndices
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    W{cUser,iBand} = complex(randn(SimParams.nRxAntenna,SimParams.maxRank),randn(SimParams.nRxAntenna,SimParams.maxRank));
                end
            end
        end

        
    case 'BF'        
        for iBand = 1:nBands
            for iBase = bsIndices
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    [U,~,~] = svd(cH{iBase,iBand}(:,:,cUser));
                    W{cUser,iBand} = U;
                end
            end
        end
        
    case 'MMSE'        
        for iBand = 1:nBands
            for iBase = bsIndices
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    for iLayer = 1:maxRank
                        if strcmpi(SimParams.sysMode,'true')
                            R = (SimParams.N + SimStructs.userStruct{cUser,1}.phyParams.restOfIF_Linear) * eye(SimParams.nRxAntenna);
                        else
                            R = SimParams.N * eye(SimParams.nRxAntenna);
                        end
                        for jBase = bsIndices
                            H = cH{jBase,iBand}(:,:,cUser);
                            for jUser = 1:usersPerCell(jBase,1)
                                R = R + H * SimStructs.baseStruct{jBase,1}.P{iBand,1}(:,:,jUser) * SimStructs.baseStruct{jBase,1}.P{iBand,1}(:,:,jUser)' * H';
                            end
                        end
                        H = cH{iBase,iBand}(:,:,cUser);
                        W{cUser,iBand}(:,iLayer) = R \ (H * SimStructs.baseStruct{iBase,1}.P{iBand,1}(:,iLayer,iUser)) + pertNoise;
                    end
                end
            end
        end
        
    case 'MMSE-X'
        
        cM = SimParams.Debug.ExchangeM{bsIndices,1};
                
        for iBand = 1:nBands
            for iBase = bsIndices
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    for iLayer = 1:maxRank
                        if strcmpi(SimParams.sysMode,'true')
                            R = (SimParams.N + SimStructs.userStruct{cUser,1}.phyParams.restOfIF_Linear) * eye(SimParams.nRxAntenna);
                        else
                            R = SimParams.N * eye(SimParams.nRxAntenna);
                        end
                        for jBase = 1:nBases
                            H = cH{jBase,iBand}(:,:,cUser);
                            if jBase ~= iBase
                                for jUser = 1:usersPerCell(jBase,1)
                                    R = R + H * SimParams.Debug.globalExchangeInfo.funcOut{1,jBase}(:,:,jUser,iBand) * SimParams.Debug.globalExchangeInfo.funcOut{1,jBase}(:,:,jUser,iBand)' * H';
                                end
                            else
                                for jUser = 1:usersPerCell(jBase,1)
                                    R = R + H * cM(:,:,jUser,iBand) * cM(:,:,jUser,iBand)' * H';
                                end
                            end
                        end
                        H = cH{iBase,iBand}(:,:,cUser);
                        SimParams.Debug.ExchangeW{cUser,iBand}(:,iLayer) = R \ (H * cM(:,iLayer,iUser,iBand)) + pertNoise;
                    end
                end
            end
        end

        
    case 'MMSE-BF'        

        for iBand = 1:nBands
            for iBase = bsIndices
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    [~,~,V] = svd(cH{iBase,iBand}(:,:,cUser));
                    M0{iBase,1}(:,:,iUser,iBand) = V(:,1:SimParams.maxRank);
                end
            end
        end
        
        for iBase = bsIndices
            totPower = norm(vec(M0{iBase,1}))^2;
            totPower = sqrt(sum(SimStructs.baseStruct{iBase,1}.sPower) / totPower);
            M0{iBase,1} = M0{iBase,1} * totPower;
        end

        for iBand = 1:nBands
            for iBase = bsIndices
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    for iLayer = 1:maxRank
                        if strcmpi(SimParams.sysMode,'true')
                            R = (SimParams.N + SimStructs.userStruct{cUser,1}.phyParams.restOfIF_Linear) * eye(SimParams.nRxAntenna);
                        else
                            R = SimParams.N * eye(SimParams.nRxAntenna);
                        end
                        for jBase = bsIndices
                            H = cH{jBase,iBand}(:,:,cUser);
                            for jUser = 1:usersPerCell(jBase,1)
                                R = R + H * M0{jBase,1}(:,:,jUser,iBand) * M0{jBase,1}(:,:,jUser,iBand)' * H';
                            end
                        end
                        H = cH{iBase,iBand}(:,:,cUser);
                        W{cUser,iBand}(:,iLayer) = R \ (H * M0{iBase,1}(:,iLayer,iUser,iBand)) + pertNoise;
                    end
                end
            end
        end
        
    case 'MMSE-Ones'        
        for iBase = bsIndices
            M0{iBase,1} = complex(ones(SimParams.nTxAntenna,maxRank,usersPerCell(iBase,1),nBands),ones(SimParams.nTxAntenna,maxRank,usersPerCell(iBase,1),nBands));
        end
        
        for iBase = bsIndices
            totPower = norm(vec(M0{iBase,1}))^2;
            totPower = sqrt(sum(SimStructs.baseStruct{iBase,1}.sPower) / totPower);
            M0{iBase,1} = M0{iBase,1} * totPower;
        end

        for iBand = 1:nBands
            for iBase = bsIndices
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    for iLayer = 1:maxRank
                        if strcmpi(SimParams.sysMode,'true')
                            R = (SimParams.N + SimStructs.userStruct{cUser,1}.phyParams.restOfIF_Linear) * eye(SimParams.nRxAntenna);
                        else
                            R = SimParams.N * eye(SimParams.nRxAntenna);
                        end
                        for jBase = bsIndices
                            H = cH{jBase,iBand}(:,:,cUser);
                            for jUser = 1:usersPerCell(jBase,1)
                                R = R + H * M0{jBase,1}(:,:,jUser,iBand) * M0{jBase,1}(:,:,jUser,iBand)' * H';
                            end
                        end
                        H = cH{iBase,iBand}(:,:,cUser);
                        W{cUser,iBand}(:,iLayer) = R \ (H * M0{iBase,1}(:,iLayer,iUser,iBand)) + pertNoise;
                    end
                end
            end
        end
        
    case 'Last'        
        for iBand = 1:nBands
            for iBase = bsIndices
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    W{cUser,iBand} = SimStructs.userStruct{cUser,1}.pW{iBand,1};
                end
            end
        end
        
    case 'FrameC'
        for iBand = 1:nBands
            for iBase = bsIndices
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    W{cUser,iBand} = SimStructs.userStruct{cUser,1}.pW{iBand,1};
                end
            end
        end
        
        for iBand = 1:nBands
            for iBase = bsIndices
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    [~,~,V] = svd(cH{iBase,iBand}(:,:,cUser));
                    M0{iBase,1}(:,:,iUser,iBand) = V(:,1:SimParams.maxRank);
                end
            end
        end
        
        for iBase = bsIndices
            totPower = norm(vec(M0{iBase,1}))^2;
            totPower = sqrt(sum(SimStructs.baseStruct{iBase,1}.sPower) / totPower);
            M0{iBase,1} = M0{iBase,1} * totPower;
        end

        for iBand = 1:nBands
            for iBase = bsIndices
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    for iLayer = 1:maxRank
                        if strcmpi(SimParams.sysMode,'true')
                            R = (SimParams.N + SimStructs.userStruct{cUser,1}.phyParams.restOfIF_Linear) * eye(SimParams.nRxAntenna);
                        else
                            R = SimParams.N * eye(SimParams.nRxAntenna);
                        end
                        for jBase = bsIndices
                            H = cH{jBase,iBand}(:,:,cUser);
                            for jUser = 1:usersPerCell(jBase,1)
                                R = R + H * M0{jBase,1}(:,:,jUser,iBand) * M0{jBase,1}(:,:,jUser,iBand)' * H';
                            end
                        end
                        H = cH{iBase,iBand}(:,:,cUser);
                        if (norm(W{cUser,iBand}(:,iLayer)) < 1e-10)
                            W{cUser,iBand}(:,iLayer) = R \ (H * M0{iBase,1}(:,iLayer,iUser,iBand)) + pertNoise;
                        end
                    end
                end
            end
        end


        
    case 'MMSE-XVAR'
        
        updateHistory = 'false';  
        for iBand = 1:nBands
            for iUser = 1:SimParams.nUsers
                SimParams.Debug.globalExchangeInfo.funcOut{6,bsIndices}{iUser,iBand} = SimStructs.userStruct{iUser,1}.pW{iBand,1};
            end
        end
        
        for iBand = 1:nBands
            for iBase = bsIndices
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    for iLayer = 1:maxRank
                        if strcmpi(SimParams.sysMode,'true')
                            R = (SimParams.N + SimStructs.userStruct{cUser,1}.phyParams.restOfIF_Linear) * eye(SimParams.nRxAntenna);
                        else
                            R = SimParams.N * eye(SimParams.nRxAntenna);
                        end
                        H = cH{iBase,iBand}(:,:,cUser);
                        
                        for jUser = 1:usersPerCell(iBase,1)
                            if jUser ~= iUser
                                R = R + H * SimStructs.baseStruct{iBase,1}.P{iBand,1}(:,:,jUser) * SimStructs.baseStruct{iBase,1}.P{iBand,1}(:,:,jUser)' * H';
                            else
                                R = R + H * SimStructs.baseStruct{iBase,1}.P{iBand,1}(:,(iLayer ~= 1:maxRank),jUser) * SimStructs.baseStruct{iBase,1}.P{iBand,1}(:,(iLayer ~= 1:maxRank),jUser)' * H';
                            end
                        end
                        
                        for jBase = 1:nBases
                            if jBase ~= iBase
                                ifLevel = (SimParams.Debug.globalExchangeInfo.gI{jBase,1}(iLayer,cUser,iBand)^2) / norm(SimStructs.userStruct{cUser,1}.pW{iBand,1}(:,iLayer))^2;
                                R = R + ifLevel * eye(SimParams.nRxAntenna) / sqrt(SimParams.nRxAntenna);
                            end
                        end                        
                        H = cH{iBase,iBand}(:,:,cUser);
                        W{cUser,iBand}(:,iLayer) = R \ (H * SimStructs.baseStruct{iBase,1}.P{iBand,1}(:,iLayer,iUser)) + pertNoise;
                    end
                end
            end
        end

        for iBand = 1:nBands
            for iUser = 1:usersPerCell(bsIndices,1)
                cUser = cellUserIndices{bsIndices,1}(iUser,1);
                SimParams.Debug.globalExchangeInfo.funcOut{6,bsIndices}{cUser,iBand} = W{cUser,iBand};
            end
        end

end

if and((~strcmpi(SimParams.SchedType,'SkipScheduling')),~(strcmpi(rxType,'MMSE-X')))
    
    mW = cell(nUsers,nBands);
    mWL = cell(nUsers,nBands);
    
    for iBand = 1:nBands
        for iBase = bsIndices
            for iUser = cellUserIndices{iBase,1}'
                mW{iUser,iBand} = ones(SimParams.nRxAntenna,maxRank) * pertNoise;
            end
        end
    end
    
    for iBand = 1:nBands
        for iBase = bsIndices
            aUsers = SimParams.Debug.schTable{iBase,iBand}.assignedUsers;
            aStream = SimParams.Debug.schTable{iBase,iBand}.assignedStreams;
            for iStream = 1:length(aUsers)
                mW{aUsers(iStream,1),iBand}(:,aStream(iStream,1)) = W{aUsers(iStream,1),iBand}(:,aStream(iStream,1));
            end
        end
    end
    
    if strcmpi(rxType,'Last')
        
        for iBand = 1:nBands
            for iBase = bsIndices
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    [~,~,V] = svd(cH{iBase,iBand}(:,:,cUser));
                    M0{iBase,1}(:,:,iUser,iBand) = V(:,1:SimParams.maxRank);
                end
            end
        end
        
        for iBase = bsIndices
            totPower = norm(vec(M0{iBase,1}))^2;
            totPower = sqrt(sum(SimStructs.baseStruct{iBase,1}.sPower) / totPower);
            M0{iBase,1} = M0{iBase,1} * totPower;
        end

        for iBand = 1:nBands
            for iBase = bsIndices
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    for iLayer = 1:maxRank
                        if strcmpi(SimParams.sysMode,'true')
                            R = (SimParams.N + SimStructs.userStruct{cUser,1}.phyParams.restOfIF_Linear) * eye(SimParams.nRxAntenna);
                        else
                            R = SimParams.N * eye(SimParams.nRxAntenna);
                        end
                        for jBase = bsIndices
                            H = cH{jBase,iBand}(:,:,cUser);
                            for jUser = 1:usersPerCell(jBase,1)
                                R = R + H * M0{jBase,1}(:,:,jUser,iBand) * M0{jBase,1}(:,:,jUser,iBand)' * H';
                            end
                        end
                        H = cH{iBase,iBand}(:,:,cUser);
                        mWL{cUser,iBand}(:,iLayer) = R \ (H * M0{iBase,1}(:,iLayer,iUser,iBand)) + pertNoise;
                    end
                end
            end
        end
        
        for iBand = 1:nBands
            for iBase = bsIndices
                aUsers = SimParams.Debug.schTable{iBase,iBand}.assignedUsers;
                aStream = SimParams.Debug.schTable{iBase,iBand}.assignedStreams;
                for iStream = 1:length(aUsers)
                    mW{aUsers(iStream,1),iBand}(:,aStream(iStream,1)) = mWL{aUsers(iStream,1),iBand}(:,aStream(iStream,1));
                end
            end
        end
      
    end
    
    W = mW;
    
end

if strcmpi(updateHistory,'true')
    for iBand = 1:nBands
        for iBase = bsIndices
            for iUser = 1:usersPerCell(iBase,1);
                cUser = cellUserIndices{iBase,1}(iUser,1);
                SimStructs.userStruct{cUser,1}.W{iBand,1} = W{cUser,iBand};
                SimStructs.userStruct{cUser,1}.pW{iBand,1} = W{cUser,iBand};
            end
        end
    end
end
