
function [SimParams,SimStructs] = initializeSCApoint(SimParams,SimStructs,bsIndex)

nBands = SimParams.nBands;
maxRank = SimParams.maxRank;

if ~ischar(bsIndex)
    
    stringCell = strsplit(SimStructs.baseStruct{bsIndex,1}.selectionType,'_');
    if length(stringCell) == 1
        cH = SimStructs.linkChan;
    else
        cH = SimStructs.prevChan;
    end
    
    selectionType = stringCell{1};
    linkedUsers = SimStructs.baseStruct{bsIndex,1}.linkedUsers;
    kUsers = length(SimStructs.baseStruct{bsIndex,1}.linkedUsers);
    
    W0 = cell(SimParams.nUsers,nBands);
    B0 = zeros(SimParams.maxRank,kUsers,nBands);
    G0 = zeros(SimParams.maxRank,kUsers,nBands);
    addNoise = ones(SimParams.maxRank,SimParams.nUsers,nBands) * 1e-40;
    M0 = zeros(SimParams.nTxAntenna,SimParams.maxRank,kUsers,nBands);
    
    for iBand = 1:nBands
        for iUser = 1:SimParams.nUsers
            W0{iUser,iBand} = SimStructs.userStruct{iUser,1}.pW{iBand,1};
        end
    end
    
    switch selectionType
        
        case 'BF'
            for iBand = 1:nBands
                for iUser = 1:kUsers
                    cUser = linkedUsers(iUser,1);
                    [~,~,V] = svd(cH{bsIndex,iBand}(:,:,cUser));
                    M0(:,:,iUser,iBand) = V(:,1:SimParams.maxRank) * diag(sum(abs(W0{cUser,iBand}),1) > 1e-20);
                end
                totPower = norm(vec(M0(:,:,:,iBand)))^2;
                totPower = sqrt(sum(SimStructs.baseStruct{bsIndex,1}.sPower) / totPower);
                M0(:,:,:,iBand) = M0(:,:,:,iBand) * totPower;
            end
            
        case 'Ones'
            for iBand = 1:nBands
                M0(:,:,:,iBand) = complex(ones(SimParams.nTxAntenna,SimParams.maxRank,kUsers),ones(SimParams.nTxAntenna,SimParams.maxRank,kUsers));
                totPower = norm(vec(M0(:,:,:,iBand)))^2;
                totPower = sqrt(sum(SimStructs.baseStruct{bsIndex,1}.sPower) / totPower);
                M0(:,:,:,iBand) = M0(:,:,:,iBand) * totPower;
            end
            
        case 'Random'
            for iBand = 1:nBands
                M0(:,:,:,iBand) = complex(randn(SimParams.nTxAntenna,SimParams.maxRank,kUsers),randn(SimParams.nTxAntenna,SimParams.maxRank,kUsers));
                totPower = norm(vec(M0(:,:,:,iBand)))^2;
                totPower = sqrt(sum(SimStructs.baseStruct{bsIndex,1}.sPower) / totPower);
                M0(:,:,:,iBand) = M0(:,:,:,iBand) * totPower;
            end
            
        case 'Last'
            for iBand = 1:nBands
                M0(:,:,:,iBand) = SimParams.Debug.globalExchangeInfo.P{bsIndex,iBand};
                for iUser = 1:kUsers
                    cUser = linkedUsers(iUser,1);
                    for jBase = 1:SimParams.nBases
                        if jBase ~= bsIndex
                            addNoise(:,cUser,iBand) = addNoise(:,cUser,iBand) + SimParams.Debug.globalExchangeInfo.gI{jBase,1}(:,cUser,iBand).^2;
                        end
                    end
                end
            end
    end
    
    for iBand = 1:nBands
        for iUser = 1:kUsers
            cUser = linkedUsers(iUser,1);
            for iLayer = 1:SimParams.maxRank
                B0(iLayer,iUser,iBand) = SimParams.N * norm(W0{cUser,iBand}(:,iLayer))^2 + addNoise(iLayer,cUser,iBand);
                for jUser = 1:kUsers
                    if iUser ~= jUser
                        B0(iLayer,iUser,iBand) = B0(iLayer,iUser,iBand) + norm(W0{cUser,iBand}(:,iLayer)' * cH{bsIndex,iBand}(:,:,cUser) * M0(:,:,jUser,iBand),2)^2;
                    else
                        B0(iLayer,iUser,iBand) = B0(iLayer,iUser,iBand) + norm(W0{cUser,iBand}(:,iLayer)' * cH{bsIndex,iBand}(:,:,cUser) * M0(:,iLayer~=1:SimParams.maxRank,iUser,iBand),2)^2;
                    end
                end
                G0(iLayer,iUser,iBand) = norm(W0{cUser,iBand}(:,iLayer)' * cH{bsIndex,iBand}(:,:,cUser) * M0(:,iLayer,iUser,iBand),2)^2;
            end
        end
    end
    
    G0 = G0 ./ B0;
    T0 = log2(1 + G0);
    
    SimParams.Debug.globalExchangeInfo.funcOut{1,bsIndex} = M0;
    SimParams.Debug.globalExchangeInfo.funcOut{2,bsIndex} = B0;
    SimParams.Debug.globalExchangeInfo.funcOut{3,bsIndex} = G0;
    SimParams.Debug.globalExchangeInfo.funcOut{4,bsIndex} = T0;
    SimParams.Debug.globalExchangeInfo.funcOut{5,bsIndex} = W0;
    
else
    
    nBases = SimParams.nBases;
    
    linkedUsers = cell(nBases,1);
    usersPerCell = zeros(nBases,1);
    W0 = cell(SimParams.nUsers,nBands);
    E0 = cell(nBases,1);M0 = cell(nBases,1);
    
    for iBand = 1:nBands
        for iUser = 1:SimParams.nUsers
            W0{iUser,iBand} = SimStructs.userStruct{iUser,1}.pW{iBand,1};
        end
    end
    
    for iBase = 1:nBases
        linkedUsers{iBase,1} = SimStructs.baseStruct{iBase,1}.linkedUsers;
        usersPerCell(iBase,1) = length(SimStructs.baseStruct{iBase,1}.linkedUsers);
    end
    
    for bsIndex = 1:nBases
        
        stringCell = strsplit(SimStructs.baseStruct{bsIndex,1}.selectionType,'_');
        if length(stringCell) == 1
            cH = SimStructs.linkChan;
        else
            cH = SimStructs.prevChan;
        end
        selectionType = stringCell{1};
        
        switch selectionType
            
            case 'BF'
                for iBand = 1:nBands
                    for iUser = 1:usersPerCell(bsIndex,1)
                        cUser = linkedUsers{bsIndex}(iUser,1);
                        [~,~,V] = svd(cH{bsIndex,iBand}(:,:,cUser));
                        M0{bsIndex,1}(:,:,iUser,iBand) = V(:,1:SimParams.maxRank) * diag(sum(abs(W0{cUser,iBand}),1) > 1e-20);
                    end
                end
                totPower = norm(vec(M0{bsIndex,1}))^2;
                totPower = sqrt(sum(SimStructs.baseStruct{bsIndex,1}.sPower) / totPower);
                M0{bsIndex,1} = M0{bsIndex,1} * totPower;
                
            case 'Ones'
                M0{bsIndex,1} = complex(ones(SimParams.nTxAntenna,SimParams.maxRank,usersPerCell(bsIndex,1),nBands),ones(SimParams.nTxAntenna,SimParams.maxRank,usersPerCell(bsIndex,1),nBands));
                totPower = norm(vec(M0{bsIndex,1}))^2;
                totPower = sqrt(sum(SimStructs.baseStruct{bsIndex,1}.sPower) / totPower);
                M0{bsIndex,1} = M0{bsIndex,1} * totPower;
                
            case 'Random'
                M0{bsIndex,1} = complex(randn(SimParams.nTxAntenna,SimParams.maxRank,usersPerCell(bsIndex,1),nBands),randn(SimParams.nTxAntenna,SimParams.maxRank,usersPerCell(bsIndex,1),nBands));
                totPower = norm(vec(M0{bsIndex,1}))^2;
                totPower = sqrt(sum(SimStructs.baseStruct{bsIndex,1}.sPower) / totPower);
                M0{bsIndex,1} = M0{bsIndex,1} * totPower;
                
            case 'Last'
                for iBand = 1:nBands
                    M0{bsIndex,1}(:,:,:,iBand) = SimParams.Debug.globalExchangeInfo.P{bsIndex,iBand};
                end    
                
            case 'FrameC'
                for iBand = 1:nBands
                    M0{bsIndex,1}(:,:,:,iBand) = SimParams.Debug.globalExchangeInfo.P{bsIndex,iBand};
                end
                
                for iBand = 1:nBands
                    for iUser = 1:usersPerCell(bsIndex,1)
                        cUser = linkedUsers{bsIndex}(iUser,1);
                        [~,~,V] = svd(cH{bsIndex,iBand}(:,:,cUser));
                        for iRank = 1:maxRank
                            if (norm(M0{bsIndex,1}(:,iRank,iUser,iBand)) < 1e-10)
                                M0{bsIndex,1}(:,iRank,iUser,iBand) = V(:,iRank) * 0.001;
                            end
                        end
                    end
                end
                
                totPower = norm(vec(M0{bsIndex,1}))^2;
                totPower = sqrt(sum(SimStructs.baseStruct{bsIndex,1}.sPower) / totPower);
                M0{bsIndex,1} = M0{bsIndex,1} * totPower;

        end
        
    end
    
    for iBase = 1:nBases
        for iBand = 1:nBands
            for iUser = 1:usersPerCell(iBase,1)
                cUser = linkedUsers{iBase,1}(iUser,1);
                for iLayer = 1:maxRank
                    currentH = cH{iBase,iBand}(:,:,cUser);
                    E0{iBase,1}(iLayer,iUser,iBand) = real(1 - W0{cUser,iBand}(:,iLayer)' * currentH * M0{iBase,1}(:,iLayer,iUser,iBand));
                end
            end
        end
    end
    
    for iBase = 1:nBases
        SimParams.Debug.globalExchangeInfo.funcOut{1,iBase} = M0{iBase,1};
        SimParams.Debug.globalExchangeInfo.funcOut{2,iBase} = E0{iBase,1};
        SimParams.Debug.globalExchangeInfo.funcOut{5,iBase} = W0;
    end
    
end

