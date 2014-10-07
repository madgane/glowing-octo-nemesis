
function [SimParams,SimStructs] = initializeSCApoint(SimParams,SimStructs,bsIndex)

cH = SimStructs.linkChan;
nBands = SimParams.nBands;
maxRank = SimParams.maxRank;
rankArray = linspace(1,maxRank,maxRank);

if ~ischar(bsIndex)
        
    linkedUsers = SimStructs.baseStruct{bsIndex,1}.linkedUsers;
    kUsers = length(SimStructs.baseStruct{bsIndex,1}.linkedUsers);
    
    W0 = cell(SimParams.nUsers,nBands);
    B0 = zeros(SimParams.maxRank,kUsers,nBands);
    G0 = zeros(SimParams.maxRank,kUsers,nBands);
    addNoise = zeros(SimParams.maxRank,SimParams.nUsers,nBands);
    M0 = zeros(SimParams.nTxAntenna,SimParams.maxRank,kUsers,nBands);
    
    for iBand = 1:nBands
        for iUser = 1:SimParams.nUsers
            W0{iUser,iBand} = SimStructs.userStruct{iUser,1}.pW{iBand,1};
        end
    end
    
    switch SimStructs.baseStruct{bsIndex,1}.selectionType
        
        case 'BF'
            for iBand = 1:nBands
                for iUser = 1:kUsers
                    cUser = linkedUsers(iUser,1);
                    [~,~,V] = svd(cH{bsIndex,iBand}(:,:,cUser));
                    M0(:,:,iUser,iBand) = V(:,1:SimParams.maxRank);
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
            
        case 'Last'
            for iBand = 1:nBands
                M0(:,:,:,iBand) = SimParams.Debug.globalExchangeInfo.P{bsIndex,iBand};
                for iUser = 1:kUsers
                    cUser = linkedUsers(iUser,1);
                    for jBase = 1:SimParams.nBases
                        if jBase ~= bsIndex
                            addNoise(:,cUser,iBand) = SimParams.Debug.globalExchangeInfo.gI{jBase,1}(:,cUser,iBand).^2;
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
        
        switch SimStructs.baseStruct{bsIndex,1}.selectionType
            
            case 'BF'
                for iBand = 1:nBands
                    for iUser = 1:usersPerCell(bsIndex,1)
                        cUser = linkedUsers{bsIndex}(iUser,1);
                        [~,~,V] = svd(cH{bsIndex,iBand}(:,:,cUser));
                        M0{bsIndex,1}(:,:,iUser,iBand) = V(:,1:SimParams.maxRank);
                    end
                    totPower = norm(vec(M0{bsIndex,1}(:,:,:,iBand)))^2;
                    totPower = sqrt(sum(SimStructs.baseStruct{bsIndex,1}.sPower) / totPower);
                    M0{bsIndex,1}(:,:,:,iBand) = M0{bsIndex,1}(:,:,:,iBand) * totPower;
                end
                SimParams.Debug.globalExchangeInfo.funcOut{3,bsIndex} = ones(maxRank,usersPerCell(bsIndex,1),nBands);
                SimParams.Debug.globalExchangeInfo.funcOut{4,bsIndex} = ones(maxRank,usersPerCell(bsIndex,1),nBands);
                
            case 'Ones'
                for iBand = 1:nBands
                    M0{bsIndex,1}(:,:,:,iBand) = complex(ones(SimParams.nTxAntenna,SimParams.maxRank,usersPerCell(bsIndex,1)),ones(SimParams.nTxAntenna,SimParams.maxRank,usersPerCell(bsIndex,1)));
                    totPower = norm(vec(M0{bsIndex,1}(:,:,:,iBand)))^2;
                    totPower = sqrt(sum(SimStructs.baseStruct{bsIndex,1}.sPower) / totPower);
                    M0{bsIndex,1}(:,:,:,iBand) = M0{bsIndex,1}(:,:,:,iBand) * totPower;
                end
                SimParams.Debug.globalExchangeInfo.funcOut{3,bsIndex} = ones(maxRank,usersPerCell(bsIndex,1),nBands);
                SimParams.Debug.globalExchangeInfo.funcOut{4,bsIndex} = ones(maxRank,usersPerCell(bsIndex,1),nBands);

                
            case 'Last'
                for iBand = 1:nBands
                    M0{bsIndex,1}(:,:,:,iBand) = SimParams.Debug.globalExchangeInfo.P{bsIndex,iBand};
                end
        end
        
    end
    
    for iBase = 1:nBases
        for iBand = 1:nBands
            for iUser = 1:usersPerCell(iBase,1)
                cUser = linkedUsers{iBase,1}(iUser,1);
                for iLayer = 1:maxRank
                    intVector = sqrt(SimParams.N) * W0{cUser,iBand}(:,iLayer)';
                    for jBase = 1:nBases
                        for jUser = 1:usersPerCell(jBase,1)
                            jxUser = linkedUsers{jBase,1}(jUser,1);
                            currentH = cH{jBase,iBand}(:,:,cUser);
                            if jxUser ~= cUser
                                intVector = [intVector, W0{cUser,iBand}(:,iLayer)' * currentH * M0{jBase,1}(:,iLayer~=rankArray,jUser,iBand)];
                            else
                                intVector = [intVector, W0{cUser,iBand}(:,iLayer)' * currentH * M0{jBase,1}(:,:,jUser,iBand)];
                            end
                        end
                    end
                    
                    currentH = cH{iBase,iBand}(:,:,cUser);
                    givenVector = (1 - W0{cUser,iBand}(:,iLayer)' * currentH * M0{iBase,1}(:,iLayer,iUser,iBand));
                    intVector = [intVector, givenVector];
                    E0{iBase,1}(iLayer,iUser,iBand) = norm(intVector)^2;
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

