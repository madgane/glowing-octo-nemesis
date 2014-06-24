function [SimParams,SimStructs] = getXScheduling(SimParams,SimStructs)

chScheduler = char(SimParams.SchedType);
uScoreIndex = find(chScheduler == '_');
if isempty(uScoreIndex)
    scheduleMethod = SimParams.SchedType;
else
    scheduleMethod = chScheduler(uScoreIndex(1,1) + 1:end);
end

for iBand = 1:SimParams.nBands
    
    switch scheduleMethod
        
        case 'EqualShare'
            
            for iBase = 1:SimParams.nBases
                
                cUsers = SimStructs.baseStruct{iBase,1}.linkedUsers;
                eH = zeros(length(cUsers),SimParams.nTxAntenna);
                if SimParams.queueWt
                    for iUser = 1:length(cUsers)
                        eH(iUser,:) = SimStructs.linkChan{iBase,iBand}(:,:,cUsers(iUser,1)) * (SimStructs.userStruct{cUsers(iUser,1),1}.weighingFactor);
                    end
                else
                    for iUser = 1:length(cUsers)
                        eH(iUser,:) = SimStructs.linkChan{iBase,iBand}(:,:,cUsers(iUser,1)) * sign(SimStructs.userStruct{cUsers(iUser,1),1}.weighingFactor);
                    end
                end
                augE = reshape(eH(:),SimParams.nTxAntenna,length(cUsers));
                
                [~,~,sortA] = qr(augE,0);
                muxIFFreeRank = SimParams.muxRank / SimParams.nBases;
                schedUsers = cUsers(sortA(1,1:muxIFFreeRank),1);
                SimStructs.baseStruct{iBase}.assignedUsers{iBand,1} = schedUsers;
                SimStructs.baseStruct{iBase}.assignedStreams{iBand,1} = ones(muxIFFreeRank,1);
                
            end
            
        case 'GroupSS'
            
            nIterations = 1;
            muxIFFreeRank = SimParams.muxRank / SimParams.nBases;
            activeUsers = cell(SimParams.nBases,1);
            
            for iIter = 1:nIterations
                
                subPspace = cell(SimParams.nBases,1);
                
                for iBase = 1:SimParams.nBases
                    uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
                    actSet = uIndices(activeUsers{iBase,1});
                    for jBase = 1:SimParams.nBases
                        if iBase ~= jBase
                            eH = SimStructs.linkChan{jBase,iBand}(:,:,actSet);
                            X = reshape(eH(:),SimParams.nTxAntenna,length(actSet));
                            subPspace{jBase,1} = [subPspace{jBase,1} X];
                        end
                    end
                end
                
                Aspace = cell(SimParams.nBases,1);
                
                for iBase = 1:SimParams.nBases
                    
                    uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
                    eH = SimStructs.linkChan{iBase,iBand}(:,:,uIndices);kUsers = length(uIndices);
                    X = reshape(eH(:),SimParams.nTxAntenna * SimParams.nRxAntenna,kUsers);
                    
                    
                    ppVolume = zeros(1,kUsers);
                    
                    
                    for iStream = 1:muxIFFreeRank
                        for iUser = 1:kUsers
                            U = [subPspace{iBase,1} , Aspace{iBase,1} , X(:,iUser)];
                            ppVolume(1,iUser) = real(det(U' * U));
                        end
                        
                        [~,sortI] = sort(ppVolume,'descend');
                        activeUsers{iBase,1}(1,iStream) = sortI(1,1);
                        Aspace{iBase,1} = X(:,activeUsers{iBase,1}(1,1:iStream));
                    end
                end
                
            end
            
            for iBase = 1:SimParams.nBases
                uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
                schedUsers = uIndices(activeUsers{iBase,1},1);
                SimStructs.baseStruct{iBase}.assignedUsers{iBand,1} = schedUsers;
                SimStructs.baseStruct{iBase}.assignedStreams{iBand,1} = ones(muxIFFreeRank,1);
            end
            
        case 'InstaSS'
            
            nIterations = 10;
            subPspace = cell(SimParams.nBases,1);
            muxIFFreeRank = SimParams.muxRank / SimParams.nBases;
            activeUsers = cell(SimParams.nBases,1);
            
            for iIter = 1:nIterations
                
                Aspace = cell(SimParams.nBases,1);
                
                for iBase = 1:SimParams.nBases
                    
                    uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
                    eH = SimStructs.linkChan{iBase,iBand}(:,:,uIndices);kUsers = length(uIndices);
                    X = reshape(eH(:),SimParams.nTxAntenna * SimParams.nRxAntenna,kUsers);
                    
                    ppVolume = zeros(1,kUsers);
                    for iStream = 1:muxIFFreeRank
                        for iUser = 1:kUsers
                            U = [subPspace{iBase,1} , Aspace{iBase,1} , X(:,iUser)];
                            ppVolume(1,iUser) = real(det(U' * U));
                        end
                        
                        [~,sortI] = sort(ppVolume,'descend');
                        activeUsers{iBase,1}(1,iStream) = sortI(1,1);
                        Aspace{iBase,1} = X(:,activeUsers{iBase,1}(1,1:iStream));
                    end
                    
                    subPspace = cell(SimParams.nBases,1);
                    for cBase = 1:SimParams.nBases
                        mIndices = SimStructs.baseStruct{cBase,1}.linkedUsers;
                        actSet = mIndices(activeUsers{cBase,1},1);
                        for jBase = 1:SimParams.nBases
                            if cBase ~= jBase
                                eH = SimStructs.linkChan{jBase,iBand}(:,:,actSet);
                                X = reshape(eH(:),SimParams.nTxAntenna,length(actSet));
                                subPspace{jBase,1} = [subPspace{jBase,1}  X];
                            end
                        end
                    end
                end
            end
            
            for iBase = 1:SimParams.nBases
                uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
                schedUsers = uIndices(activeUsers{iBase,1},1);
                SimStructs.baseStruct{iBase}.assignedUsers{iBand,1} = schedUsers;
                SimStructs.baseStruct{iBase}.assignedStreams{iBand,1} = ones(muxIFFreeRank,1);
            end
            
        case 'OrthoSS'
            
            nIterations = 10;
            subPspace = cell(SimParams.nBases,1);
            subDspace = cell(SimParams.nBases,1);
            muxIFFreeRank = SimParams.muxRank / SimParams.nBases;
            activeUsers = cell(SimParams.nBases,1);
            
            for iIter = 1:nIterations
                
                Aspace = cell(SimParams.nBases,1);
                Dspace = cell(SimParams.nBases,1);
                
                for mBase = 1:SimParams.nBases
                    if isempty(subPspace{mBase,1})
                        subPspace{mBase,1} = zeros(SimParams.nTxAntenna);
                    end
                    if isempty(subDspace{mBase,1})
                        subDspace{mBase,1} = zeros(SimParams.nTxAntenna);
                    end
                    if isempty(Aspace{mBase,1})
                        Aspace{mBase,1} = zeros(SimParams.nTxAntenna);
                    end
                    if isempty(Dspace{mBase,1})
                        Dspace{mBase,1} = zeros(SimParams.nTxAntenna);
                    end
                end
                
                
                for iBase = 1:SimParams.nBases
                    
                    uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
                    eH = SimStructs.linkChan{iBase,iBand}(:,:,uIndices);kUsers = length(uIndices);
                    X = reshape(eH(:),SimParams.nTxAntenna * SimParams.nRxAntenna,kUsers);
                    
                    ppVolumeN = zeros(1,kUsers);
                    ppVolumeD = zeros(1,kUsers);
                    for iStream = 1:muxIFFreeRank
                        for iUser = 1:kUsers
                            U = subPspace{iBase,1} + Aspace{iBase,1};
                            K = subDspace{iBase,1} + Dspace{iBase,1};
                            if ~isempty(U)
                                ppVolumeN(1,iUser) = norm(U' * X(:,iUser));
                            else
                                ppVolumeN(1,iUser) = norm(X(:,iUser));
                            end
                            if ~isempty(K)
                                ppVolumeD(1,iUser) = norm(K' * X(:,iUser));
                            else
                                ppVolumeD(1,iUser) = 1;
                            end
                        end
                        
                        if iStream ~= 1
                            ppVolume = atan(ppVolumeN ./ ppVolumeD);
                        else
                            ppVolume = ppVolumeN;
                        end
                        
                        [~,sortI] = sort(ppVolume,'descend');
                        activeUsers{iBase,1}(1,iStream) = sortI(1,1);
                        M = X(:,activeUsers{iBase,1}(1,1:iStream));
                        Aspace{iBase,1} = eye(SimParams.nTxAntenna) - M * inv(M' * M) * M';
                        Dspace{iBase,1} = M * inv(M' * M) * M';
                    end
                    
                    subPspace = cell(SimParams.nBases,1);
                    for mBase = 1:SimParams.nBases
                        if isempty(subPspace{mBase,1})
                            subPspace{mBase,1} = zeros(SimParams.nTxAntenna);
                        end
                        if isempty(subDspace{mBase,1})
                            subDspace{mBase,1} = zeros(SimParams.nTxAntenna);
                        end
                        if isempty(Aspace{mBase,1})
                            Aspace{mBase,1} = zeros(SimParams.nTxAntenna);
                        end
                        if isempty(Dspace{mBase,1})
                            Dspace{mBase,1} = zeros(SimParams.nTxAntenna);
                        end
                    end
                    
                    
                    for cBase = 1:SimParams.nBases
                        mIndices = SimStructs.baseStruct{cBase,1}.linkedUsers;
                        actSet = mIndices(activeUsers{cBase,1},1);
                        for jBase = 1:SimParams.nBases
                            if cBase ~= jBase
                                eH = SimStructs.linkChan{jBase,iBand}(:,:,actSet);
                                X = reshape(eH(:),SimParams.nTxAntenna,length(actSet));
                                subPspace{jBase,1} = subPspace{jBase,1} + (eye(SimParams.nTxAntenna) - X * inv(X' * X) * X');
                                subDspace{jBase,1} = subPspace{jBase,1} + X * inv(X' * X) * X';
                            end
                        end
                    end
                end
                
            end
            
            for iBase = 1:SimParams.nBases
                uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
                schedUsers = uIndices(activeUsers{iBase,1},1);
                SimStructs.baseStruct{iBase}.assignedUsers{iBand,1} = schedUsers;
                SimStructs.baseStruct{iBase}.assignedStreams{iBand,1} = ones(muxIFFreeRank,1);
            end
            
        case 'StreamSearch'
            
            nIterations = 10;
            muxIFFreeRank = SimParams.muxRank / SimParams.nBases;
            
            Nspace = cell(SimParams.nBases,1);
            activeUsers = cell(SimParams.nBases,1);
            activeStreams = cell(SimParams.nBases,1);
            completeH = cell(SimParams.nUsers,SimParams.nBases);
            
            for iUser = 1:SimParams.nUsers
                
                bNode = SimStructs.userStruct{iUser,1}.baseNode;
                Hd = SimStructs.linkChan{bNode,iBand}(:,:,iUser);
                [U,~,~] = svd(Hd);
                
                for iBase = 1:SimParams.nBases
                    Hk = SimStructs.linkChan{iBase,iBand}(:,:,iUser);M = U' * Hk;
                    if SimParams.queueWt
                        M = M(1:SimParams.maxRank,:).' * (SimStructs.userStruct{iUser,1}.weighingFactor);
                    else
                        M = M(1:SimParams.maxRank,:).' * sign(SimStructs.userStruct{iUser,1}.weighingFactor);
                    end
                    
                    completeH{iUser,iBase} = [completeH{iUser,iBase} M];
                end
                
            end
            
            uLocs = cell(SimParams.nBases,1);
            for iBase = 1:SimParams.nBases
                uLocs{iBase,1} = zeros(SimParams.nUsers * SimParams.maxRank,2);
            end
            
            for iIter = 1:nIterations
                
                Dspace = cell(SimParams.nBases,1);
                for iBase = 1:SimParams.nBases
                    lkUsers = SimStructs.baseStruct{iBase,1}.linkedUsers;
                    ppVolume = zeros(1,length(lkUsers) * SimParams.maxRank);
                    for iStream = 1:muxIFFreeRank
                        iIndex = 0;
                        for iUser = 1:length(lkUsers)
                            for jStream = 1:SimParams.maxRank
                                iIndex = iIndex + 1;
                                U = [Nspace{iBase,1} Dspace{iBase,1} completeH{lkUsers(iUser,1),iBase}(:,jStream)];
                                ppVolume(1,iIndex) = real(det(U' * U));
                                uLocs{iBase,1}(iIndex,:) = [lkUsers(iUser,1),jStream];
                            end
                        end
                        
                        [~,sortI] = sort(ppVolume,'descend');maxI = sortI(1,1);
                        activeUsers{iBase,1}(1,iStream) = uLocs{iBase,1}(maxI,1);
                        activeStreams{iBase,1}(1,iStream) = uLocs{iBase,1}(maxI,2);
                        Dspace{iBase,1} = [Dspace{iBase,1} completeH{activeUsers{iBase,1}(1,iStream),iBase}(:,activeStreams{iBase,1}(1,iStream))];
                    end
                    
                    Nspace = cell(SimParams.nBases,1);
                    for kBase = 1:SimParams.nBases
                        for jBase = 1:SimParams.nBases
                            if kBase ~= jBase
                                for iUser = 1:length(activeUsers{kBase,1})
                                    cUser = activeUsers{kBase,1}(1,iUser);cStream = activeStreams{kBase,1}(1,iUser);
                                    Nspace{jBase,1} = [Nspace{jBase,1} completeH{cUser,jBase}(:,cStream)];
                                end
                            end
                        end
                    end
                end
            end
            
            for iBase = 1:SimParams.nBases
                SimStructs.baseStruct{iBase}.assignedUsers{iBand,1} = activeUsers{iBase,1}';
                SimStructs.baseStruct{iBase}.assignedStreams{iBand,1} = activeStreams{iBase,1}';
            end
            
        case 'IterativeSearch'
            
            Haug = [];iIndex = 0;
            Ucell = cell(SimParams.nUsers,SimParams.nBases);
            uLocs = zeros(SimParams.nUsers * SimParams.maxRank * SimParams.nBases,3);
            
            linkedUsers = cell(SimParams.nBases,1);
            for iBase = 1:SimParams.nBases
                linkedUsers{iBase,1} = SimStructs.baseStruct{iBase,1}.linkedUsers;
            end
            
            for iUser = 1:SimParams.nUsers
                for iBase = 1:SimParams.nBases
                    
                    H = SimStructs.linkChan{iBase,iBand}(:,:,iUser);
                    [U,~,~] = svd(H);M = U' * H;
                    Ucell{iUser,iBase} = U;
                    
                    if SimParams.queueWt
                        M = M.' * (SimStructs.userStruct{iUser,1}.weighingFactor);
                    else
                        M = M.' * sign(SimStructs.userStruct{iUser,1}.weighingFactor);
                    end
                    
                    if ~sum(linkedUsers{iBase,1} == iUser)
                        M = M * 0;
                    end
                    
                    Haug = [Haug M(:,1:SimParams.maxRank)];
                    
                    for iStream = 1:SimParams.maxRank
                        iIndex = iIndex + 1;
                        uLocs(iIndex,:) = [iUser iStream iBase];
                    end
                    
                end
            end
            
            Hm = Haug;
            nIterations = 100;
            xIndex = SimParams.nUsers * SimParams.maxRank * SimParams.nBases;
            
            ppVolume = zeros(xIndex,1);
            Nspace = cell(SimParams.nBases,1);
            assIndex = cell(SimParams.nBases,SimParams.nBands);
            assUsers = cell(SimParams.nBases,SimParams.nBands);
            assStreams = cell(SimParams.nBases,SimParams.nBands);
            
            for iIteration = 1:nIterations
                
                Aspace = cell(SimParams.nBases,1);
                userIterIndex = zeros(SimParams.nBases,1);
                
                for iRank = 1:SimParams.muxRank
                    for iUser = 1:xIndex
                        cBase = uLocs(iUser,3);
                        M = [Aspace{cBase,1} Nspace{cBase,1} Hm(:,iUser)];
                        ppVolume(iUser,1) = real(det(M' * M));
                    end
                    
                    rIndex = 1;continueAgain = 1;
                    [~,sortI] = sort(ppVolume,'descend');maxI = sortI(1,1);
                    
                    while continueAgain
                        maxI = sortI(rIndex,1);
                        cUser = uLocs(maxI,1);cStream = uLocs(maxI,2);cBase = uLocs(maxI,3);
                        [~,totalLength] = size(Aspace{cBase,1});totalLength = totalLength + 1;
                        if totalLength > (SimParams.muxRank / SimParams.nBases)
                            continueAgain = 1;
                            rIndex = rIndex + 1;
                        else
                            continueAgain = 0;
                            userIterIndex(cBase,1) = userIterIndex(cBase,1) + 1;
                        end
                    end
                    
                    Aspace{cBase,1} = [Aspace{cBase,1} Hm(:,maxI)];
                    assIndex{cBase,iBand}(userIterIndex(cBase,1),1) = maxI;
                    assUsers{cBase,iBand}(userIterIndex(cBase,1),1) = cUser;
                    assStreams{cBase,iBand}(userIterIndex(cBase,1),1) = cStream;
                    
                    Nspace = cell(SimParams.nBases,1);
                    for jBase = 1:SimParams.nBases
                        for iBase = 1:SimParams.nBases
                            if iBase ~= jBase
                                for iUser = 1:length(assUsers{jBase,1})
                                    kUser = assUsers{jBase,1}(iUser,1);kStream = assStreams{jBase,1}(iUser,1);
                                    M = Ucell{kUser,jBase}(:,kStream)' * SimStructs.linkChan{iBase,iBand}(:,:,kUser);
                                    Nspace{iBase,1} = [Nspace{iBase,1} M.'];
                                end
                            end
                        end
                    end
                end
            end
            
            for iBase = 1:SimParams.nBases
                SimStructs.baseStruct{iBase}.assignedUsers{iBand,1} = assUsers{iBase,iBand};
                SimStructs.baseStruct{iBase}.assignedStreams{iBand,1} = assStreams{iBase,iBand};
            end
            
        case 'GroupProjection'
            
            nIterations = 10;
            muxIFFreeRank = SimParams.muxRank;
            
            Nspace = cell(SimParams.nBases,1);
            activeUsers = cell(SimParams.nBases,1);
            activeStreams = cell(SimParams.nBases,1);
            completeH = cell(SimParams.nUsers,SimParams.nBases);
            
            for iUser = 1:SimParams.nUsers
                
                bNode = SimStructs.userStruct{iUser,1}.baseNode;
                Hd = SimStructs.linkChan{bNode,iBand}(:,:,iUser);
                [U,~,~] = svd(Hd);
                
                for iBase = 1:SimParams.nBases
                    Hk = SimStructs.linkChan{iBase,iBand}(:,:,iUser);M = U' * Hk;
                    if SimParams.queueWt
                        M = M(1:SimParams.maxRank,:).' * (SimStructs.userStruct{iUser,1}.weighingFactor);
                    else
                        M = M(1:SimParams.maxRank,:).' * sign(SimStructs.userStruct{iUser,1}.weighingFactor);
                    end
                    
                    completeH{iUser,iBase} = [completeH{iUser,iBase} M];
                end
                
            end
            
            uLocs = cell(SimParams.nBases,1);
            for iBase = 1:SimParams.nBases
                uLocs{iBase,1} = zeros(SimParams.nUsers * SimParams.maxRank,2);
            end
            
            for iIter = 1:nIterations
                
                Dspace = cell(SimParams.nBases,1);
                for iBase = 1:SimParams.nBases
                    lkUsers = SimStructs.baseStruct{iBase,1}.linkedUsers;
                    ppVolume = zeros(1,length(lkUsers) * SimParams.maxRank);
                    for iStream = 1:muxIFFreeRank
                        iIndex = 0;
                        for iUser = 1:length(lkUsers)
                            for jStream = 1:SimParams.maxRank
                                iIndex = iIndex + 1;
                                U = [Nspace{iBase,1} Dspace{iBase,1} completeH{lkUsers(iUser,1),iBase}(:,jStream)];
                                ppVolume(1,iIndex) = real(det(U' * U));
                                uLocs{iBase,1}(iIndex,:) = [lkUsers(iUser,1),jStream];
                            end
                        end
                        
                        [~,sortI] = sort(ppVolume,'descend');maxI = sortI(1,1);
                        activeUsers{iBase,1}(1,iStream) = uLocs{iBase,1}(maxI,1);
                        activeStreams{iBase,1}(1,iStream) = uLocs{iBase,1}(maxI,2);
                        Dspace{iBase,1} = [Dspace{iBase,1} completeH{activeUsers{iBase,1}(1,iStream),iBase}(:,activeStreams{iBase,1}(1,iStream))];
                    end
                    
                    Nspace = cell(SimParams.nBases,1);
                    for kBase = 1:SimParams.nBases
                        for jBase = 1:SimParams.nBases
                            if kBase ~= jBase
                                for iUser = 1:length(activeUsers{kBase,1})
                                    cUser = activeUsers{kBase,1}(1,iUser);cStream = activeStreams{kBase,1}(1,iUser);
                                    Nspace{jBase,1} = [Nspace{jBase,1} completeH{cUser,jBase}(:,cStream)];
                                end
                            end
                        end
                    end
                end
            end
            
            for iBase = 1:SimParams.nBases
                SimStructs.baseStruct{iBase}.assignedUsers{iBand,1} = activeUsers{iBase,1}';
                SimStructs.baseStruct{iBase}.assignedStreams{iBand,1} = activeStreams{iBase,1}';
            end
           
            
    end
    
    
end
