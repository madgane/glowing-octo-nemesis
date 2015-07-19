function [SimParams,SimStructs] = getRealTimeQWSRM(SimParams,SimStructs)

SimParams.Debug.reSchedule = 'true';

proLogue;
if SimParams.iDrop == 1
    for iUser = 1:nUsers
        SimStructs.userStruct{iUser,1}.pW = cell(nBands,1);
    end
    
    SimParams.Debug.globalExchangeInfo.D = cell(nBases,1);
    SimParams.Debug.globalExchangeInfo.I = cell(nBases,1);
    SimParams.Debug.globalExchangeInfo.gI = cell(nBases,1);
    SimParams.Debug.globalExchangeInfo.P = cell(nBases,nBands);
    SimParams.Debug.globalExchangeInfo.funcOut = cell(5,nBases);
end

switch selectionMethod
    
    case 'perBSAlloc'
        
        maxIterations = 25;
        for iBase = 1:nBases
            SimParams.Debug.globalExchangeInfo.gI{iBase,1} = zeros(maxRank,nUsers,nBands);
        end
        
        isTDM = SimParams.additionalParams;
        
        for iExchangeOTA = 1:1
            
            for iBase = 1:nBases
                
                if strcmpi(isTDM,'TDM')
                    if (mod(SimParams.iDrop - 1,nBases) ~= (iBase - 1))
                        for iBand = 1:nBands
                            SimStructs.baseStruct{iBase,1}.P{iBand,1} = zeros(SimParams.nTxAntenna,maxRank,usersPerCell(iBase,1));
                            SimParams.Debug.globalExchangeInfo.P{iBase,iBand} = zeros(SimParams.nTxAntenna,maxRank,usersPerCell(iBase,1));
                        end
                        [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE',iBase);
                        continue;
                    end
                end
                
                cvx_objective = -50;
                for inLoop = 1:maxIterations
                    
                    kUsers = usersPerCell(iBase,1);
                    if inLoop == 1
                        SimStructs.baseStruct{iBase,1}.selectionType = 'BF';
                        [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'BF',iBase);
                        [SimParams, SimStructs] = initializeSCApoint(SimParams,SimStructs,iBase);
                    else
                        SimStructs.baseStruct{iBase,1}.selectionType = 'Last';
                        [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE',iBase);
                        [SimParams, SimStructs] = initializeSCApoint(SimParams,SimStructs,iBase);
                    end
                    
                    M0 = SimParams.Debug.globalExchangeInfo.funcOut{1,iBase};B0 = SimParams.Debug.globalExchangeInfo.funcOut{2,iBase};
                    W0 = SimParams.Debug.globalExchangeInfo.funcOut{5,iBase};
                    
                    cvx_begin
                    
                    variables epiObjective userObjective(kUsers,1)
                    variables T(maxRank,kUsers,nBands) B(maxRank,kUsers,nBands) G(maxRank,kUsers,nBands)
                    variable M(SimParams.nTxAntenna,maxRank,kUsers,nBands) complex
                    
                    for iUser = 1:kUsers
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        userWts(cUser,1) * abs(QueuedPkts(cUser,1) - sum(vec(T(:,iUser,:)))) <= userObjective(iUser,1);
                    end
                    
                    epiObjective >= norm(userObjective,qExponent);
                    minimize(epiObjective);
                    
                    for iBand = 1:nBands
                        for iUser = 1:kUsers
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for iLayer = 1:maxRank
                                intVector = sqrt(SimParams.N) * W0{cUser,iBand}(:,iLayer)';
                                for jUser = 1:kUsers
                                    if jUser ~= iUser
                                        intVector = [intVector, W0{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * M(:,:,jUser,iBand)];
                                    else
                                        intVector = [intVector, W0{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * M(:,iLayer ~= rankArray,iUser,iBand)];
                                    end
                                end
                                
                                intVector * intVector' <= B(iLayer,iUser,iBand);
                                log(1 + G(iLayer,iUser,iBand)) >= T(iLayer,iUser,iBand) * log(2);
                                
                                P = real(W0{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * M(:,iLayer,iUser,iBand));
                                Q = imag(W0{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * M(:,iLayer,iUser,iBand));
                                P0 = real(W0{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * M0(:,iLayer,iUser,iBand));
                                Q0 = imag(W0{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * M0(:,iLayer,iUser,iBand));
                                (P0^2 + Q0^2) / B0(iLayer,iUser,iBand) + (2 / B0(iLayer,iUser,iBand)) * (P0 * (P - P0) + Q0 * (Q - Q0)) ...
                                    - ((P0^2 + Q0^2) / (B0(iLayer,iUser,iBand)^2)) * (B(iLayer,iUser,iBand) - B0(iLayer,iUser,iBand)) >= G(iLayer,iUser,iBand);
                            end
                        end
                    end
                    
                    if strcmpi(isTDM,'TDM')
                        norm(vec(M))  <= sqrt(sum(SimStructs.baseStruct{iBase,1}.sPower(1,:)) * SimParams.nBases);
                    else
                        norm(vec(M))  <= sqrt(sum(SimStructs.baseStruct{iBase,1}.sPower(1,:)));
                    end
                    
                    cvx_end
                    
                    if strfind(cvx_status,'Solved')
                        M0 = full(M);
                    else
                        display(cvx_status);
                        for iBand = 1:nBands
                            M0(:,:,:,iBand) = SimParams.Debug.globalExchangeInfo.P{iBase,iBand} / sqrt(2);
                        end
                    end
                    
                    if norm(cvx_objective - cvx_optval,1) <= epsilonT
                        break;
                    else
                        cvx_objective = cvx_optval;
                    end
                    
                    for iBand = 1:nBands
                        SimStructs.baseStruct{iBase,1}.P{iBand,1} = M0(:,:,:,iBand);
                        SimParams.Debug.globalExchangeInfo.P{iBase,iBand} = M0(:,:,:,iBand);
                    end
                    
                end
            end
            
            [SimParams,SimStructs] = updateIteratePerformance(SimParams,SimStructs);
        end
        
    case 'distBSAlloc'
        
        stInstant = 0;
        SimParams.currentQueue = 100;
        
        if SimParams.nExchangesOTA ~= 0
            SimParams.BITFactor = 1 - (SimParams.nExchangesOTA / SimParams.nSymbolsBIT);
        end
        
        for iExchangeOTA = stInstant:SimParams.nExchangesOTA
            
            stepFactor = 2;
            switch iExchangeOTA
                
                case -1
                    for iBase = 1:nBases
                        if or((SimParams.distIteration - 1) == 0,mod((SimParams.iDrop - 1),SimParams.exchangeResetInterval) == 0)
                            fprintf('Resetting History for BS - %d \n',iBase);
                            SimStructs.baseStruct{iBase,1}.selectionType = 'BF_Prev';
                            [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE-BF_Prev',iBase);
                            SimParams.Debug.globalExchangeInfo.gI{iBase,1} = zeros(maxRank,nUsers,nBands);
                            SimParams.Debug.globalExchangeInfo.D{iBase,1} = ones(maxRank,nUsers,nBands,nBases);
                        else
                            if iBase == 1
                                fprintf('Reusing History \n');
                            end
                            SimStructs.baseStruct{iBase,1}.selectionType = 'Last_Prev';
                            [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'Last_Prev',iBase);
                        end
                    end
                    cH = SimStructs.prevChan;
                    maxBackHaulExchanges = SimParams.nExchangesOBH;
                case 0
                    if stInstant == 0
                        for iBase = 1:nBases
                            if or((SimParams.distIteration - 1) == 0,mod((SimParams.iDrop - 1),SimParams.exchangeResetInterval) == 0)
                                fprintf('Resetting History for BS - %d \n',iBase);
                                SimStructs.baseStruct{iBase,1}.selectionType = 'BF';
                                [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE-BF',iBase);
                                SimParams.Debug.globalExchangeInfo.gI{iBase,1} = zeros(maxRank,nUsers,nBands);
                                SimParams.Debug.globalExchangeInfo.D{iBase,1} = zeros(maxRank,nUsers,nBands,nBases);
                            else
                                if iBase == 1
                                    fprintf('Reusing History \n');
                                end
                                SimStructs.baseStruct{iBase,1}.selectionType = 'Last';
                                [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'Last',iBase);
                            end
                        end
                    else
                        for iBase = 1:nBases
                            SimStructs.baseStruct{iBase,1}.selectionType = 'Last';
                            [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'Last',iBase);
                        end
                    end
                    cH = SimStructs.linkChan;
                    maxBackHaulExchanges = SimParams.nExchangesOBH;
                    fprintf('OTA Performed - %d \n',iExchangeOTA);
                otherwise
                    cH = SimStructs.linkChan;
                    maxBackHaulExchanges = SimParams.nExchangesOBH;
                    fprintf('OTA Performed - %d \n',iExchangeOTA);
                    [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE');
            end
            
            for iExchangeBH = 1:maxBackHaulExchanges
                
                for iBase = 1:nBases
                    
                    if iExchangeBH ~= 1
                        SimStructs.baseStruct{iBase,1}.selectionType = 'Last';
                    end
                    
                    if and(strcmpi(SimParams.additionalParams,'H-MMSE'),(iExchangeBH ~= 1))
                        [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE-XVAR',iBase);
                        W0 = SimParams.Debug.globalExchangeInfo.funcOut{6,iBase};
                    else
                        for iBand = 1:nBands
                            for iUser = 1:SimParams.nUsers
                                W0{iUser,iBand} = SimStructs.userStruct{iUser,1}.pW{iBand,1};
                            end
                        end
                    end
                    
                    kUsers = usersPerCell(iBase,1);
                    SimParams.Debug.exchangeIndex = iExchangeBH + iExchangeOTA;
                    [SimParams, SimStructs] = initializeSCApoint(SimParams,SimStructs,iBase);
                    M0 = SimParams.Debug.globalExchangeInfo.funcOut{1,iBase};B0 = SimParams.Debug.globalExchangeInfo.funcOut{2,iBase};
                    
                    cvx_begin
                    
                    expression T(maxRank,kUsers,nBands)
                    
                    variable M(SimParams.nTxAntenna,maxRank,kUsers,nBands) complex
                    variables Tx(maxRank,kUsers,nBands) B(maxRank,kUsers,nBands) G(maxRank,kUsers,nBands)
                    variables I(maxRank,nUsers,nBands,nBases) userObjective(kUsers,1) epiObjective
                    
                    T = SimParams.BITFactor * Tx;
                    
                    for iUser = 1:kUsers
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        userWts(cUser,1) * abs(QueuedPkts(cUser,1) - sum(vec(T(:,iUser,:)))) <= userObjective(iUser,1);
                    end
                    
                    augmentedTerms = 0;
                    for jBase = 1:nBases
                        if jBase ~= iBase
                            vecA = SimParams.Debug.globalExchangeInfo.gI{jBase,1}(:,cellUserIndices{iBase,1},:) - I(:,cellUserIndices{iBase,1},:,jBase);
                            vecB = vecA .* SimParams.Debug.globalExchangeInfo.D{iBase,1}(:,cellUserIndices{iBase,1},:,jBase);
                            augmentedTerms = augmentedTerms + sum(vecB(:)) + stepFactor * 0.5 * sum(pow_abs(vecA(:),2));
                            
                            vecA = SimParams.Debug.globalExchangeInfo.gI{iBase,1}(:,cellUserIndices{jBase,1},:) - I(:,cellUserIndices{jBase,1},:,iBase);
                            vecB = vecA .* SimParams.Debug.globalExchangeInfo.D{iBase,1}(:,cellUserIndices{jBase,1},:,iBase);
                            augmentedTerms = augmentedTerms + sum(vecB(:)) + stepFactor * 0.5 * sum(pow_abs(vecA(:),2));
                        end
                    end
                    
                    epiObjective >= norm(userObjective,qExponent) + augmentedTerms;
                    minimize(epiObjective);
                    
                    for iBand = 1:nBands
                        for iUser = 1:kUsers
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for iLayer = 1:maxRank
                                
                                intVector = sqrt(SimParams.N) * W0{cUser,iBand}(:,iLayer)';
                                for jUser = 1:kUsers
                                    if jUser ~= iUser
                                        intVector = [intVector, W0{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * M(:,:,jUser,iBand)];
                                    else
                                        intVector = [intVector, W0{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * M(:,iLayer ~= rankArray,iUser,iBand)];
                                    end
                                end
                                
                                for jBase = 1:nBases
                                    if jBase ~= iBase
                                        intVector = [intVector, I(iLayer,cUser,iBand,jBase)];
                                    end
                                end
                                
                                norm(intVector) <= sqrt(B(iLayer,iUser,iBand));
                                log(1 + G(iLayer,iUser,iBand)) >= T(iLayer,iUser,iBand) * log(2);
                                
                                for jUser = 1:nUsers
                                    nCellIndex = SimStructs.userStruct{jUser,1}.baseNode;
                                    if nCellIndex ~= iBase
                                        for jLayer = 1:maxRank
                                            intVector = [];
                                            for inUser = 1:kUsers
                                                intVector = [intVector, W0{jUser,iBand}(:,jLayer)' * cH{iBase,iBand}(:,:,jUser) * M(:,:,inUser,iBand)];
                                            end
                                            norm(intVector,2) <= I(:,jUser,iBand,iBase);
                                        end
                                    end
                                end
                                
                                P = real(W0{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * M(:,iLayer,iUser,iBand));
                                Q = imag(W0{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * M(:,iLayer,iUser,iBand));
                                P0 = real(W0{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * M0(:,iLayer,iUser,iBand));
                                Q0 = imag(W0{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * M0(:,iLayer,iUser,iBand));
                                
                                (P0^2 + Q0^2) / B0(iLayer,iUser,iBand) + (2 / B0(iLayer,iUser,iBand)) * (P0 * (P - P0) + Q0 * (Q - Q0)) ...
                                    - ((P0^2 + Q0^2) / (B0(iLayer,iUser,iBand)^2)) * (B(iLayer,iUser,iBand) - B0(iLayer,iUser,iBand)) >= G(iLayer,iUser,iBand);
                                
                            end
                        end
                    end
                    
                    vec(M)' * vec(M) <= sum(SimStructs.baseStruct{iBase,1}.sPower(1,:));
                    
                    cvx_end
                    
                    if strfind(cvx_status,'Solved')
                        M0 = full(M);
                    else
                        display(cvx_status);
                        for iBand = 1:nBands
                            M0 = SimParams.Debug.globalExchangeInfo.P{iBase,iBand} / sqrt(2);
                        end
                    end
                    
                    for iBand = 1:nBands
                        SimStructs.baseStruct{iBase,1}.P{iBand,1} = M0(:,:,:,iBand);
                        SimParams.Debug.globalExchangeInfo.P{iBase,iBand} = M0(:,:,:,iBand);
                    end
                    SimParams.Debug.globalExchangeInfo.I{iBase,1} = I;
                    
                    SimParams.Debug.globalExchangeInfo.funcOut{1,iBase} = M0;
                    SimParams.Debug.globalExchangeInfo.funcOut{2,iBase} = B0;
                    SimParams.Debug.globalExchangeInfo.funcOut{5,iBase} = W0;
                end
                
                tempTensor = zeros(maxRank,nUsers,nBands,nBases);
                for iBase = 1:nBases
                    tempTensor = tempTensor + SimParams.Debug.globalExchangeInfo.I{iBase,1};
                end
                for iBase = 1:nBases
                    SimParams.Debug.globalExchangeInfo.gI{iBase,1} = tempTensor(:,:,:,iBase) / 2;
                end
                
                for iBase = 1:nBases
                    for jBase = 1:nBases
                        if iBase ~= jBase
                            vecA = SimParams.Debug.globalExchangeInfo.gI{jBase,1}(:,cellUserIndices{iBase,1},:) - SimParams.Debug.globalExchangeInfo.I{iBase,1}(:,cellUserIndices{iBase,1},:,jBase);
                            SimParams.Debug.globalExchangeInfo.D{iBase,1}(:,cellUserIndices{iBase,1},:,jBase) = SimParams.Debug.globalExchangeInfo.D{iBase,1}(:,cellUserIndices{iBase,1},:,jBase) ...
                                + stepFactor * vecA;
                            
                            vecA = SimParams.Debug.globalExchangeInfo.gI{iBase,1}(:,cellUserIndices{jBase,1},:) - SimParams.Debug.globalExchangeInfo.I{iBase,1}(:,cellUserIndices{jBase,1},:,iBase);
                            SimParams.Debug.globalExchangeInfo.D{iBase,1}(:,cellUserIndices{jBase,1},:,iBase) = SimParams.Debug.globalExchangeInfo.D{iBase,1}(:,cellUserIndices{jBase,1},:,iBase) ...
                                + stepFactor * vecA;
                        end
                    end
                end
                
                [SimParams,SimStructs] = updateIteratePerformance(SimParams,SimStructs);
            end
            
            if SimParams.currentQueue < epsilonT
                break;
            end
            
        end
        
    case 'distMSEAllocA'
        
        stInstant = -1;
        stepIndex = 0.1;
        SimParams.currentQueue = 100;
        M0 = cell(nBases,1);E0 = cell(nBases,1);
        alphaLKN = cell(nBases,1);lambdaLKN = cell(nBases,1);
        
        for iExchangeOTA = stInstant:SimParams.nExchangesOTA
            
            switch iExchangeOTA
                case -1
                    for iBase = 1:nBases
                        if or((SimParams.distIteration - 1) == 0,mod((SimParams.iDrop - 1),SimParams.exchangeResetInterval) == 0)
                            fprintf('Resetting History for BS - %d \n',iBase);
                            SimStructs.baseStruct{iBase,1}.selectionType = 'BF_Prev';
                            [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE-BF_Prev',iBase);
                            SimParams.Debug.globalExchangeInfo.funcOut{3,iBase} = ones(maxRank,usersPerCell(iBase,1),nBands);
                            SimParams.Debug.globalExchangeInfo.funcOut{4,iBase} = ones(maxRank,usersPerCell(iBase,1),nBands);
                        else
                            if iBase == 1
                                fprintf('Reusing History \n');
                            end
                            SimStructs.baseStruct{iBase,1}.selectionType = 'Last_Prev';
                            [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'Last',iBase);
                        end
                    end
                    cH = SimStructs.prevChan;
                    maxBackHaulExchanges = SimParams.nExchangesOBH;
                case 0
                    if stInstant == 0
                        for iBase = 1:nBases
                            if or(SimParams.distIteration == 1,mod(SimParams.iDrop,SimParams.exchangeResetInterval) == 0)
                                fprintf('Resetting History for BS - %d \n',iBase);
                                SimStructs.baseStruct{iBase,1}.selectionType = 'BF';
                                [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE-BF',iBase);
                                SimParams.Debug.globalExchangeInfo.funcOut{3,iBase} = ones(maxRank,usersPerCell(iBase,1),nBands);
                                SimParams.Debug.globalExchangeInfo.funcOut{4,iBase} = ones(maxRank,usersPerCell(iBase,1),nBands);
                            else
                                if iBase == 1
                                    fprintf('Reusing History \n');
                                end
                                SimStructs.baseStruct{iBase,1}.selectionType = 'Last';
                                [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'Last',iBase);
                            end
                        end
                    else
                        for iBase = 1:nBases
                            SimStructs.baseStruct{iBase,1}.selectionType = 'Last';
                            [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'Last',iBase);
                        end
                    end
                    maxBackHaulExchanges = 1;
                    cH = SimStructs.linkChan;
                    fprintf('OTA Performed - %d \n',iExchangeOTA);
                otherwise
                    for iBase = 1:nBases
                        SimStructs.baseStruct{iBase,1}.selectionType = 'Last';
                    end
                    cH = SimStructs.linkChan;
                    maxBackHaulExchanges = 1;
                    fprintf('OTA Performed - %d \n',iExchangeOTA);
                    [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE');
            end
            
            [SimParams, SimStructs] = initializeSCApoint(SimParams,SimStructs,'MSE');
            
            for iBase = 1:nBases
                M0{iBase,1} = SimParams.Debug.globalExchangeInfo.funcOut{1,iBase};
                E0{iBase,1} = SimParams.Debug.globalExchangeInfo.funcOut{2,iBase};
                alphaLKN{iBase,1} = SimParams.Debug.globalExchangeInfo.funcOut{3,iBase};
                lambdaLKN{iBase,1} = SimParams.Debug.globalExchangeInfo.funcOut{4,iBase};
            end
            
            W = SimParams.Debug.globalExchangeInfo.funcOut{5,1};
            
            for iExchangeBH = 1:maxBackHaulExchanges
                
                if iExchangeBH ~= 1
                    for iBase = 1:nBases
                        SimStructs.baseStruct{iBase,1}.selectionType = 'Last';
                    end
                end
                
                for iBase = 1:nBases
                    
                    kUsers = usersPerCell(iBase,1);
                    for iBand = 1:nBands
                        for iUser = 1:kUsers
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for iLayer = 1:maxRank
                                intVector = sqrt(SimParams.N) * W{cUser,iBand}(:,iLayer)';
                                for jBase = 1:nBases
                                    P = M0{jBase,1};
                                    for jUser = 1:usersPerCell(jBase,1)
                                        jxUser = cellUserIndices{jBase,1}(jUser,1);
                                        currentH = cH{jBase,iBand}(:,:,cUser);
                                        if jxUser == cUser
                                            intVector = [intVector, W{cUser,iBand}(:,iLayer)' * currentH * P(:,iLayer~=rankArray,jUser,iBand)];
                                        else
                                            intVector = [intVector, W{cUser,iBand}(:,iLayer)' * currentH * P(:,:,jUser,iBand)];
                                        end
                                    end
                                end
                                
                                currentH = cH{iBase,iBand}(:,:,cUser);
                                givenVector = (1 - W{cUser,iBand}(:,iLayer)' * currentH * M0{iBase,1}(:,iLayer,iUser,iBand));
                                intVector = [intVector, givenVector];
                                E(iLayer,iUser,iBand) = norm(intVector,2)^2;
                            end
                        end
                    end
                    
                    for iBand = 1:nBands
                        for iUser = 1:kUsers
                            t(:,iUser,iBand) = -log2(E0{iBase,1}(:,iUser,iBand)) - (E(:,iUser,iBand) - E0{iBase,1}(:,iUser,iBand)) ./ (E0{iBase,1}(:,iUser,iBand) * log(2));
                        end
                    end
                    
                    for iBand = 1:nBands
                        for iUser = 1:kUsers
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for iRank = 1:maxRank
                                tempT = t(:,iUser,:);
                                lambdaLKN{iBase,1}(iRank,iUser,iBand) = qExponent * (QueuedPkts(cUser,1) - sum(tempT(:)))^(qExponent - 1) / log(2);
                                if QueuedPkts(cUser,1) < sum(tempT(:))
                                    if mod(qExponent,2) ~= 0
                                        lambdaLKN{iBase,1}(iRank,iUser,iBand) = 0;
                                    end
                                end
                                alphaLKN{iBase,1}(iRank,iUser,iBand) = alphaLKN{iBase,1}(iRank,iUser,iBand) + stepIndex * (max(lambdaLKN{iBase,1}(iRank,iUser,iBand) / E(iRank,iUser,iBand),0) - alphaLKN{iBase,1}(iRank,iUser,iBand));
                            end
                        end
                    end
                    
                    SimParams.Debug.globalExchangeInfo.funcOut{2,iBase} = E;
                    SimParams.Debug.globalExchangeInfo.funcOut{3,iBase} = alphaLKN{iBase,1};
                    SimParams.Debug.globalExchangeInfo.funcOut{4,iBase} = lambdaLKN{iBase,1};
                    
                end
                
                for iBase = 1:nBases
                    
                    kUsers = usersPerCell(iBase,1);
                    
                    for iBand = 1:nBands
                        for iUser = 1:kUsers
                            for iLayer = 1:maxRank
                                I = zeros(SimParams.nTxAntenna);
                                for jBase = 1:nBases
                                    for jUser = 1:usersPerCell(jBase,1)
                                        jxUser = cellUserIndices{jBase,1}(jUser,1);
                                        for jLayer = 1:maxRank
                                            I = I + cH{iBase,iBand}(:,:,jxUser)' * W{jxUser,iBand}(:,jLayer) * W{jxUser,iBand}(:,jLayer)' * cH{iBase,iBand}(:,:,jxUser) * SimParams.Debug.globalExchangeInfo.funcOut{3,jBase}(jLayer,jUser,iBand);
                                        end
                                    end
                                end
                                R(:,:,iLayer,iUser,iBand) = I;
                            end
                        end
                    end
                    
                    muMax = 100000;
                    muMin = 0;
                    iterateAgain = 1;
                    while iterateAgain
                        totalPower = 0;
                        currentMu = (muMax + muMin) / 2;
                        for iBand = 1:nBands
                            for iUser = 1:kUsers
                                cUser = cellUserIndices{iBase,1}(iUser,1);
                                for iLayer = 1:maxRank
                                    M(:,iLayer,iUser,iBand) = (currentMu * eye(SimParams.nTxAntenna) + R(:,:,iLayer,iUser,iBand)) \ (SimParams.Debug.globalExchangeInfo.funcOut{3,iBase}(iLayer,iUser,iBand) * cH{iBase,iBand}(:,:,cUser)' * W{cUser,iBand}(:,iLayer));
                                    totalPower = totalPower + real(trace(M(:,iLayer,iUser,iBand) * M(:,iLayer,iUser,iBand)'));
                                end
                            end
                        end
                        
                        if totalPower > sum(SimStructs.baseStruct{iBase,1}.sPower)
                            muMin = currentMu;
                        else
                            muMax = currentMu;
                        end
                        
                        if abs(muMin - muMax) <= 1e-6
                            iterateAgain = 0;
                        end
                    end
                    
                    for iBand = 1:nBands
                        SimStructs.baseStruct{iBase,1}.P{iBand,1} = M(:,:,:,iBand);
                        SimParams.Debug.globalExchangeInfo.P{iBase,iBand} = M(:,:,:,iBand);
                    end
                    
                    SimParams.Debug.globalExchangeInfo.funcOut{1,iBase} = M;
                    
                end
                
                
                [SimParams,SimStructs] = updateIteratePerformance(SimParams,SimStructs);
            end
            
            if SimParams.currentQueue < epsilonT
                break;
            end
            
        end
        
    case 'distMSEAllocB'
        
        stepIndex = 0.25;
        E0 = cell(nBases,1);
        SimParams.currentQueue = 100;
        alphaLKN = cell(nBases,1);lambdaLKN = cell(nBases,1);
        
        if SimParams.nExchangesOTA ~= 0
            SimParams.BITFactor = 1 - (SimParams.nExchangesOTA / SimParams.nSymbolsBIT);
        end
        
        for iExchangeOTA = 0:SimParams.nExchangesOTA
            
            switch iExchangeOTA
                case 0
                    for iBase = 1:nBases
                        if or((SimParams.distIteration - 1) == 0,mod((SimParams.iDrop - 1),SimParams.exchangeResetInterval) == 0)
                            fprintf('Resetting History for BS - %d \n',iBase);
                            SimStructs.baseStruct{iBase,1}.selectionType = 'BF';
                            [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE-BF',iBase);
                            SimParams.Debug.globalExchangeInfo.funcOut{3,iBase} = ones(maxRank,usersPerCell(iBase,1),nBands);
                            SimParams.Debug.globalExchangeInfo.funcOut{4,iBase} = ones(maxRank,usersPerCell(iBase,1),nBands);
                        else
                            if iBase == 1
                                fprintf('Reusing History \n');
                            end
                            SimStructs.baseStruct{iBase,1}.selectionType = 'Last';
                            [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'Last',iBase);
                        end
                    end
                    cH = SimStructs.linkChan;
                    fprintf('OTA Performed - %d \n',iExchangeOTA);
                otherwise
                    for iBase = 1:nBases
                        if iBase == 1
                            fprintf('OTA Performed - %d \n',iExchangeOTA);
                        end
                        SimStructs.baseStruct{iBase,1}.selectionType = 'Last';
                        [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE');
                    end
                    cH = SimStructs.linkChan;
            end
            
            [SimParams, SimStructs] = initializeSCApoint(SimParams,SimStructs,'MSE');
            
            for iBase = 1:nBases
                E0{iBase,1} = SimParams.Debug.globalExchangeInfo.funcOut{2,iBase};
                alphaLKN{iBase,1} = SimParams.Debug.globalExchangeInfo.funcOut{3,iBase};
                lambdaLKN{iBase,1} = SimParams.Debug.globalExchangeInfo.funcOut{4,iBase};
            end
            
            W = SimParams.Debug.globalExchangeInfo.funcOut{5,1};
            
            for iBase = 1:nBases
                
                kUsers = usersPerCell(iBase,1);
                
                for iBand = 1:nBands
                    for iUser = 1:kUsers
                        for iLayer = 1:maxRank
                            I = zeros(SimParams.nTxAntenna);
                            for jBase = 1:nBases
                                for jUser = 1:usersPerCell(jBase,1)
                                    jxUser = cellUserIndices{jBase,1}(jUser,1);
                                    for jLayer = 1:maxRank
                                        I = I + cH{iBase,iBand}(:,:,jxUser)' * W{jxUser,iBand}(:,jLayer) * W{jxUser,iBand}(:,jLayer)' * cH{iBase,iBand}(:,:,jxUser) * alphaLKN{jBase,1}(jLayer,jUser,iBand);
                                    end
                                end
                            end
                            R(:,:,iLayer,iUser,iBand) = I;
                        end
                    end
                end
                
                muMax = 100000;
                muMin = 0;
                iterateAgain = 1;
                while iterateAgain
                    totalPower = 0;
                    currentMu = (muMax + muMin) / 2;
                    for iBand = 1:nBands
                        for iUser = 1:kUsers
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for iLayer = 1:maxRank
                                M(:,iLayer,iUser,iBand) = (currentMu * eye(SimParams.nTxAntenna) + R(:,:,iLayer,iUser,iBand)) \ (alphaLKN{iBase,1}(iLayer,iUser,iBand) * cH{iBase,iBand}(:,:,cUser)' * W{cUser,iBand}(:,iLayer));
                                totalPower = totalPower + real(trace(M(:,iLayer,iUser,iBand) * M(:,iLayer,iUser,iBand)'));
                            end
                        end
                    end
                    
                    if totalPower > sum(SimStructs.baseStruct{iBase,1}.sPower)
                        muMin = currentMu;
                    else
                        muMax = currentMu;
                    end
                    
                    if abs(muMin - muMax) <= 1e-6
                        iterateAgain = 0;
                    end
                end
                
                for iBand = 1:nBands
                    SimStructs.baseStruct{iBase,1}.P{iBand,1} = M(:,:,:,iBand);
                    SimParams.Debug.globalExchangeInfo.P{iBase,iBand} = M(:,:,:,iBand);
                end
                
                SimParams.Debug.globalExchangeInfo.funcOut{1,iBase} = M;
                
            end
            
            [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE');
            for iBand = 1:nBands
                for iUser = 1:SimParams.nUsers
                    W{iUser,iBand} = SimStructs.userStruct{iUser,1}.pW{iBand,1};
                end
            end
            
            for iBase = 1:nBases
                
                P = SimParams.Debug.globalExchangeInfo.funcOut{1,iBase};
                for iBand = 1:nBands
                    for iUser = 1:kUsers
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        currentH = cH{iBase,iBand}(:,:,cUser);
                        for iLayer = 1:maxRank
                            E(iLayer,iUser,iBand) = (1 - W{cUser,iBand}(:,iLayer)' * currentH * P(:,iLayer,iUser,iBand));
                        end
                    end
                end
                
                for iBand = 1:nBands
                    for iUser = 1:kUsers
                        t(:,iUser,iBand) = -log2(E0{iBase,1}(:,iUser,iBand)) - (E(:,iUser,iBand) - E0{iBase,1}(:,iUser,iBand)) ./ (E0{iBase,1}(:,iUser,iBand) * log(2));
                    end
                end
                
                for iBand = 1:nBands
                    for iUser = 1:kUsers
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        for iRank = 1:maxRank
                            tempT = t(:,iUser,:) * SimParams.BITFactor;
                            lambdaLKN{iBase,1}(iRank,iUser,iBand) = qExponent * (QueuedPkts(cUser,1) - sum(tempT(:)))^(qExponent - 1) / log(2);
                            if QueuedPkts(cUser,1) < sum(tempT(:))
                                if mod(qExponent,2) ~= 0
                                    lambdaLKN{iBase,1}(iRank,iUser,iBand) = 0;
                                end
                            end
                            alphaLKN{iBase,1}(iRank,iUser,iBand) = alphaLKN{iBase,1}(iRank,iUser,iBand) + stepIndex * (max(lambdaLKN{iBase,1}(iRank,iUser,iBand) / E(iRank,iUser,iBand),0) - alphaLKN{iBase,1}(iRank,iUser,iBand));
                        end
                    end
                end
                
                SimParams.Debug.globalExchangeInfo.funcOut{2,iBase} = E;
                SimParams.Debug.globalExchangeInfo.funcOut{3,iBase} = alphaLKN{iBase,1};
                SimParams.Debug.globalExchangeInfo.funcOut{4,iBase} = lambdaLKN{iBase,1};
                
            end
            
            [SimParams,SimStructs] = updateIteratePerformance(SimParams,SimStructs);
            
            if SimParams.currentQueue < epsilonT
                break;
            end
            
        end
                
    case 'distMSEAllocC'
        
        stepIndex = 0.25;
        E0 = cell(nBases,1);
        SimParams.currentQueue = 100;
        alphaLKN = cell(nBases,1);lambdaLKN = cell(nBases,1);cW = cell(nBases,1);alphaLKN_O = cell(nBases,1);
        
        if SimParams.nExchangesOTA ~= 0
            SimParams.BITFactor = 1 - (SimParams.nExchangesOTA / SimParams.nSymbolsBIT);
        end
        
        for iExchangeOTA = 0:SimParams.nExchangesOTA
            
            switch iExchangeOTA
                case 0
                    for iBase = 1:nBases
                        if or((SimParams.distIteration - 1) == 0,mod((SimParams.iDrop - 1),SimParams.exchangeResetInterval) == 0)
                            fprintf('Resetting History for BS - %d \n',iBase);
                            SimStructs.baseStruct{iBase,1}.selectionType = 'BF';
                            [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE-BF',iBase);
                            SimParams.Debug.globalExchangeInfo.funcOut{3,iBase} = ones(maxRank,usersPerCell(iBase,1),nBands);
                            SimParams.Debug.globalExchangeInfo.funcOut{4,iBase} = ones(maxRank,usersPerCell(iBase,1),nBands);
                        else
                            if iBase == 1
                                fprintf('Reusing History \n');
                            end
                            SimStructs.baseStruct{iBase,1}.selectionType = 'FrameC';
                            [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'FrameC',iBase);
                        end
                    end
                    cH = SimStructs.linkChan;
                    maxBHIterations = max(20,SimParams.nExchangesOBH);
                    fprintf('OTA Performed - %d \n',iExchangeOTA);
                otherwise
                    for iBase = 1:nBases
                        if iBase == 1
                            fprintf('OTA Performed - %d \n',iExchangeOTA);
                        end
                        SimStructs.baseStruct{iBase,1}.selectionType = 'Last';
                        [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE');
                    end
                    cH = SimStructs.linkChan;
                    maxBHIterations = SimParams.nExchangesOBH;
            end
            
            [SimParams, SimStructs] = initializeSCApoint(SimParams,SimStructs,'MSE');
            W = SimParams.Debug.globalExchangeInfo.funcOut{5,1};
            
            for iBase = 1:nBases
                cW{iBase,1} = W;
                E0{iBase,1} = SimParams.Debug.globalExchangeInfo.funcOut{2,iBase};
                alphaLKN{iBase,1} = SimParams.Debug.globalExchangeInfo.funcOut{3,iBase};
                alphaLKN_O{iBase,1} = SimParams.Debug.globalExchangeInfo.funcOut{3,iBase};
                lambdaLKN{iBase,1} = SimParams.Debug.globalExchangeInfo.funcOut{4,iBase};
                
                
                switch (SimParams.additionalParams)                    
                    case 'Optimal'
                        if ~iExchangeOTA
                            for iBand = 1:nBands
                                for iUser = 1:SimParams.nUsers
                                    if isempty(find(iUser == cellUserIndices{iBase,1}))
                                        cW{iBase,1}{iUser,iBand} = zeros(size(W{iUser,iBand}));
                                    end
                                end
                            end
                        end
                    case 'E-Optimal'
                        maxBHIterations = iExchangeOTA + 1;
                        if maxBHIterations >= SimParams.nExchangesOTA
                            maxBHIterations = SimParams.nExchangesOBH;
                        end     
                end
                
            end
            
            for iBase = 1:nBases
                
                if strcmpi(SimParams.additionalParams,'S-Optimal')
                    switch (iExchangeOTA)
                        case {0,SimParams.nExchangesOTA}
                            maxBHIterations = SimParams.nExchangesOBH;
                        otherwise
                            if (mod((iExchangeOTA + iBase),nBases) == 0)
                                maxBHIterations = SimParams.nExchangesOBH;
                            else
                                maxBHIterations = 1;
                            end
                    end                    
                end
                
                kUsers = usersPerCell(iBase,1);
                for iExchangeOBH = 1:maxBHIterations
                    
                    for iBand = 1:nBands
                        for iUser = 1:kUsers
                            for iLayer = 1:maxRank
                                I = zeros(SimParams.nTxAntenna);
                                for jBase = 1:nBases
                                    if jBase ~= iBase
                                        for jUser = 1:usersPerCell(jBase,1)
                                            jxUser = cellUserIndices{jBase,1}(jUser,1);
                                            for jLayer = 1:maxRank
                                                I = I + cH{iBase,iBand}(:,:,jxUser)' * cW{iBase,1}{jxUser,iBand}(:,jLayer) * cW{iBase,1}{jxUser,iBand}(:,jLayer)' * cH{iBase,iBand}(:,:,jxUser) * alphaLKN_O{jBase,1}(jLayer,jUser,iBand);
                                            end
                                        end
                                    else
                                        for jUser = 1:usersPerCell(jBase,1)
                                            jxUser = cellUserIndices{jBase,1}(jUser,1);
                                            for jLayer = 1:maxRank
                                                I = I + cH{iBase,iBand}(:,:,jxUser)' * cW{iBase,1}{jxUser,iBand}(:,jLayer) * cW{iBase,1}{jxUser,iBand}(:,jLayer)' * cH{iBase,iBand}(:,:,jxUser) * alphaLKN_O{jBase,1}(jLayer,jUser,iBand);
                                            end
                                        end
                                    end
                                end
                                R(:,:,iLayer,iUser,iBand) = I;
                            end
                        end
                    end
                    
                    muMax = 100000;
                    muMin = 0;
                    iterateAgain = 1;
                    while iterateAgain
                        totalPower = 0;
                        currentMu = (muMax + muMin) / 2;
                        for iBand = 1:nBands
                            for iUser = 1:kUsers
                                cUser = cellUserIndices{iBase,1}(iUser,1);
                                for iLayer = 1:maxRank
                                    M(:,iLayer,iUser,iBand) = (currentMu * eye(SimParams.nTxAntenna) + R(:,:,iLayer,iUser,iBand)) \ (alphaLKN_O{iBase,1}(iLayer,iUser,iBand) * cH{iBase,iBand}(:,:,cUser)' * cW{iBase,1}{cUser,iBand}(:,iLayer));
                                    totalPower = totalPower + real(trace(M(:,iLayer,iUser,iBand) * M(:,iLayer,iUser,iBand)'));
                                end
                            end
                        end
                        
                        if totalPower > sum(SimStructs.baseStruct{iBase,1}.sPower)
                            muMin = currentMu;
                        else
                            muMax = currentMu;
                        end
                        
                        if abs(muMin - muMax) <= 1e-6
                            iterateAgain = 0;
                        end
                    end
                    
                    SimParams.Debug.ExchangeM{iBase,1} = M;
                    SimParams.Debug.ExchangeW = cW{iBase,1};
                    [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE-X',iBase);
                    cW{iBase,1} = SimParams.Debug.ExchangeW;
                    
                    for iBand = 1:nBands
                        for iUser = 1:kUsers
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            currentH = cH{iBase,iBand}(:,:,cUser);
                            for iLayer = 1:maxRank
                                E(iLayer,iUser,iBand) = (1 - cW{iBase,1}{cUser,iBand}(:,iLayer)' * currentH * M(:,iLayer,iUser,iBand));
                            end
                        end
                    end
                    
                    for iBand = 1:nBands
                        for iUser = 1:kUsers
                            t(:,iUser,iBand) = -log2(E0{iBase,1}(:,iUser,iBand)) - (E(:,iUser,iBand) - E0{iBase,1}(:,iUser,iBand)) ./ (E0{iBase,1}(:,iUser,iBand) * log(2));
                        end
                    end
                    
                    for iBand = 1:nBands
                        for iUser = 1:kUsers
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for iRank = 1:maxRank
                                tempT = t(:,iUser,:) * SimParams.BITFactor;
                                lambdaLKN{iBase,1}(iRank,iUser,iBand) = qExponent * (QueuedPkts(cUser,1) - sum(tempT(:)))^(qExponent - 1) / log(2);
                                if QueuedPkts(cUser,1) < sum(tempT(:))
                                    if mod(qExponent,2) ~= 0
                                        lambdaLKN{iBase,1}(iRank,iUser,iBand) = 0;
                                    end
                                end
                                alphaLKN_O{iBase,1}(iRank,iUser,iBand) = alphaLKN_O{iBase,1}(iRank,iUser,iBand) + stepIndex * (max(lambdaLKN{iBase,1}(iRank,iUser,iBand) / E(iRank,iUser,iBand),0) - alphaLKN_O{iBase,1}(iRank,iUser,iBand));
                            end
                        end
                    end
                    
                end
                
            end
            
            for iBase = 1:nBases
                for iBand = 1:nBands
                    SimStructs.baseStruct{iBase,1}.P{iBand,1} = SimParams.Debug.ExchangeM{iBase,1}(:,:,:,iBand);
                    SimParams.Debug.globalExchangeInfo.P{iBase,iBand} = SimParams.Debug.ExchangeM{iBase,1}(:,:,:,iBand);
                end
                
                alphaLKN{iBase,1} = alphaLKN_O{iBase,1};
                SimParams.Debug.globalExchangeInfo.funcOut{1,iBase} = SimParams.Debug.ExchangeM{iBase,1};
            end
            
            [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE');
            
            for iBand = 1:nBands
                for iUser = 1:SimParams.nUsers
                    W{iUser,iBand} = SimStructs.userStruct{iUser,1}.pW{iBand,1};
                end
            end
            
            for iBase = 1:nBases
                
                P = SimParams.Debug.globalExchangeInfo.funcOut{1,iBase};
                for iBand = 1:nBands
                    for iUser = 1:kUsers
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        currentH = cH{iBase,iBand}(:,:,cUser);
                        for iLayer = 1:maxRank
                            E(iLayer,iUser,iBand) = (1 - W{cUser,iBand}(:,iLayer)' * currentH * P(:,iLayer,iUser,iBand));
                        end
                    end
                end
                
                for iBand = 1:nBands
                    for iUser = 1:kUsers
                        t(:,iUser,iBand) = -log2(E0{iBase,1}(:,iUser,iBand)) - (E(:,iUser,iBand) - E0{iBase,1}(:,iUser,iBand)) ./ (E0{iBase,1}(:,iUser,iBand) * log(2));
                    end
                end
                
                for iBand = 1:nBands
                    for iUser = 1:kUsers
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        for iRank = 1:maxRank
                            tempT = t(:,iUser,:) * SimParams.BITFactor;
                            lambdaLKN{iBase,1}(iRank,iUser,iBand) = qExponent * (QueuedPkts(cUser,1) - sum(tempT(:)))^(qExponent - 1) / log(2);
                            if QueuedPkts(cUser,1) < sum(tempT(:))
                                if mod(qExponent,2) ~= 0
                                    lambdaLKN{iBase,1}(iRank,iUser,iBand) = 0;
                                end
                            end
                            alphaLKN{iBase,1}(iRank,iUser,iBand) = alphaLKN{iBase,1}(iRank,iUser,iBand) + stepIndex * (max(lambdaLKN{iBase,1}(iRank,iUser,iBand) / E(iRank,iUser,iBand),0) - alphaLKN{iBase,1}(iRank,iUser,iBand));
                        end
                    end
                end
                
                SimParams.Debug.globalExchangeInfo.funcOut{2,iBase} = E;
                SimParams.Debug.globalExchangeInfo.funcOut{3,iBase} = alphaLKN{iBase,1};
                SimParams.Debug.globalExchangeInfo.funcOut{4,iBase} = lambdaLKN{iBase,1};
                
            end
            
            [SimParams,SimStructs] = updateIteratePerformance(SimParams,SimStructs);
            
            if SimParams.currentQueue < epsilonT
                break;
            end
            
        end
        
    case 'distMSEAllocD'
        
        stepIndex = 0.25;
        E0 = cell(nBases,1);
        SimParams.currentQueue = 100;
        alphaLKN = cell(nBases,1);lambdaLKN = cell(nBases,1);cW = cell(nBases,1);alphaLKN_O = cell(nBases,1);
        
        if SimParams.nExchangesOTA ~= 0
            SimParams.BITFactor = 1 - (SimParams.nExchangesOTA / SimParams.nSymbolsBIT);
        end
        
        userGroupingScript;
        
        for iExchangeOTA = 0:SimParams.nExchangesOTA
            
            switch iExchangeOTA
                case 0
                    for iBase = 1:nBases
                        if or((SimParams.distIteration - 1) == 0,mod((SimParams.iDrop - 1),SimParams.exchangeResetInterval) == 0)
                            fprintf('Resetting History for BS - %d \n',iBase);
                            SimStructs.baseStruct{iBase,1}.selectionType = 'BF';
                            [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE-BF',iBase);
                            SimParams.Debug.globalExchangeInfo.funcOut{3,iBase} = ones(maxRank,usersPerCell(iBase,1),nBands);
                            SimParams.Debug.globalExchangeInfo.funcOut{4,iBase} = ones(maxRank,usersPerCell(iBase,1),nBands);
                        else
                            if iBase == 1
                                fprintf('Reusing History \n');
                            end
                            SimStructs.baseStruct{iBase,1}.selectionType = 'FrameC';
                            [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'FrameC',iBase);
                        end
                    end
                    cH = SimStructs.linkChan;
                    maxBHIterations = max(20,SimParams.nExchangesOBH);
                    fprintf('OTA Performed - %d \n',iExchangeOTA);
                otherwise
                    for iBase = 1:nBases
                        if iBase == 1
                            fprintf('OTA Performed - %d \n',iExchangeOTA);
                        end
                        SimStructs.baseStruct{iBase,1}.selectionType = 'Last';
                        [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE');
                    end
                    cH = SimStructs.linkChan;
                    maxBHIterations = SimParams.nExchangesOBH;
            end
            
            [SimParams, SimStructs] = initializeSCApoint(SimParams,SimStructs,'MSE');
            W = SimParams.Debug.globalExchangeInfo.funcOut{5,1};
            
            for iBase = 1:nBases
                cW{iBase,1} = W;
                E0{iBase,1} = SimParams.Debug.globalExchangeInfo.funcOut{2,iBase};
                alphaLKN{iBase,1} = SimParams.Debug.globalExchangeInfo.funcOut{3,iBase};
                alphaLKN_O{iBase,1} = SimParams.Debug.globalExchangeInfo.funcOut{3,iBase};
                lambdaLKN{iBase,1} = SimParams.Debug.globalExchangeInfo.funcOut{4,iBase};
                
                
                switch (SimParams.additionalParams)                    
                    case 'Optimal'
                        if ~iExchangeOTA
                            for iBand = 1:nBands
                                for iUser = 1:SimParams.nUsers
                                    if isempty(find(iUser == cellUserIndices{iBase,1}))
                                        cW{iBase,1}{iUser,iBand} = zeros(size(W{iUser,iBand}));
                                    end
                                end
                            end
                        end
                    case 'E-Optimal'
                        maxBHIterations = iExchangeOTA + 1;
                        if maxBHIterations >= SimParams.nExchangesOTA
                            maxBHIterations = SimParams.nExchangesOBH;
                        end     
                end
                
            end
            
            for iBase = 1:nBases
                
                if strcmpi(SimParams.additionalParams,'S-Optimal')
                    switch (iExchangeOTA)
                        case {0,SimParams.nExchangesOTA}
                            maxBHIterations = SimParams.nExchangesOBH;
                        otherwise
                            if (mod((iExchangeOTA + iBase),nBases) == 0)
                                maxBHIterations = SimParams.nExchangesOBH;
                            else
                                maxBHIterations = 1;
                            end
                    end                    
                end
                
                kUsers = usersPerCell(iBase,1);
                for iExchangeOBH = 1:maxBHIterations
                    
                    for iBand = 1:nBands
                        for iUser = 1:kUsers
                            for iLayer = 1:maxRank
                                I = zeros(SimParams.nTxAntenna);
                                for iNeighbor = 1:length(ifUsersOfEachCell{iBase,1})
                                    xNeighbor = ifUsersOfEachCell{iBase,1}(1,iNeighbor);
                                    for jLayer = 1:maxRank
                                        I = I + cH{iBase,iBand}(:,:,xNeighbor)' * cW{iBase,1}{xNeighbor,iBand}(:,jLayer) * cW{iBase,1}{xNeighbor,iBand}(:,jLayer)' * cH{iBase,iBand}(:,:,xNeighbor) * alphaLKN_O{ifUsersOfEachCell{iBase,3}(1,iNeighbor),1}(jLayer,ifUsersOfEachCell{iBase,2}(1,iNeighbor),iBand);
                                    end
                                    
                                end
                                
                                for jUser = 1:usersPerCell(iBase,1)
                                    jxUser = cellUserIndices{iBase,1}(jUser,1);
                                    for jLayer = 1:maxRank
                                        I = I + cH{iBase,iBand}(:,:,jxUser)' * cW{iBase,1}{jxUser,iBand}(:,jLayer) * cW{iBase,1}{jxUser,iBand}(:,jLayer)' * cH{iBase,iBand}(:,:,jxUser) * alphaLKN_O{iBase,1}(jLayer,jUser,iBand);
                                    end
                                end
                                R(:,:,iLayer,iUser,iBand) = I;
                            end
                        end
                    end
                    
                    muMax = 1e10;
                    muMin = 0;
                    iterateAgain = 1;
                    while iterateAgain
                        totalPower = 0;
                        currentMu = (muMax + muMin) / 2;
                        for iBand = 1:nBands
                            for iUser = 1:kUsers
                                cUser = cellUserIndices{iBase,1}(iUser,1);
                                for iLayer = 1:maxRank
                                    M(:,iLayer,iUser,iBand) = (currentMu * eye(SimParams.nTxAntenna) + R(:,:,iLayer,iUser,iBand)) \ (alphaLKN_O{iBase,1}(iLayer,iUser,iBand) * cH{iBase,iBand}(:,:,cUser)' * cW{iBase,1}{cUser,iBand}(:,iLayer));
                                    totalPower = totalPower + real(trace(M(:,iLayer,iUser,iBand) * M(:,iLayer,iUser,iBand)'));
                                end
                            end
                        end
                        
                        if totalPower > sum(SimStructs.baseStruct{iBase,1}.sPower)
                            muMin = currentMu;
                        else
                            muMax = currentMu;
                        end
                        
                        if abs(muMin - muMax) <= 1e-6
                            iterateAgain = 0;
                        end
                    end
                    
                    SimParams.Debug.ExchangeM{iBase,1} = M;
                    SimParams.Debug.ExchangeW = cW{iBase,1};
                    [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE-X',iBase);
                    cW{iBase,1} = SimParams.Debug.ExchangeW;
                    
                    for iBand = 1:nBands
                        for iUser = 1:kUsers
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            currentH = cH{iBase,iBand}(:,:,cUser);
                            for iLayer = 1:maxRank
                                E(iLayer,iUser,iBand) = (1 - cW{iBase,1}{cUser,iBand}(:,iLayer)' * currentH * M(:,iLayer,iUser,iBand));
                            end
                        end
                    end
                    
                    for iBand = 1:nBands
                        for iUser = 1:kUsers
                            t(:,iUser,iBand) = -log2(E0{iBase,1}(:,iUser,iBand)) - (E(:,iUser,iBand) - E0{iBase,1}(:,iUser,iBand)) ./ (E0{iBase,1}(:,iUser,iBand) * log(2));
                        end
                    end
                    
                    for iBand = 1:nBands
                        for iUser = 1:kUsers
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for iRank = 1:maxRank
                                tempT = t(:,iUser,:) * SimParams.BITFactor;
                                lambdaLKN{iBase,1}(iRank,iUser,iBand) = qExponent * (QueuedPkts(cUser,1) - sum(tempT(:)))^(qExponent - 1) / log(2);
                                if QueuedPkts(cUser,1) < sum(tempT(:))
                                    if mod(qExponent,2) ~= 0
                                        lambdaLKN{iBase,1}(iRank,iUser,iBand) = 0;
                                    end
                                end
                                alphaLKN_O{iBase,1}(iRank,iUser,iBand) = alphaLKN_O{iBase,1}(iRank,iUser,iBand) + stepIndex * (max(lambdaLKN{iBase,1}(iRank,iUser,iBand) / E(iRank,iUser,iBand),0) - alphaLKN_O{iBase,1}(iRank,iUser,iBand));
                            end
                        end
                    end
                    
                end
                
            end
            
            for iBase = 1:nBases
                for iBand = 1:nBands
                    SimStructs.baseStruct{iBase,1}.P{iBand,1} = SimParams.Debug.ExchangeM{iBase,1}(:,:,:,iBand);
                    SimParams.Debug.globalExchangeInfo.P{iBase,iBand} = SimParams.Debug.ExchangeM{iBase,1}(:,:,:,iBand);
                end
                
                alphaLKN{iBase,1} = alphaLKN_O{iBase,1};
                SimParams.Debug.globalExchangeInfo.funcOut{1,iBase} = SimParams.Debug.ExchangeM{iBase,1};
            end
            
            [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE');
            
            for iBand = 1:nBands
                for iUser = 1:SimParams.nUsers
                    W{iUser,iBand} = SimStructs.userStruct{iUser,1}.pW{iBand,1};
                end
            end
            
            for iBase = 1:nBases
                
                P = SimParams.Debug.globalExchangeInfo.funcOut{1,iBase};
                for iBand = 1:nBands
                    for iUser = 1:kUsers
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        currentH = cH{iBase,iBand}(:,:,cUser);
                        for iLayer = 1:maxRank
                            E(iLayer,iUser,iBand) = (1 - W{cUser,iBand}(:,iLayer)' * currentH * P(:,iLayer,iUser,iBand));
                        end
                    end
                end
                
                for iBand = 1:nBands
                    for iUser = 1:kUsers
                        t(:,iUser,iBand) = -log2(E0{iBase,1}(:,iUser,iBand)) - (E(:,iUser,iBand) - E0{iBase,1}(:,iUser,iBand)) ./ (E0{iBase,1}(:,iUser,iBand) * log(2));
                    end
                end
                
                for iBand = 1:nBands
                    for iUser = 1:kUsers
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        for iRank = 1:maxRank
                            tempT = t(:,iUser,:) * SimParams.BITFactor;
                            lambdaLKN{iBase,1}(iRank,iUser,iBand) = qExponent * (QueuedPkts(cUser,1) - sum(tempT(:)))^(qExponent - 1) / log(2);
                            if QueuedPkts(cUser,1) < sum(tempT(:))
                                if mod(qExponent,2) ~= 0
                                    lambdaLKN{iBase,1}(iRank,iUser,iBand) = 0;
                                end
                            end
                            alphaLKN{iBase,1}(iRank,iUser,iBand) = alphaLKN{iBase,1}(iRank,iUser,iBand) + stepIndex * (max(lambdaLKN{iBase,1}(iRank,iUser,iBand) / E(iRank,iUser,iBand),0) - alphaLKN{iBase,1}(iRank,iUser,iBand));
                        end
                    end
                end
                
                SimParams.Debug.globalExchangeInfo.funcOut{2,iBase} = E;
                SimParams.Debug.globalExchangeInfo.funcOut{3,iBase} = alphaLKN{iBase,1};
                SimParams.Debug.globalExchangeInfo.funcOut{4,iBase} = lambdaLKN{iBase,1};
                
            end
            
            [SimParams,SimStructs] = updateIteratePerformance(SimParams,SimStructs);
            
            if SimParams.currentQueue < epsilonT
                break;
            end
            
        end
        
end

% Update the receivers upon reception by all users (combined detection)
[SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE');


