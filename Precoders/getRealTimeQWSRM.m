function [SimParams,SimStructs] = getRealTimeQWSRM(SimParams,SimStructs)

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
        
        for iExchangeOTA = 1:1
            
            for iBase = 1:nBases
                
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
                    
                    vec(M)' * vec(M) <= sum(SimStructs.baseStruct{iBase,1}.sPower(1,:));
                    
                    cvx_end
                    
                    if strfind(cvx_status,'Solved')
                        M0 = full(M);
                    else
                        display(cvx_status);
                        break;
                    end
                    
                    fprintf('[%d,%f] \t',iBase,cvx_optval);
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
                
                fprintf('\n');
            end
            
            [SimParams,SimStructs] = updateIteratePerformance(SimParams,SimStructs);
        end
        
    case 'distBSAlloc'
        
        stInstant = -1;
        stepIndex = 10;
        for iExchangeOTA = stInstant:SimParams.nExchangesOTA
            
            switch iExchangeOTA
                
                case -1
                    for iBase = 1:nBases
                        resetPoint = (iBase - 1) * floor(SimParams.exchangeResetInterval / nBases);
                        if or(SimParams.distIteration == 1,mod(SimParams.iDrop,SimParams.exchangeResetInterval) == resetPoint)
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
                    cH = SimStructs.prevChan;
                    maxBackHaulExchanges = SimParams.nExchangesOBH;
                case 0
                    if stInstant == 0
                        for iBase = 1:nBases
                            resetPoint = (iBase - 1) * floor(SimParams.exchangeResetInterval / nBases);
                            if or(SimParams.distIteration == 1,mod(SimParams.iDrop,SimParams.exchangeResetInterval) == resetPoint)
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
                    maxBackHaulExchanges = 1;
                    fprintf('OTA Performed - %d \n',iExchangeOTA);
                otherwise
                    cH = SimStructs.linkChan;
                    maxBackHaulExchanges = 1;
                    fprintf('OTA Performed - %d \n',iExchangeOTA);
                    [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE');
            end
            
            for iExchangeBH = 1:maxBackHaulExchanges
                
                for iBase = 1:nBases
                    
                    if iExchangeBH ~= 1
                        SimStructs.baseStruct{iBase,1}.selectionType = 'Last';
                    end
                    
                    kUsers = usersPerCell(iBase,1);
                    [SimParams, SimStructs] = initializeSCApoint(SimParams,SimStructs,iBase);
                    M0 = SimParams.Debug.globalExchangeInfo.funcOut{1,iBase};B0 = SimParams.Debug.globalExchangeInfo.funcOut{2,iBase};
                    W0 = SimParams.Debug.globalExchangeInfo.funcOut{5,iBase};
                    
                    cvx_begin
                    
                    variable M(SimParams.nTxAntenna,maxRank,kUsers,nBands) complex
                    variables T(maxRank,kUsers,nBands) B(maxRank,kUsers,nBands) G(maxRank,kUsers,nBands)
                    variables I(maxRank,nUsers,nBands,nBases) userObjective(kUsers,1) epiObjective
                    
                    for iUser = 1:kUsers
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        userWts(cUser,1) * abs(QueuedPkts(cUser,1) - sum(vec(T(:,iUser,:)))) <= userObjective(iUser,1);
                    end
                    
                    augmentedTerms = 0;
                    for jBand = 1:nBands
                        for jUser = 1:nUsers
                            nCellIndex = SimStructs.userStruct{jUser,1}.baseNode;
                            if nCellIndex ~= iBase
                                augmentedTerms = augmentedTerms + sum((SimParams.Debug.globalExchangeInfo.gI{iBase,1}(:,jUser,jBand) - I(:,jUser,jBand,iBase)) .* SimParams.Debug.globalExchangeInfo.D{iBase,1}(:,jUser,jBand,iBase)) ...
                                    + (stepIndex / 2) * sum(vec(SimParams.Debug.globalExchangeInfo.gI{iBase,1}(:,jUser,jBand) - I(:,jUser,jBand,iBase)).^2);
                            end
                        end
                        
                        for jUser = 1:kUsers
                            cUser = cellUserIndices{iBase,1}(jUser,1);
                            for jBase = 1:nBases
                                if jBase ~= iBase
                                    augmentedTerms = augmentedTerms + sum((SimParams.Debug.globalExchangeInfo.gI{jBase,1}(:,cUser,jBand) - I(:,cUser,jBand,jBase)) .* SimParams.Debug.globalExchangeInfo.D{iBase,1}(:,cUser,jBand,jBase)) ...
                                        + (stepIndex / 2) * sum(vec(SimParams.Debug.globalExchangeInfo.gI{jBase,1}(:,cUser,jBand) - I(:,cUser,jBand,jBase)).^2);
                                end
                            end
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
                        break;
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
                
                for jBand = 1:nBands
                    for iUser = 1:nUsers
                        dCell = SimStructs.userStruct{iUser,1}.baseNode;
                        for jBase = 1:nBases
                            if jBase ~= dCell
                                SimParams.Debug.globalExchangeInfo.gI{jBase,1}(:,iUser,jBand) = 0.5 * (SimParams.Debug.globalExchangeInfo.I{jBase,1}(:,iUser,jBand,jBase) + SimParams.Debug.globalExchangeInfo.I{dCell,1}(:,iUser,jBand,jBase));
                                SimParams.Debug.globalExchangeInfo.D{dCell,1}(:,iUser,jBand,jBase) = (SimParams.Debug.globalExchangeInfo.gI{jBase,1}(:,iUser,jBand) - SimParams.Debug.globalExchangeInfo.I{dCell,1}(:,iUser,jBand,jBase)) * stepIndex + SimParams.Debug.globalExchangeInfo.D{dCell,1}(:,iUser,jBand,jBase);
                                SimParams.Debug.globalExchangeInfo.D{jBase,1}(:,iUser,jBand,jBase) = (SimParams.Debug.globalExchangeInfo.gI{jBase,1}(:,iUser,jBand) - SimParams.Debug.globalExchangeInfo.I{jBase,1}(:,iUser,jBand,jBase)) * stepIndex + SimParams.Debug.globalExchangeInfo.D{jBase,1}(:,iUser,jBand,jBase);
                            end
                        end
                    end
                end
                
                [SimParams,SimStructs] = updateIteratePerformance(SimParams,SimStructs);
            end
        end
        
    case 'distMSEAlloc'
        
        stInstant = 0;
        stepIndex = 0.5;
        M0 = cell(nBases,1);E0 = cell(nBases,1);
        alphaLKN = cell(nBases,1);lambdaLKN = cell(nBases,1);
        
        for iExchangeOTA = stInstant:SimParams.nExchangesOTA
            
            switch iExchangeOTA
                case -1
                    for iBase = 1:nBases
                        resetPoint = (iBase - 1) * floor(SimParams.exchangeResetInterval / nBases);
                        if or(SimParams.distIteration == 1,mod(SimParams.iDrop,SimParams.exchangeResetInterval) == resetPoint)
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
                    cH = SimStructs.prevChan;
                    maxBackHaulExchanges = SimParams.nExchangesOBH;
                case 0
                    if stInstant == 0
                        for iBase = 1:nBases
                            resetPoint = (iBase - 1) * floor(SimParams.exchangeResetInterval / nBases);
                            if or(SimParams.distIteration == 1,mod(SimParams.iDrop,SimParams.exchangeResetInterval) == resetPoint)
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
                        if iBase == 1
                            fprintf('OTA Performed - %d \n',iExchangeOTA);
                        end
                        [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE');
                    end
                    cH = SimStructs.linkChan;
                    maxBackHaulExchanges = 1;
            end
            
            for iExchangeBH = 1:maxBackHaulExchanges
                
                if iExchangeBH ~= 1
                    for iBase = 1:nBases
                        SimStructs.baseStruct{iBase,1}.selectionType = 'Last';
                    end
                end
                
                [SimParams, SimStructs] = initializeSCApoint(SimParams,SimStructs,'MSE');
                
                for iBase = 1:nBases
                    M0{iBase,1} = SimParams.Debug.globalExchangeInfo.funcOut{1,iBase};
                    if iExchangeBH == 1
                        E0{iBase,1} = SimParams.Debug.globalExchangeInfo.funcOut{2,iBase};
                    end
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
                
                if iExchangeBH == 1                                   
                    [SimParams, SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE');                    
                    for iBand = 1:nBands
                        for iUser = 1:SimParams.nUsers
                            W{iUser,iBand} = SimStructs.userStruct{iUser,1}.pW{iBand,1};
                        end
                    end
                end
        
                for iBase = 1:nBases
                    
                    for iBand = 1:nBands
                        for iUser = 1:kUsers
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for iLayer = 1:maxRank
                                intVector = sqrt(SimParams.N) * W{cUser,iBand}(:,iLayer)';
                                for jBase = 1:nBases
                                    P = SimParams.Debug.globalExchangeInfo.P{iBase,iBand};
                                    for jUser = 1:usersPerCell(jBase,1)
                                        jxUser = cellUserIndices{jBase,1}(jUser,1);
                                        currentH = cH{jBase,iBand}(:,:,cUser);
                                        if jxUser == cUser
                                            intVector = [intVector, W{cUser,iBand}(:,iLayer)' * currentH * P(:,iLayer~=rankArray,jUser)];
                                        else
                                            intVector = [intVector, W{cUser,iBand}(:,iLayer)' * currentH * P(:,:,jUser)];
                                        end
                                    end
                                end
                                
                                currentH = cH{iBase,iBand}(:,:,cUser);
                                givenVector = (1 - W{cUser,iBand}(:,iLayer)' * currentH * P(:,iLayer,iUser));
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
                                if mod(qExponent,2) ~= 0
                                    if QueuedPkts(cUser,1) < sum(tempT(:))
                                        lambdaLKN{iBase,1}(iRank,iUser,iBand) = 0;
                                    end
                                end
                                alphaLKN{iBase,1}(iRank,iUser,iBand) = alphaLKN{iBase,1}(iRank,iUser,iBand) + 1 * (max(lambdaLKN{iBase,1}(iRank,iUser,iBand) / E(iRank,iUser,iBand),0) - alphaLKN{iBase,1}(iRank,iUser,iBand));
                            end
                        end
                    end
                    
                    SimParams.Debug.globalExchangeInfo.funcOut{2,iBase} = E;
                    SimParams.Debug.globalExchangeInfo.funcOut{3,iBase} = alphaLKN{iBase,1};
                    SimParams.Debug.globalExchangeInfo.funcOut{4,iBase} = lambdaLKN{iBase,1};
                    
                end
                
                [SimParams,SimStructs] = updateIteratePerformance(SimParams,SimStructs);
            end
            
        end
        
        
    case 'distYBSAlloc'
        
        stepIndex = 10;
        for iExchangeOTA = -1:SimParams.nExchangesOTA
            
            switch iExchangeOTA
                
                case -1
                    for iBase = 1:nBases
                        resetPoint = (iBase - 1) * floor(SimParams.exchangeResetInterval / nBases);
                        if or(SimParams.distIteration == 1,mod(SimParams.iDrop,SimParams.exchangeResetInterval) == resetPoint)
                            fprintf('Resetting History for BS - %d \n',iBase);
                            SimStructs.baseStruct{iBase,1}.selectionType = 'BF';
                            [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'BF',iBase);
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
                    cH = SimStructs.prevChan;
                    maxBackHaulExchanges = SimParams.nExchangesOBH;
                case 0
                    for iBase = 1:nBases
                        SimStructs.baseStruct{iBase,1}.selectionType = 'Last';
                        [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'Last',iBase);
                    end
                    cH = SimStructs.linkChan;
                    maxBackHaulExchanges = 1;
                    fprintf('OTA Performed - %d \n',iExchangeOTA);
                otherwise
                    cH = SimStructs.linkChan;
                    maxBackHaulExchanges = 1;
                    fprintf('OTA Performed - %d \n',iExchangeOTA);
                    [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE');
            end
            
            for iExchangeBH = 1:maxBackHaulExchanges
                
                for iBase = 1:nBases
                    
                    if iExchangeBH ~= 1
                        SimStructs.baseStruct{iBase,1}.selectionType = 'Last';
                    end
                    
                    constVector = [];
                    kUsers = usersPerCell(iBase,1);
                    [SimParams, SimStructs] = initializeSCApoint(SimParams,SimStructs,iBase);
                    M0 = SimParams.Debug.globalExchangeInfo.funcOut{1,iBase};B0 = SimParams.Debug.globalExchangeInfo.funcOut{2,iBase};
                    T0 = SimParams.Debug.globalExchangeInfo.funcOut{4,iBase};W0 = SimParams.Debug.globalExchangeInfo.funcOut{5,iBase};
                    
                    M = sdpvar(SimParams.nTxAntenna,maxRank,kUsers,nBands,'full','complex');
                    T = sdpvar(maxRank,kUsers,nBands,'full');B = sdpvar(maxRank,kUsers,nBands,'full');G = sdpvar(maxRank,kUsers,nBands,'full');
                    I = sdpvar(maxRank,nUsers,nBands,nBases,'full');userObjective = sdpvar(kUsers,1);epiObjective = sdpvar(1,1);
                    
                    for iUser = 1:kUsers
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        constVector = [constVector, userWts(cUser,1) * abs(QueuedPkts(cUser,1) - sum(vector(T(:,iUser,:)))) <= userObjective(iUser,1)];
                    end
                    
                    augmentedTerms = 0;
                    for jBand = 1:nBands
                        for jUser = 1:nUsers
                            nCellIndex = SimStructs.userStruct{jUser,1}.baseNode;
                            if nCellIndex ~= iBase
                                augmentedTerms = augmentedTerms + sum((SimParams.Debug.globalExchangeInfo.gI{iBase,1}(:,jUser,jBand) - I(:,jUser,jBand,iBase)) .* SimParams.Debug.globalExchangeInfo.D{iBase,1}(:,jUser,jBand,iBase)) ...
                                    + (stepIndex / 2) * sum(vector(SimParams.Debug.globalExchangeInfo.gI{iBase,1}(:,jUser,jBand) - I(:,jUser,jBand,iBase)).^2);
                            end
                        end
                        
                        for jUser = 1:kUsers
                            cUser = cellUserIndices{iBase,1}(jUser,1);
                            for jBase = 1:nBases
                                if jBase ~= iBase
                                    augmentedTerms = augmentedTerms + sum((SimParams.Debug.globalExchangeInfo.gI{jBase,1}(:,cUser,jBand) - I(:,cUser,jBand,jBase)) .* SimParams.Debug.globalExchangeInfo.D{iBase,1}(:,cUser,jBand,jBase)) ...
                                        + (stepIndex / 2) * sum(vector(SimParams.Debug.globalExchangeInfo.gI{jBase,1}(:,cUser,jBand) - I(:,cUser,jBand,jBase)).^2);
                                end
                            end
                        end
                    end
                    
                    constVector = [constVector, epiObjective >= norm(userObjective,qExponent) + augmentedTerms];
                    
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
                                
                                constVector = [constVector, norm(intVector) <= sqrt(B(iLayer,iUser,iBand))];
                                
                                for jUser = 1:nUsers
                                    nCellIndex = SimStructs.userStruct{jUser,1}.baseNode;
                                    if nCellIndex ~= iBase
                                        for jLayer = 1:maxRank
                                            intVector = [];
                                            for inUser = 1:kUsers
                                                intVector = [intVector, W0{jUser,iBand}(:,jLayer)' * cH{iBase,iBand}(:,:,jUser) * M(:,:,inUser,iBand)];
                                            end
                                            constVector = [constVector, norm(intVector,2) <= I(:,jUser,iBand,iBase)];
                                        end
                                    end
                                end
                                
                                %                                 constVector = [constVector, T(iLayer,iUser,iBand) >= T0(iLayer,iUser,iBand) * (1 - epsilon)];
                                %                                 constVector = [constVector, T(iLayer,iUser,iBand) <= T0(iLayer,iUser,iBand) * (1 + epsilon)];
                                constVector = [constVector, 1 + G(iLayer,iUser,iBand) >= ...
                                    exp(T0(iLayer,iUser,iBand) * log(2)) * (1 + log(2) * (T(iLayer,iUser,iBand) - T0(iLayer,iUser,iBand)) + 0.5 * log(2)^2 * (T(iLayer,iUser,iBand) - T0(iLayer,iUser,iBand))^2)];
                                
                                P = real(W0{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * M(:,iLayer,iUser,iBand));
                                Q = imag(W0{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * M(:,iLayer,iUser,iBand));
                                P0 = real(W0{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * M0(:,iLayer,iUser,iBand));
                                Q0 = imag(W0{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * M0(:,iLayer,iUser,iBand));
                                constVector = [constVector, (P0^2 + Q0^2) / B0(iLayer,iUser,iBand) + (2 / B0(iLayer,iUser,iBand)) * (P0 * (P - P0) + Q0 * (Q - Q0)) ...
                                    - ((P0^2 + Q0^2) / (B0(iLayer,iUser,iBand)^2)) * (B(iLayer,iUser,iBand) - B0(iLayer,iUser,iBand)) >= G(iLayer,iUser,iBand)];
                                
                            end
                        end
                    end
                    
                    solverSettings = sdpsettings('verbose',0,'solver','Mosek');
                    constVector = [constVector, vector(M)' * vector(M) <= sum(SimStructs.baseStruct{iBase,1}.sPower(1,:))];
                    solutionOut = solvesdp(constVector,epiObjective,solverSettings);
                    
                    if solutionOut.problem == 0
                        M0 = full(double(M));
                    else
                        display(solutionOut);
                        break;
                    end
                    
                    for iBand = 1:nBands
                        SimStructs.baseStruct{iBase,1}.P{iBand,1} = M0(:,:,:,iBand);
                        SimParams.Debug.globalExchangeInfo.P{iBase,iBand} = M0(:,:,:,iBand);
                    end
                    SimParams.Debug.globalExchangeInfo.I{iBase,1} = double(I);
                    
                    SimParams.Debug.globalExchangeInfo.funcOut{1,iBase} = M0;
                    SimParams.Debug.globalExchangeInfo.funcOut{2,iBase} = double(B);
                    SimParams.Debug.globalExchangeInfo.funcOut{4,iBase} = double(T);
                    SimParams.Debug.globalExchangeInfo.funcOut{5,iBase} = W0;
                end
                
                for jBand = 1:nBands
                    for iUser = 1:nUsers
                        dCell = SimStructs.userStruct{iUser,1}.baseNode;
                        for jBase = 1:nBases
                            if jBase ~= dCell
                                SimParams.Debug.globalExchangeInfo.gI{jBase,1}(:,iUser,jBand) = 0.5 * (SimParams.Debug.globalExchangeInfo.I{jBase,1}(:,iUser,jBand,jBase) + SimParams.Debug.globalExchangeInfo.I{dCell,1}(:,iUser,jBand,jBase));
                                SimParams.Debug.globalExchangeInfo.D{dCell,1}(:,iUser,jBand,jBase) = (SimParams.Debug.globalExchangeInfo.gI{jBase,1}(:,iUser,jBand) - SimParams.Debug.globalExchangeInfo.I{dCell,1}(:,iUser,jBand,jBase)) * stepIndex + SimParams.Debug.globalExchangeInfo.D{dCell,1}(:,iUser,jBand,jBase);
                                SimParams.Debug.globalExchangeInfo.D{jBase,1}(:,iUser,jBand,jBase) = (SimParams.Debug.globalExchangeInfo.gI{jBase,1}(:,iUser,jBand) - SimParams.Debug.globalExchangeInfo.I{jBase,1}(:,iUser,jBand,jBase)) * stepIndex + SimParams.Debug.globalExchangeInfo.D{jBase,1}(:,iUser,jBand,jBase);
                            end
                        end
                    end
                end
                
                [SimParams,SimStructs] = updateIteratePerformance(SimParams,SimStructs);
            end
        end
        
        
end
