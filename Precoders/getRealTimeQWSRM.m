function [SimParams,SimStructs] = getRealTimeQWSRM(SimParams,SimStructs)

if 0
    
    preLogue;
    pNoise = 1e-30;
    sdpOptions = sdpsettings('verbose',0,'solver','gurobi');
    SimParams.Debug.globalExchangeInfo.funcOut = cell(5,SimParams.nBases);
    
    if SimParams.iDrop == 1
        for iUser = 1:nUsers
            SimStructs.userStruct{iUser,1}.pW = cell(nBands,1);
        end
    end
    
    for iBase = 1:nBases
        SimStructs.baseStruct{iBase,1}.selectionType = 'BF';
    end
    
    switch selectionMethod
        
        case 'perBSAlloc'
            
            for iBase = 1:nBases
                
                epsilon = 0.01;
                maxRank = SimParams.maxRank;
                kUsers = usersPerCell(iBase,1);
                taylorApprox = SimParams.additionalParams;
                
                currentObj = -50;
                currentIteration = 0;
                totalRanks = 1:maxRank;
                
                [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'BF',iBase);
                [SimParams, SimStructs] = initializeSCApoint(SimParams,SimStructs,iBase);
                M0 = SimParams.Debug.globalExchangeInfo.funcOut{1,iBase};B0 = SimParams.Debug.globalExchangeInfo.funcOut{2,iBase};
                G0 = SimParams.Debug.globalExchangeInfo.funcOut{3,iBase};T0 = SimParams.Debug.globalExchangeInfo.funcOut{4,iBase};
                W0 = SimParams.Debug.globalExchangeInfo.funcOut{5,iBase};
                
                while currentIteration < maxIterations
                    
                    M = cell(nBands,1);
                    epiObjective = sdpvar(1);
                    userObjective = sdpvar(kUsers,1,'full');
                    T = sdpvar(maxRank,kUsers,nBands,'full');
                    B = sdpvar(maxRank,kUsers,nBands,'full');
                    G = sdpvar(maxRank,kUsers,nBands,'full');
                    
                    for iBand = 1:nBands
                        M{iBand,1} = sdpvar(SimParams.nTxAntenna,maxRank,kUsers,'full','complex');
                    end
                    
                    gConstraints = [];
                    for iUser = 1:kUsers
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        gConstraints = [gConstraints, userWts(cUser,1) * abs(QueuedPkts(cUser,1) - sum(vec(T(:,iUser,:)))) <= userObjective(iUser,1)];
                    end
                    
                    gConstraints = [gConstraints, epiObjective >= norm(userObjective,qExponent)];
                    for iBand = 1:nBands
                        for iUser = 1:kUsers
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for iLayer = 1:maxRank
                                intVector = sqrt(SimParams.N) * W0{cUser,iBand}(:,iLayer);
                                for jUser = 1:kUsers
                                    if jUser ~= iUser
                                        intVector = [intVector; W0{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * M{iBand,1}(:,:,jUser)];
                                    else
                                        intVector = [intVector; W0{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * M{iBand,1}(:,iLayer ~= totalRanks,iUser)];
                                    end
                                end
                                
                                gConstraints = [gConstraints, intVector' * intVector <= B(iLayer,iUser,iBand)];
                                if strcmpi(taylorApprox,'LOG')
                                    gConstraints = [gConstraints, G(iLayer,iUser,iBand) >= G0(iLayer,iUser,iBand) * (1 - epsilon)];
                                    gConstraints = [gConstraints, G(iLayer,iUser,iBand) <= G0(iLayer,iUser,iBand) * (1 + epsilon)];
                                    gConstraints = [gConstraints, log(1 + G0(iLayer,iUser,iBand)) + (G(iLayer,iUser,iBand) - G0(iLayer,iUser,iBand)) / (1 + G0(iLayer,iUser,iBand)) - (G(iLayer,iUser,iBand) - G0(iLayer,iUser,iBand))^2 / (2 * (1 + G0(iLayer,iUser,iBand))^2) ...
                                        >= T(iLayer,iUser,iBand) * log(2)];
                                else
                                    gConstraints = [gConstraints, T(iLayer,iUser,iBand) >= T0(iLayer,iUser,iBand) * (1 - epsilon)];
                                    gConstraints = [gConstraints, T(iLayer,iUser,iBand) <= T0(iLayer,iUser,iBand) * (1 + epsilon)];
                                    gConstraints = [gConstraints, 1 + G(iLayer,iUser,iBand) >= ...
                                        exp(T0(iLayer,iUser,iBand) * log(2)) * (1 + log(2) * (T(iLayer,iUser,iBand) - T0(iLayer,iUser,iBand)) + 0.5 * log(2)^2 * (T(iLayer,iUser,iBand) - T0(iLayer,iUser,iBand))^2)];
                                end
                                
                                H = cH{iBase,iBand}(:,:,cUser);
                                gConstraints = [gConstraints, abs(W0{cUser,iBand}(:,iLayer)' * H * M0{iBand,1}(:,iLayer,iUser))^2 / B0(iLayer,iUser,iBand) ...
                                    + 2 * (M0{iBand,1}(:,iLayer,iUser)' * H' * W0{cUser,iBand}(:,iLayer) * W0{cUser,iBand}(:,iLayer)' * H * (M{iBand,1}(:,iLayer,iUser) - M0{iBand,1}(:,iLayer,iUser))) / B0(iLayer,iUser,iBand) ...
                                    - (abs(W0{cUser,iBand}(:,iLayer)' * H * M0{iBand,1}(:,iLayer,iUser))^2 * (B(iLayer,iUser,iBand) - B0(iLayer,iUser,iBand))) / B0(iLayer,iUser,iBand)^2 >= G(iLayer,iUser,iBand)];
                            end
                        end
                    end
                    
                    totPower = 0;
                    for iBand = 1:nBands
                        totPower = totPower + M{iBand,1}(:)' * M{iBand,1}(:);
                    end
                    gConstraints = [gConstraints, totPower < sum(SimStructs.baseStruct{iBase,1}.sPower(1,:))];
                    sdpSol = solvesdp(gConstraints,epiObjective,sdpOptions);
                    
                    if sdpSol.problem == 0
                        B0 = full(double(B));
                        for iBand = 1:nBands
                            M0{iBand,1} = full(double(M{iBand,1})) + pNoise;
                        end
                        if strcmpi(taylorApprox,'LOG')
                            G0 = full(double(G));
                        else
                            T0 = full(double(T));
                        end
                        currentIteration = currentIteration + 1;
                    else
                        display(sdpSol);
                        break;
                    end
                    
                    epiObjective = double(epiObjective);
                    if norm(epiObjective - currentObj,2) < 1e-5
                        break;
                    else
                        currentObj = epiObjective;
                    end
                    
                    for iBand = 1:nBands
                        SimStructs.baseStruct{iBase,1}.P{iBand,1} = M0{iBand,1};
                    end
                    
                    [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE',iBase);
                    
                end
                
            end
            
            [SimParams,SimStructs] = updateIteratePerformance(SimParams,SimStructs);
            
        case 'distBSAlloc'
            
            epsilon = 0.25;
            stepIndex = 10;
            IGlobal = cell(nBases,1);
            DGlobal = cell(nBases,1);
            fW = cell(nUsers,nBands);
            maxRank = SimParams.maxRank;
            
            for iExchange = 1:SimParams.nExchangesOverBH
                
                for iBase = 1:nBases
                    
                    if and((SimParams.distIteration == 1),(iExchange == 1))
                        IGlobal{iBase,1} = ones(maxRank,nUsers,nBands) * sqrt(10);
                        DGlobal{iBase,1} = zeros(maxRank,nUsers,nBands,nBases);
                        SimParams.Debug.globalExchangeInfo.D{iBase,1} = zeros(maxRank,nUsers,nBands,nBases);
                        for iUser = 1:usersPerCell(iBase,1)
                            for iBand = 1:nBands
                                cUser = cellUserIndices{iBase,1}(iUser,1);
                                fW{cUser,iBand} = ones(SimParams.nRxAntenna,1);
                            end
                        end
                        [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'Ones',iBase);
                    else
                        IGlobal{iBase,1} = SimParams.Debug.globalExchangeInfo.gI{iBase,1};
                        DGlobal{iBase,1} = SimParams.Debug.globalExchangeInfo.D{iBase,1};
                        for iUser = 1:usersPerCell(iBase,1)
                            for iBand = 1:nBands
                                cUser = cellUserIndices{iBase,1}(iUser,1);
                                fW{cUser,iBand} = SimParams.Debug.globalExchangeInfo.W{iBand,1};
                            end
                        end
                    end
                    
                    renewSearch = (iExchange == 1);
                    if renewSearch == 1
                        DGlobal{iBase,1} = zeros(maxRank,nUsers,nBands,nBases);
                        for iUser = 1:usersPerCell(iBase,1)
                            for iBand = 1:nBands
                                cUser = cellUserIndices{iBase,1}(iUser,1);
                                fW{cUser,iBand} = ones(SimParams.nRxAntenna,1);
                            end
                        end
                        [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'Ones',iBase);
                    end
                    
                    resetPoint = 0 * (SimParams.exchangeResetInterval / SimParams.nBases) * (iBase - 1);
                    if and(mod(SimParams.distIteration,SimParams.exchangeResetInterval) == resetPoint,iExchange == 1)
                        DGlobal{iBase,1} = zeros(maxRank,nUsers,nBands,nBases);
                        IGlobal{iBase,1} = ones(maxRank,nUsers,nBands) * sqrt(10);
                        SimParams.Debug.globalExchangeInfo.gI{iBase,1} = ones(maxRank,nUsers,nBands) * sqrt(10);
                        fprintf('Resetting ADMM Variables for BS - %d \n',iBase);
                        for iUser = 1:usersPerCell(iBase,1)
                            for iBand = 1:nBands
                                cUser = cellUserIndices{iBase,1}(iUser,1);
                                fW{cUser,iBand} = ones(SimParams.nRxAntenna,1);
                            end
                        end
                        [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'Ones',iBase);
                    end
                end
                
                for iBase = 1:nBases
                    
                    kUsers = usersPerCell(iBase,1);
                    taylorApprox = SimParams.additionalParams;
                    
                    currentObj = -50;
                    currentIteration = 0;
                    totalRanks = 1:maxRank;
                    
                    if renewSearch
                        SimStructs.baseStruct{iBase,1}.selectionType = 'BF';
                    else
                        SimStructs.baseStruct{iBase,1}.selectionType = 'Last';
                    end
                    
                    [SimParams, SimStructs] = initializeSCApoint(SimParams,SimStructs,iBase);
                    M0 = SimParams.Debug.globalExchangeInfo.funcOut{1,iBase};B0 = SimParams.Debug.globalExchangeInfo.funcOut{2,iBase};
                    G0 = SimParams.Debug.globalExchangeInfo.funcOut{3,iBase};T0 = SimParams.Debug.globalExchangeInfo.funcOut{4,iBase};
                    
                    while currentIteration < maxIterations
                        
                        M = cell(nBands,1);
                        epiObjective = sdpvar(1);
                        userObjective = sdpvar(kUsers,1,'full');
                        T = sdpvar(maxRank,kUsers,nBands,'full');
                        B = sdpvar(maxRank,kUsers,nBands,'full');
                        G = sdpvar(maxRank,kUsers,nBands,'full');
                        I = sdpvar(maxRank,nUsers,nBands,nBases,'full');
                        
                        for iBand = 1:nBands
                            M{iBand,1} = sdpvar(SimParams.nTxAntenna,maxRank,kUsers,'full','complex');
                        end
                        
                        gConstraints = [];
                        for iUser = 1:kUsers
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            gConstraints = [gConstraints, userWts(cUser,1) * abs(QueuedPkts(cUser,1) - sum(vec(T(:,iUser,:)))) <= userObjective(iUser,1)];
                        end
                        
                        augmentedTerms = 0;
                        for jBand = 1:nBands
                            for jUser = 1:nUsers
                                nCellIndex = SimStructs.userStruct{jUser,1}.baseNode;
                                if nCellIndex ~= iBase
                                    augmentedTerms = augmentedTerms + sum((IGlobal{iBase,1}(:,jUser,jBand) - I(:,jUser,jBand,iBase)) .* DGlobal{iBase,1}(:,jUser,jBand,iBase)) ...
                                        + (stepIndex / 2) * norm(IGlobal{iBase,1}(:,jUser,jBand) - I(:,jUser,jBand,iBase))^2;
                                end
                            end
                            
                            for jUser = 1:kUsers
                                cUser = cellUserIndices{iBase,1}(jUser,1);
                                for jBase = 1:nBases
                                    if jBase ~= iBase
                                        augmentedTerms = augmentedTerms + sum((IGlobal{jBase,1}(:,cUser,jBand) - I(:,cUser,jBand,jBase)) .* DGlobal{iBase,1}(:,cUser,jBand,jBase)) ...
                                            + (stepIndex / 2) * norm(IGlobal{jBase,1}(:,cUser,jBand) - I(:,cUser,jBand,jBase))^2;
                                    end
                                end
                            end
                        end
                        
                        gConstraints = [gConstraints, epiObjective >= norm(userObjective,qExponent) + augmentedTerms];
                        
                        for iBand = 1:nBands
                            for iUser = 1:kUsers
                                cUser = cellUserIndices{iBase,1}(iUser,1);
                                for iLayer = 1:maxRank
                                    intVector = sqrt(SimParams.N) * fW{cUser,iBand}(:,iLayer)';
                                    for jUser = 1:kUsers
                                        if jUser ~= iUser
                                            intVector = [intVector, fW{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * M{iBand,1}(:,:,jUser)];
                                        else
                                            intVector = [intVector, fW{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * M{iBand,1}(:,iLayer ~= totalRanks,iUser)];
                                        end
                                    end
                                    
                                    for jBase = 1:nBases
                                        if jBase ~= iBase
                                            intVector = [intVector, I(iLayer,cUser,iBand,jBase)];
                                        end
                                    end
                                    
                                    gConstraints = [gConstraints, intVector * intVector' <= B(iLayer,iUser,iBand)];
                                    
                                    for jUser = 1:nUsers
                                        nCellIndex = SimStructs.userStruct{jUser,1}.baseNode;
                                        if nCellIndex ~= iBase
                                            for jLayer = 1:maxRank
                                                intVector = [];
                                                for inUser = 1:kUsers
                                                    intVector = [intVector, fW{jUser,iBand}(:,jLayer)' * cH{iBase,iBand}(:,:,jUser) * M{iBand,1}(:,:,inUser)];
                                                end
                                                gConstraints = [gConstraints, norm(intVector,2) <= I(:,jUser,iBand,iBase)];
                                            end
                                        end
                                    end
                                    
                                    if strcmpi(taylorApprox,'LOG')
                                        gConstraints = [gConstraints, G(iLayer,iUser,iBand) >= G0(iLayer,iUser,iBand) * (1 - epsilon)];
                                        gConstraints = [gConstraints, G(iLayer,iUser,iBand) <= G0(iLayer,iUser,iBand) * (1 + epsilon)];
                                        gConstraints = [gConstraints, log(1 + G0(iLayer,iUser,iBand)) + (G(iLayer,iUser,iBand) - G0(iLayer,iUser,iBand)) / (1 + G0(iLayer,iUser,iBand)) - (G(iLayer,iUser,iBand) - G0(iLayer,iUser,iBand))^2 / (2 * (1 + G0(iLayer,iUser,iBand))^2) ...
                                            >= T(iLayer,iUser,iBand) * log(2)];
                                    else
                                        gConstraints = [gConstraints, T(iLayer,iUser,iBand) >= T0(iLayer,iUser,iBand) * (1 - epsilon)];
                                        gConstraints = [gConstraints, T(iLayer,iUser,iBand) <= T0(iLayer,iUser,iBand) * (1 + epsilon)];
                                        gConstraints = [gConstraints, 1 + G(iLayer,iUser,iBand) >= ...
                                            exp(T0(iLayer,iUser,iBand) * log(2)) * (1 + log(2) * (T(iLayer,iUser,iBand) - T0(iLayer,iUser,iBand)) + 0.5 * log(2)^2 * (T(iLayer,iUser,iBand) - T0(iLayer,iUser,iBand))^2)];
                                    end
                                    
                                    H = cH{iBase,iBand}(:,:,cUser);
                                    gConstraints = [gConstraints, abs(fW{cUser,iBand}(:,iLayer)' * H * M0{iBand,1}(:,iLayer,iUser))^2 / B0(iLayer,iUser,iBand) ...
                                        + 2 * (M0{iBand,1}(:,iLayer,iUser)' * H' * fW{cUser,iBand}(:,iLayer) * fW{cUser,iBand}(:,iLayer)' * H * (M{iBand,1}(:,iLayer,iUser) - M0{iBand,1}(:,iLayer,iUser))) / B0(iLayer,iUser,iBand) ...
                                        - (abs(fW{cUser,iBand}(:,iLayer)' * H * M0{iBand,1}(:,iLayer,iUser))^2 * (B(iLayer,iUser,iBand) - B0(iLayer,iUser,iBand))) / B0(iLayer,iUser,iBand)^2 >= G(iLayer,iUser,iBand)];
                                end
                            end
                        end
                        
                        totPower = 0;
                        for iBand = 1:nBands
                            totPower = totPower + norm(M{iBand,1}(:))^2;
                        end
                        gConstraints = [gConstraints, totPower <= sum(SimStructs.baseStruct{iBase,1}.sPower(1,:))];
                        
                        if isnan(gConstraints)
                            keyboard;
                        end
                        sdpSol = solvesdp(gConstraints,epiObjective,sdpOptions);
                        
                        if sum(sdpSol.problem == [0 4])
                            B0 = full(double(B));
                            for iBand = 1:nBands
                                M0{iBand,1} = full(double(M{iBand,1})) + pNoise;
                            end
                            if strcmpi(taylorApprox,'LOG')
                                G0 = full(double(G));
                            else
                                T0 = full(double(T));
                            end
                            currentIteration = currentIteration + 1;
                        else
                            display(sdpSol);
                            break;
                        end
                        
                        epiObjective = double(epiObjective);
                        fprintf('%.2f, \t',epiObjective);
                        if norm(epiObjective - currentObj,2) < 1e-5
                            break;
                        else
                            currentObj = epiObjective;
                        end
                        
                        for iBand = 1:nBands
                            SimStructs.baseStruct{iBase,1}.P{iBand,1} = M0{iBand,1};
                        end
                        
                    end
                    
                    for iBand = 1:nBands
                        SimParams.Debug.globalExchangeInfo.P{iBase,iBand} = M0{iBand,1};
                        SimStructs.baseStruct{iBase,1}.P{iBand,1} = M0{iBand,1};
                    end
                    
                    fprintf('\n');
                    SimParams.Debug.globalExchangeInfo.I{iBase,1} = double(I);
                    
                end
                
                for jBand = 1:nBands
                    for iUser = 1:nUsers
                        dCell = SimStructs.userStruct{iUser,1}.baseNode;
                        for jBase = 1:nBases
                            if jBase ~= dCell
                                IGlobal{jBase,1}(:,iUser,jBand) = 0.5 * (SimParams.Debug.globalExchangeInfo.I{jBase,1}(:,iUser,jBand,jBase) + SimParams.Debug.globalExchangeInfo.I{dCell,1}(:,iUser,jBand,jBase));
                                SimParams.Debug.globalExchangeInfo.D{dCell,1}(:,iUser,jBand,jBase) = (IGlobal{jBase,1}(:,iUser,jBand) - SimParams.Debug.globalExchangeInfo.I{dCell,1}(:,iUser,jBand,jBase)) * stepIndex + SimParams.Debug.globalExchangeInfo.D{dCell,1}(:,iUser,jBand,jBase);
                                SimParams.Debug.globalExchangeInfo.D{jBase,1}(:,iUser,jBand,jBase) = (IGlobal{jBase,1}(:,iUser,jBand) - SimParams.Debug.globalExchangeInfo.I{jBase,1}(:,iUser,jBand,jBase)) * stepIndex + SimParams.Debug.globalExchangeInfo.D{jBase,1}(:,iUser,jBand,jBase);
                            end
                        end
                    end
                end
                
                SimParams.Debug.globalExchangeInfo.gI = IGlobal;
                [SimParams,SimStructs] = updateIteratePerformance(SimParams,SimStructs);
                [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'Ones');
                
                for iBand = 1:nBands
                    for iUser = 1:nUsers
                        fW{iUser,iBand} = SimStructs.userStruct{iUser,1}.W{iBand,1};
                        SimParams.Debug.globalExchangeInfo.W{iUser,iBand} = SimStructs.userStruct{iUser,1}.W{iBand,1};
                    end
                end
                
                fprintf('Exchange Happened - %d \n',iExchange);
                
            end
            
    end
    
else
    
    preLogue;
    pNoise = 1e-30;
    SimParams.Debug.globalExchangeInfo.funcOut = cell(5,SimParams.nBases);
    
    if SimParams.iDrop == 1
        for iUser = 1:nUsers
            SimStructs.userStruct{iUser,1}.pW = cell(nBands,1);
        end
    end
    
    for iBase = 1:nBases
        SimStructs.baseStruct{iBase,1}.selectionType = 'BF';
    end
    
    switch selectionMethod
        
        case 'perBSAlloc'
            
            for iBase = 1:nBases
                
                maxRank = SimParams.maxRank;
                kUsers = usersPerCell(iBase,1);
                
                currentObj = -50;
                currentIteration = 0;
                totalRanks = 1:maxRank;
                
                [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'BF',iBase);
                [SimParams, SimStructs] = initializeSCApoint(SimParams,SimStructs,iBase);
                M0 = SimParams.Debug.globalExchangeInfo.funcOut{1,iBase};B0 = SimParams.Debug.globalExchangeInfo.funcOut{2,iBase};
                W0 = SimParams.Debug.globalExchangeInfo.funcOut{5,iBase};
                
                while currentIteration < maxIterations
                    
                    cvx_begin
                    
                    variable M(SimParams.nTxAntenna,maxRank,kUsers,nBands) complex
                    variables T(maxRank,kUsers,nBands) B(maxRank,kUsers,nBands) G(maxRank,kUsers,nBands)
                    variable epiObjective userObjective(kUsers,1)
                    
                    minimize(epiObjective);
                    
                    subject to
                    
                    for iUser = 1:kUsers
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        userWts(cUser,1) * abs(QueuedPkts(cUser,1) - sum(vec(T(:,iUser,:)))) <= userObjective(iUser,1);
                    end
                    
                    epiObjective >= norm(userObjective,qExponent);
                    for iBand = 1:nBands
                        for iUser = 1:kUsers
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for iLayer = 1:maxRank
                                intVector = sqrt(SimParams.N) * W0{cUser,iBand}(:,iLayer)';
                                for jUser = 1:kUsers
                                    if jUser ~= iUser
                                        intVector = [intVector, W0{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * M(:,:,jUser,iBand)];
                                    else
                                        intVector = [intVector, W0{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * M(:,iLayer ~= totalRanks,iUser,iBand)];
                                    end
                                end
                                
                                intVector * intVector' <= B(iLayer,iUser,iBand);
                                log(1 + G(iLayer,iUser,iBand)) >= T(iLayer,iUser,iBand) * log(2);
                                
                                H = cH{iBase,iBand}(:,:,cUser);
                                abs(W0{cUser,iBand}(:,iLayer)' * H * M0(:,iLayer,iUser,iBand))^2 / B0(iLayer,iUser,iBand) ...
                                    + 2 * (M0(:,iLayer,iUser,iBand)' * H' * W0{cUser,iBand}(:,iLayer) * W0{cUser,iBand}(:,iLayer)' * H * (M(:,iLayer,iUser,iBand) - M0(:,iLayer,iUser,iBand))) / B0(iLayer,iUser,iBand) ...
                                    - (abs(W0{cUser,iBand}(:,iLayer)' * H * M0(:,iLayer,iUser,iBand))^2 * (B(iLayer,iUser,iBand) - B0(iLayer,iUser,iBand))) / B0(iLayer,iUser,iBand)^2 >= G(iLayer,iUser,iBand);
                            end
                        end
                    end
                    
                    norm(vec(M),2) <= sqrt(sum(SimStructs.baseStruct{iBase,1}.sPower(1,:)));
                    
                    cvx_end
                    
                    if strfind(cvx_status,'Solved')
                        B0 = full(B);
                        M0 = full(M) + pNoise;
                        currentIteration = currentIteration + 1;
                    else
                        display(cvx_status);
                        break;
                    end
                    
                    epiObjective = cvx_optval;
                    if norm(epiObjective - currentObj,2) < 1e-5
                        break;
                    else
                        currentObj = epiObjective;
                    end
                    
                    for iBand = 1:nBands
                        SimStructs.baseStruct{iBase,1}.P{iBand,1} = M0(:,:,:,iBand);
                    end
                    
                    [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE',iBase);
                    
                end
                
            end
            
            [SimParams,SimStructs] = updateIteratePerformance(SimParams,SimStructs);
            
        case 'distBSAlloc'
            
            stepIndex = 5;
            IGlobal = cell(nBases,1);
            DGlobal = cell(nBases,1);
            fW = cell(nUsers,nBands);
            maxRank = SimParams.maxRank;
            
            for iExchange = 1:SimParams.nExchangesOverBH
                
                for iBase = 1:nBases
                    
                    if and((SimParams.distIteration == 1),(iExchange == 1))
                        IGlobal{iBase,1} = ones(maxRank,nUsers,nBands) * sqrt(10);
                        DGlobal{iBase,1} = zeros(maxRank,nUsers,nBands,nBases);
                        SimParams.Debug.globalExchangeInfo.D{iBase,1} = zeros(maxRank,nUsers,nBands,nBases);
                        for iUser = 1:usersPerCell(iBase,1)
                            for iBand = 1:nBands
                                cUser = cellUserIndices{iBase,1}(iUser,1);
                                fW{cUser,iBand} = ones(SimParams.nRxAntenna,1);
                            end
                        end
                        [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'Ones',iBase);
                    else
                        IGlobal{iBase,1} = SimParams.Debug.globalExchangeInfo.gI{iBase,1};
                        DGlobal{iBase,1} = SimParams.Debug.globalExchangeInfo.D{iBase,1};
                        for iUser = 1:usersPerCell(iBase,1)
                            for iBand = 1:nBands
                                cUser = cellUserIndices{iBase,1}(iUser,1);
                                fW{cUser,iBand} = SimParams.Debug.globalExchangeInfo.W{iBand,1};
                            end
                        end
                    end
                    
                    renewSearch = (iExchange == 1);
                    if renewSearch == 1
                        DGlobal{iBase,1} = zeros(maxRank,nUsers,nBands,nBases);
                        for iUser = 1:usersPerCell(iBase,1)
                            for iBand = 1:nBands
                                cUser = cellUserIndices{iBase,1}(iUser,1);
                                fW{cUser,iBand} = ones(SimParams.nRxAntenna,1);
                            end
                        end
                        [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'Ones',iBase);
                    end
                    
                    resetPoint = 0 * (SimParams.exchangeResetInterval / SimParams.nBases) * (iBase - 1);
                    if and(mod(SimParams.distIteration,SimParams.exchangeResetInterval) == resetPoint,iExchange == 1)
                        DGlobal{iBase,1} = zeros(maxRank,nUsers,nBands,nBases);
                        IGlobal{iBase,1} = ones(maxRank,nUsers,nBands) * sqrt(10);
                        SimParams.Debug.globalExchangeInfo.gI{iBase,1} = ones(maxRank,nUsers,nBands) * sqrt(10);
                        fprintf('Resetting ADMM Variables for BS - %d \n',iBase);
                        for iUser = 1:usersPerCell(iBase,1)
                            for iBand = 1:nBands
                                cUser = cellUserIndices{iBase,1}(iUser,1);
                                fW{cUser,iBand} = ones(SimParams.nRxAntenna,1);
                            end
                        end
                        [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'Ones',iBase);
                    end
                end
                
                for iBase = 1:nBases
                    
                    currentObj = -50;
                    currentIteration = 0;
                    totalRanks = 1:maxRank;
                    kUsers = usersPerCell(iBase,1);
                                        
                    if renewSearch
                        SimStructs.baseStruct{iBase,1}.selectionType = 'BF';
                    else
                        SimStructs.baseStruct{iBase,1}.selectionType = 'Last';
                    end
                    
                    [SimParams, SimStructs] = initializeSCApoint(SimParams,SimStructs,iBase);
                    M0 = SimParams.Debug.globalExchangeInfo.funcOut{1,iBase};B0 = SimParams.Debug.globalExchangeInfo.funcOut{2,iBase};
                    
                    while currentIteration < maxIterations
                        
                        cvx_begin
                        
                        variable M(SimParams.nTxAntenna,maxRank,kUsers,nBands) complex
                        variables T(maxRank,kUsers,nBands) B(maxRank,kUsers,nBands) G(maxRank,kUsers,nBands)
                        variables I(maxRank,nUsers,nBands,nBases) epiObjective userObjective(kUsers,1)
                        
                        for iUser = 1:kUsers
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            userWts(cUser,1) * abs(QueuedPkts(cUser,1) - sum(vec(T(:,iUser,:)))) <= userObjective(iUser,1);
                        end
                        
                        augmentedTerms = 0;
                        for jBand = 1:nBands
                            for jUser = 1:nUsers
                                nCellIndex = SimStructs.userStruct{jUser,1}.baseNode;
                                if nCellIndex ~= iBase
                                    augmentedTerms = augmentedTerms + sum((IGlobal{iBase,1}(:,jUser,jBand) - I(:,jUser,jBand,iBase)) .* DGlobal{iBase,1}(:,jUser,jBand,iBase)) ...
                                        + (stepIndex / 2) * sum(vec(IGlobal{iBase,1}(:,jUser,jBand) - I(:,jUser,jBand,iBase)).^2);
                                end
                            end
                            
                            for jUser = 1:kUsers
                                cUser = cellUserIndices{iBase,1}(jUser,1);
                                for jBase = 1:nBases
                                    if jBase ~= iBase
                                        augmentedTerms = augmentedTerms + sum((IGlobal{jBase,1}(:,cUser,jBand) - I(:,cUser,jBand,jBase)) .* DGlobal{iBase,1}(:,cUser,jBand,jBase)) ...
                                            + (stepIndex / 2) * sum(vec(IGlobal{jBase,1}(:,cUser,jBand) - I(:,cUser,jBand,jBase)).^2);
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
                                    intVector = sqrt(SimParams.N) * fW{cUser,iBand}(:,iLayer)';
                                    for jUser = 1:kUsers
                                        if jUser ~= iUser
                                            intVector = [intVector, fW{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * M(:,:,jUser,iBand)];
                                        else
                                            intVector = [intVector, fW{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * M(:,iLayer ~= totalRanks,iUser,iBand)];
                                        end
                                    end
                                    
                                    for jBase = 1:nBases
                                        if jBase ~= iBase
                                            intVector = [intVector, I(iLayer,cUser,iBand,jBase)];
                                        end
                                    end
                                    
                                    norm(intVector) <= sqrt(B(iLayer,iUser,iBand));
                                    
                                    for jUser = 1:nUsers
                                        nCellIndex = SimStructs.userStruct{jUser,1}.baseNode;
                                        if nCellIndex ~= iBase
                                            for jLayer = 1:maxRank
                                                intVector = [];
                                                for inUser = 1:kUsers
                                                    intVector = [intVector, fW{jUser,iBand}(:,jLayer)' * cH{iBase,iBand}(:,:,jUser) * M(:,:,inUser,iBand)];
                                                end
                                                norm(intVector,2) <= I(:,jUser,iBand,iBase);
                                            end
                                        end
                                    end
                                    
                                    log(1 + G(iLayer,iUser,iBand)) >= T(iLayer,iUser,iBand) * log(2);
                                    
                                    H = cH{iBase,iBand}(:,:,cUser);
                                    abs(fW{cUser,iBand}(:,iLayer)' * H * M0(:,iLayer,iUser,iBand))^2 / B0(iLayer,iUser,iBand) ...
                                        + 2 * real(M0(:,iLayer,iUser,iBand)' * H' * fW{cUser,iBand}(:,iLayer) * fW{cUser,iBand}(:,iLayer)' * H * (M(:,iLayer,iUser,iBand) - M0(:,iLayer,iUser,iBand))) / B0(iLayer,iUser,iBand) ...
                                        - (abs(fW{cUser,iBand}(:,iLayer)' * H * M0(:,iLayer,iUser,iBand))^2 * (B(iLayer,iUser,iBand) - B0(iLayer,iUser,iBand))) / B0(iLayer,iUser,iBand)^2 >= G(iLayer,iUser,iBand);
                                    imag(M0(:,iLayer,iUser,iBand)' * H' * fW{cUser,iBand}(:,iLayer) * fW{cUser,iBand}(:,iLayer)' * H * (M(:,iLayer,iUser,iBand) - M0(:,iLayer,iUser,iBand))) == 0;
                                end
                            end
                        end
                        
                        norm(vec(M),2) <= sqrt(sum(SimStructs.baseStruct{iBase,1}.sPower(1,:)));
                        
                        cvx_end
                        
                        if strfind(cvx_status,'Solved')
                            B0 = full(B);M0 = full(M);
                            currentIteration = currentIteration + 1;
                        else
                            display(cvx_status);
                            break;
                        end
                        
                        epiObjective = cvx_optval;
                        fprintf('%.2f, \t',epiObjective);
                        if norm(epiObjective - currentObj,2) < 1e-5
                            break;
                        else
                            currentObj = epiObjective;
                        end
                        
                        for iBand = 1:nBands
                            SimStructs.baseStruct{iBase,1}.P{iBand,1} = M0(:,:,:,iBand);
                        end
                        
                    end
                    
                    for iBand = 1:nBands
                        SimParams.Debug.globalExchangeInfo.P{iBase,iBand} = M0(:,:,:,iBand);
                        SimStructs.baseStruct{iBase,1}.P{iBand,1} = M0(:,:,:,iBand);
                    end
                    
                    fprintf('\n');
                    SimParams.Debug.globalExchangeInfo.I{iBase,1} = double(I);
                    
                end
                
                for jBand = 1:nBands
                    for iUser = 1:nUsers
                        dCell = SimStructs.userStruct{iUser,1}.baseNode;
                        for jBase = 1:nBases
                            if jBase ~= dCell
                                IGlobal{jBase,1}(:,iUser,jBand) = 0.5 * (SimParams.Debug.globalExchangeInfo.I{jBase,1}(:,iUser,jBand,jBase) + SimParams.Debug.globalExchangeInfo.I{dCell,1}(:,iUser,jBand,jBase));
                                SimParams.Debug.globalExchangeInfo.D{dCell,1}(:,iUser,jBand,jBase) = (IGlobal{jBase,1}(:,iUser,jBand) - SimParams.Debug.globalExchangeInfo.I{dCell,1}(:,iUser,jBand,jBase)) * stepIndex + SimParams.Debug.globalExchangeInfo.D{dCell,1}(:,iUser,jBand,jBase);
                                SimParams.Debug.globalExchangeInfo.D{jBase,1}(:,iUser,jBand,jBase) = (IGlobal{jBase,1}(:,iUser,jBand) - SimParams.Debug.globalExchangeInfo.I{jBase,1}(:,iUser,jBand,jBase)) * stepIndex + SimParams.Debug.globalExchangeInfo.D{jBase,1}(:,iUser,jBand,jBase);
                            end
                        end
                    end
                end
                
                SimParams.Debug.globalExchangeInfo.gI = IGlobal;
                [SimParams,SimStructs] = updateIteratePerformance(SimParams,SimStructs);
                [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'Ones');
                
                for iBand = 1:nBands
                    for iUser = 1:nUsers
                        fW{iUser,iBand} = SimStructs.userStruct{iUser,1}.W{iBand,1};
                        SimParams.Debug.globalExchangeInfo.W{iUser,iBand} = SimStructs.userStruct{iUser,1}.W{iBand,1};
                    end
                end
                
                fprintf('Exchange Happened - %d \n',iExchange);
                
            end
            
    end
    
end

end