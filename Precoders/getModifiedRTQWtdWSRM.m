function [SimParams,SimStructs] = getModifiedRTQWtdWSRM(SimParams,SimStructs)

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

nPreExchanges = 1;
nSlots = SimParams.exchangeResetInterval;
mdpFactor = 1 - (SimParams.userDoppler / norm(SimParams.userDoppler)^2);
if or((SimParams.distIteration - 1) == 0,mod((SimParams.iDrop - 1),SimParams.exchangeResetInterval) == 0)
    
    xIndex = 0;
    reIterate = 1;
    currentIteration = 0;
    cvx_hist = -500 * ones(2,1);
    maxRank = SimParams.maxRank;
    
    avgArrival = zeros(nUsers,1);
    [pt_o,qt_o,bt_o,Wx] = randomizeInitialSCApoint(SimParams,SimStructs);
    
    vW = cell(nUsers,nBands,nSlots);
    p_o = zeros(maxRank,nUsers,nBands,nSlots);
    q_o = zeros(maxRank,nUsers,nBands,nSlots);
    b_o = zeros(maxRank,nUsers,nBands,nSlots);
    
    for iSlot = 1:nSlots
        for iUser = 1:nUsers
            for iBand = 1:nBands
                vW{iUser,iBand,iSlot} = Wx{iUser,iBand};
            end
        end
        p_o(:,:,:,iSlot) = pt_o;
        q_o(:,:,:,iSlot) = qt_o;
        b_o(:,:,:,iSlot) = bt_o;
    end
    
    for iUser = 1:nUsers
        avgArrival(iUser,1) = SimStructs.userStruct{iUser,1}.trafficConfig.avgArrRate;
    end
    
    while reIterate
        
        cvx_begin
        
        expressions p(maxRank,nUsers,nBands,nSlots) q(maxRank,nUsers,nBands,nSlots)
        variable M(SimParams.nTxAntenna,maxRank,nUsers,nBands,nSlots) complex
        variables t(maxRank,nUsers,nBands,nSlots) b(maxRank,nUsers,nBands,nSlots) g(maxRank,nUsers,nBands,nSlots) B(nUsers,nSlots)
        variables epiObjective
        
        minimize(epiObjective)
        
        subject to
        
        epiObjective >= norm(B(:,nSlots),qExponent);
        
        for iUser = 1:nUsers
            ((QueuedPkts(iUser,1) - sum(vec(t(:,iUser,:,1))))) <= B(iUser,1);
            for iSlot = 2:nSlots
                ((B(iUser,iSlot - 1) - (mdpFactor(iUser,1)^(iSlot - 1)) * sum(vec(t(:,iUser,:,iSlot))))) <= B(iUser,iSlot);
            end
        end
        
        for iSlot = 1:nSlots
            
            for iBase = 1:nBases
                for iBand = 1:nBands
                    for iUser = 1:usersPerCell(iBase,1)
                        
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        for iLayer = 1:maxRank
                            intVector = sqrt(SimParams.N) * vW{cUser,iBand,iSlot}(:,iLayer);
                            for jBase = 1:nBases
                                currentH = cH{jBase,iBand}(:,:,cUser);
                                for jUser = 1:usersPerCell(jBase,1)
                                    rUser = cellUserIndices{jBase,1}(jUser,1);
                                    if rUser ~= cUser
                                        for jLayer = 1:maxRank
                                            intVector = [intVector ; vW{cUser,iBand,iSlot}(:,iLayer)' * currentH * M(:,jLayer,rUser,iBand,iSlot)];
                                        end
                                    else
                                        for jLayer = 1:maxRank
                                            if jLayer ~= iLayer
                                                intVector = [intVector ; vW{cUser,iBand,iSlot}(:,iLayer)' * currentH * M(:,jLayer,rUser,iBand,iSlot)];
                                            end
                                        end
                                    end
                                end
                            end
                            
                            norm(intVector,2) <= sqrt(b(iLayer,cUser,iBand,iSlot));
                            log(1 + g(iLayer,cUser,iBand,iSlot)) >= t(iLayer,cUser,iBand,iSlot) * log(2);
                            
                            currentH = cH{iBase,iBand}(:,:,cUser);
                            p(iLayer,cUser,iBand,iSlot) = real(vW{cUser,iBand,iSlot}(:,iLayer)' * currentH * M(:,iLayer,cUser,iBand,iSlot));
                            q(iLayer,cUser,iBand,iSlot) = imag(vW{cUser,iBand,iSlot}(:,iLayer)' * currentH * M(:,iLayer,cUser,iBand,iSlot));
                            
                            (p_o(iLayer,cUser,iBand,iSlot)^2 + q_o(iLayer,cUser,iBand,iSlot)^2) / (b_o(iLayer,cUser,iBand,iSlot)) + ...
                                (2 / b_o(iLayer,cUser,iBand,iSlot)) * (p_o(iLayer,cUser,iBand,iSlot) * (p(iLayer,cUser,iBand,iSlot) - p_o(iLayer,cUser,iBand,iSlot))) + ...
                                (2 / b_o(iLayer,cUser,iBand,iSlot)) * (q_o(iLayer,cUser,iBand,iSlot) * (q(iLayer,cUser,iBand,iSlot) - q_o(iLayer,cUser,iBand,iSlot))) - ...
                                (p_o(iLayer,cUser,iBand,iSlot)^2 + q_o(iLayer,cUser,iBand,iSlot)^2) / (b_o(iLayer,cUser,iBand,iSlot)^2) * ...
                                (b(iLayer,cUser,iBand,iSlot) - b_o(iLayer,cUser,iBand,iSlot)) >= g(iLayer,cUser,iBand,iSlot);
                            
                        end
                    end
                end
                
                if strcmp(globalMode,'false')
                    for iBand = 1:nBands
                        norm(vec(M(:,:,cellUserIndices{iBase,1},iBand,iSlot)),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
                    end
                else
                    norm(vec(M(:,:,cellUserIndices{iBase,1},:,iSlot)),2) <= sqrt(sum(SimStructs.baseStruct{iBase,1}.sPower(1,:)));
                end
                
            end
        end
        
        cvx_end
        
        display([B',sum(B)']);
        display(epiObjective);
        if strfind(cvx_status,'Solved')
            
            M = full(M);b_o = full(b);
            p_o = full(p);q_o = full(q);
            if min(abs(cvx_optval - cvx_hist)) <= epsilonT
                reIterate = 0;
            else
                xIndex = xIndex + 1;
                cvx_hist(mod(xIndex,2) + 1,1) = cvx_optval;
            end
            
            for iSlot = 1:nSlots
                for iBand = 1:nBands
                    for iBase = 1:nBases
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for iLayer = 1:maxRank
                                R = SimParams.N * eye(SimParams.nRxAntenna);
                                for jBase = 1:nBases
                                    for jUser = 1:usersPerCell(jBase,1)
                                        rUser = cellUserIndices{jBase,1}(jUser,1);
                                        H = cH{jBase,iBand}(:,:,cUser);
                                        R = R + H * M(:,:,rUser,iBand,iSlot) * M(:,:,rUser,iBand,iSlot)' * H';
                                    end
                                end
                                H = cH{iBase,iBand}(:,:,cUser);
                                vW{cUser,iBand,iSlot}(:,iLayer) = R \ (H * M(:,iLayer,cUser,iBand,iSlot));
                            end
                        end
                    end
                end
                for iBand = 1:nBands
                    for iUser = 1:nUsers
                        baseNode = SimStructs.userStruct{iUser,1}.baseNode;
                        channelH = cH{baseNode,iBand}(:,:,iUser);
                        for iLayer = 1:maxRank
                            p_o(iLayer,iUser,iBand,iSlot) = real(vW{iUser,iBand,iSlot}(:,iLayer)' * channelH * M(:,iLayer,iUser,iBand,iSlot));
                            q_o(iLayer,iUser,iBand,iSlot) = imag(vW{iUser,iBand,iSlot}(:,iLayer)' * channelH * M(:,iLayer,iUser,iBand,iSlot));
                        end
                    end
                end
            end
        else
            b_o = b_o * 2;
            display('Failed CVX !');
        end
        
        currentIteration = currentIteration + 1;
        if currentIteration >= nPreExchanges
            reIterate = 0;
        end
    end
    
    SimParams.Debug.resetCounter = 0;
    SimParams.Debug.dataExchange{1,1} = M;
    SimParams.Debug.dataExchange{2,1} = vW;
    SimParams.Debug.dataExchange{3,1} = p_o;
    SimParams.Debug.dataExchange{4,1} = b_o;
    SimParams.Debug.dataExchange{5,1} = cell(nBases,nSlots);
    
    for iSlot = 1:nSlots
        for iBase = 1:nBases
            SimParams.Debug.dataExchange{5,1}{iBase,iSlot} = zeros(maxRank,nUsers,nBands,nSlots);
            for iBand = 1:nBands
                for iUser = 1:length(cellNeighbourIndices{iBase,1})
                    cUser = cellNeighbourIndices{iBase,1}(iUser,1);
                    for iRank = 1:maxRank
                        noiseEnergy = 0;
                        for jUser = 1:length(cellUserIndices{iBase,1})
                            noiseEnergy = noiseEnergy + norm(vW{cUser,iSlot}(:,iRank)' * cH{iBase,iBand}(:,:,cUser) * M(:,:,cellUserIndices{iBase,1}(jUser,1),iBand,iSlot))^2;
                        end
                        SimParams.Debug.dataExchange{5,1}{iBase,iSlot}(iRank,cUser,iBand) = noiseEnergy;
                    end
                end
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = M(:,:,cellUserIndices{iBase,1},iBand,1);
            end
        end
    end
    
    clear R;
    [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE');
    updateIteratePerformance(SimParams,SimStructs);

end

SimParams.Debug.resetCounter = SimParams.Debug.resetCounter + 1;

switch selectionMethod
    
    case 'distMSEAllocB'
        
        stepIndex = 0.25;
        E0 = cell(nBases,1);
        SimParams.currentQueue = 100;
        alphaLKN = cell(nBases,1);lambdaLKN = cell(nBases,1);
        
        cH = SimStructs.linkChan;
        for iExchangeOTA = 0:SimParams.nExchangesOTA
            
            switch iExchangeOTA
                case 0
                    for iBase = 1:nBases
                        for iBand = 1:nBands
                            SimStructs.baseStruct{iBase,1}.P{iBand,1} = SimParams.Debug.dataExchange{1,1}(:,:,cellUserIndices{iBase,1},iBand,SimParams.Debug.resetCounter);
                            SimParams.Debug.globalExchangeInfo.P{iBase,iBand} = SimParams.Debug.dataExchange{1,1}(:,:,cellUserIndices{iBase,1},iBand,SimParams.Debug.resetCounter);
                        end
                        SimParams.Debug.globalExchangeInfo.funcOut{3,iBase} = abs(SimParams.Debug.dataExchange{3,1}(:,cellUserIndices{iBase,1},:,SimParams.Debug.resetCounter));
                        SimParams.Debug.globalExchangeInfo.funcOut{4,iBase} = ones(maxRank,usersPerCell(iBase,1),nBands);
                    end
                    fprintf('OTA Performed - %d \n',iExchangeOTA);
            end
            
            for iBase = 1:nBases
                SimStructs.baseStruct{iBase,1}.selectionType = 'Last';
            end
            [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE');
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
            
            [SimParams,SimStructs] = updateIteratePerformance(SimParams,SimStructs);
            
            if SimParams.currentQueue < epsilonT
                break;
            end
            
        end
        
    case 'distBSAlloc'
        
        cH = SimStructs.linkChan;
        SimParams.currentQueue = 100;
        for iExchangeOTA = 0:SimParams.nExchangesOTA
            
            switch iExchangeOTA
                case 0
                    for iBase = 1:nBases
                        SimStructs.baseStruct{iBase,1}.selectionType = 'Last';
                        for iBand = 1:nBands
                            SimStructs.baseStruct{iBase,1}.P{iBand,1} = SimParams.Debug.dataExchange{1,1}(:,:,cellUserIndices{iBase,1},iBand,SimParams.Debug.resetCounter);
                            SimParams.Debug.globalExchangeInfo.P{iBase,iBand} = SimParams.Debug.dataExchange{1,1}(:,:,cellUserIndices{iBase,1},iBand,SimParams.Debug.resetCounter);
                            SimParams.Debug.globalExchangeInfo.gI{iBase,1}(:,:,iBand) = sqrt(SimParams.Debug.dataExchange{5,1}{iBase,SimParams.Debug.resetCounter}(:,:,iBand));
                            SimParams.Debug.globalExchangeInfo.D{iBase,1} = ones(maxRank,nUsers,nBands,nBases);
                        end
                    end
                    
                    maxBackHaulExchanges = SimParams.nExchangesOBH;
                    if SimParams.Debug.resetCounter == 1
                        break;
                    end
                otherwise
                    [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE');
            end

            stepFactor = 10;
            for iExchangeBH = 1:maxBackHaulExchanges  
                
                for iBase = 1:nBases
                    
                    for iBand = 1:nBands
                        for iUser = 1:SimParams.nUsers
                            W0{iUser,iBand} = SimStructs.userStruct{iUser,1}.pW{iBand,1};
                        end
                    end

                    kUsers = usersPerCell(iBase,1);
                    SimParams.Debug.exchangeIndex = iExchangeBH + iExchangeOTA;
                    [SimParams, SimStructs] = initializeSCApoint(SimParams,SimStructs,iBase);
                    M0 = SimParams.Debug.globalExchangeInfo.funcOut{1,iBase};B0 = SimParams.Debug.globalExchangeInfo.funcOut{2,iBase};
                    
                    cvx_begin
                    
                    variable M(SimParams.nTxAntenna,maxRank,kUsers,nBands) complex
                    variables T(maxRank,kUsers,nBands) B(maxRank,kUsers,nBands) G(maxRank,kUsers,nBands)
                    variables I(maxRank,nUsers,nBands,nBases) userObjective(kUsers,1) epiObjective
                    
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
                        break;
                    end
                    
                    for iBand = 1:nBands
                        SimStructs.baseStruct{iBase,1}.P{iBand,1} = M0(:,:,:,iBand);
                        SimParams.Debug.globalExchangeInfo.P{iBase,iBand} = M0(:,:,:,iBand);
                    end
                    SimParams.Debug.globalExchangeInfo.I{iBase,1} = full(I);
                    
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
        
    case 'centAlloc'
        
        for iBase = 1:nBases
            for iBand = 1:nBands
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = SimParams.Debug.dataExchange{1,1}(:,:,cellUserIndices{iBase,1},iBand,SimParams.Debug.resetCounter);
            end
        end
        
end

% Update the receivers upon reception by all users (combined detection)
[SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE');

