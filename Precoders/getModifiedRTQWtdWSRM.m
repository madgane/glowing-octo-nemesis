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

nPreExchanges = 10;
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
            abs((QueuedPkts(iUser,1) - mdpFactor(iUser,1) * sum(vec(t(:,iUser,:,1))))) <= B(iUser,1);
            for iSlot = 2:nSlots
                abs((B(iUser,iSlot - 1) - (mdpFactor(iUser,1)^iSlot) * sum(vec(t(:,iUser,:,iSlot))))) <= B(iUser,iSlot);
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
    SimParams.Debug.dataExchange{4,1} = q_o;
    SimParams.Debug.dataExchange{5,1} = b_o;
    
    clear R;
    
end

switch selectionMethod
    
    case 'distMSEAllocB'
        
        stepIndex = 0.25;
        E0 = cell(nBases,1);
        SimParams.currentQueue = 100;
        alphaLKN = cell(nBases,1);lambdaLKN = cell(nBases,1);
        
        cH = SimStructs.linkChan;
        SimParams.Debug.resetCounter = SimParams.Debug.resetCounter + 1;
        
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
        
end

% Update the receivers upon reception by all users (combined detection)
[SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE');

