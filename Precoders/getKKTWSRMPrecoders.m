function [SimParams,SimStructs] = getKKTWSRMPrecoders(SimParams,SimStructs)

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
    
    case 'ADMMMethod'
        
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
        
    case 'AP-GP-Approx'
        
        stInstant = 1;
        stepInitialize = 1;
        SimParams.currentQueue = 100;
        M = zeros(SimParams.nTxAntenna,nUsers);
        gammaLKN = zeros(nUsers,1);betaLKN = zeros(nUsers,1);
                
        R = zeros(SimParams.nTxAntenna,SimParams.nTxAntenna,nUsers);
        targetRate = zeros(SimParams.nExchangesOTA,SimParams.nExchangesOBH);

        for iUser = 1:nUsers
            xBase = SimStructs.userStruct{iUser,1}.baseNode;
            M(:,iUser) = cH{xBase,1}(:,:,iUser)' / norm(cH{xBase,1}(:,:,iUser)');
            M(:,iUser) = M(:,iUser) * sqrt(sum(SimStructs.baseStruct{xBase,1}.sPower) / usersPerCell(xBase,1));
        end   

        for iExchangeOTA = stInstant:SimParams.nExchangesOTA
            
            stepSize = stepInitialize;

            for iUser = 1:nUsers
                betaLKN(iUser,1) = SimParams.N;
                for jUser = 1:nUsers
                    if iUser ~= jUser
                        xBase = SimStructs.userStruct{jUser,1}.baseNode;
                        betaLKN(iUser,1) = betaLKN(iUser,1) + abs(cH{xBase,1}(:,:,iUser) * M(:,jUser))^2;
                    end
                end
            end
            
            for iUser = 1:nUsers
                xBase = SimStructs.userStruct{iUser,1}.baseNode;
                gammaLKN(iUser,1) = abs(cH{xBase,1}(:,:,iUser) * M(:,iUser))^2 / betaLKN(iUser,1);
            end

            phiLKN = sqrt(gammaLKN ./ betaLKN);            
            alphaLKN = ones(nUsers,1) * 0.1;
            deltaLKN = 0.5 * sqrt(gammaLKN ./ betaLKN) .* alphaLKN;
            
            for iExchangeBH = 1:SimParams.nExchangesOBH
                
                stepSize = stepSize * 0.9;
                                                
                for iUser = 1:nUsers
                    R(:,:,iUser) = zeros(SimParams.nTxAntenna,SimParams.nTxAntenna);
                    xBase = SimStructs.userStruct{iUser,1}.baseNode;
                    for jUser = 1:nUsers
                        if iUser ~= jUser
                            R(:,:,iUser) = R(:,:,iUser) + cH{xBase,1}(:,:,jUser)' * cH{xBase,1}(:,:,jUser) * deltaLKN(jUser,1);
                        end
                    end
                end
                
                for iBase = 1:nBases
                    muMax = 100000;
                    muMin = 0;
                    iterateAgain = 1;
                    while iterateAgain
                        totalPower = 0;
                        currentMu = (muMax + muMin) / 2;
                        for iBand = 1:nBands
                            for iUser = 1:usersPerCell(iBase,1)
                                cUser = cellUserIndices{iBase,1}(iUser,1);
                                M(:,cUser) = 0.5 * alphaLKN(cUser,1) * (pinv(currentMu * eye(SimParams.nTxAntenna) + R(:,:,iUser)) * (cH{iBase,1}(:,:,cUser)'));
                                totalPower = totalPower + real(trace(M(:,iUser) * M(:,iUser)'));
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
                end
                
                for iUser = 1:nUsers
                    betaLKN(iUser,1) = SimParams.N;
                    for jUser = 1:nUsers
                        if iUser ~= jUser
                            xBase = SimStructs.userStruct{jUser,1}.baseNode;
                            betaLKN(iUser,1) = betaLKN(iUser,1) + abs(cH{xBase,1}(:,:,iUser) * M(:,jUser))^2;
                        end
                    end
                end

                for iUser = 1:nUsers
                    xBase = SimStructs.userStruct{iUser,1}.baseNode;
                    bracketTerm = (1 / (2 * phiLKN(iUser,1))) * gammaLKN(iUser,1) + 0.5 * phiLKN(iUser,1) * betaLKN(iUser,1) - abs(cH{xBase,1}(:,:,iUser) * M(:,iUser));
                    alphaLKN(iUser,1) = alphaLKN(iUser,1) + stepSize * bracketTerm;
                    deltaLKN(iUser,1) = alphaLKN(iUser,1) * phiLKN(iUser,1) * 0.5;
                end
                
                alphaLKN = (alphaLKN > 0) .* alphaLKN;
                deltaLKN = (deltaLKN > 0) .* deltaLKN;
                
            end
            
            gammaLKN = (2 * (phiLKN./alphaLKN)) - 1;
            gammaLKN = (gammaLKN > 0) .* gammaLKN;
            
            for iUser = 1:nUsers
                targetRate(iExchangeOTA,iExchangeBH) = targetRate(iExchangeOTA,iExchangeBH) + (log2(1 + gammaLKN(iUser,1)));
            end
            
        end
        

end

% Update the receivers upon reception by all users (combined detection)
[SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE');


