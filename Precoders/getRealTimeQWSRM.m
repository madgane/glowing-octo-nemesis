function [SimParams,SimStructs] = getRealTimeQWSRM(SimParams,SimStructs)

preLogue;

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
        
        maxIterations = 50;  
        for iBase = 1:nBases
            SimParams.Debug.globalExchangeInfo.gI{iBase,1} = zeros(maxRank,nUsers,nBands);
        end
        
        for exchangeOTA = 1:1

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
                    
                    [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE',iBase);
                end
                
                fprintf('\n');
            end
            
            [SimParams,SimStructs] = updateIteratePerformance(SimParams,SimStructs);
        end
        
    case 'distBSAlloc'
        
        stepIndex = 10;
        for iExchangeOTA = 1:SimParams.nExchangesOTA
            
            if iExchangeOTA == 1
                for iBase = 1:nBases
                    SimStructs.baseStruct{iBase,1}.selectionType = 'BF';
                    [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'BF',iBase);
                    SimParams.Debug.globalExchangeInfo.gI{iBase,1} = zeros(maxRank,nUsers,nBands);
                    SimParams.Debug.globalExchangeInfo.D{iBase,1} = zeros(maxRank,nUsers,nBands,nBases);                                       
                end
                cH = SimStructs.prevChan;
                maxBackHaulExchanges = SimParams.nExchangesOBH;
            else
                cH = SimStructs.linkChan;
                maxBackHaulExchanges = 1;
            end
            
            for iExchangeBH = 1:maxBackHaulExchanges
            
                for iBase = 1:nBases
                    
                    if ~and(iExchangeBH == 1,iExchangeOTA == 1)
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
            
            fprintf('OTA Performed - %d \n',iExchangeOTA);
            [SimParams,SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE');
            
        end
        
end

end