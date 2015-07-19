
function [SimParams,SimStructs] = getQWtdWSRMDist(SimParams,SimStructs)

proLogue;
alpha = 10;
sumDeviationH = -50;

mIterationsSG = SimParams.nExchangesOBH;
mIterationsSCA = SimParams.nExchangesOTA;
SimParams.distDecompSteps = mIterationsSG;

switch selectionMethod
    
    case 'PrimalMethod'
        
        nLayers = SimParams.maxRank;
        cellP = cell(nBases,1);cellQ = cell(nBases,1);cellB = cell(nBases,1);
        cellM = cell(nBases,1);cellD = cell(nBases,1);cellBH = cell(nBases,1);
        
        xIteration = 0;
        scaContinue = 1;
        currentIF = zeros(nLayers,nUsers,nBases,nBands);
        [p_o,q_o,b_o,W] = randomizeInitialSCApoint(SimParams,SimStructs);
        
        for iBase = 1:nBases
            cellP{iBase,1} = p_o(:,cellUserIndices{iBase,1},:);
            cellQ{iBase,1} = q_o(:,cellUserIndices{iBase,1},:);
            cellB{iBase,1} = b_o(:,cellUserIndices{iBase,1},:);
        end
        
        if SimParams.distIteration == 1
            for iBase = 1:nBases
                currentIF(:,:,iBase,:) = b_o / nBases;
            end
        else
            currentIF = SimParams.Debug.DataExchange{2,1};
        end
        
        while scaContinue
            
            yIteration = 0;
            masterContinue = 1;
            xIteration = xIteration + 1;
            if xIteration >= mIterationsSCA
                scaContinue = 0;
            end
            
            while masterContinue
                
                yIteration = yIteration + 1;
                if yIteration >= mIterationsSG
                    masterContinue = 0;
                end
                
                for iBase = 1:nBases
                    
                    kUsers = usersPerCell(iBase,1);
                    
                    cvx_begin
                    
                    dual variables dualD{nLayers,nUsers,nBands}
                    expressions p(nLayers,kUsers,nBands) q(nLayers,kUsers,nBands)
                    variable M(SimParams.nTxAntenna,nLayers,kUsers,nBands) complex
                    variables t(nLayers,kUsers,nBands) b(nLayers,kUsers,nBands) g(nLayers,kUsers,nBands)
                    variables userObjective(kUsers,1) epiObjective
                    
                    minimize(epiObjective)
                    
                    subject to
                    
                    for iUser = 1:kUsers
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        abs(QueuedPkts(cUser,1) - sum(vec(t(:,iUser,:)))) <= userObjective(iUser,1);
                    end
                    
                    epiObjective >= norm(userObjective,qExponent);
                    
                    for iBand = 1:nBands
                        
                        for iUser = 1:kUsers
                            
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            H = cH{iBase,iBand}(:,:,cUser);
                            
                            for iLayer = 1:nLayers
                                
                                intVector = sqrt(SimParams.N) * W{cUser,iBand}(:,iLayer);
                                
                                for jBase = 1:nBases
                                    if jBase ~= iBase
                                        intVector = [intVector ; sqrt(currentIF(iLayer,cUser,jBase,iBand))];
                                    end
                                end
                                
                                for jUser = 1:kUsers
                                    if jUser ~= iUser
                                        for jLayer = 1:nLayers
                                            intVector = [intVector ; W{cUser,iBand}(:,iLayer)' * H * M(:,jLayer,jUser,iBand)];
                                        end
                                    else
                                        for jLayer = 1:nLayers
                                            if jLayer ~= iLayer
                                                intVector = [intVector ; W{cUser,iBand}(:,iLayer)' * H * M(:,jLayer,jUser,iBand)];
                                            end
                                        end
                                    end
                                end
                                
                                dualD{iLayer,cUser,iBand} : norm(intVector,2) <= sqrt(b(iLayer,iUser,iBand));
                                log(1 + g(iLayer,iUser,iBand)) >= t(iLayer,iUser,iBand) * log(2);
                                
                                p(iLayer,iUser,iBand) = real(W{cUser,iBand}(:,iLayer)' * H * M(:,iLayer,iUser,iBand));
                                q(iLayer,iUser,iBand) = imag(W{cUser,iBand}(:,iLayer)' * H * M(:,iLayer,iUser,iBand));
                                
                                (cellP{iBase,1}(iLayer,iUser,iBand)^2 + cellQ{iBase,1}(iLayer,iUser,iBand)^2) / (cellB{iBase,1}(iLayer,iUser,iBand)) + ...
                                    (2 / cellB{iBase,1}(iLayer,iUser,iBand)) * (cellP{iBase,1}(iLayer,iUser,iBand) * (p(iLayer,iUser,iBand) - cellP{iBase,1}(iLayer,iUser,iBand))) + ...
                                    (2 / cellB{iBase,1}(iLayer,iUser,iBand)) * (cellQ{iBase,1}(iLayer,iUser,iBand) * (q(iLayer,iUser,iBand) - cellQ{iBase,1}(iLayer,iUser,iBand))) - ...
                                    (cellP{iBase,1}(iLayer,iUser,iBand)^2 + cellQ{iBase,1}(iLayer,iUser,iBand)^2) / (cellB{iBase,1}(iLayer,iUser,iBand)^2) * ...
                                    (b(iLayer,iUser,iBand) - cellB{iBase,1}(iLayer,iUser,iBand)) >= g(iLayer,iUser,iBand);
                                
                            end
                            
                        end
                        
                        for iUser = 1:length(cellNeighbourIndices{iBase,1})
                            cUser = cellNeighbourIndices{iBase,1}(iUser,1);
                            for iLayer = 1:nLayers
                                intVector = [];
                                H = cH{iBase,iBand}(:,:,cUser);
                                for jUser = 1:kUsers
                                    for jLayer = 1:nLayers
                                        intVector = [intVector ; W{cUser,iBand}(:,iLayer)' * H * M(:,jLayer,jUser,iBand)];
                                    end
                                end
                                dualD{iLayer,cUser,iBand} : norm(intVector,2) <= sqrt(currentIF(iLayer,cUser,iBase,iBand));
                            end
                        end
                        
                    end
                    
                    if strcmp(globalMode,'false')
                        for iBand = 1:nBands
                            norm(vec(M(:,:,:,iBand)),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
                        end
                    else
                        norm(vec(M),2) <= sqrt(sum(SimStructs.baseStruct{iBase,1}.sPower(1,:)));
                    end
                    
                    cvx_end
                    
                    if iBase == 1
                        status = cvx_status;
                    else
                        status = strcat(status,'-',cvx_status);
                    end
                    
                    if strfind(cvx_status,'Solved')
                        cellM{iBase,1} = M;
                        cellBH{iBase,1} = b;
                        cellD{iBase,1} = dualD;
                    else
                        cellM{iBase,1} = zeros(size(M));
                    end
                    
                end
                
                currentIFH = currentIF;
                currentIF = zeros(nLayers,nUsers,nBases,nBands);
                [SimParams,SimStructs] = updateIteratePerformance(SimParams,SimStructs,cellM,W);
                
                for iBand = 1:nBands
                    for iBase = 1:nBases
                        for iUser = 1:length(cellNeighbourIndices{iBase,1})
                            for iLayer = 1:nLayers
                                cUser = cellNeighbourIndices{iBase,1}(iUser,1);
                                baseNode = SimStructs.userStruct{cUser,1}.baseNode;
                                currentIF(iLayer,cUser,iBase,iBand) = currentIFH(iLayer,cUser,iBase,iBand) - alpha * (cellD{baseNode,1}{iLayer,cUser,iBand} - cellD{iBase,1}{iLayer,cUser,iBand});
                            end
                        end
                    end
                end
                
                if strcmp(SimParams.DebugMode,'true')
                    status = strcat(status,'-',sprintf('%d',yIteration),'-',sprintf('%d',xIteration));
                    display(status);
                    %display(currentIF);
                end
                
                currentIF = max(currentIF,0);
                if norm(vec(currentIF - currentIFH),2) <= epsilonT
                    masterContinue = 0;
                end
                
            end
            
            sumDeviation = sum(cell2mat(SimParams.Debug.tempResource{3,SimParams.iDrop}));
            if abs(sumDeviation(1,end) - sumDeviationH) < epsilonT
                scaContinue = 0;
            else
                sumDeviationH = sumDeviation(1,end);
            end
            
            for iBase = 1:nBases
                cellB{iBase,1} = cellBH{iBase,1};
                for iBand = 1:nBands
                    for iUser = 1:usersPerCell(iBase,1)
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        for iLayer = 1:nLayers
                            cellP{iBase,1}(iLayer,iUser,iBand) = real(W{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * cellM{iBase,1}(:,iLayer,iUser,iBand));
                            cellQ{iBase,1}(iLayer,iUser,iBand) = imag(W{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * cellM{iBase,1}(:,iLayer,iUser,iBand));
                        end
                    end
                end
            end
            
            for iBand = 1:nBands
                for iUser = 1:nUsers
                    baseNode = SimStructs.userStruct{iUser,1}.baseNode;
                    for iLayer = 1:nLayers
                        R = SimParams.N * eye(SimParams.nRxAntenna);
                        for iBase = 1:nBases
                            for jUser = 1:usersPerCell(iBase,1)
                                H = cH{iBase,iBand}(:,:,iUser);
                                M = cellM{iBase,1}(:,:,jUser,iBand);
                                R = R + H * (M * M') * H';
                            end
                        end
                        H = cH{baseNode,iBand}(:,:,iUser);
                        xUser = (iUser == cellUserIndices{baseNode,1});
                        W{iUser,iBand}(:,iLayer) = R \ (H * cellM{baseNode,1}(:,iLayer,xUser,iBand));
                    end
                end
            end
            
        end
        
        for iBand = 1:nBands
            for iBase = 1:nBases
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = cellM{iBase,1}(:,:,:,iBand);
            end
        end
        
        SimParams.Debug.DataExchange{1,1} = cellM;
        SimParams.Debug.DataExchange{2,1} = currentIF;
        
    case 'ADMMMethod'
        
        nLayers = SimParams.maxRank;
        cellP = cell(nBases,1);cellQ = cell(nBases,1);cellB = cell(nBases,1);
        cellM = cell(nBases,1);cellX = cell(nBases,1);cellBH = cell(nBases,1);
        
        xIteration = 0;
        scaContinue = 1;
        currentDual = cell(nBases,1);
        cellXGlobal = zeros(nLayers,nUsers,nBases,nBands);
        [p_o,q_o,b_o,W] = randomizeInitialSCApoint(SimParams,SimStructs);
        
        for iBase = 1:nBases
            cellP{iBase,1} = p_o(:,cellUserIndices{iBase,1},:);
            cellQ{iBase,1} = q_o(:,cellUserIndices{iBase,1},:);
            cellB{iBase,1} = b_o(:,cellUserIndices{iBase,1},:);
            cellX{iBase,1} = zeros(nLayers,nUsers,nBases,nBands);
        end
        
        if SimParams.distIteration == 1
            for iBase = 1:nBases
                cellXGlobal(:,:,iBase,:) = zeros(nLayers,nUsers,nBands);
                currentDual{iBase,1} = zeros(nLayers,nUsers,nBases,nBands);
            end
        else
            cellXGlobal = SimParams.Debug.DataExchange{2,1};
            currentDual = SimParams.Debug.DataExchange{3,1};
        end
        
        while scaContinue
            
            
            yIteration = 0;
            masterContinue = 1;
            
            xIteration = xIteration + 1;
            if xIteration >= mIterationsSCA
                scaContinue = 0;
            end
            
            stepSize = alpha;
            while masterContinue
                
                stepSize = stepSize * 0.9;
                yIteration = yIteration + 1;
                if yIteration >= mIterationsSG
                    masterContinue = 0;
                end
                
                for iBase = 1:nBases
                    
                    kUsers = usersPerCell(iBase,1);
                    
                    cvx_begin
                    
                    expressions p(nLayers,kUsers,nBands) q(nLayers,kUsers,nBands) tempFirst tempSecond tempADMM
                    variable M(SimParams.nTxAntenna,nLayers,kUsers,nBands) complex
                    variables t(nLayers,kUsers,nBands) b(nLayers,kUsers,nBands) g(nLayers,kUsers,nBands) x(nLayers,nUsers,nBases,nBands)
                    variables userObjective(kUsers,1) epiObjective
                    
                    minimize(epiObjective)
                    
                    subject to
                    
                    for iUser = 1:kUsers
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        abs(QueuedPkts(cUser,1) - sum(vec(t(:,iUser,:)))) <= userObjective(iUser,1);
                    end
                    
                    tempFirst = 0;tempSecond = 0;tempADMM = 0;
                    
                    for iBand = 1:nBands
                        
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for jBase = 1:nBases
                                if jBase ~= iBase
                                    tempFirst = tempFirst + sum(currentDual{iBase,1}(:,cUser,jBase,iBand) .* x(:,cUser,jBase,iBand));
                                    tempADMM = tempADMM + sum(vec(x(:,cUser,jBase,iBand) - cellXGlobal(:,cUser,jBase,iBand)).^2);
                                end
                            end
                        end
                        
                        for iUser = 1:length(cellNeighbourIndices{iBase,1})
                            cUser = cellNeighbourIndices{iBase,1}(iUser,1);
                            tempSecond = tempSecond + sum(currentDual{iBase,1}(:,cUser,iBase,iBand) .* x(:,cUser,iBase,iBand));
                            tempADMM = tempADMM + sum(vec(x(:,cUser,iBase,iBand) - cellXGlobal(:,cUser,iBase,iBand)).^2);
                        end
                        
                    end
                    
                    epiObjective >= norm(userObjective,qExponent) + tempFirst + tempSecond + tempADMM * (stepSize / 2);
                    
                    for iBand = 1:nBands
                        
                        for iUser = 1:kUsers
                            
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            H = cH{iBase,iBand}(:,:,cUser);
                            
                            for iLayer = 1:nLayers
                                
                                intVector = sqrt(SimParams.N) * W{cUser,iBand}(:,iLayer);
                                
                                for jBase = 1:nBases
                                    if jBase ~= iBase
                                        intVector = [intVector ; x(iLayer,cUser,jBase,iBand)];
                                    end
                                end
                                
                                for jUser = 1:kUsers
                                    if jUser ~= iUser
                                        for jLayer = 1:nLayers
                                            intVector = [intVector ; W{cUser,iBand}(:,iLayer)' * H * M(:,jLayer,jUser,iBand)];
                                        end
                                    else
                                        for jLayer = 1:nLayers
                                            if jLayer ~= iLayer
                                                intVector = [intVector ; W{cUser,iBand}(:,iLayer)' * H * M(:,jLayer,jUser,iBand)];
                                            end
                                        end
                                    end
                                end
                                
                                norm(intVector,2) <= sqrt(b(iLayer,iUser,iBand));
                                log(1 + g(iLayer,iUser,iBand)) >= t(iLayer,iUser,iBand) * log(2);
                                
                                p(iLayer,iUser,iBand) = real(W{cUser,iBand}(:,iLayer)' * H * M(:,iLayer,iUser,iBand));
                                q(iLayer,iUser,iBand) = imag(W{cUser,iBand}(:,iLayer)' * H * M(:,iLayer,iUser,iBand));
                                
                                (cellP{iBase,1}(iLayer,iUser,iBand)^2 + cellQ{iBase,1}(iLayer,iUser,iBand)^2) / (cellB{iBase,1}(iLayer,iUser,iBand)) + ...
                                    (2 / cellB{iBase,1}(iLayer,iUser,iBand)) * (cellP{iBase,1}(iLayer,iUser,iBand) * (p(iLayer,iUser,iBand) - cellP{iBase,1}(iLayer,iUser,iBand))) + ...
                                    (2 / cellB{iBase,1}(iLayer,iUser,iBand)) * (cellQ{iBase,1}(iLayer,iUser,iBand) * (q(iLayer,iUser,iBand) - cellQ{iBase,1}(iLayer,iUser,iBand))) - ...
                                    (cellP{iBase,1}(iLayer,iUser,iBand)^2 + cellQ{iBase,1}(iLayer,iUser,iBand)^2) / (cellB{iBase,1}(iLayer,iUser,iBand)^2) * ...
                                    (b(iLayer,iUser,iBand) - cellB{iBase,1}(iLayer,iUser,iBand)) >= g(iLayer,iUser,iBand);
                                
                            end
                            
                        end
                        
                        for iUser = 1:length(cellNeighbourIndices{iBase,1})
                            cUser = cellNeighbourIndices{iBase,1}(iUser,1);
                            for iLayer = 1:nLayers
                                intVector = [];
                                H = cH{iBase,iBand}(:,:,cUser);
                                for jUser = 1:kUsers
                                    for jLayer = 1:nLayers
                                        intVector = [intVector ; W{cUser,iBand}(:,iLayer)' * H * M(:,jLayer,jUser,iBand)];
                                    end
                                end
                                norm(intVector,2) <= x(iLayer,cUser,iBase,iBand);
                            end
                        end
                        
                    end
                    
                    if strcmp(globalMode,'false')
                        for iBand = 1:nBands
                            norm(vec(M(:,:,:,iBand)),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
                        end
                    else
                        norm(vec(M),2) <= sqrt(sum(SimStructs.baseStruct{iBase,1}.sPower(1,:)));
                    end
                    
                    cvx_end
                    
                    if iBase == 1
                        status = cvx_status;
                    else
                        status = strcat(status,'-',cvx_status);
                    end
                    
                    if strfind(cvx_status,'Solved')
                        cellM{iBase,1} = M;
                        cellBH{iBase,1} = b;
                        cellX{iBase,1} = x;
                    else
                        display('Not solved !');
                        cellM{iBase,1} = zeros(size(M));
                    end
                end
                
                [SimParams,SimStructs] = updateIteratePerformance(SimParams,SimStructs,cellM,W);
                
                for iBand = 1:nBands
                    for iBase = 1:nBases
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for jBase = 1:nBases
                                if jBase ~= iBase
                                    cellXGlobal(:,cUser,jBase,iBand) = (cellX{iBase,1}(:,cUser,jBase,iBand) + cellX{jBase,1}(:,cUser,jBase,iBand)) / 2;
                                end
                            end
                        end
                    end
                end
                
                for iBand = 1:nBands
                    for iBase = 1:nBases
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for jBase = 1:nBases
                                if jBase ~= iBase
                                    currentDual{iBase,1}(:,cUser,jBase,iBand) = currentDual{iBase,1}(:,cUser,jBase,iBand) + stepSize * (cellX{iBase,1}(:,cUser,jBase,iBand) - cellXGlobal(:,cUser,jBase,iBand));
                                end
                            end
                        end
                        for iUser = 1:length(cellNeighbourIndices{iBase,1})
                            cUser = cellNeighbourIndices{iBase,1}(iUser,1);
                            currentDual{iBase,1}(:,cUser,iBase,iBand) = currentDual{iBase,1}(:,cUser,iBase,iBand) + stepSize * (cellX{iBase,1}(:,cUser,iBase,iBand) - cellXGlobal(:,cUser,iBase,iBand));
                        end
                    end
                end
                
                if strcmp(SimParams.DebugMode,'true')
                    status = strcat(status,'-',sprintf('%d',yIteration),'-',sprintf('%d',xIteration));
                    display(status);
                end
            end
            
            sumDeviation = sum(cell2mat(SimParams.Debug.tempResource{3,SimParams.iDrop}));
            if abs(sumDeviation(1,end) - sumDeviationH) < epsilonT
                scaContinue = 0;
            else
                sumDeviationH = sumDeviation(1,end);
            end
            
            for iBase = 1:nBases
                cellB{iBase,1} = cellBH{iBase,1};
                for iBand = 1:nBands
                    for iUser = 1:usersPerCell(iBase,1)
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        for iLayer = 1:nLayers
                            cellP{iBase,1}(iLayer,iUser,iBand) = real(W{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * cellM{iBase,1}(:,iLayer,iUser,iBand));
                            cellQ{iBase,1}(iLayer,iUser,iBand) = imag(W{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * cellM{iBase,1}(:,iLayer,iUser,iBand));
                        end
                    end
                end
            end
            
            for iBand = 1:nBands
                for iUser = 1:nUsers
                    baseNode = SimStructs.userStruct{iUser,1}.baseNode;
                    for iLayer = 1:nLayers
                        R = SimParams.N * eye(SimParams.nRxAntenna);
                        for iBase = 1:nBases
                            for jUser = 1:usersPerCell(iBase,1)
                                H = cH{iBase,iBand}(:,:,iUser);
                                M = cellM{iBase,1}(:,:,jUser,iBand);
                                R = R + H * (M * M') * H';
                            end
                        end
                        H = cH{baseNode,iBand}(:,:,iUser);
                        xUser = (iUser == cellUserIndices{baseNode,1});
                        W{iUser,iBand}(:,iLayer) = R \ (H * cellM{baseNode,1}(:,iLayer,xUser,iBand));
                    end
                end
            end
            
        end
        
        SimParams.Debug.DataExchange{1,1} = cellM;
        SimParams.Debug.DataExchange{2,1} = cellXGlobal;
        SimParams.Debug.DataExchange{3,1} = currentDual;
        
        for iBand = 1:nBands
            for iBase = 1:nBases
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = cellM{iBase,1}(:,:,:,iBand);
            end
        end
        
    case 'PrimalMSEMethod'
        
        nLayers = SimParams.maxRank;
        cellD = cell(nBases,1);cellM = cell(nBases,1);cellTH = cell(nBases,1);
        
        xIteration = 0;
        scaContinue = 1;
        currentIF = zeros(nLayers,nUsers,nBases,nBands);
        [initialMSE,W,currentF] = randomizeInitialMSESCApoint(SimParams,SimStructs);
        
        for iBase = 1:nBases
            cellTH{iBase,1} = initialMSE(:,cellUserIndices{iBase,1},:);
        end
        
        if SimParams.distIteration == 1
            for iBase = 1:nBases
                currentIF(:,:,iBase,:) = currentF / nBases;
            end
        else
            currentIF = SimParams.Debug.DataExchange{2,1};
        end
        
        while scaContinue
            
            yIteration = 0;
            masterContinue = 1;
            xIteration = xIteration + 1;
            if xIteration >= mIterationsSCA
                scaContinue = 0;
            end
            
            cellT = cellTH;
            while masterContinue
                
                yIteration = yIteration + 1;
                if yIteration >= mIterationsSG
                    masterContinue = 0;
                end
                
                for iBase = 1:nBases
                    
                    mseError_o = cellT{iBase,1};
                    kUsers = usersPerCell(iBase,1);
                    
                    cvx_begin
                    
                    dual variables dualD{nLayers,nUsers,nBands}
                    expressions p(nLayers,kUsers,nBands) q(nLayers,kUsers,nBands)
                    variable M(SimParams.nTxAntenna,nLayers,kUsers,nBands) complex
                    variables t(nLayers,kUsers,nBands) b(nLayers,kUsers,nBands) g(nLayers,kUsers,nBands) mseError(nLayers,kUsers,nBands)
                    variables userObjective(kUsers,1) epiObjective
                    
                    minimize(epiObjective)
                    
                    subject to
                    
                    for iUser = 1:kUsers
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        abs(QueuedPkts(cUser,1) - sum(vec(t(:,iUser,:)))) <= userObjective(iUser,1);
                    end
                    
                    epiObjective >= norm(userObjective,qExponent);
                    
                    for iBand = 1:nBands
                        
                        for iUser = 1:kUsers
                            
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            H = cH{iBase,iBand}(:,:,cUser);
                            
                            for iLayer = 1:nLayers
                                
                                intVector = sqrt(SimParams.N) * W{cUser,iBand}(:,iLayer);
                                
                                for jBase = 1:nBases
                                    if jBase ~= iBase
                                        intVector = [intVector ; sqrt(currentIF(iLayer,cUser,jBase,iBand))];
                                    end
                                end
                                
                                for jUser = 1:kUsers
                                    if jUser ~= iUser
                                        for jLayer = 1:nLayers
                                            intVector = [intVector ; W{cUser,iBand}(:,iLayer)' * H * M(:,jLayer,jUser,iBand)];
                                        end
                                    else
                                        for jLayer = 1:nLayers
                                            if jLayer ~= iLayer
                                                intVector = [intVector ; W{cUser,iBand}(:,iLayer)' * H * M(:,jLayer,jUser,iBand)];
                                            end
                                        end
                                    end
                                end
                                
                                givenVector = (1 - W{cUser,iBand}(:,iLayer)' * H * M(:,iLayer,iUser,iBand));
                                intVector = [intVector ; givenVector];
                                
                                dualD{iLayer,cUser,iBand} : norm(intVector,2) <= sqrt(mseError(iLayer,iUser,iBand));
                                (mseError(iLayer,iUser,iBand) - mseError_o(iLayer,iUser,iBand)) / mseError_o(iLayer,iUser,iBand) + log(mseError_o(iLayer,iUser,iBand)) <= -t(iLayer,iUser,iBand) * log(2);
                                
                            end
                            
                        end
                        
                        for iUser = 1:length(cellNeighbourIndices{iBase,1})
                            cUser = cellNeighbourIndices{iBase,1}(iUser,1);
                            for iLayer = 1:nLayers
                                intVector = [];
                                H = cH{iBase,iBand}(:,:,cUser);
                                for jUser = 1:kUsers
                                    for jLayer = 1:nLayers
                                        intVector = [intVector ; W{cUser,iBand}(:,iLayer)' * H * M(:,jLayer,jUser,iBand)];
                                    end
                                end
                                dualD{iLayer,cUser,iBand} : norm(intVector,2) <= sqrt(currentIF(iLayer,cUser,iBase,iBand));
                            end
                        end
                        
                    end
                    
                    if strcmp(globalMode,'false')
                        for iBand = 1:nBands
                            norm(vec(M(:,:,:,iBand)),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
                        end
                    else
                        norm(vec(M(:,:,:,:)),2) <= sqrt(sum(SimStructs.baseStruct{iBase,1}.sPower(1,:)));
                    end
                    
                    cvx_end
                    
                    if iBase == 1
                        status = cvx_status;
                    else
                        status = strcat(status,'-',cvx_status);
                    end
                    
                    if strfind(cvx_status,'Solved')
                        
                        cellM{iBase,1} = M;
                        cellD{iBase,1} = dualD;
                        cellTH{iBase,1} = mseError;
                        
                    else
                        
                        cellM{iBase,1} = zeros(size(M));
                        
                    end
                    
                end
                
                currentIFH = currentIF;
                currentIF = zeros(nLayers,nUsers,nBases,nBands);
                [SimParams,SimStructs] = updateIteratePerformance(SimParams,SimStructs,cellM,W);
                
                for iBand = 1:nBands
                    for iBase = 1:nBases
                        for iUser = 1:length(cellNeighbourIndices{iBase,1})
                            for iLayer = 1:nLayers
                                cUser = cellNeighbourIndices{iBase,1}(iUser,1);
                                baseNode = SimStructs.userStruct{cUser,1}.baseNode;
                                currentIF(iLayer,cUser,iBase,iBand) = currentIFH(iLayer,cUser,iBase,iBand) - alpha * (cellD{baseNode,1}{iLayer,cUser,iBand} - cellD{iBase,1}{iLayer,cUser,iBand});
                            end
                        end
                    end
                end
                
                if strcmp(SimParams.DebugMode,'true')
                    status = strcat(status,'-',sprintf('%d',yIteration),'-',sprintf('%d',xIteration));
                    display(status);
                    %display(currentIF);
                end
                currentIF = max(currentIF,0);
                if norm(vec(currentIF - currentIFH),2) <= epsilonT
                    masterContinue = 0;
                end
                
            end
            
            sumDeviation = sum(cell2mat(SimParams.Debug.tempResource{3,SimParams.iDrop}));
            if abs(sumDeviation(1,end) - sumDeviationH) < epsilonT
                scaContinue = 0;
            else
                sumDeviationH = sumDeviation(1,end);
            end
            
            for iBand = 1:nBands
                for iUser = 1:nUsers
                    baseNode = SimStructs.userStruct{iUser,1}.baseNode;
                    for iLayer = 1:nLayers
                        R = SimParams.N * eye(SimParams.nRxAntenna);
                        for iBase = 1:nBases
                            for jUser = 1:usersPerCell(iBase,1)
                                H = cH{iBase,iBand}(:,:,iUser);
                                M = cellM{iBase,1}(:,:,jUser,iBand);
                                R = R + H * (M * M') * H';
                            end
                        end
                        H = cH{baseNode,iBand}(:,:,iUser);
                        xUser = (iUser == cellUserIndices{baseNode,1});
                        W{iUser,iBand}(:,iLayer) = R \ (H * cellM{baseNode,1}(:,iLayer,xUser,iBand));
                    end
                end
            end
            
        end
        
        SimParams.Debug.DataExchange{1,1} = cellM;
        SimParams.Debug.DataExchange{2,1} = currentIF;
        
        for iBand = 1:nBands
            for iBase = 1:nBases
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = cellM{iBase,1}(:,:,:,iBand);
            end
        end
        
    case 'ADMMMSEMethod'
        
        nLayers = SimParams.maxRank;
        cellM = cell(nBases,1);cellX = cell(nBases,1);cellBH = cell(nBases,1);
        
        xIteration = 0;
        scaContinue = 1;
        currentDual = cell(nBases,1);
        [initialMSE,W,~] = randomizeInitialMSESCApoint(SimParams,SimStructs);
        
        for iBase = 1:nBases
            cellX{iBase,1} = ones(nLayers,nUsers,nBases,nBands);
            cellBH{iBase,1} = initialMSE(:,cellUserIndices{iBase,1},:);
        end
        
        if SimParams.distIteration == 1
            for iBase = 1:nBases
                cellXGlobal = zeros(nLayers,nUsers,nBases,nBands);
                currentDual{iBase,1} = zeros(nLayers,nUsers,nBases,nBands);
            end
        else
            cellXGlobal = SimParams.Debug.DataExchange{2,1};
            currentDual = SimParams.Debug.DataExchange{3,1};
        end
        
        while scaContinue
            
            yIteration = 0;
            masterContinue = 1;
            xIteration = xIteration + 1;
            if xIteration >= mIterationsSCA
                scaContinue = 0;
            end
            
            cellMSE = cellBH;
            while masterContinue
                
                yIteration = yIteration + 1;
                if yIteration >= mIterationsSG
                    masterContinue = 0;
                end
                
                for iBase = 1:nBases
                    
                    mseError_o = cellMSE{iBase,1};
                    kUsers = usersPerCell(iBase,1);
                    
                    cvx_begin
                    
                    expressions p(nLayers,kUsers,nBands) q(nLayers,kUsers,nBands) tempFirst tempSecond tempADMM
                    variable M(SimParams.nTxAntenna,nLayers,kUsers,nBands) complex
                    variables t(nLayers,kUsers,nBands) b(nLayers,kUsers,nBands) g(nLayers,kUsers,nBands) x(nLayers,nUsers,nBases,nBands) mseError(nLayers,nUsers,nBands)
                    variables userObjective(kUsers,1) epiObjective
                    
                    minimize(epiObjective)
                    
                    subject to
                    
                    for iUser = 1:kUsers
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        abs(QueuedPkts(cUser,1) - sum(vec(t(:,iUser,:)))) <= userObjective(iUser,1);
                    end
                    
                    tempFirst = 0;tempSecond = 0;tempADMM = 0;
                    
                    for iBand = 1:nBands
                        
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for jBase = 1:nBases
                                if jBase ~= iBase
                                    tempFirst = tempFirst + sum(currentDual{iBase,1}(:,cUser,jBase,iBand) .* x(:,cUser,jBase,iBand));
                                    tempADMM = tempADMM + vec(x(:,cUser,jBase,iBand) - cellXGlobal(:,cUser,jBase,iBand)).^2;
                                end
                            end
                        end
                        
                        for iUser = 1:length(cellNeighbourIndices{iBase,1})
                            cUser = cellNeighbourIndices{iBase,1}(iUser,1);
                            tempSecond = tempSecond + sum(currentDual{iBase,1}(:,cUser,iBase,iBand) .* x(:,cUser,iBase,iBand));
                            tempADMM = tempADMM + vec(x(:,cUser,iBase,iBand) - cellXGlobal(:,cUser,iBase,iBand)).^2;
                        end
                        
                    end
                    
                    epiObjective >= norm(userObjective,qExponent) + tempFirst + tempSecond + tempADMM * (alpha / 2);
                    
                    for iBand = 1:nBands
                        
                        for iUser = 1:kUsers
                            
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            H = cH{iBase,iBand}(:,:,cUser);
                            
                            for iLayer = 1:nLayers
                                
                                intVector = sqrt(SimParams.N) * W{cUser,iBand}(:,iLayer);
                                
                                for jBase = 1:nBases
                                    if jBase ~= iBase
                                        intVector = [intVector ; x(iLayer,cUser,jBase,iBand)];
                                    end
                                end
                                
                                for jUser = 1:kUsers
                                    if jUser ~= iUser
                                        for jLayer = 1:nLayers
                                            intVector = [intVector ; W{cUser,iBand}(:,iLayer)' * H * M(:,jLayer,jUser,iBand)];
                                        end
                                    else
                                        for jLayer = 1:nLayers
                                            if jLayer ~= iLayer
                                                intVector = [intVector ; W{cUser,iBand}(:,iLayer)' * H * M(:,jLayer,jUser,iBand)];
                                            end
                                        end
                                    end
                                end
                                
                                givenVector = (1 - W{cUser,iBand}(:,iLayer)' * H * M(:,iLayer,iUser,iBand));
                                intVector = [intVector ; givenVector];
                                
                                norm(intVector,2) <= sqrt(mseError(iLayer,iUser,iBand));
                                (mseError(iLayer,iUser,iBand) - mseError_o(iLayer,iUser,iBand)) / mseError_o(iLayer,iUser,iBand) + log(mseError_o(iLayer,iUser,iBand)) <= -t(iLayer,iUser,iBand) * log(2);
                                
                            end
                            
                        end
                        
                        for iUser = 1:length(cellNeighbourIndices{iBase,1})
                            cUser = cellNeighbourIndices{iBase,1}(iUser,1);
                            for iLayer = 1:nLayers
                                intVector = [];
                                H = cH{iBase,iBand}(:,:,cUser);
                                for jUser = 1:kUsers
                                    for jLayer = 1:nLayers
                                        intVector = [intVector ; W{cUser,iBand}(:,iLayer)' * H * M(:,jLayer,jUser,iBand)];
                                    end
                                end
                                norm(intVector,2) <= x(iLayer,cUser,iBase,iBand);
                            end
                        end
                        
                    end
                    
                    if strcmp(globalMode,'false')
                        for iBand = 1:nBands
                            norm(vec(M(:,:,:,iBand)),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
                        end
                    else
                        norm(vec(M(:,:,:,:)),2) <= sqrt(sum(SimStructs.baseStruct{iBase,1}.sPower(1,:)));
                    end
                    
                    cvx_end
                    
                    if iBase == 1
                        status = cvx_status;
                    else
                        status = strcat(status,'-',cvx_status);
                    end
                    
                    if strfind(cvx_status,'Solved')
                        cellM{iBase,1} = M;
                        cellX{iBase,1} = x;
                        cellBH{iBase,1} = mseError;
                    else
                        cellM{iBase,1} = zeros(size(M));
                    end
                end
                
                [SimParams,SimStructs] = updateIteratePerformance(SimParams,SimStructs,cellM,W);
                
                for iBand = 1:nBands
                    for iBase = 1:nBases
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for jBase = 1:nBases
                                if jBase ~= iBase
                                    cellXGlobal(:,cUser,jBase,iBand) = (cellX{iBase,1}(:,cUser,jBase,iBand) + cellX{jBase,1}(:,cUser,jBase,iBand)) / 2;
                                end
                            end
                        end
                    end
                end
                
                for iBand = 1:nBands
                    for iBase = 1:nBases
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for jBase = 1:nBases
                                if jBase ~= iBase
                                    currentDual{iBase,1}(:,cUser,jBase,iBand) = currentDual{iBase,1}(:,cUser,jBase,iBand) + alpha * (cellX{iBase,1}(:,cUser,jBase,iBand) - cellXGlobal(:,cUser,jBase,iBand));
                                end
                            end
                        end
                        for iUser = 1:length(cellNeighbourIndices{iBase,1})
                            cUser = cellNeighbourIndices{iBase,1}(iUser,1);
                            currentDual{iBase,1}(:,cUser,iBase,iBand) = currentDual{iBase,1}(:,cUser,iBase,iBand) + alpha * (cellX{iBase,1}(:,cUser,iBase,iBand) - cellXGlobal(:,cUser,iBase,iBand));
                        end
                    end
                end
                
                if strcmp(SimParams.DebugMode,'true')
                    status = strcat(status,'-',sprintf('%d',yIteration),'-',sprintf('%d',xIteration));
                    display(status);
                    %display(currentDual);
                    %display([squeeze(cellX{1}) squeeze(cellX{2})]);
                end
            end
            
            sumDeviation = sum(cell2mat(SimParams.Debug.tempResource{3,SimParams.iDrop}));
            if abs(sumDeviation(1,end) - sumDeviationH) < epsilonT
                scaContinue = 0;
            else
                sumDeviationH = sumDeviation(1,end);
            end
            
            for iBand = 1:nBands
                for iUser = 1:nUsers
                    baseNode = SimStructs.userStruct{iUser,1}.baseNode;
                    for iLayer = 1:nLayers
                        R = SimParams.N * eye(SimParams.nRxAntenna);
                        for iBase = 1:nBases
                            for jUser = 1:usersPerCell(iBase,1)
                                H = cH{iBase,iBand}(:,:,iUser);
                                M = cellM{iBase,1}(:,:,jUser,iBand);
                                R = R + H * (M * M') * H';
                            end
                        end
                        H = cH{baseNode,iBand}(:,:,iUser);
                        xUser = (iUser == cellUserIndices{baseNode,1});
                        W{iUser,iBand}(:,iLayer) = R \ (H * cellM{baseNode,1}(:,iLayer,xUser,iBand));
                    end
                end
            end
        end
        
        SimParams.Debug.DataExchange{1,1} = cellM;
        SimParams.Debug.DataExchange{2,1} = cellXGlobal;
        SimParams.Debug.DataExchange{3,1} = currentDual;
        
        for iBand = 1:nBands
            for iBase = 1:nBases
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = cellM{iBase,1}(:,:,:,iBand);
            end
        end
        
    case 'MSEKKTMethod'
        
        maxRank = SimParams.maxRank;
        
        xIndex = 1;
        currentIteration = 0;
        cvx_hist = -500 * ones(2,1);
        SimParams.distDecompSteps = 1;
        
        M = cell(nUsers,nBands);
        R = cell(maxRank,nUsers,nBands);
        [mseError_o,W] = randomizeInitialMSESCApoint(SimParams,SimStructs);
        alphaLKN = ones(maxRank,nUsers,nBands);
        lambdaLKN = ones(maxRank,nUsers,nBands);
        
        while currentIteration < mIterationsSCA
            
            for iBand = 1:nBands
                for iUser = 1:nUsers
                    xNode = SimStructs.userStruct{iUser,1}.baseNode;
                    for iLayer = 1:maxRank
                        I = zeros(SimParams.nTxAntenna);
                        for jUser = 1:nUsers
                            for jLayer = 1:maxRank
                                I = I + cH{xNode,iBand}(:,:,jUser)' * W{jUser,iBand}(:,jLayer) * W{jUser,iBand}(:,jLayer)' * cH{xNode,iBand}(:,:,jUser) * alphaLKN(jLayer,jUser,iBand);
                            end
                        end
                        R{iLayer,iUser,iBand} = I;
                    end
                end
            end
            
            if strcmpi(globalMode,'false')
                
                for iBase = 1:nBases
                    for iBand = 1:nBands
                        muMin = 0;
                        muMax = 100000;
                        iterateAgain = 1;
                        while iterateAgain
                            totalPower = 0;
                            currentMu = (muMax + muMin) / 2;
                            for iUser = 1:usersPerCell(iBase,1)
                                cUser = cellUserIndices{iBase,1}(iUser,1);
                                for iLayer = 1:maxRank
                                    M{cUser,iBand}(:,iLayer) = (currentMu * eye(SimParams.nTxAntenna) + R{iLayer,cUser,iBand}) \ (alphaLKN(iLayer,cUser,iBand) * cH{iBase,iBand}(:,:,cUser)' * W{cUser,iBand}(:,iLayer));
                                    totalPower = totalPower + real(trace(M{cUser,iBand}(:,iLayer) * M{cUser,iBand}(:,iLayer)'));
                                end
                            end
                            
                            if totalPower > (SimStructs.baseStruct{iBase,1}.sPower(1,iBand))
                                muMin = currentMu;
                            else
                                muMax = currentMu;
                            end
                            if abs(muMin - muMax) <= 1e-6
                                iterateAgain = 0;
                            end
                        end
                    end
                end
            else
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
                                for iLayer = 1:maxRank
                                    M{cUser,iBand}(:,iLayer) = (currentMu * eye(SimParams.nTxAntenna) + R{iLayer,cUser,iBand}) \ (alphaLKN(iLayer,cUser,iBand) * cH{iBase,iBand}(:,:,cUser)' * W{cUser,iBand}(:,iLayer));
                                    totalPower = totalPower + real(trace(M{cUser,iBand}(:,iLayer) * M{cUser,iBand}(:,iLayer)'));
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
                end
            end
            
            for iBand = 1:nBands
                for iBase = 1:nBases
                    for iUser = 1:usersPerCell(iBase,1)
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        for iLayer = 1:maxRank
                            X = SimParams.N * eye(SimParams.nRxAntenna);
                            for jBase = 1:nBases
                                for jUser = 1:usersPerCell(jBase,1)
                                    rUser = cellUserIndices{jBase,1}(jUser,1);
                                    H = cH{jBase,iBand}(:,:,cUser);
                                    X = X + H * M{rUser,iBand} * M{rUser,iBand}' * H';
                                end
                            end
                            H = cH{iBase,iBand}(:,:,cUser);
                            W{cUser,iBand}(:,iLayer) = X \ (H * M{cUser,iBand}(:,iLayer));
                        end
                    end
                end
            end
            
            for iBand = 1:nBands
                for iUser = 1:nUsers
                    baseNode = SimStructs.userStruct{iUser,1}.baseNode;
                    for iLayer = 1:maxRank
                        intVector = sqrt(SimParams.N) * W{iUser,iBand}(:,iLayer);
                        for jUser = 1:nUsers
                            ifNode = SimStructs.userStruct{jUser,1}.baseNode;
                            currentH = cH{ifNode,iBand}(:,:,iUser);
                            if jUser == iUser
                                for jLayer = 1:maxRank
                                    if jLayer ~= iLayer
                                        intVector = [intVector ; W{iUser,iBand}(:,iLayer)' * currentH * M{jUser,iBand}(:,jLayer)];
                                    end
                                end
                            else
                                for jLayer = 1:maxRank
                                    intVector = [intVector ; W{iUser,iBand}(:,iLayer)' * currentH * M{jUser,iBand}(:,jLayer)];
                                end
                            end
                        end
                        
                        currentH = cH{baseNode,iBand}(:,:,iUser);
                        givenVector = (1 - W{iUser,iBand}(:,iLayer)' * currentH * M{iUser,iBand}(:,iLayer));
                        intVector = [intVector ; givenVector];
                        mseError(iLayer,iUser,iBand) = norm(intVector,2)^2;
                    end
                end
            end
            
            for iBand = 1:nBands
                for iUser = 1:nUsers
                    for iRank = 1:maxRank
                        t(iRank,iUser,iBand) = -log2(mseError_o(iRank,iUser,iBand)) - (mseError(iRank,iUser,iBand) - mseError_o(iRank,iUser,iBand)) / (mseError_o(iRank,iUser,iBand) * log(2));
                    end
                end
            end
            
            for iBand = 1:nBands
                for iUser = 1:nUsers
                    for iRank = 1:maxRank
                        tempT = t(:,iUser,:);
                        lambdaLKN(iRank,iUser,iBand) = qExponent * (QueuedPkts(iUser,1) - sum(tempT(:)))^(qExponent - 1) / log(2);
                        if mod(qExponent,2) ~= 0
                            deviation = QueuedPkts(iUser,1) - sum(tempT(:));
                            if deviation < 0
                                lambdaLKN(iRank,iUser,iBand) = -qExponent * (QueuedPkts(iUser,1) - sum(tempT(:)))^(qExponent - 1) / log(2);
                            elseif deviation > 0
                                lambdaLKN(iRank,iUser,iBand) = qExponent * (QueuedPkts(iUser,1) - sum(tempT(:)))^(qExponent - 1) / log(2);
                            else
                                lambdaLKN(iRank,iUser,iBand) = 0;
                            end
                        end
                        alphaLKN(iRank,iUser,iBand) = alphaLKN(iRank,iUser,iBand) + 0.1 * (max(lambdaLKN(iRank,iUser,iBand) / mseError(iRank,iUser,iBand),0) - alphaLKN(iRank,iUser,iBand));
                    end
                end
            end
            
            cvx_optval = 0;
            for iUser = 1:nUsers
                tempT = t(:,iUser,:);
                cvx_optval = cvx_optval + abs(QueuedPkts(iUser,1) - sum(tempT(:)));
            end
            
            mseError_o = mseError;
            cellP = cell(nBases,1);
            for iBase = 1:nBases
                cellP{iBase,1} = zeros(SimParams.nTxAntenna,maxRank,usersPerCell(iBase,1),nBands);
                for iBand = 1:nBands
                    for iUser = 1:usersPerCell(iBase,1)
                        for iLayer = 1:maxRank
                            cellP{iBase,1}(:,iLayer,iUser,iBand) = M{cellUserIndices{iBase,1}(iUser,1),iBand}(:,iLayer);
                        end
                    end
                end
            end
            
            [SimParams,SimStructs] = updateIteratePerformance(SimParams,SimStructs,cellP,W);
            
            if (SimParams.currentQueue <= epsilonT)
                reIterate = 0;
            end
            
            tBandUser = SimParams.Debug.privateExchanges.resAllocation;
            qDeviation = QueuedPkts(iUser,1) - sum(tBandUser(:,iUser));
            totalDeviation = norm(qDeviation,qExponent);
            
            if totalDeviation(end) < epsilonT
                reIterate = 0;
            end
            
            if min(abs(cvx_optval - cvx_hist)) <= epsilonT
                reIterate = 0;
            else
                xIndex = xIndex + 1;
                cvx_hist(mod(xIndex,2) + 1,1) = cvx_optval;
            end
            
            currentIteration = currentIteration + 1;
            
        end
        
        SimParams.Debug.DataExchange{1,1} = cellP;
        SimParams.Debug.DataExchange{2,1} = alphaLKN;
        SimParams.Debug.DataExchange{3,1} = lambdaLKN;
        
        for iBase = 1:nBases
            for iBand = 1:nBands
                P = [];
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = zeros(SimParams.nTxAntenna,usersPerCell(iBase,1));
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    P = [P M{cUser,iBand}];
                end
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = P;
            end
        end
        
    case 'ScutariMethod'
        
        prevQ = 1e50;
        prevRate = 1e50;
        xIteration = 0;
        scaContinue = 1;
        cellMH = cell(nBases,1);
        nLayers = SimParams.maxRank;
        cellExchanges = cell(nBases,3);
        [cellP,cellQ,cellB,W,cellM] = randomizeInitialSCApoint(SimParams,SimStructs);
        
        while scaContinue
            
            for iBase = 1:nBases
                
                kUsers = usersPerCell(iBase,1);
                
                cvx_begin
                
                dual variables dualD{nLayers,nUsers,nBands}
                expressions p(nLayers,nUsers,nBands) q(nLayers,nUsers,nBands)
                variable M(SimParams.nTxAntenna,nLayers,nUsers,nBands) complex
                variables t(nLayers,nUsers,nBands) b(nLayers,nUsers,nBands) g(nLayers,nUsers,nBands)
                variables userObjective(nUsers,1) epiObjective
                
                minimize(epiObjective)
                
                subject to
                
                for cUser = 1:nUsers
                    abs(QueuedPkts(cUser,1) - sum(vec(t(:,cUser,:)))) <= userObjective(cUser,1);
                end
                
                epiObjective >= norm(userObjective,qExponent);
                
                for iBand = 1:nBands
                    for iUser = 1:kUsers
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        H = cH{iBase,iBand}(:,:,cUser);
                        
                        for iLayer = 1:nLayers
                            intVector = sqrt(SimParams.N) * W{cUser,iBand}(:,iLayer)';
                            for jUser = 1:nUsers
                                baseNode = SimStructs.userStruct{jUser,1}.baseNode;
                                iH = cH{baseNode,iBand}(:,:,cUser);
                                if jUser ~= cUser
                                    if baseNode ~= iBase
                                        intVector = [intVector W{cUser,iBand}(:,iLayer)' * iH * cellM(:,:,jUser,iBand)];
                                    else
                                        intVector = [intVector W{cUser,iBand}(:,iLayer)' * iH * M(:,:,jUser,iBand)];
                                    end
                                else
                                    for jLayer = 1:nLayers
                                        intVector = [intVector W{cUser,iBand}(:,iLayer)' * iH * M(:,iLayer ~= (1:nLayers),jUser,iBand)];
                                    end
                                end
                            end
                            
                            norm(intVector,2) <= sqrt(b(iLayer,cUser,iBand));
                            log(1 + g(iLayer,cUser,iBand)) >= t(iLayer,cUser,iBand) * log(2);
                            p(iLayer,cUser,iBand) = real(W{cUser,iBand}(:,iLayer)' * H * M(:,iLayer,cUser,iBand));
                            q(iLayer,cUser,iBand) = imag(W{cUser,iBand}(:,iLayer)' * H * M(:,iLayer,cUser,iBand));
                            (cellP(iLayer,cUser,iBand)^2 + cellQ(iLayer,cUser,iBand)^2) / (cellB(iLayer,cUser,iBand)) + ...
                                (2 / cellB(iLayer,cUser,iBand)) * (cellP(iLayer,cUser,iBand) * (p(iLayer,cUser,iBand) - cellP(iLayer,cUser,iBand))) + ...
                                (2 / cellB(iLayer,cUser,iBand)) * (cellQ(iLayer,cUser,iBand) * (q(iLayer,cUser,iBand) - cellQ(iLayer,cUser,iBand))) - ...
                                (cellP(iLayer,cUser,iBand)^2 + cellQ(iLayer,cUser,iBand)^2) / (cellB(iLayer,cUser,iBand)^2) * ...
                                (b(iLayer,cUser,iBand) - cellB(iLayer,cUser,iBand)) >= g(iLayer,cUser,iBand);
                        end
                    end
                    
                    for iUser = 1:length(cellNeighbourIndices{iBase,1})
                        
                        cUser = cellNeighbourIndices{iBase,1}(iUser,1);
                        gBaseNode = SimStructs.userStruct{cUser,1}.baseNode;
                        nH = cH{gBaseNode,iBand}(:,:,cUser);
                        
                        for iLayer = 1:nLayers
                            intVector = sqrt(SimParams.N) * W{cUser,iBand}(:,iLayer)';
                            for jUser = 1:nUsers
                                baseNode = SimStructs.userStruct{jUser,1}.baseNode;
                                iH = cH{baseNode,iBand}(:,:,cUser);
                                
                                if jUser ~= cUser
                                    if baseNode ~= iBase
                                        intVector = [intVector W{cUser,iBand}(:,iLayer)' * iH * cellM(:,:,jUser,iBand)];
                                    else
                                        intVector = [intVector W{cUser,iBand}(:,iLayer)' * iH * M(:,:,jUser,iBand)];
                                    end
                                else
                                    for jLayer = 1:nLayers
                                        intVector = [intVector W{cUser,iBand}(:,iLayer)' * iH * cellM(:,iLayer ~= (1:nLayers),jUser,iBand)];
                                    end
                                end
                            end
                            
                            norm(intVector,2) <= sqrt(b(iLayer,cUser,iBand));
                            log(1 + g(iLayer,cUser,iBand)) >= t(iLayer,cUser,iBand) * log(2);
                            p(iLayer,cUser,iBand) = cellP(iLayer,cUser,iBand);
                            q(iLayer,cUser,iBand) = cellQ(iLayer,cUser,iBand);
                            (cellP(iLayer,cUser,iBand)^2 + cellQ(iLayer,cUser,iBand)^2) / (cellB(iLayer,cUser,iBand)) - ...
                                (cellP(iLayer,cUser,iBand)^2 + cellQ(iLayer,cUser,iBand)^2) / (cellB(iLayer,cUser,iBand)^2) * ...
                                (b(iLayer,cUser,iBand) - cellB(iLayer,cUser,iBand)) >= g(iLayer,cUser,iBand);
                        end
                    end
                end
                
                if strcmp(globalMode,'false')
                    for iBand = 1:nBands
                        norm(vec(M(:,:,cellUserIndices{iBase,1},iBand)),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
                    end
                else
                    norm(vec(M),2) <= sqrt(sum(SimStructs.baseStruct{iBase,1}.sPower(1,:)));
                end
                
                cvx_end
                
                if iBase == 1
                    status = cvx_status;
                else
                    status = strcat(status,'-',cvx_status);
                end
                
                if strfind(cvx_status,'Solved')
                    cellMH{iBase,1} = M;
                    cellExchanges{iBase,1} = p;
                    cellExchanges{iBase,2} = q;
                    cellExchanges{iBase,3} = b;
                else
                    display('Failed');
                    cellMH{iBase,1} = zeros(size(M));
                end
                
            end
            
            cellMM = cell(nBases,1);
            for iBase = 1:nBases
                cellMM{iBase,1} = cellMH{iBase,1}(:,:,cellUserIndices{iBase,1},:);
            end
            
            [SimParams,SimStructs] = updateIteratePerformance(SimParams,SimStructs,cellMM,W);
            
            for iBase = 1:nBases
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    cellM(:,:,cUser,:) = cellM(:,:,cUser,:) + alpha * (cellMH{iBase,1}(:,:,cUser,:) - cellM(:,:,cUser,:));
                end
            end
            
            for iBand = 1:nBands
                for iUser = 1:nUsers
                    baseNode = SimStructs.userStruct{iUser,1}.baseNode;
                    for iLayer = 1:nLayers
                        R = SimParams.N * eye(SimParams.nRxAntenna);
                        for iBase = 1:nBases
                            for jUser = 1:usersPerCell(iBase,1)
                                jcUser = cellUserIndices{iBase,1}(jUser,1);
                                H = cH{iBase,iBand}(:,:,iUser);
                                M = cellM(:,:,jcUser,iBand);
                                R = R + H * (M * M') * H';
                            end
                        end
                        H = cH{baseNode,iBand}(:,:,iUser);
                        W{iUser,iBand}(:,iLayer) = R \ (H * cellM(:,iLayer,iUser,iBand));
                    end
                end
            end
            
            %             for iBase = 1:nBases
            %                 for iBand = 1:nBands
            %                     for iUser = 1:usersPerCell(iBase,1)
            %                         cUser = cellUserIndices{iBase,1}(iUser,1);
            %                         for iLayer = 1:nLayers
            %                             cellP(iLayer,cUser,iBand) = real(W{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * cellM(:,iLayer,cUser,iBand));
            %                             cellQ(iLayer,cUser,iBand) = imag(W{cUser,iBand}(:,iLayer)' * cH{iBase,iBand}(:,:,cUser) * cellM(:,iLayer,cUser,iBand));
            %                         end
            %                     end
            %                 end
            %             end
            
            %             for iBand = 1:nBands
            %                 for iUser = 1:nUsers
            %                     for iLayer = 1:nLayers
            %                         tempVector = W{iUser,iBand}(:,iLayer)' * sqrt(SimParams.N);
            %                         for jUser = 1:nUsers
            %                             nBaseNode = SimStructs.userStruct{jUser,1}.baseNode;
            %                             if jUser ~= iUser
            %                                 tempVector = [tempVector W{iUser,iBand}(:,iLayer)' * cH{nBaseNode,iBand}(:,:,iUser) * cellM(:,:,jUser,iBand)];
            %                             else
            %                                 tempVector = [tempVector W{iUser,iBand}(:,iLayer)' * cH{nBaseNode,iBand}(:,:,iUser) * cellM(:,iLayer~=(1:nLayers),jUser,iBand)];
            %                             end
            %                         end
            %                         cellB(iLayer,iUser,iBand) = norm(tempVector,2)^2;
            %                     end
            %                 end
            %             end
            
            for iBase = 1:nBases
                cellP(:,cellUserIndices{iBase,1},:) = cellExchanges{iBase,1}(:,cellUserIndices{iBase,1},:);
                cellQ(:,cellUserIndices{iBase,1},:) = cellExchanges{iBase,2}(:,cellUserIndices{iBase,1},:);
                cellB(:,cellUserIndices{iBase,1},:) = cellExchanges{iBase,3}(:,cellUserIndices{iBase,1},:);
            end
            
            if xIteration < mIterationsSCA
                xIteration = xIteration + 1;
            else
                scaContinue = 0;
            end
            
            cRate = sum(cell2mat(SimParams.Debug.tempResource{2,SimParams.iDrop}));
            if abs(prevRate - cRate(1,end)) < epsilonT
                scaContinue = 0;
            else
                prevRate = cRate(1,end);
            end
            
            cQueues = sum(cell2mat(SimParams.Debug.tempResource{3,SimParams.iDrop}));
            if abs(prevQ - cQueues(1,end)) < epsilonT
                scaContinue = 0;
            else
                prevQ = cQueues(1,end);
            end
            
        end
        
        for iBand = 1:nBands
            for iBase = 1:nBases
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = cellM(:,:,cellUserIndices{iBase,1},iBand);
            end
        end
        
    case 'KKTMethod-1'
        
        maxRank = SimParams.maxRank;
        
        reIterate = 1;
        maxIterations = 100;
        currentIteration = 0;
        SimParams.distDecompSteps = 1;
        
        M = cell(nUsers,nBands);
        xM = cell(nUsers,nBands);
        W = cell(nUsers,nBands);
        R = cell(maxRank,nUsers,nBands);
        
        betaLKN_o = zeros(maxRank,nUsers,nBands);
        gammaLKN_o = zeros(maxRank,nUsers,nBands);
        
        for iUser = 1:nUsers
            for iBand = 1:nBands
                xM{iUser,iBand} = complex(randn(SimParams.nTxAntenna,maxRank),randn(SimParams.nTxAntenna,maxRank)) / sqrt(SimParams.nTxAntenna);
                W{iUser,iBand} = ones(SimParams.nRxAntenna);
            end
        end
        
        rho = 1e-4;
        alphaLKN_1 = ones(maxRank,nUsers,nBands);
        alphaLKN_2 = ones(maxRank,nUsers,nBands);
        
        while reIterate
            
            for iBand = 1:nBands
                for iUser = 1:nUsers
                    for iLayer = 1:maxRank
                        I = SimParams.N;
                        for jUser = 1:nUsers
                            if jUser ~= iUser
                                for jLayer = 1:maxRank
                                    xNode = SimStructs.userStruct{jUser,1}.baseNode;
                                    I = I + norm(cH{xNode,iBand}(:,:,iUser) * xM{jUser,iBand}(:,jLayer))^2;
                                end
                            end
                        end
                        betaLKN_o(iLayer,iUser,iBand) = I;
                    end
                    xNode = SimStructs.userStruct{iUser,1}.baseNode;
                    gammaLKN_o(iLayer,iUser,iBand) = norm(cH{xNode,iBand}(:,:,iUser) * xM{iUser,iBand}(:,iLayer))^2 / betaLKN_o(iLayer,iUser,iBand);
                end
            end
            
            phiLKN = gammaLKN_o ./ betaLKN_o;
            
            for dualIterate = 1:1e3
                
                for iBand = 1:nBands
                    for iUser = 1:nUsers
                        xNode = SimStructs.userStruct{iUser,1}.baseNode;
                        for iLayer = 1:maxRank
                            I = zeros(SimParams.nTxAntenna);
                            for jUser = 1:nUsers
                                if jUser ~= iUser
                                    for jLayer = 1:maxRank
                                        I = I + cH{xNode,iBand}(:,:,jUser)' * cH{xNode,iBand}(:,:,jUser) * alphaLKN_2(jLayer,jUser,iBand);
                                    end
                                end
                            end
                            R{iLayer,iUser,iBand} = I;
                        end
                    end
                end
                
                if strcmpi(globalMode,'false')
                    
                    for iBase = 1:nBases
                        for iBand = 1:nBands
                            muMin = 0;
                            muMax = 100000;
                            iterateAgain = 1;
                            while iterateAgain
                                totalPower = 0;
                                currentMu = (muMax + muMin) / 2;
                                for iUser = 1:usersPerCell(iBase,1)
                                    cUser = cellUserIndices{iBase,1}(iUser,1);
                                    for iLayer = 1:maxRank
                                        M{cUser,iBand}(:,iLayer) = (currentMu * eye(SimParams.nTxAntenna) + R{iLayer,cUser,iBand}) \ (alphaLKN_1(iLayer,cUser,iBand) * cH{iBase,iBand}(:,:,cUser)' * cH{iBase,iBand}(:,:,cUser) * xM{cUser,iBand}(:,iLayer));
                                        totalPower = totalPower + real(trace(M{cUser,iBand}(:,iLayer) * M{cUser,iBand}(:,iLayer)'));
                                    end
                                end
                                
                                if totalPower > (SimStructs.baseStruct{iBase,1}.sPower(1,iBand))
                                    muMin = currentMu;
                                else
                                    muMax = currentMu;
                                end
                                if abs(muMin - muMax) <= 1e-6
                                    iterateAgain = 0;
                                end
                            end
                        end
                    end
                else
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
                                    for iLayer = 1:maxRank
                                        M{cUser,iBand}(:,iLayer) = (currentMu * eye(SimParams.nTxAntenna) + R{iLayer,cUser,iBand}) \ (alphaLKN_1(iLayer,cUser,iBand) * cH{iBase,iBand}(:,:,cUser)' * cH{iBase,iBand}(:,:,cUser) * xM{cUser,iBand}(:,iLayer));
                                        totalPower = totalPower + real(trace(M{cUser,iBand}(:,iLayer) * M{cUser,iBand}(:,iLayer)'));
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
                    end
                end
                
                phiLKN(phiLKN < 1e-3) = 1e-3;
                phiLKN(phiLKN == Inf) = 1e3;
                betaLKN = alphaLKN_2 ./ (alphaLKN_1 .* phiLKN);
                gammaLKN = -0.5 + 0.5 * sqrt(1 + (4 * phiLKN ./ alphaLKN_1));
                
                for iBand = 1:nBands
                    for iBase = 1:nBases
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for iLayer = 1:maxRank
                                X = SimParams.N * eye(SimParams.nRxAntenna);
                                for jBase = 1:nBases
                                    for jUser = 1:usersPerCell(jBase,1)
                                        rUser = cellUserIndices{jBase,1}(jUser,1);
                                        if rUser ~= cUser
                                            H = cH{jBase,iBand}(:,:,cUser);
                                            X = X + H * M{rUser,iBand} * M{rUser,iBand}' * H';
                                        end
                                    end
                                end
                                alphaLKN_2(iLayer,cUser,iBand) = alphaLKN_2(iLayer,cUser,iBand) + rho * (real(trace(X)) - betaLKN(iLayer,cUser,iBand));
                            end
                        end
                    end
                end
                
                for iBand = 1:nBands
                    for iBase = 1:nBases
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for iLayer = 1:maxRank
                                H = cH{iBase,iBand}(:,:,cUser);
                                firstTerm = - norm(H * xM{cUser,iBand}(:,iLayer))^2;
                                secTerm = - 2 * xM{cUser,iBand}(:,iLayer)' * (H' * H) * (M{cUser,iBand}(:,iLayer) - xM{cUser,iBand}(:,iLayer));
                                thTerm = phiLKN(iLayer,cUser,iBand) * betaLKN(iLayer,cUser,iBand)^2 / 2;
                                fourthTerm = gammaLKN(iLayer,cUser,iBand)^2 / (2 * phiLKN(iLayer,cUser,iBand));
                                if phiLKN(iLayer,cUser,iBand) == 0
                                    alphaLKN_1(iLayer,cUser,iBand) = 0;
                                else
                                    alphaLKN_1(iLayer,cUser,iBand) = alphaLKN_1(iLayer,cUser,iBand) + rho * (firstTerm + secTerm + thTerm + fourthTerm);
                                end
                            end
                        end
                    end
                end
                
                alphaLKN_1 = max(real(alphaLKN_1),0);
                alphaLKN_2 = max(real(alphaLKN_2),0);
                
            end
            
            xM = M;
            cellP = cell(nBases,1);
            for iBase = 1:nBases
                cellP{iBase,1} = zeros(SimParams.nTxAntenna,maxRank,usersPerCell(iBase,1),nBands);
                for iBand = 1:nBands
                    for iUser = 1:usersPerCell(iBase,1)
                        for iLayer = 1:maxRank
                            cellP{iBase,1}(:,iLayer,iUser,iBand) = M{cellUserIndices{iBase,1}(iUser,1),iBand}(:,iLayer);
                        end
                    end
                end
            end
            
            [SimParams,SimStructs] = updateIteratePerformance(SimParams,SimStructs,cellP,W);
            
            tBandUser = SimParams.Debug.privateExchanges.resAllocation;
            qDeviation = QueuedPkts(iUser,1) - sum(tBandUser(:,iUser));
            totalDeviation = norm(qDeviation,qExponent);
            
            if totalDeviation(end) < epsilonT
                reIterate = 0;
            end
            
            currentIteration = currentIteration + 1;
            if currentIteration > maxIterations
                reIterate = 0;
            end
            
        end
        
        for iBase = 1:nBases
            for iBand = 1:nBands
                P = [];
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = zeros(SimParams.nTxAntenna,usersPerCell(iBase,1));
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    P = [P M{cUser,iBand}];
                end
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = P;
            end
        end
        
    case 'KKTMethod-2'
        
        maxRank = SimParams.maxRank;
        
        reIterate = 1;
        maxIterations = 50;
        currentIteration = 0;
        SimParams.distDecompSteps = 1;
        
        M = cell(nUsers,nBands);
        xM = cell(nUsers,nBands);
        W = cell(nUsers,nBands);
        R = cell(maxRank,nUsers,nBands);
        
        betaLKN = zeros(maxRank,nUsers,nBands);
        betaLKN_o = zeros(maxRank,nUsers,nBands);
        gammaLKN_o = zeros(maxRank,nUsers,nBands);
        
        for iUser = 1:nUsers
            for iBand = 1:nBands
                xM{iUser,iBand} = complex(randn(SimParams.nTxAntenna,maxRank),randn(SimParams.nTxAntenna,maxRank)) / sqrt(SimParams.nTxAntenna);
                W{iUser,iBand} = ones(SimParams.nRxAntenna);
            end
        end
        
        rho = 1e-3;
        alphaLKN_1 = ones(maxRank,nUsers,nBands);
        alphaLKN_2 = ones(maxRank,nUsers,nBands);
        
        while reIterate
            
            for iBand = 1:nBands
                for iUser = 1:nUsers
                    for iLayer = 1:maxRank
                        I = SimParams.N;
                        for jUser = 1:nUsers
                            if jUser ~= iUser
                                for jLayer = 1:maxRank
                                    xNode = SimStructs.userStruct{jUser,1}.baseNode;
                                    I = I + norm(cH{xNode,iBand}(:,:,iUser) * xM{jUser,iBand}(:,jLayer))^2;
                                end
                            end
                        end
                        betaLKN_o(iLayer,iUser,iBand) = I;
                    end
                    xNode = SimStructs.userStruct{iUser,1}.baseNode;
                    gammaLKN_o(iLayer,iUser,iBand) = norm(cH{xNode,iBand}(:,:,iUser) * xM{iUser,iBand}(:,iLayer))^2 / betaLKN_o(iLayer,iUser,iBand);
                end
            end
            
            phiLKN = sqrt(gammaLKN_o ./ betaLKN_o);
            
            for dualIterate = 1:100
                
                for iBand = 1:nBands
                    for iUser = 1:nUsers
                        xNode = SimStructs.userStruct{iUser,1}.baseNode;
                        for iLayer = 1:maxRank
                            I = zeros(SimParams.nTxAntenna);
                            for jUser = 1:nUsers
                                if jUser ~= iUser
                                    for jLayer = 1:maxRank
                                        I = I + cH{xNode,iBand}(:,:,jUser)' * cH{xNode,iBand}(:,:,jUser) * alphaLKN_2(jLayer,jUser,iBand);
                                    end
                                end
                            end
                            R{iLayer,iUser,iBand} = I;
                        end
                    end
                end
                
                if strcmpi(globalMode,'false')
                    
                    for iBase = 1:nBases
                        for iBand = 1:nBands
                            muMin = 0;
                            muMax = 100000;
                            iterateAgain = 1;
                            while iterateAgain
                                totalPower = 0;
                                currentMu = (muMax + muMin) / 2;
                                for iUser = 1:usersPerCell(iBase,1)
                                    cUser = cellUserIndices{iBase,1}(iUser,1);
                                    for iLayer = 1:maxRank
                                        M{cUser,iBand}(:,iLayer) = (currentMu * eye(SimParams.nTxAntenna) + R{iLayer,cUser,iBand}) \ (0.5 * alphaLKN_1(iLayer,cUser,iBand) * cH{iBase,iBand}(:,:,cUser)');
                                        totalPower = totalPower + real(trace(M{cUser,iBand}(:,iLayer) * M{cUser,iBand}(:,iLayer)'));
                                    end
                                end
                                
                                if totalPower > (SimStructs.baseStruct{iBase,1}.sPower(1,iBand))
                                    muMin = currentMu;
                                else
                                    muMax = currentMu;
                                end
                                if abs(muMin - muMax) <= 1e-6
                                    iterateAgain = 0;
                                end
                            end
                        end
                    end
                else
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
                                    for iLayer = 1:maxRank
                                        M{cUser,iBand}(:,iLayer) = (currentMu * eye(SimParams.nTxAntenna) + R{iLayer,cUser,iBand}) \ (0.5 * alphaLKN_1(iLayer,cUser,iBand) * cH{iBase,iBand}(:,:,cUser)');
                                        totalPower = totalPower + real(trace(M{cUser,iBand}(:,iLayer) * M{cUser,iBand}(:,iLayer)'));
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
                    end
                end
                
                gammaLKN = (2 * phiLKN ./ alphaLKN_1) - 1;
                
                for iBand = 1:nBands
                    for iBase = 1:nBases
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for iLayer = 1:maxRank
                                X = SimParams.N * eye(SimParams.nRxAntenna);
                                for jBase = 1:nBases
                                    for jUser = 1:usersPerCell(jBase,1)
                                        rUser = cellUserIndices{jBase,1}(jUser,1);
                                        if rUser ~= cUser
                                            H = cH{jBase,iBand}(:,:,cUser);
                                            X = X + H * M{rUser,iBand} * M{rUser,iBand}' * H';
                                        end
                                    end
                                end
                                betaLKN(iLayer,cUser,iBand) = real(trace(X));
                            end
                        end
                    end
                end
                
                for iBand = 1:nBands
                    for iBase = 1:nBases
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for iLayer = 1:maxRank
                                H = cH{iBase,iBand}(:,:,cUser);
                                firstTerm = - H * M{cUser,iBand}(:,iLayer);
                                thTerm = phiLKN(iLayer,cUser,iBand) * betaLKN(iLayer,cUser,iBand) / 2;
                                fourthTerm = gammaLKN(iLayer,cUser,iBand) / (2 * phiLKN(iLayer,cUser,iBand));
                                if phiLKN(iLayer,cUser,iBand) == 0
                                    alphaLKN_1(iLayer,cUser,iBand) = 0;
                                else
                                    alphaLKN_1(iLayer,cUser,iBand) = alphaLKN_1(iLayer,cUser,iBand) + rho * (firstTerm + thTerm + fourthTerm);
                                end
                            end
                        end
                    end
                end
                
                alphaLKN_2 = 0.5 * alphaLKN_1 .* phiLKN;
                
                alphaLKN_1 = max(real(alphaLKN_1),0);
                alphaLKN_2 = max(real(alphaLKN_2),0);
                
            end
            
            xM = M;
            cellP = cell(nBases,1);
            for iBase = 1:nBases
                cellP{iBase,1} = zeros(SimParams.nTxAntenna,maxRank,usersPerCell(iBase,1),nBands);
                for iBand = 1:nBands
                    for iUser = 1:usersPerCell(iBase,1)
                        for iLayer = 1:maxRank
                            cellP{iBase,1}(:,iLayer,iUser,iBand) = M{cellUserIndices{iBase,1}(iUser,1),iBand}(:,iLayer);
                        end
                    end
                end
            end
            
            [SimParams,SimStructs] = updateIteratePerformance(SimParams,SimStructs,cellP,W);
            
            tBandUser = SimParams.Debug.privateExchanges.resAllocation;
            qDeviation = QueuedPkts(iUser,1) - sum(tBandUser(:,iUser));
            totalDeviation = norm(qDeviation,qExponent);
            
            if totalDeviation(end) < epsilonT
                reIterate = 0;
            end
            
            currentIteration = currentIteration + 1;
            if currentIteration > maxIterations
                reIterate = 0;
            end
            
        end
        
        for iBase = 1:nBases
            for iBand = 1:nBands
                P = [];
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = zeros(SimParams.nTxAntenna,usersPerCell(iBase,1));
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    P = [P M{cUser,iBand}];
                end
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = P;
            end
        end
        
end

for iUser = 1:nUsers
    SimParams.Debug.tempResource{2,SimParams.iDrop}{iUser,1} = SimParams.Debug.tempResource{2,SimParams.iDrop}{iUser,1};
    for iBand = 1:nBands
        SimParams.Debug.tempResource{4,SimParams.iDrop}{iUser,iBand} = SimParams.Debug.tempResource{4,SimParams.iDrop}{iUser,iBand};
        SimStructs.userStruct{iUser,1}.W{iBand,1} = W{iUser,iBand};
    end
end
