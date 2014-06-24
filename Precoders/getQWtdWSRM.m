
function [SimParams,SimStructs] = getQWtdWSRM(SimParams,SimStructs)

epsilonT = 1e-5;
maxIterations = 50;
cH = SimStructs.linkChan;
nBases = SimParams.nBases;
nBands = SimParams.nBands;
globalMode = SimParams.totalPwrDistOverSC;

updatePrecoders = 'true';
usersPerCell = zeros(nBases,1);
cellUserIndices = cell(nBases,1);

% Debug Buffers initialization

SimParams.Debug.tempResource{2,SimParams.iDrop} = cell(SimParams.nUsers,1);
SimParams.Debug.tempResource{3,SimParams.iDrop} = cell(SimParams.nUsers,1);
SimParams.Debug.tempResource{4,SimParams.iDrop} = cell(SimParams.nUsers,SimParams.nBands);

for iBase = 1:nBases
    for iBand = 1:nBands
        cellUserIndices{iBase,1} = [cellUserIndices{iBase,1} ; SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1}];
    end
    cellUserIndices{iBase,1} = unique(cellUserIndices{iBase,1});
    usersPerCell(iBase,1) = length(cellUserIndices{iBase,1});
end

nUsers = sum(usersPerCell);
QueuedPkts = zeros(nUsers,1);
bandRateMax = zeros(SimParams.nUsers,nBands);

for iBase = 1:nBases
    for iUser = 1:usersPerCell(iBase,1)
        cUser = cellUserIndices{iBase,1}(iUser,1);
        QueuedPkts(cUser,1) = SimStructs.userStruct{cUser,1}.trafficStats.backLogPkt;
    end
end

userWts = ones(nUsers,1);
underscore_location = strfind(SimParams.weightedSumRateMethod,'_');
if isempty(underscore_location)
    qExponent = 1;
    selectionMethod = SimParams.weightedSumRateMethod;
else
    qExponent = str2double(SimParams.weightedSumRateMethod(underscore_location + 1:end));
    selectionMethod = SimParams.weightedSumRateMethod(1:underscore_location-1);
end

switch selectionMethod
    
    case 'BandAlloc'
        
        currentDesign = SimParams.PrecodingMethod;
        currentApproach = SimParams.weightedSumRateMethod;
        
        SimParams.PrecodingMethod = 'Best_WMMSE_Method';
        SimParams.weightedSumRateMethod = 'PerformScheduling';
        
        for iBand = 1:nBands
            for iBase = 1:nBases
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    SimStructs.userStruct{cUser,1}.weighingFactor = QueuedPkts(cUser,1);
                end
            end
            
            SimParams.Debug.privateExchanges.takeOverBand = iBand;
            [SimParams, SimStructs] = getWeightedMMSEDesign(SimParams,SimStructs);
            
            [SimParams,SimStructs] = performDummyReception(SimParams,SimStructs,iBand);
            QueuedPkts = max(QueuedPkts - SimParams.Debug.privateExchanges.resAllocation(iBand,:)',0);
            
        end
        
        for iBase = 1:nBases
            for iUser = 1:usersPerCell(iBase,1)
                cUser = cellUserIndices{iBase,1}(iUser,1);
                sumRateOverBand = [];bandMaxRate = 0;
                for iBand = 1:nBands
                    sumRateOverBand = [sumRateOverBand (SimParams.Debug.tempResource{4,SimParams.iDrop}{cUser,iBand} + bandMaxRate)];
                    bandMaxRate = sumRateOverBand(1,end);
                end
                SimParams.Debug.tempResource{2,SimParams.iDrop}{cUser,1} = sumRateOverBand;
            end
        end
        
        SimParams.PrecodingMethod = currentDesign;
        SimParams.weightedSumRateMethod = currentApproach;
        SimParams.privateExchanges = rmfield(SimParams.Debug.privateExchanges,'takeOverBand');
        
        updatePrecoders = 'false';
        
    case 'GenAlloc'
        
        xIndex = 0;
        reIterate = 1;
        currentIteration = 0;
        cvx_hist = -500 * ones(2,1);
        maxRank = SimParams.maxRank;
        [p_o,q_o,b_o,vW] = randomizeInitialSCApoint(SimParams,SimStructs);
        
        while reIterate
            
            cvx_begin
            
            expressions p(maxRank,nUsers,nBands) q(maxRank,nUsers,nBands)
            variable M(SimParams.nTxAntenna,maxRank,nUsers,nBands) complex
            variables t(maxRank,nUsers,nBands) b(maxRank,nUsers,nBands) g(maxRank,nUsers,nBands)
            variables userObjective(nUsers,1) epiObjective
            
            minimize(epiObjective)
            
            subject to
            
            for iUser = 1:nUsers
                userWts(iUser,1) * abs(QueuedPkts(iUser,1) - sum(vec(t(:,iUser,:)))) <= userObjective(iUser,1);
            end
            
            epiObjective >= norm(userObjective,qExponent);
            
            for iBase = 1:nBases
                for iBand = 1:nBands
                    for iUser = 1:usersPerCell(iBase,1)
                        
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        for iLayer = 1:maxRank
                            intVector = sqrt(SimParams.N) * vW{cUser,iBand}(:,iLayer);
                            
                            for jBase = 1:nBases
                                currentH = cH{jBase,iBand}(:,:,cUser);
                                for jUser = 1:usersPerCell(jBase,1)
                                    rUser = cellUserIndices{jBase,1}(jUser,1);
                                    if rUser ~= cUser
                                        for jLayer = 1:maxRank
                                            intVector = [intVector ; vW{cUser,iBand}(:,iLayer)' * currentH * M(:,jLayer,rUser,iBand)];
                                        end
                                    else
                                        for jLayer = 1:maxRank
                                            if jLayer ~= iLayer
                                                intVector = [intVector ; vW{cUser,iBand}(:,iLayer)' * currentH * M(:,jLayer,rUser,iBand)];
                                            end
                                        end
                                    end
                                end
                            end
                            
                            norm(intVector,2) <= sqrt(b(iLayer,cUser,iBand));
                            log(1 + g(iLayer,cUser,iBand)) >= t(iLayer,cUser,iBand) * log(2);
                            
                            currentH = cH{iBase,iBand}(:,:,cUser);
                            p(iLayer,cUser,iBand) = real(vW{cUser,iBand}(:,iLayer)' * currentH * M(:,iLayer,cUser,iBand));
                            q(iLayer,cUser,iBand) = imag(vW{cUser,iBand}(:,iLayer)' * currentH * M(:,iLayer,cUser,iBand));
                            
                            (p_o(iLayer,cUser,iBand)^2 + q_o(iLayer,cUser,iBand)^2) / (b_o(iLayer,cUser,iBand)) + ...
                                (2 / b_o(iLayer,cUser,iBand)) * (p_o(iLayer,cUser,iBand) * (p(iLayer,cUser,iBand) - p_o(iLayer,cUser,iBand))) + ...
                                (2 / b_o(iLayer,cUser,iBand)) * (q_o(iLayer,cUser,iBand) * (q(iLayer,cUser,iBand) - q_o(iLayer,cUser,iBand))) - ...
                                (p_o(iLayer,cUser,iBand)^2 + q_o(iLayer,cUser,iBand)^2) / (b_o(iLayer,cUser,iBand)^2) * ...
                                (b(iLayer,cUser,iBand) - b_o(iLayer,cUser,iBand)) >= g(iLayer,cUser,iBand);
                            
                        end
                    end
                end
                
                if strcmp(globalMode,'false')
                    for iBand = 1:nBands
                        norm(vec(M(:,:,cellUserIndices{iBase,1},iBand)),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
                    end
                else
                    norm(vec(M(:,:,cellUserIndices{iBase,1},:)),2) <= sqrt(sum(SimStructs.baseStruct{iBase,1}.sPower(1,:)));
                end
                
            end
            
            cvx_end
            
            if strfind(cvx_status,'Solved')
                
                b_o = b;
                for iBand = 1:nBands
                    for iBase = 1:nBases
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            currentH = cH{iBase,iBand}(:,:,cUser);
                            for iLayer = 1:maxRank
                                p_o(iLayer,cUser,iBand) = real(vW{cUser,iBand}(:,iLayer)' * currentH * M(:,iLayer,cUser,iBand));
                                q_o(iLayer,cUser,iBand) = imag(vW{cUser,iBand}(:,iLayer)' * currentH * M(:,iLayer,cUser,iBand));
                            end
                        end
                    end
                end
                
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
                                        R = R + H * M(:,:,rUser,iBand) * M(:,:,rUser,iBand)' * H';
                                    end
                                end
                                H = cH{iBase,iBand}(:,:,cUser);
                                vW{cUser,iBand}(:,iLayer) = R \ (H * M(:,iLayer,cUser,iBand));
                            end
                        end
                    end
                end
                
                if min(abs(cvx_optval - cvx_hist)) <= epsilonT
                    reIterate = 0;
                else
                    xIndex = xIndex + 1;
                    cvx_hist(mod(xIndex,2) + 1,1) = cvx_optval;
                end
                
            else
                b_o = b_o * 2;
                display('Failed CVX !');
            end
            
            currentIteration = currentIteration + 1;
            if currentIteration >= maxIterations
                reIterate = 0;
            end
            
            [SimParams,SimStructs] = updateIteratePerformance(SimParams,SimStructs,M,vW);
            
        end
        
        for iBase = 1:nBases
            for iBand = 1:nBands
                P = [];
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = zeros(SimParams.nTxAntenna,usersPerCell(iBase,1));
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    P = [P M(:,:,cUser,iBand)];
                end
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = P;
            end
        end
        
    case 'GenBandAlloc'
        
        vW = cell(nUsers,nBands);
        for iBand = 1:nBands
            
            xIndex = 0;
            reIterate = 1;
            currentIteration = 0;
            maxRank = SimParams.maxRank;
            cvx_hist = -500 * ones(2,1);
            [p_o,q_o,b_o,xW] = randomizeInitialSCApoint(SimParams,SimStructs,iBand);
            
            for iUser = 1:nUsers
                vW{iUser,iBand} = xW{iUser,1};
            end
            
            while reIterate
                
                cvx_begin
                
                expressions p(maxRank,nUsers) q(maxRank,nUsers)
                variable M(SimParams.nTxAntenna,maxRank,nUsers) complex
                variables t(maxRank,nUsers) b(maxRank,nUsers) g(maxRank,nUsers)
                variables userObjective(nUsers,1) epiObjective
                
                minimize(epiObjective)
                
                subject to
                
                for iUser = 1:nUsers
                    userWts(iUser,1) * abs(QueuedPkts(iUser,1) - sum(vec(t(:,iUser)))) <= userObjective(iUser,1);
                end
                
                epiObjective >= norm(userObjective,qExponent);
                
                for iBase = 1:nBases
                    for iUser = 1:usersPerCell(iBase,1)
                        
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        for iLayer = 1:maxRank
                            intVector = sqrt(SimParams.N) * vW{cUser,iBand}(:,iLayer);
                            
                            for jBase = 1:nBases
                                currentH = cH{jBase,iBand}(:,:,cUser);
                                for jUser = 1:usersPerCell(jBase,1)
                                    rUser = cellUserIndices{jBase,1}(jUser,1);
                                    if rUser ~= cUser
                                        for jLayer = 1:maxRank
                                            intVector = [intVector ; vW{cUser,iBand}(:,iLayer)' * currentH * M(:,jLayer,rUser)];
                                        end
                                    else
                                        for jLayer = 1:maxRank
                                            if jLayer ~= iLayer
                                                intVector = [intVector ; vW{cUser,iBand}(:,iLayer)' * currentH * M(:,jLayer,rUser)];
                                            end
                                        end
                                    end
                                end
                            end
                            
                            norm(intVector,2) <= sqrt(b(iLayer,cUser));
                            log(1 + g(iLayer,cUser)) >= t(iLayer,cUser) * log(2);
                            
                            currentH = cH{iBase,iBand}(:,:,cUser);
                            p(iLayer,cUser) = real(vW{cUser,iBand}(:,iLayer)' * currentH * M(:,iLayer,cUser));
                            q(iLayer,cUser) = imag(vW{cUser,iBand}(:,iLayer)' * currentH * M(:,iLayer,cUser));
                            
                            (p_o(iLayer,cUser)^2 + q_o(iLayer,cUser)^2) / (b_o(iLayer,cUser)) + ...
                                (2 / b_o(iLayer,cUser)) * (p_o(iLayer,cUser) * (p(iLayer,cUser) - p_o(iLayer,cUser))) + ...
                                (2 / b_o(iLayer,cUser)) * (q_o(iLayer,cUser) * (q(iLayer,cUser) - q_o(iLayer,cUser))) - ...
                                (p_o(iLayer,cUser)^2 + q_o(iLayer,cUser)^2) / (b_o(iLayer,cUser)^2) * ...
                                (b(iLayer,cUser) - b_o(iLayer,cUser)) >= g(iLayer,cUser);
                            
                        end
                        
                    end
                    
                    norm(vec(M(:,:,cellUserIndices{iBase,1})),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand));                    
                end
                
                cvx_end
                
                if strfind(cvx_status,'Solved')
                    
                    b_o = b;
                    for iBase = 1:nBases
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            currentH = cH{iBase,iBand}(:,:,cUser);
                            for iLayer = 1:maxRank
                                p_o(iLayer,cUser) = real(vW{cUser,iBand}(:,iLayer)' * currentH * M(:,iLayer,cUser));
                                q_o(iLayer,cUser) = imag(vW{cUser,iBand}(:,iLayer)' * currentH * M(:,iLayer,cUser));
                            end
                            
                            qDeviation = max(QueuedPkts(cUser,1) - sum(vec(t(:,cUser,:))),0);
                            SimParams.Debug.tempResource{2,SimParams.iDrop}{cUser,1} = [SimParams.Debug.tempResource{2,SimParams.iDrop}{cUser,1} sum(vec(t(:,cUser,:))) + sum(bandRateMax(cUser,1:(iBand - 1)))];
                            SimParams.Debug.tempResource{3,SimParams.iDrop}{cUser,1} = [SimParams.Debug.tempResource{3,SimParams.iDrop}{cUser,1} qDeviation];
                            SimParams.Debug.tempResource{4,SimParams.iDrop}{cUser,iBand} = [SimParams.Debug.tempResource{4,SimParams.iDrop}{cUser,iBand} sum(vec(t(:,cUser,1)))];
                        end
                    end
                    
                    for iBase = 1:nBases
                        for iUser = 1:usersPerCell(iBase,1)
                            cUser = cellUserIndices{iBase,1}(iUser,1);
                            for iLayer = 1:maxRank
                                R = SimParams.N * eye(SimParams.nRxAntenna);
                                for jBase = 1:nBases
                                    for jUser = 1:usersPerCell(jBase,1)
                                        rUser = cellUserIndices{jBase,1}(jUser,1);
                                        H = cH{jBase,iBand}(:,:,cUser);
                                        R = R + H * M(:,:,rUser) * M(:,:,rUser)' * H';
                                    end
                                end
                                H = cH{iBase,iBand}(:,:,cUser);
                                vW{cUser,iBand}(:,iLayer) = R \ (H * M(:,iLayer,cUser));
                            end
                        end
                    end
                    
                    if min(abs(cvx_optval - cvx_hist)) <= epsilonT
                        reIterate = 0;
                    else
                        xIndex = xIndex + 1;
                        cvx_hist(mod(xIndex,2) + 1,1) = cvx_optval;
                    end
                else
                    b_o = b_o * 2;
                    display('Failed CVX !');
                end
                
                currentIteration = currentIteration + 1;
                if currentIteration >= maxIterations
                    reIterate = 0;
                end
                
            end
            
            for iBase = 1:nBases
                P = [];
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = zeros(SimParams.nTxAntenna,usersPerCell(iBase,1));
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    P = [P M(:,:,cUser)];
                end
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = P;
            end
            
            for iUser = 1:nUsers
                SimStructs.userStruct{iUser,1}.W{iBand,1} = vW{iUser,iBand};
            end

            [SimParams,SimStructs] = performDummyReception(SimParams,SimStructs,iBand);
            QueuedPkts = max(QueuedPkts - SimParams.Debug.privateExchanges.resAllocation(iBand,:)',0);
            bandRateMax(:,iBand) = SimParams.Debug.privateExchanges.resAllocation(iBand,:)';
            
        end
        
    case 'GenMSEAlloc'
        
        maxRank = SimParams.maxRank;
        
        xIndex = 0;
        reIterate = 1;
        currentIteration = 0;
        cvx_hist = -500 * ones(2,1);
        [mseError_o,vW] = randomizeInitialMSESCApoint(SimParams,SimStructs);
        
        while reIterate
            
            cvx_begin
            
            expression givenVector
            variable M(SimParams.nTxAntenna,maxRank,nUsers,nBands) complex
            variables t(maxRank,nUsers,nBands) e(maxRank,nUsers,nBands) mseError(maxRank,nUsers,nBands)
            variables userObjective(nUsers,1) epiObjective
            
            minimize(epiObjective)
            
            subject to
            
            for iUser = 1:nUsers
                userWts(iUser,1) * abs(QueuedPkts(iUser,1) - sum(vec(t(:,iUser,:)))) <= userObjective(iUser,1);
            end
            
            epiObjective >= norm(userObjective,qExponent);
            
            for iBand = 1:nBands
                for iUser = 1:nUsers
                    baseNode = SimStructs.userStruct{iUser,1}.baseNode;
                    for iLayer = 1:maxRank
                        intVector = sqrt(SimParams.N) * vW{iUser,iBand}(:,iLayer);
                        for jUser = 1:nUsers
                            ifNode = SimStructs.userStruct{jUser,1}.baseNode;
                            currentH = cH{ifNode,iBand}(:,:,iUser);
                            if jUser == iUser
                                for jLayer = 1:maxRank
                                    if jLayer ~= iLayer
                                        intVector = [intVector ; vW{iUser,iBand}(:,iLayer)' * currentH * M(:,jLayer,jUser,iBand)];
                                    end
                                end
                            else
                                for jLayer = 1:maxRank
                                    intVector = [intVector ; vW{iUser,iBand}(:,iLayer)' * currentH * M(:,jLayer,jUser,iBand)];
                                end
                            end
                        end
                        
                        currentH = cH{baseNode,iBand}(:,:,iUser);
                        givenVector = (1 - vW{iUser,iBand}(:,iLayer)' * currentH * M(:,iLayer,iUser,iBand));
                        intVector = [intVector ; givenVector];
                        norm(intVector,2) <= sqrt(mseError(iLayer,iUser,iBand));
                        (mseError(iLayer,iUser,iBand) - mseError_o(iLayer,iUser,iBand)) / mseError_o(iLayer,iUser,iBand) + log(mseError_o(iLayer,iUser,iBand)) <= -t(iLayer,iUser,iBand) * log(2);
                        
                    end
                end
                
            end
            
            for iBase = 1:nBases
                if strcmp(globalMode,'false')
                    for iBand = 1:nBands
                        norm(vec(M(:,:,cellUserIndices{iBase,1},iBand)),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
                    end
                else
                    norm(vec(M(:,:,cellUserIndices{iBase,1},:)),2) <= sqrt(sum(SimStructs.baseStruct{iBase,1}.sPower(1,:)));
                end
            end
            
            cvx_end
            
            if strfind(cvx_status,'Solved')
                
                mseError_o = mseError;
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
                                        R = R + H * M(:,:,rUser,iBand) * M(:,:,rUser,iBand)' * H';
                                    end
                                end
                                H = cH{iBase,iBand}(:,:,cUser);
                                vW{cUser,iBand}(:,iLayer) = R \ (H * M(:,iLayer,cUser,iBand));
                            end
                        end
                    end
                end
                
                if min(abs(cvx_optval - cvx_hist)) <= epsilonT
                    reIterate = 0;
                else
                    xIndex = xIndex + 1;
                    cvx_hist(mod(xIndex,2) + 1,1) = cvx_optval;
                end
                
            else
                display('Failed CVX !');
            end
            
            currentIteration = currentIteration + 1;
            if currentIteration >= maxIterations
                reIterate = 0;
            end
            
            [SimParams,SimStructs] = updateIteratePerformance(SimParams,SimStructs,M,vW);
            
        end
        
        for iBase = 1:nBases
            for iBand = 1:nBands
                P = [];
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = zeros(SimParams.nTxAntenna,usersPerCell(iBase,1));
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    P = [P M(:,:,cUser,iBand)];
                end
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = P;
            end
        end
        
    case 'JointWSRMAlloc'
        
        maxRank = SimParams.maxRank;
        
        xIndex = 0;
        reIterate = 1;
        currentIteration = 0;
        cvx_hist = -500 * ones(2,1);
        [mseError_o,vW] = randomizeInitialMSESCApoint(SimParams,SimStructs);
        
        while reIterate
            
            cvx_begin
            
            expression givenVector
            variable M(SimParams.nTxAntenna,maxRank,nUsers,nBands) complex
            variables t(maxRank,nUsers,nBands) e(maxRank,nUsers,nBands) mseError(maxRank,nUsers,nBands)
            variables userObjective(nUsers,1) epiObjective
            
            maximize(epiObjective)
            
            subject to
            
            for iUser = 1:nUsers
                userWts(iUser,1) * QueuedPkts(iUser,1) * (sum(vec(t(:,iUser,:)))) >= userObjective(iUser,1);
            end
            
            epiObjective <= sum(userObjective);
            
            for iBand = 1:nBands
                for iUser = 1:nUsers
                    baseNode = SimStructs.userStruct{iUser,1}.baseNode;
                    for iLayer = 1:maxRank
                        intVector = sqrt(SimParams.N) * vW{iUser,iBand}(:,iLayer);
                        for jUser = 1:nUsers
                            ifNode = SimStructs.userStruct{jUser,1}.baseNode;
                            currentH = cH{ifNode,iBand}(:,:,iUser);
                            if jUser == iUser
                                for jLayer = 1:maxRank
                                    if jLayer ~= iLayer
                                        intVector = [intVector ; vW{iUser,iBand}(:,iLayer)' * currentH * M(:,jLayer,jUser,iBand)];
                                    end
                                end
                            else
                                for jLayer = 1:maxRank
                                    intVector = [intVector ; vW{iUser,iBand}(:,iLayer)' * currentH * M(:,jLayer,jUser,iBand)];
                                end
                            end
                        end
                        
                        currentH = cH{baseNode,iBand}(:,:,iUser);
                        givenVector = (1 - vW{iUser,iBand}(:,iLayer)' * currentH * M(:,iLayer,iUser,iBand));
                        intVector = [intVector ; givenVector];
                        norm(intVector,2) <= sqrt(mseError(iLayer,iUser,iBand));
                        (mseError(iLayer,iUser,iBand) - mseError_o(iLayer,iUser,iBand)) / mseError_o(iLayer,iUser,iBand) + log(mseError_o(iLayer,iUser,iBand)) <= -t(iLayer,iUser,iBand) * log(2);
                        
                    end
                end
                
            end
            
            for iUser = 1:nUsers
                sum(vec(t(:,iUser,:))) <= QueuedPkts(iUser,1);
            end
            
            for iBase = 1:nBases
                if strcmp(globalMode,'false')
                    for iBand = 1:nBands
                        norm(vec(M(:,:,cellUserIndices{iBase,1},iBand)),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
                    end
                else
                    norm(vec(M(:,:,cellUserIndices{iBase,1},:)),2) <= sqrt(sum(SimStructs.baseStruct{iBase,1}.sPower(1,:)));
                end
            end
            
            cvx_end
            
            if strfind(cvx_status,'Solved')
                
                mseError_o = mseError;
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
                                        R = R + H * M(:,:,rUser,iBand) * M(:,:,rUser,iBand)' * H';
                                    end
                                end
                                H = cH{iBase,iBand}(:,:,cUser);
                                vW{cUser,iBand}(:,iLayer) = R \ (H * M(:,iLayer,cUser,iBand));
                            end
                        end
                    end
                end
                
                if min(abs(cvx_optval - cvx_hist)) <= epsilonT
                    reIterate = 0;
                else
                    xIndex = xIndex + 1;
                    cvx_hist(mod(xIndex,2) + 1,1) = cvx_optval;
                end
                
            else
                display('Failed CVX !');
            end
            
            currentIteration = currentIteration + 1;
            if currentIteration >= maxIterations
                reIterate = 0;
            end
            
            [SimParams,SimStructs] = updateIteratePerformance(SimParams,SimStructs,M,vW);
            
        end
        
        for iBase = 1:nBases
            for iBand = 1:nBands
                P = [];
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = zeros(SimParams.nTxAntenna,usersPerCell(iBase,1));
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    P = [P M(:,:,cUser,iBand)];
                end
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = P;
            end
        end
        
        
end

for iUser = 1:nUsers
    SimParams.Debug.tempResource{2,SimParams.iDrop}{iUser,1} = SimParams.Debug.tempResource{2,SimParams.iDrop}{iUser,1};
    for iBand = 1:nBands
        SimParams.Debug.tempResource{4,SimParams.iDrop}{iUser,iBand} = SimParams.Debug.tempResource{4,SimParams.iDrop}{iUser,iBand};
        if strcmp(updatePrecoders,'true')
            SimStructs.userStruct{iUser,1}.W{iBand,1} = vW{iUser,iBand};
        end
    end
end

end

% if 0
%
%         case 'JointAlloc'
%
%         xIndex = 0;
%         reIterate = 1;
%         currentIteration = 0;
%         cvx_hist = -500 * ones(2,1);
%         [p_o,q_o,b_o] = randomizeInitialSCApoint(SimParams,SimStructs);
%         p_o = reshape(p_o,nUsers,nBands);q_o = reshape(q_o,nUsers,nBands);b_o = reshape(b_o,nUsers,nBands);
%
%         while reIterate
%
%             cvx_begin
%
%             expressions p(nUsers,nBands) q(nUsers,nBands)
%             variable M(SimParams.nTxAntenna,nUsers,nBands) complex
%             variables t(nUsers,nBands) b(nUsers,nBands) g(nUsers,nBands)
%             variables userObjective(nUsers,1) epiObjective
%
%             minimize(epiObjective)
%
%             subject to
%
%             for iUser = 1:nUsers
%                 userWts(iUser,1) * abs(QueuedPkts(iUser,1) - sum(t(iUser,:))) <= userObjective(iUser,1);
%             end
%
%             epiObjective >= norm(userObjective,qExponent);
%
%             for iBase = 1:nBases
%                 for iBand = 1:nBands
%                     for iUser = 1:usersPerCell(iBase,1)
%                         cUser = cellUserIndices{iBase,1}(iUser,1);
%                         intVector = sqrt(SimParams.N);
%
%                         for jBase = 1:nBases
%                             currentH = cH{jBase,iBand}(:,:,cUser);
%                             for jUser = 1:usersPerCell(jBase,1)
%                                 rUser = cellUserIndices{jBase,1}(jUser,1);
%                                 if rUser ~= cUser
%                                     intVector = [intVector ; currentH * M(:,rUser,iBand)];
%                                 end
%                             end
%                         end
%
%                         norm(intVector,2) <= sqrt(b(cUser,iBand));
%                         log(1 + g(cUser,iBand)) >= t(cUser,iBand) * log(2);
%
%                         currentH = cH{iBase,iBand}(:,:,cUser);
%                         p(cUser,iBand) = real(currentH * M(:,cUser,iBand));
%                         q(cUser,iBand) = imag(currentH * M(:,cUser,iBand));
%
%                         (p_o(cUser,iBand)^2 + q_o(cUser,iBand)^2) / (b_o(cUser,iBand)) + ...
%                             (2 / b_o(cUser,iBand)) * (p_o(cUser,iBand) * (p(cUser,iBand) - p_o(cUser,iBand))) + ...
%                             (2 / b_o(cUser,iBand)) * (q_o(cUser,iBand) * (q(cUser,iBand) - q_o(cUser,iBand))) - ...
%                             (p_o(cUser,iBand)^2 + q_o(cUser,iBand)^2) / (2 * b_o(cUser,iBand)^2) * ...
%                             (b(cUser,iBand) - b_o(cUser,iBand)) >= g(cUser,iBand);
%
%                     end
%
%                     norm(vec(M(:,cellUserIndices{iBase,1},iBand)),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
%
%                 end
%
%             end
%
%             cvx_end
%
%             if strfind(cvx_status,'Solved')
%
%                 b_o = b;
%                 M = full(M);
%                 for iBand = 1:nBands
%                     for iBase = 1:nBases
%                         for iUser = 1:usersPerCell(iBase,1)
%                             cUser = cellUserIndices{iBase,1}(iUser,1);
%                             currentH = cH{iBase,iBand}(:,:,cUser);
%                             p_o(cUser,iBand) = real(currentH * M(:,cUser,iBand));
%                             q_o(cUser,iBand) = imag(currentH * M(:,cUser,iBand));
%
%                             if iBand == 1
%                                 qDeviation = max(QueuedPkts(cUser,1) - sum(vec(t(cUser,:))),0);
%                                 SimParams.Debug.tempResource{2,SimParams.iDrop}{cUser,1} = [SimParams.Debug.tempResource{2,SimParams.iDrop}{cUser,1} sum(vec(t(cUser,:)))];
%                                 SimParams.Debug.tempResource{3,SimParams.iDrop}{cUser,1} = [SimParams.Debug.tempResource{3,SimParams.iDrop}{cUser,1} qDeviation];
%                             end
%                             SimParams.Debug.tempResource{4,SimParams.iDrop}{cUser,iBand} = [SimParams.Debug.tempResource{4,SimParams.iDrop}{cUser,iBand} t(cUser,iBand)];
%                         end
%                     end
%                 end
%
%                 if min(abs(cvx_optval - cvx_hist)) <= epsilonT
%                     reIterate = 0;
%                 else
%                     xIndex = xIndex + 1;
%                     cvx_hist(mod(xIndex,2) + 1,1) = cvx_optval;
%                 end
%             else
%                 b_o = b_o * 2;
%                 display('Failed CVX !');
%             end
%
%             currentIteration = currentIteration + 1;
%             if currentIteration > maxIterations
%                 reIterate = 0;
%             end
%
%         end
%
%         updatePrecoders = 'false';
%         for iBase = 1:nBases
%             for iBand = 1:nBands
%                 SimStructs.baseStruct{iBase,1}.P{iBand,1} = zeros(SimParams.nTxAntenna,usersPerCell(iBase,1));
%                 for iUser = 1:usersPerCell(iBase,1)
%                     cUser = cellUserIndices{iBase,1}(iUser,1);
%                     SimStructs.baseStruct{iBase,1}.P{iBand,1}(:,iUser) = M(:,cUser,iBand);
%                 end
%             end
%         end
%
%
%     case 'GMApproxAlloc'
%
%         reIterate = 1;
%         cvx_hist = -500 * ones(2,1);
%
%         xIndex = 0;
%         currentIteration = 0;
%         p_o = ones(nUsers,nBands) / nUsers;
%         q_o = ones(nUsers,nBands) / nUsers;
%         phi = ones(nUsers,nBands) * 10 + rand(nUsers,nBands);
%
%         while reIterate
%
%             cvx_begin
%
%             expressions p(nUsers,nBands) q(nUsers,nBands)
%             variable M(SimParams.nTxAntenna,nUsers,nBands) complex
%             variables t(nUsers,nBands) b(nUsers,nBands) g(nUsers,nBands)
%             variables userObjective(nUsers,1) epiObjective
%
%             minimize(epiObjective)
%
%             subject to
%
%             for iUser = 1:nUsers
%                 userWts(iUser,1) * abs(QueuedPkts(iUser,1) - sum(t(iUser,:))) <= userObjective(iUser,1);
%             end
%
%             epiObjective >= norm(userObjective,qExponent);
%
%             for iBase = 1:nBases
%                 for iBand = 1:nBands
%                     for iUser = 1:usersPerCell(iBase,1)
%                         cUser = cellUserIndices{iBase,1}(iUser,1);
%                         intVector = sqrt(SimParams.N);
%
%                         for jBase = 1:nBases
%                             currentH = cH{jBase,iBand}(:,:,cUser);
%                             for jUser = 1:usersPerCell(jBase,1)
%                                 rUser = cellUserIndices{jBase,1}(jUser,1);
%                                 if rUser ~= cUser
%                                     intVector = [intVector ; currentH * M(:,rUser,iBand)];
%                                 end
%                             end
%                         end
%
%                         norm(intVector,2) <= sqrt(b(cUser,iBand));
%                         log(g(cUser,iBand)) >= t(cUser,iBand) * log(2);
%
%                         currentH = cH{iBase,iBand}(:,:,cUser);
%                         p(cUser,iBand) = real(currentH * M(:,cUser,iBand));
%                         q(cUser,iBand) = imag(currentH * M(:,cUser,iBand));
%
%                         p_o(cUser,iBand)^2 + q_o(cUser,iBand)^2 + ...
%                             2 * (p_o(cUser,iBand) * (p(cUser,iBand) - p_o(cUser,iBand))) + ...
%                             2 * (q_o(cUser,iBand) * (q(cUser,iBand) - q_o(cUser,iBand))) + b(cUser,iBand) >= ...
%                             0.5 * (phi(cUser,iBand) * g(cUser,iBand)^2 + (1 / phi(cUser,iBand)) * b(cUser,iBand)^2);
%
%                     end
%
%
%
%                     norm(vec(M(:,cellUserIndices{iBase,1},iBand)),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
%
%                 end
%
%                 if strcmp(globalMode,'false')
%                     for iBand = 1:nBands
%                         norm(vec(M(:,cellUserIndices{iBase,1},iBand)),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
%                     end
%                 else
%                     norm(vec(M(:,cellUserIndices{iBase,1},:)),2) <= sqrt(sum(SimStructs.baseStruct{iBase,1}.sPower(1,:)));
%                 end
%
%             end
%
%             cvx_end
%
%             if strfind(cvx_status,'Solved')
%
%                 phi = b ./ g;
%                 for iBand = 1:nBands
%                     for iBase = 1:nBases
%                         for iUser = 1:usersPerCell(iBase,1)
%                             cUser = cellUserIndices{iBase,1}(iUser,1);
%                             currentH = cH{iBase,iBand}(:,:,cUser);
%                             p_o(cUser,iBand) = real(currentH * M(:,cUser,iBand));
%                             q_o(cUser,iBand) = imag(currentH * M(:,cUser,iBand));
%
%                             if iBand == 1
%                                 qDeviation = max(QueuedPkts(cUser,1) - sum(vec(t(cUser,:))),0);
%                                 SimParams.Debug.tempResource{2,SimParams.iDrop}{cUser,1} = [SimParams.Debug.tempResource{2,SimParams.iDrop}{cUser,1} sum(vec(t(cUser,:)))];
%                                 SimParams.Debug.tempResource{3,SimParams.iDrop}{cUser,1} = [SimParams.Debug.tempResource{3,SimParams.iDrop}{cUser,1} qDeviation];
%                             end
%
%                             SimParams.Debug.tempResource{4,SimParams.iDrop}{cUser,iBand} = [SimParams.Debug.tempResource{4,SimParams.iDrop}{cUser,iBand} t(cUser,iBand)];
%
%                         end
%                     end
%                 end
%
%                 if min(abs(cvx_optval - cvx_hist)) <= epsilonT
%                     reIterate = 0;
%                 else
%                     xIndex = xIndex + 1;
%                     cvx_hist(mod(xIndex,2) + 1,1) = cvx_optval;
%                 end
%             else
%                 display('Failed CVX !');
%                 phi = ones(nUsers,nBands) * 10 + rand(nUsers,nBands);
%             end
%
%             currentIteration = currentIteration + 1;
%             if currentIteration > maxIterations
%                 reIterate = 0;
%             end
%
%         end
%
%         updatePrecoders = 'false';
%         for iBase = 1:nBases
%             for iBand = 1:nBands
%                 SimStructs.baseStruct{iBase,1}.P{iBand,1} = zeros(SimParams.nTxAntenna,usersPerCell(iBase,1));
%                 for iUser = 1:usersPerCell(iBase,1)
%                     cUser = cellUserIndices{iBase,1}(iUser,1);
%                     SimStructs.baseStruct{iBase,1}.P{iBand,1}(:,iUser) = M(:,cUser,iBand);
%                 end
%             end
%         end
%
%     case 'JointBandAlloc'
%
%         for iBand = 1:nBands
%
%             xIndex = 0;
%             reIterate = 1;
%             currentIteration = 0;
%             cvx_hist = -500 * ones(2,1);
%             [p_o,q_o,b_o] = randomizeInitialSCApoint(SimParams,SimStructs,iBand);
%             p_o = squeeze(p_o)';q_o = squeeze(q_o)';b_o = squeeze(b_o)';
%
%             while reIterate
%
%                 cvx_begin
%
%                 expressions p(nUsers,1) q(nUsers,1)
%                 variable M(SimParams.nTxAntenna,nUsers) complex
%                 variables t(nUsers,1) b(nUsers,1) g(nUsers,1)
%                 variables userObjective(nUsers,1) epiObjective
%
%                 minimize(epiObjective)
%
%                 subject to
%
%                 for iUser = 1:nUsers
%                     userWts(iUser,1) * abs(QueuedPkts(iUser,1) - t(iUser,1)) <= userObjective(iUser,1);
%                 end
%
%                 epiObjective >= norm(userObjective,qExponent);
%
%                 for iBase = 1:nBases
%                     for iUser = 1:usersPerCell(iBase,1)
%                         cUser = cellUserIndices{iBase,1}(iUser,1);
%                         intVector = sqrt(SimParams.N);
%
%                         for jBase = 1:nBases
%                             currentH = cH{jBase,iBand}(:,:,cUser);
%                             for jUser = 1:usersPerCell(jBase,1)
%                                 rUser = cellUserIndices{jBase,1}(jUser,1);
%                                 if rUser ~= cUser
%                                     intVector = [intVector ; currentH * M(:,rUser)];
%                                 end
%                             end
%                         end
%
%                         norm(intVector,2) <= sqrt(b(cUser,1));
%                         log(1 + g(cUser,1)) >= t(cUser,1) * log(2);
%
%                         currentH = cH{iBase,iBand}(:,:,cUser);
%                         p(cUser,1) = real(currentH * M(:,cUser));
%                         q(cUser,1) = imag(currentH * M(:,cUser));
%
%                         (p_o(cUser,1)^2 + q_o(cUser,1)^2) / (b_o(cUser,1)) + ...
%                             (2 / b_o(cUser,1)) * (p_o(cUser,1) * (p(cUser,1) - p_o(cUser,1))) + ...
%                             (2 / b_o(cUser,1)) * (q_o(cUser,1) * (q(cUser,1) - q_o(cUser,1))) - ...
%                             (p_o(cUser,1)^2 + q_o(cUser,1)^2) / (2 * b_o(cUser,1)^2) * ...
%                             (b(cUser,1) - b_o(cUser,1)) >= g(cUser,1);
%
%                     end
%
%                     norm(vec(M(:,cellUserIndices{iBase,1})),2) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
%
%                 end
%
%                 cvx_end
%
%                 if strfind(cvx_status,'Solved')
%
%                     b_o = b;
%                     for iBase = 1:nBases
%                         for iUser = 1:usersPerCell(iBase,1)
%                             cUser = cellUserIndices{iBase,1}(iUser,1);
%                             currentH = cH{iBase,iBand}(:,:,cUser);
%                             p_o(cUser,1) = real(currentH * M(:,cUser));
%                             q_o(cUser,1) = imag(currentH * M(:,cUser));
%
%                             qDeviation = max(QueuedPkts(cUser,1) - sum(vec(t(cUser,:))),0);
%                             SimParams.Debug.tempResource{2,SimParams.iDrop}{cUser,1} = [SimParams.Debug.tempResource{2,SimParams.iDrop}{cUser,1} sum(vec(t(cUser,:))) + sum(bandRateMax(cUser,1:(iBand-1)))];
%                             SimParams.Debug.tempResource{3,SimParams.iDrop}{cUser,1} = [SimParams.Debug.tempResource{3,SimParams.iDrop}{cUser,1} qDeviation];
%                             SimParams.Debug.tempResource{4,SimParams.iDrop}{cUser,iBand} = [SimParams.Debug.tempResource{4,SimParams.iDrop}{cUser,iBand} t(cUser,1)];
%                         end
%                     end
%
%                     if min(abs(cvx_optval - cvx_hist)) <= epsilonT
%                         reIterate = 0;
%                     else
%                         xIndex = xIndex + 1;
%                         cvx_hist(mod(xIndex,2) + 1,1) = cvx_optval;
%                     end
%                 else
%                     b_o = b_o * 2;
%                     display('Failed CVX !');
%                 end
%
%                 currentIteration = currentIteration + 1;
%                 if currentIteration > maxIterations
%                     reIterate = 0;
%                 end
%
%             end
%
%             for iBase = 1:nBases
%                 SimStructs.baseStruct{iBase,1}.P{iBand,1} = zeros(SimParams.nTxAntenna,usersPerCell(iBase,1));
%                 for iUser = 1:usersPerCell(iBase,1)
%                     cUser = cellUserIndices{iBase,1}(iUser,1);
%                     SimStructs.baseStruct{iBase,1}.P{iBand,1}(:,iUser) = M(:,cUser);
%                 end
%             end
%
%             [SimParams,SimStructs] = performDummyReception(SimParams,SimStructs,iBand);
%             QueuedPkts = max(QueuedPkts - SimParams.Debug.privateExchanges.resAllocation(iBand,:)',0);
%             bandRateMax(:,iBand) = SimParams.Debug.privateExchanges.resAllocation(iBand,:)';
%         end
%
%         updatePrecoders = 'false';
%
% end
