
function [SimParams,SimStructs] = getQWtdWSRM(SimParams,SimStructs)

proLogue;
updatePrecoders = 'true';
maxIterations = SimParams.nExchangesOTA;
bandRateMax = zeros(SimParams.nUsers,nBands);

switch selectionMethod
    
    case 'BandAlloc'
        
        currentDesign = SimParams.PrecodingMethod;
        currentApproach = SimParams.weightedSumRateMethod;
        
        bandChannels = zeros(SimParams.nUsers,nBands);
        for iBand = 1:nBands
            for iUser = 1:SimParams.nUsers
                baseID = SimStructs.userStruct{iUser,1}.baseNode;
                bandChannels(iUser,iBand) = norm(cH{baseID,iBand}(:,:,iUser));% * QueuedPkts(iUser,1);
            end
        end
        
        xSort = sort(bandChannels,'descend');
        [~,bandSort] = sort(xSort(1,:));

        bandSort = [4 1 3 2];        
        SimParams.PrecodingMethod = 'Best_WMMSE_Method';
        SimParams.weightedSumRateMethod = 'PerformScheduling';
        
        for iBand = 1:nBands
            for iBase = 1:nBases
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    SimStructs.userStruct{cUser,1}.weighingFactor = QueuedPkts(cUser,1);
                end
            end
            
            xBand = bandSort(1,iBand);
            SimParams.Debug.privateExchanges.takeOverBand = xBand;
            [SimParams, SimStructs] = getWeightedMMSEDesign(SimParams,SimStructs);
            
            [SimParams,SimStructs] = performDummyReception(SimParams,SimStructs,xBand);
            QueuedPkts = max(QueuedPkts - SimParams.Debug.privateExchanges.resAllocation(xBand,:)',0);            
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
        receiverMode = SimParams.additionalParams;
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
                
                M = full(M);b_o = full(b);  
                p_o = full(p);q_o = full(q);
                if min(abs(cvx_optval - cvx_hist)) <= epsilonT
                    reIterate = 0;
                else
                    xIndex = xIndex + 1;
                    cvx_hist(mod(xIndex,2) + 1,1) = cvx_optval;
                end
                
                switch receiverMode
                    
                    case 'MMSE'
                    
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
                                    vW{cUser,iBand}(:,iLayer) = R \ (H * M(:,iLayer,cUser,iBand)) + epsilonT;
                                end
                            end
                        end
                    end                    
                    
                    case 'Optimal'
                        [~,~,~,vW] = findOptimalW(SimParams,SimStructs,M,vW,p_o,q_o,b_o);
                    otherwise
                        updatePrecoders = 'false';
                        if reIterate == 1                            
                            display('No receiver selected and therefore not updating !');
                        end
                end
                
                
            else
                b_o = b_o * 2;
                display('Failed CVX !');
            end
            
%             for iBand = 1:nBands
%                 for iUser = 1:nUsers
%                     baseNode = SimStructs.userStruct{iUser,1}.baseNode;
%                     channelH = cH{baseNode,iBand}(:,:,iUser);
%                     for iLayer = 1:maxRank
%                         p_o(iLayer,iUser,iBand) = real(vW{iUser,iBand}(:,iLayer)' * channelH * M(:,iLayer,iUser,iBand));
%                         q_o(iLayer,iUser,iBand) = imag(vW{iUser,iBand}(:,iLayer)' * channelH * M(:,iLayer,iUser,iBand));
%                     end
%                 end
%             end
            
            currentIteration = currentIteration + 1;
            if currentIteration >= maxIterations
                reIterate = 0;
            end
            
            [SimParams,SimStructs] = updateIteratePerformance(SimParams,SimStructs,M,vW);
            
            if (SimParams.currentQueue <= epsilonT)
                reIterate = 0;
            end

        end
        
        for iBase = 1:nBases
            for iBand = 1:nBands
                P = zeros(SimParams.nTxAntenna,maxRank,usersPerCell(iBase,1));
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = zeros(SimParams.nTxAntenna,usersPerCell(iBase,1));
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    P(:,:,iUser) = M(:,:,cUser,iBand);
                end
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = P;
            end
        end


        [SimParams, SimStructs] = getReceiveEqualizer(SimParams,SimStructs,'MMSE');
                
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
            
            t >= 0;
            
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
            
            if (SimParams.currentQueue <= epsilonT)
                reIterate = 0;
            end
            
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
            
            if (strcmpi(SimParams.additionalParams,'Optimal'))
                for iUser = 1:nUsers
                    sum(vec(t(:,iUser,:))) <= QueuedPkts(iUser,1);
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
            
            if (SimParams.currentQueue <= epsilonT)
                reIterate = 0;
            end
            
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
        
    case 'GenAllocY'
        
        epsilon = 0.5;
        pNoise = 1e-30;
        taylorApprox = 'EXP';
        
        xIndex = 0;
        reIterate = 1;
        currentIteration = 0;
        cvx_hist = -500 * ones(2,1);
        maxRank = SimParams.maxRank;
        [p_o,q_o,b_o,vW] = randomizeInitialSCApoint(SimParams,SimStructs);
        
        if strcmpi(taylorApprox,'LOG')
            g_o = (p_o.^2 + q_o.^2)./b_o;
        else
            t_o = log(ones(size(p_o)) + (p_o.^2 + q_o.^2)./b_o);
        end
        
        while reIterate
            
            M = sdpvar(SimParams.nTxAntenna,maxRank,nUsers,nBands,'full','complex');
            t = sdpvar(maxRank,nUsers,nBands,'full');b = sdpvar(maxRank,nUsers,nBands,'full');g = sdpvar(maxRank,nUsers,nBands,'full');
            userObjective = sdpvar(nUsers,1,'full');epiObjective = sdpvar(1);

            G = [];
            for iUser = 1:nUsers
                G = [G, userWts(iUser,1) * abs(QueuedPkts(iUser,1) - sum(vec(t(:,iUser,:)))) <= userObjective(iUser,1)];
            end
            
            G = [G, epiObjective >= norm(userObjective,qExponent)];
            
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
                            
                            G = [G, intVector' * intVector <= b(iLayer,cUser,iBand)];
                            
                            if strcmpi(taylorApprox,'LOG')
                                G = [G, g(iLayer,cUser,iBand) >= g_o(iLayer,cUser,iBand) * (1 - epsilon)];
                                G = [G, g(iLayer,cUser,iBand) <= g_o(iLayer,cUser,iBand) * (1 + epsilon)];
                                G = [G, log(1 + g_o(iLayer,cUser,iBand)) + (g(iLayer,cUser,iBand) - g_o(iLayer,cUser,iBand)) / (1 + g_o(iLayer,cUser,iBand)) - (g(iLayer,cUser,iBand) - g_o(iLayer,cUser,iBand))^2 / (2 * (1 + g_o(iLayer,cUser,iBand))^2) ...
                                    >= t(iLayer,cUser,iBand) * log(2)];
                            else                                
                                G = [G, t(iLayer,cUser,iBand) >= t_o(iLayer,cUser,iBand) * (1 - epsilon)];
                                G = [G, t(iLayer,cUser,iBand) <= t_o(iLayer,cUser,iBand) * (1 + epsilon)];
                                G = [G, 1 + g(iLayer,cUser,iBand) >= ...
                                    exp(t_o(iLayer,cUser,iBand) * log(2)) * (1 + log(2) * (t(iLayer,cUser,iBand) - t_o(iLayer,cUser,iBand)) + 0.5 * log(2)^2 * (t(iLayer,cUser,iBand) - t_o(iLayer,cUser,iBand))^2)];
                            end
                            
                            currentH = cH{iBase,iBand}(:,:,cUser);
                            p = real(vW{cUser,iBand}(:,iLayer)' * currentH * M(:,iLayer,cUser,iBand));
                            q = imag(vW{cUser,iBand}(:,iLayer)' * currentH * M(:,iLayer,cUser,iBand));                            
                            G = [G, (p_o(iLayer,cUser,iBand)^2 + q_o(iLayer,cUser,iBand)^2) / (b_o(iLayer,cUser,iBand)) + ...
                                (2 / b_o(iLayer,cUser,iBand)) * (p_o(iLayer,cUser,iBand) * (p - p_o(iLayer,cUser,iBand))) + ...
                                (2 / b_o(iLayer,cUser,iBand)) * (q_o(iLayer,cUser,iBand) * (q - q_o(iLayer,cUser,iBand))) - ...
                                (p_o(iLayer,cUser,iBand)^2 + q_o(iLayer,cUser,iBand)^2) / (b_o(iLayer,cUser,iBand)^2) * ...
                                (b(iLayer,cUser,iBand) - b_o(iLayer,cUser,iBand)) >= g(iLayer,cUser,iBand)];
                            
                        end
                    end
                end
                
                if strcmp(globalMode,'false')
                    for iBand = 1:nBands
                        tempM = vec(M(:,:,cellUserIndices{iBase,1},iBand));
                        G = [G, tempM' * tempM <= (SimStructs.baseStruct{iBase,1}.sPower(1,iBand))];
                    end
                else
                    tempM = vec(M(:,:,cellUserIndices{iBase,1},:));
                    G = [G, tempM' * tempM <= sum(SimStructs.baseStruct{iBase,1}.sPower(1,:))];
                end
                
            end
            
            sdpOptions = sdpsettings('verbose',0,'solver','Gurobi');
            sdpSol = solvesdp(G,epiObjective,sdpOptions);

            if sdpSol.problem == 0
                M = full(double(M)) + pNoise;
                b_o = full(double(b));
                if strcmpi(taylorApprox,'LOG')
                    g_o = full(double(g));
                else
                    t_o = full(double(t));
                end
                for iBand = 1:nBands
                    for iUser = 1:nUsers
                        currentH = cH{SimStructs.userStruct{iUser,1}.baseNode,iBand}(:,:,iUser);
                        for iLayer = 1:maxRank
                            p_o(iLayer,iUser,iBand) = real(vW{iUser,iBand}(:,iLayer)' * currentH * M(:,iLayer,iUser,iBand));
                            q_o(iLayer,iUser,iBand) = imag(vW{iUser,iBand}(:,iLayer)' * currentH * M(:,iLayer,iUser,iBand));
                        end
                    end
                end
                
                cvx_optval = double(epiObjective);
                if min(abs(cvx_optval - cvx_hist)) <= epsilonT
                    reIterate = 0;
                else
                    xIndex = xIndex + 1;
                    cvx_hist(mod(xIndex,2) + 1,1) = cvx_optval;
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
                
            else
                b_o = b_o * 2;
                display('Failed to converge !');
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
