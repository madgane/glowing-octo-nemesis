
function [p_o,q_o,b_o,vW] = findOptimalW(SimParams,SimStructs,M,vW,p_o,q_o,b_o)

epsilonT = 1e-5;
maxIterations = 100;
cH = SimStructs.linkChan;
nBases = SimParams.nBases;
nBands = SimParams.nBands;
usersPerCell = zeros(nBases,1);
cellUserIndices = cell(nBases,1);

closedFormSolution = 1;

for iBase = 1:nBases
    for iBand = 1:nBands
        cellUserIndices{iBase,1} = [cellUserIndices{iBase,1} ; SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1}];
    end
    cellUserIndices{iBase,1} = unique(cellUserIndices{iBase,1});
    usersPerCell(iBase,1) = length(cellUserIndices{iBase,1});
end

nUsers = sum(usersPerCell);
QueuedPkts = zeros(nUsers,1);

for iBase = 1:nBases
    for iUser = 1:usersPerCell(iBase,1)
        cUser = cellUserIndices{iBase,1}(iUser,1);
        QueuedPkts(cUser,1) = SimStructs.userStruct{cUser,1}.trafficStats.backLogPkt;
    end
end

reIterate = 1;
currentIteration = 0;
maxRank = SimParams.maxRank;

if closedFormSolution
    
    cW = cell(nUsers,nBands);
    
    while reIterate        
        for iBand = 1:nBands
            for iUser = 1:nUsers
                baseNode = SimStructs.userStruct{iUser,1}.baseNode;
                for iLayer = 1:maxRank
                    oldB = SimParams.N * trace(vW{iUser,iBand}(:,iLayer) * vW{iUser,iBand}(:,iLayer)');
                    IF = SimParams.N * eye(SimParams.nRxAntenna);
                    for jUser = 1:nUsers
                        ifNode = SimStructs.userStruct{jUser,1}.baseNode;
                        currentH = cH{ifNode,iBand}(:,:,iUser);
                        if jUser ~= iUser
                            IF = IF + currentH * M(:,:,jUser,iBand) * M(:,:,jUser,iBand)' * currentH';
                            oldB = oldB + trace(vW{iUser,iBand}(:,iLayer)' * currentH * M(:,:,jUser,iBand) * M(:,:,jUser,iBand)' * currentH' * vW{iUser,iBand}(:,iLayer));
                        else
                            IF = IF + currentH * M(:,iLayer ~= (1:maxRank),jUser,iBand) * M(:,iLayer ~= (1:maxRank),jUser,iBand)' * currentH';
                            oldB = oldB + trace(vW{iUser,iBand}(:,iLayer)' * currentH * M(:,iLayer ~= (1:maxRank),jUser,iBand) * M(:,iLayer ~= (1:maxRank),jUser,iBand)' * currentH' * vW{iUser,iBand}(:,iLayer));
                        end
                    end
                    currentH = cH{baseNode,iBand}(:,:,iUser);
                    oldS = trace(vW{iUser,iBand}(:,iLayer)' * currentH * M(:,iLayer,iUser,iBand) * M(:,iLayer,iUser,iBand)' * currentH' * vW{iUser,iBand}(:,iLayer));
                    FF = (real(oldB) / real(oldS)) * (currentH * M(:,iLayer,iUser,iBand) * M(:,iLayer,iUser,iBand)' * currentH') * vW{iUser,iBand}(:,iLayer);
                    if real(oldS) == 0
                        cW{iUser,iBand}(:,iLayer) = zeros(size(vW{iUser,iBand}(:,iLayer)));
                    else
                        cW{iUser,iBand}(:,iLayer) = IF \ FF;
                    end
                end
            end
        end
        
        cwmat = cell2mat(cW);
        vwmat = cell2mat(vW);
        
        if norm(cwmat(:) - vwmat(:),2) < epsilonT
            reIterate = 0;
        end
        
        vW = cW;
        currentIteration = currentIteration + 1;
        if currentIteration >= maxIterations
            reIterate = 0;
        end
    end
    
    
else
    
    xIndex = 0;
    cvx_hist = -500 * ones(2,1);
    
    userWts = ones(nUsers,1);
    underscore_location = strfind(SimParams.weightedSumRateMethod,'_');
    if isempty(underscore_location)
        qExponent = 1;
    else
        qExponent = str2double(SimParams.weightedSumRateMethod(underscore_location + 1:end));
    end
    
    while reIterate
        
        cvx_begin
        
        expressions p(maxRank,nUsers,nBands) q(maxRank,nUsers,nBands)
        variable W(SimParams.nRxAntenna,maxRank,nUsers,nBands) complex
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
                        intVector = sqrt(SimParams.N) * W(:,iLayer,cUser,iBand);
                        
                        for jBase = 1:nBases
                            currentH = cH{jBase,iBand}(:,:,cUser);
                            for jUser = 1:usersPerCell(jBase,1)
                                rUser = cellUserIndices{jBase,1}(jUser,1);
                                if rUser ~= cUser
                                    for jLayer = 1:maxRank
                                        intVector = [intVector ; W(:,iLayer,cUser,iBand)' * currentH * M(:,jLayer,rUser,iBand)];
                                    end
                                else
                                    for jLayer = 1:maxRank
                                        if jLayer ~= iLayer
                                            intVector = [intVector ; W(:,iLayer,cUser,iBand)' * currentH * M(:,jLayer,rUser,iBand)];
                                        end
                                    end
                                end
                            end
                        end
                        
                        norm(intVector,2) <= sqrt(b(iLayer,cUser,iBand));
                        log(1 + g(iLayer,cUser,iBand)) >= t(iLayer,cUser,iBand) * log(2);
                        
                        currentH = cH{iBase,iBand}(:,:,cUser);
                        p(iLayer,cUser,iBand) = real(W(:,iLayer,cUser,iBand)' * currentH * M(:,iLayer,cUser,iBand));
                        q(iLayer,cUser,iBand) = imag(W(:,iLayer,cUser,iBand)' * currentH * M(:,iLayer,cUser,iBand));
                        
                        (p_o(iLayer,cUser,iBand)^2 + q_o(iLayer,cUser,iBand)^2) / (b_o(iLayer,cUser,iBand)) + ...
                            (2 / b_o(iLayer,cUser,iBand)) * (p_o(iLayer,cUser,iBand) * (p(iLayer,cUser,iBand) - p_o(iLayer,cUser,iBand))) + ...
                            (2 / b_o(iLayer,cUser,iBand)) * (q_o(iLayer,cUser,iBand) * (q(iLayer,cUser,iBand) - q_o(iLayer,cUser,iBand))) - ...
                            (p_o(iLayer,cUser,iBand)^2 + q_o(iLayer,cUser,iBand)^2) / (b_o(iLayer,cUser,iBand)^2) * ...
                            (b(iLayer,cUser,iBand) - b_o(iLayer,cUser,iBand)) >= g(iLayer,cUser,iBand);
                    end
                end
            end
        end
        
        cvx_end
        
        if strfind(cvx_status,'Solved')
            
            W = full(W);
            b_o = full(b);p_o = full(p);q_o = full(q);
            for iBand = 1:nBands
                for iBase = 1:nBases
                    for iUser = 1:usersPerCell(iBase,1)
                        cUser = cellUserIndices{iBase,1}(iUser,1);
                        for iLayer = 1:maxRank
                            vW{cUser,iBand}(:,iLayer) = W(:,iLayer,cUser,iBand);
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
            break;
        end
        
        currentIteration = currentIteration + 1;
        if currentIteration >= maxIterations
            reIterate = 0;
        end
        
        %[SimParams,SimStructs] = updateIteratePerformance(SimParams,SimStructs,M,vW);
        
    end
    
end

end