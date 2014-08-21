
function [SimParams,SimStructs] = getMultiCastSDPAS(SimParams,SimStructs,nIterations)

initMultiCastVariables;
U = cell(nBases,nBands);
for iBase = 1:nBases
    for iBand = 1:nBands
        U{iBase,iBand} = ones(SimParams.nTxAntenna);
    end
end

weightedDualInf = 1;
prevObjective = -100;
enabledAntenna = cell(nBases,nBands);
nEnabledAntenna = zeros(nBases,nBands);

while weightedDualInf
    
    gConstraints = [];
    X = cell(nBases,nBands);
    Xtilde = cell(nBases,nBands);
    for iBand = 1:nBands
        for iBase = 1:nBases
            Xtilde{iBase,iBand} = sdpvar(SimParams.nTxAntenna,SimParams.nTxAntenna,'full');
            X{iBase,iBand} = sdpvar(SimParams.nTxAntenna,SimParams.nTxAntenna,nGroupsPerCell(iBase,1),'hermitian','complex');
        end
    end
    
    for iBand = 1:nBands
        for iBase = 1:nBases
            for iGroup = 1:nGroupsPerCell(iBase,1)
                groupUsers = SimStructs.baseStruct{iBase,1}.mcGroup{iGroup,1};
                for iUser = 1:length(groupUsers)
                    cUser = groupUsers(iUser,1);
                    tempSum = -SimParams.N * reqSINRPerUser(cUser,1);
                    for jBase = 1:nBases
                        for jGroup = 1:nGroupsPerCell(jBase,1)
                            Hsdp = cH{jBase,iBand}(:,:,cUser)' * cH{jBase,iBand}(:,:,cUser);
                            if and((iBase == jBase),(iGroup == jGroup))
                                tempSum = tempSum + trace(Hsdp * X{jBase,iBand}(:,:,jGroup));
                            else
                                tempSum = tempSum - reqSINRPerUser(cUser,1) * trace(Hsdp * X{jBase,iBand}(:,:,jGroup));
                            end
                        end
                    end
                    gConstraints = [gConstraints, tempSum >= 0];
                end
                gConstraints = [gConstraints, X{iBase,iBand}(:,:,iGroup) >= 0];
            end
        end
    end
    
    for iBase = 1:nBases
        for iBand = 1:nBands
            for iGroup = 1:nGroupsPerCell(iBase,1)
                gConstraints = [gConstraints, Xtilde{iBase,iBand} >= abs(X{iBase,iBand}(:,:,iGroup))];
            end
        end
    end
    
    objective = 0;
    for iBand = 1:nBands
        for iBase = 1:nBases
            objective = objective + trace(Xtilde{iBase,iBand} * U{iBase,iBand});
        end
    end
    
    options = sdpsettings('verbose',0,'solver','DSDP');
    solverOut = solvesdp(gConstraints,objective,options);
    SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray,SimParams.iGroupArray) = solverOut.solvertime + SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray,SimParams.iGroupArray);
    
    if solverOut.problem
        display(solverOut);
        display('SDP Failed !');
        for iBase = 1:nBases
            for iBand = 1:nBands
                SimStructs.baseStruct{iBase,1}.P_SDP{iBand,1} = zeros(SimParams.nTxAntenna,SimParams.nTxAntenna,nGroupsPerCell(iBase,1));
                SimStructs.baseStruct{iBase,1}.PG{iBand,1} = zeros(SimParams.nTxAntenna,nGroupsPerCell(iBase,1));
            end
        end
        return;
    else
        for iBase = 1:nBases
            for iBand = 1:nBands
                nEnabledAntenna(iBase,iBand) = sum(diag(double(Xtilde{iBase,iBand})) > epsilonT);
                if nEnabledAntenna(iBase,iBand) >= SimParams.nTxAntennaEnabled
                    U{iBase,iBand} = 1./(double(Xtilde{iBase,iBand}) + epsilonT);
                end
            end
        end
    end
    
    display(nEnabledAntenna);
    objective = double(objective);
    if abs(objective - prevObjective) < epsilonT
        weightedDualInf = 0;
    else
        prevObjective = objective;
    end
    
end

bisectionBasedDualSearch = 0;
for iBase = 1:nBases
    for iBand = 1:nBands
        if nEnabledAntenna(iBase,iBand) ~= SimParams.nTxAntennaEnabled
            bisectionBasedDualSearch = 1;
        end
    end
end

if ~bisectionBasedDualSearch
    reducedSDPMultiCasting = 1;
end

if bisectionBasedDualSearch
    
    bisectionSearch = 1;
    dualLambdaMin = zeros(nBases,nBands);
    dualLambdaMax = 1e5 * ones(nBases,nBands);
    
    while bisectionSearch        
               
        gConstraints = [];
        X = cell(nBases,nBands);
        Xtilde = cell(nBases,nBands);
        dualLambda = (dualLambdaMax + dualLambdaMin) / 2;
        
        for iBand = 1:nBands
            for iBase = 1:nBases
                Xtilde{iBase,iBand} = sdpvar(SimParams.nTxAntenna,SimParams.nTxAntenna,'full');
                X{iBase,iBand} = sdpvar(SimParams.nTxAntenna,SimParams.nTxAntenna,nGroupsPerCell(iBase,1),'hermitian','complex');
            end
        end
        
        for iBand = 1:nBands
            for iBase = 1:nBases
                for iGroup = 1:nGroupsPerCell(iBase,1)
                    groupUsers = SimStructs.baseStruct{iBase,1}.mcGroup{iGroup,1};
                    for iUser = 1:length(groupUsers)
                        cUser = groupUsers(iUser,1);
                        tempSum = -SimParams.N * reqSINRPerUser(cUser,1);
                        for jBase = 1:nBases
                            for jGroup = 1:nGroupsPerCell(jBase,1)
                                Hsdp = cH{jBase,iBand}(:,:,cUser)' * cH{jBase,iBand}(:,:,cUser);
                                if and((iBase == jBase),(iGroup == jGroup))
                                    tempSum = tempSum + trace(Hsdp * X{jBase,iBand}(:,:,jGroup));
                                else
                                    tempSum = tempSum - reqSINRPerUser(cUser,1) * trace(Hsdp * X{jBase,iBand}(:,:,jGroup));
                                end
                            end
                        end
                        gConstraints = [gConstraints, tempSum >= 0];
                    end
                    gConstraints = [gConstraints, X{iBase,iBand}(:,:,iGroup) >= 0];
                end
            end
        end
        
        for iBase = 1:nBases
            for iBand = 1:nBands
                for iGroup = 1:nGroupsPerCell(iBase,1)
                    gConstraints = [gConstraints, Xtilde{iBase,iBand} >= abs(X{iBase,iBand}(:,:,iGroup))];
                end
            end
        end
        
        objective = 0;
        for iBand = 1:nBands
            for iBase = 1:nBases
                for iGroup = 1:nGroupsPerCell(iBase,1)
                    objective = objective + real(trace(X{iBase,iBand}(:,:,iGroup)));
                end
                objective = objective + trace(Xtilde{iBase,iBand} * U{iBase,iBand}) * dualLambda(iBase,iBand);
            end
        end
        
        options = sdpsettings('verbose',0,'solver','Mosek');
        solverOut = solvesdp(gConstraints,objective,options);
        SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray,SimParams.iGroupArray) = solverOut.solvertime + SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray,SimParams.iGroupArray);
        
        if solverOut.problem == 0
            for iBase = 1:nBases
                for iBand = 1:nBands
                    nEnabledAntenna(iBase,iBand) = sum(diag(double(Xtilde{iBase,iBand})) > epsilonT);
                    if nEnabledAntenna(iBase,iBand) ~= SimParams.nTxAntennaEnabled
                        if nEnabledAntenna(iBase,iBand) >= SimParams.nTxAntennaEnabled
                            dualLambdaMin(iBase,iBand) = dualLambda(iBase,iBand);
                        else
                            dualLambdaMax(iBase,iBand) = dualLambda(iBase,iBand);
                        end
                    end
                end
            end
        else
            for iBase = 1:nBases
                for iBand = 1:nBands
                    dualLambdaMax(iBase,iBand) = dualLambdaMax(iBase,iBand) * 2;
                end
            end
        end
        
        conditionMet = 0;
        display(nEnabledAntenna);
        for iBand = 1:nBands
            for iBase = 1:nBases
                if nEnabledAntenna(iBase,iBand) == SimParams.nTxAntennaEnabled
                    conditionMet = conditionMet + 1;
                end
            end
        end
        
        if conditionMet == (nBases * nBands)
            bisectionSearch = 0;
            reducedSDPMultiCasting = 1;
        end
    end
end

for iBase = 1:nBases
    for iBand = 1:nBands
        enabledAntenna{iBase,iBand} = find(diag(double(Xtilde{iBase,iBand})) > epsilonT);
    end
end

if reducedSDPMultiCasting
    SimParams.Debug.MultiCastSDPExchange = enabledAntenna;
    [SimParams, SimStructs] = getMultiCastSDP(SimParams,SimStructs,nIterations);    
end

end




