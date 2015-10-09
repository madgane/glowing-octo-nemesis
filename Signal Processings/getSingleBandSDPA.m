
function [SimParams,SimStructs] = getSingleBandSDPA(SimParams,SimStructs,nIterations)

initMultiCastVariables;

weightedDualInf = 0;
if SimParams.iAntennaArray == 1
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
    
    if (SimParams.nTxAntennaEnabled == SimParams.nTxAntenna)
        for iBase = 1:nBases
            for iBand = 1:nBands
                enabledAntenna{iBase,iBand} = linspace(1,SimParams.nTxAntenna,SimParams.nTxAntenna)';
            end
        end
        
        SimParams.Debug.MultiCastSDPExchange = enabledAntenna;
        [SimParams, SimStructs] = getMultiCastSDP(SimParams,SimStructs,nIterations);
    end    
end

while weightedDualInf
    
    iBand = 1;
    gConstraints = [];
    X = cell(nBases,1);
    Xtilde = cell(nBases,1);
    for iBase = 1:nBases
        Xtilde{iBase,iBand} = sdpvar(SimParams.nTxAntenna,SimParams.nTxAntenna,'full');
        X{iBase,iBand} = sdpvar(SimParams.nTxAntenna,SimParams.nTxAntenna,nGroupsPerCell(iBase,1),'hermitian','complex');
    end
    
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
    
    for iBase = 1:nBases
        for iGroup = 1:nGroupsPerCell(iBase,1)
            gConstraints = [gConstraints, Xtilde{iBase,iBand} >= abs(X{iBase,iBand}(:,:,iGroup))];
        end
    end
    
    objective = 0;
    for iBand = 1:nBands
        for iBase = 1:nBases
            objective = objective + trace(Xtilde{iBase,iBand} * U{iBase,iBand});
        end
    end
    
    options = sdpsettings('verbose',0,'solver','Mosek');
    solverOut = solvesdp(gConstraints,objective,options);
    SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray) = solverOut.solvertime + SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray);
    
    if ~((solverOut.problem == 0) || (solverOut.problem == 3) || (solverOut.problem == 4))
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
                U{iBase,iBand} = 1./(double(Xtilde{iBase,iBand}) + 1e-3);
            end
        end
    end
    
    objective = double(objective);
    if abs(objective - prevObjective) / abs(prevObjective) < epsilonT
        weightedDualInf = 0;
    else
        prevObjective = objective;
    end
    
    fprintf('Enabled Antenna Pattern - \n');
    fprintf('%3.4f \t',diag(value(Xtilde{1,1})));
    fprintf('\n');
    
end

if SimParams.iAntennaArray == 1
    SimParams.Debug.SDPDebug.nEnabledAntenna = nEnabledAntenna;
    SimParams.Debug.SDPDebug.X = X;
    SimParams.Debug.SDPDebug.Xtilde = Xtilde;
    SimParams.Debug.SDPDebug.U = U;
else
    nEnabledAntenna = SimParams.Debug.SDPDebug.nEnabledAntenna;
    X = SimParams.Debug.SDPDebug.X;
    Xtilde = SimParams.Debug.SDPDebug.Xtilde;
    U = SimParams.Debug.SDPDebug.U;
end

bisectionBasedDualSearch = 0;
fprintf('Minumum number of active antenna elements required - %d \n',nEnabledAntenna);

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
    
    iBand = 1;
    bisectionSearch = 1;
    
    dualLambdaMin = 0;
    dualLambdaMax = maxObj;
    
    while bisectionSearch
        
        gConstraints = [];
        X = cell(nBases,nBands);
        Xtilde = cell(nBases,nBands);
        dualLambda = (dualLambdaMax + dualLambdaMin) / 2;
        
        for iBase = 1:nBases
            Xtilde{iBase,iBand} = sdpvar(SimParams.nTxAntenna,SimParams.nTxAntenna,'full');
            X{iBase,iBand} = sdpvar(SimParams.nTxAntenna,SimParams.nTxAntenna,nGroupsPerCell(iBase,1),'hermitian','complex');
        end
        
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
        
        for iBase = 1:nBases
            for iGroup = 1:nGroupsPerCell(iBase,1)
                gConstraints = [gConstraints, Xtilde{iBase,iBand} >= abs(X{iBase,iBand}(:,:,iGroup))];
            end
        end
        
        aObjective = 0;
        objective = 0;
        for iBand = 1:nBands
            for iBase = 1:nBases
                for iGroup = 1:nGroupsPerCell(iBase,1)
                    objective = objective + real(trace(X{iBase,iBand}(:,:,iGroup)));
                end
                aObjective = aObjective + trace(Xtilde{iBase,iBand} * U{iBase,iBand});
            end
        end
        
        objective = objective + aObjective * dualLambda;
        options = sdpsettings('verbose',0,'solver','Mosek');
        solverOut = solvesdp(gConstraints,objective,options);
        SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray) = solverOut.solvertime + SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray);
        
        if ((solverOut.problem == 0) || (solverOut.problem == 3) || (solverOut.problem == 4))
            nEnabledAntenna = sum(diag(value(Xtilde{1,1})) > epsilonT);
            if nEnabledAntenna < SimParams.nTxAntennaEnabled
                dualLambdaMax = dualLambda;
            else
                dualLambdaMin = dualLambda;
            end
        else
            fprintf('Solver Problem ----- Quitting ! \n');
            yalmiperror(solverOut.problem);
        end
        
        fprintf('Enabled Antenna Pattern - \n');
        fprintf('%3.4f \t',diag(value(Xtilde{1,1})));
        fprintf('\nCurrent number of active antenna elements - %d with DualLambda of - %4.2f \n',nEnabledAntenna,dualLambda);
        
        if nEnabledAntenna(iBase,iBand) == SimParams.nTxAntennaEnabled
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
    [SimParams, SimStructs] = getSingleBandSDP(SimParams,SimStructs,nIterations);
end

end




