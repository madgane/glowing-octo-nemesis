
function [SimParams,SimStructs] = getMultiCastSDP(SimParams,SimStructs,nIterations)

initMultiCastVariables;
if isfield(SimParams.Debug,'MultiCastSDPExchange')
    nTxAntenna = SimParams.nTxAntennaEnabled;
    enabledAntenna = SimParams.Debug.MultiCastSDPExchange;
else
    enabledAntenna = cell(nBases,nBands);
    for iBase = 1:nBases
        for iBand = 1:nBands
            enabledAntenna{iBase,iBand} = linspace(1,SimParams.nTxAntenna,SimParams.nTxAntenna);
        end
    end
    nTxAntenna = SimParams.nTxAntenna;
end


gConstraints = [];
X = cell(nBases,nBands);
for iBand = 1:nBands
    for iBase = 1:nBases
        X{iBase,iBand} = sdpvar(nTxAntenna,nTxAntenna,nGroupsPerCell(iBase,1),'hermitian','complex');
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
                        Hsdp = cH{jBase,iBand}(:,enabledAntenna{jBase,iBand},cUser)' * cH{jBase,iBand}(:,enabledAntenna{jBase,iBand},cUser);
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

objective = 0;
for iBand = 1:nBands
    for iBase = 1:nBases
        for iGroup = 1:nGroupsPerCell(iBase,1)
            objective = objective + real(trace(X{iBase,iBand}(:,:,iGroup)));
        end
    end
end

options = sdpsettings('verbose',0,'solver','DSDP');
solverOut = solvesdp(gConstraints,objective,options);
SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray,SimParams.iGroupArray) = solverOut.solvertime + SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray,SimParams.iGroupArray);

if solverOut.problem ~= 0
    display(solverOut);
    display('SDP Failed !');
    for iBase = 1:nBases
        for iBand = 1:nBands
            SimStructs.baseStruct{iBase,1}.P_SDP{iBand,1} = double(X{iBase,iBand});
            SimStructs.baseStruct{iBase,1}.PG{iBand,1} = zeros(nTxAntenna,nGroupsPerCell(iBase,1));
        end
    end
    return;
end    

maxPower = 1e20;
Y = cell(nBases,nBands);
basePower = double(objective);

for iBase = 1:nBases
    for iBand = 1:nBands
        SimStructs.baseStruct{iBase,1}.P_SDP{iBand,1} = double(X{iBase,iBand});
    end
end

for iIterate = 1:nIterations
    
    reIterate = 1;
    while reIterate
        for iBase = 1:nBases
            for iBand = 1:nBands
                dX = double(X{iBase,iBand});
                Y{iBase,iBand} = zeros(nTxAntenna,nGroupsPerCell(iBase,1));
                for iGroup = 1:nGroupsPerCell(iBase,1)
                    [P, D] = eig(dX(:,:,iGroup)); 
                    P = P(:,(diag(D) >= epsilonT));
                    D = diag(D); D = D(D >= epsilonT);
                    randSelection = complex(randn(length(D),1),randn(length(D),1));
                    Y{iBase,iBand}(:,iGroup) = P * diag(sqrt(D)) * randSelection;
                end
            end
        end
        
        gConstraints = [];
        grpPwr = cell(nBases,1);
        for iBase = 1:nBases
            grpPwr{iBase,1} = sdpvar(nGroupsPerCell(iBase,1),nBands);
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
                                Hsdp = cH{jBase,iBand}(:,enabledAntenna{jBase,iBand},cUser);
                                if and((iBase == jBase),(iGroup == jGroup))
                                    tempSum = tempSum + grpPwr{jBase,1}(jGroup,iBand) * abs(Hsdp * Y{jBase,iBand}(:,jGroup))^2;
                                else
                                    tempSum = tempSum - reqSINRPerUser(cUser,1) * grpPwr{jBase,1}(jGroup,iBand) * abs(Hsdp * Y{jBase,iBand}(:,jGroup))^2;
                                end
                            end
                        end
                        gConstraints = [gConstraints, tempSum >= 0];
                    end
                    gConstraints = [gConstraints, grpPwr{iBase,1}(:,iBand) >= 0];
                end
            end
        end
        
        objective = 0;
        for iBase = 1:nBases
            objective = objective + sum(grpPwr{iBase,1}(:));
        end
        
        options = sdpsettings('verbose',0,'solver','linprog');
        solverOut = solvesdp(gConstraints,objective,options);
        SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray,SimParams.iGroupArray) = solverOut.solvertime + SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray,SimParams.iGroupArray);
        
        if ~solverOut.problem
            reIterate = 0;
        else
            fprintf('Linear Programming Failed ! - %d \n',iIterate);
        end
    end
    
    totalPower = double(objective);
    if totalPower < maxPower
        for iBase = 1:nBases
            grpP = double(grpPwr{iBase,1});
            for iBand = 1:nBands
                tempPrecoder = Y{iBase,iBand} * diag(sqrt(grpP(:,iBand)));
                SimStructs.baseStruct{iBase,1}.PG{iBand,1} = zeros(SimParams.nTxAntenna,nGroupsPerCell(iBase,1));
                SimStructs.baseStruct{iBase,1}.PG{iBand,1}(enabledAntenna{iBase,iBand},:) = tempPrecoder;
            end
        end
        maxPower = totalPower;
    end
    
    if (totalPower - basePower) < epsilonT
        break;
    end
    
end

end
