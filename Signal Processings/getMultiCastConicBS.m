
function [SimParams,SimStructs] = getMultiCastConicBS(SimParams,SimStructs)

initMultiCastVariables;
rX = SimParams.Debug.tempResource{2,1}{1,1};
iX = SimParams.Debug.tempResource{3,1}{1,1};
bX = SimParams.Debug.tempResource{4,1}{1,1};
gX = SimParams.Debug.tempResource{4,1}{2,1};

iterateSCA = 1;
iIterateSCA = 0;
minPower = 1e20;
iterateSCAMax = 50;

while iterateSCA
    
    gConstraints = [];
    X = cell(nBases,nBands);
    binVar = cell(nBases,1);
    aPowVar = cell(nBases,nBands);
    for iBase = 1:nBases
        binVar{iBase,1} = binvar(SimParams.nTxAntenna,1,'full');
        for iBand = 1:nBands
            aPowVar{iBase,iBand} = sdpvar(SimParams.nTxAntenna,1,'full');
            X{iBase,iBand} = sdpvar(SimParams.nTxAntenna,nGroupsPerCell(iBase,1),'full','complex');
        end
    end
    
    Beta = sdpvar(nUsers,nBands,'full');
    Gamma = sdpvar(nUsers,nBands,'full');
    
    for iBase = 1:nBases
        for iGroup = 1:nGroupsPerCell(iBase,1)
            groupUsers = SimStructs.baseStruct{iBase,1}.mcGroup{iGroup,1};
            for iUser = 1:length(groupUsers)
                cUser = groupUsers(iUser,1);
                feasConstraint = 0;
                for iBand = 1:nBands
                    Hsdp = cH{iBase,iBand}(:,:,cUser);
                    
                    tempVec = [sqrt(SimParams.N)];
                    riX = real(Hsdp * X{iBase,iBand}(:,iGroup));
                    iiX = imag(Hsdp * X{iBase,iBand}(:,iGroup));
                    fixedPoint = (rX(cUser,iBand)^2 + iX(cUser,iBand)^2) / bX(cUser,iBand);
                    tempExpression = fixedPoint + 2 * (rX(cUser,iBand) * (riX - rX(cUser,iBand)) + iX(cUser,iBand) * (iiX - iX(cUser,iBand))) ...
                        - (fixedPoint / bX(cUser,iBand)) * (Beta(cUser,iBand) - bX(cUser,iBand));
                    
                    for jBase = 1:nBases
                        for jGroup = 1:nGroupsPerCell(jBase,1)
                            Hsdp = cH{jBase,iBand}(:,:,cUser);
                            if ~and((iBase == jBase),(iGroup == jGroup))
                                riX = real(Hsdp * X{jBase,iBand}(:,jGroup));
                                iiX = imag(Hsdp * X{jBase,iBand}(:,jGroup));
                                tempVec = [tempVec; riX; iiX];
                            end
                        end
                    end
                    
                    fixedPoint = 1 + gX(cUser,iBand);
                    feasConstraint = feasConstraint + log(fixedPoint) + (1 / fixedPoint) * (Gamma(cUser,iBand) - gX(cUser,iBand))...
                        - (1 / (2 * fixedPoint^2)) * (Gamma(cUser,iBand) - gX(cUser,iBand))^2;
                    gConstraints = [gConstraints, Gamma(cUser,iBand) - tempExpression <= 0];
                    gConstraints = [gConstraints, (tempVec' * tempVec) - Beta(cUser,iBand) <= 0];
                end
                gConstraints = [gConstraints, QueuedNats(cUser,1) - feasConstraint <= 0];
            end
        end
    end
    
    for iBase = 1:nBases
        for iBand = 1:nBands
            for iAntenna = 1:SimParams.nTxAntenna
                tempVector = [2 * X{iBase,iBand}(iAntenna,:), (aPowVar{iBase,iBand}(iAntenna,1) - binVar{iBase,1}(iAntenna,1))];
                gConstraints = [gConstraints, cone(tempVector,(aPowVar{iBase,iBand}(iAntenna,1) + binVar{iBase,1}(iAntenna,1)))];
            end
        end
    end
    
    objective = 0;
    for iBand = 1:nBands
        for iBase = 1:nBases
            objective = objective + aPowVar{iBase,iBand}' * aPowVar{iBase,iBand};
        end
        objective = objective + max((sum(binVar{iBase,1}) - SimParams.nTxAntennaEnabled),0) * 1e3;
    end
    
    options = sdpsettings('verbose',0,'solver','Mosek');
    solverOut = solvesdp(gConstraints,objective,options);
    SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray,SimParams.iGroupArray) = solverOut.solvertime + SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray,SimParams.iGroupArray);
    
    if solverOut.problem == 0
        
        for iBase = 1:nBases
            for iGroup = 1:nGroupsPerCell(iBase,1)
                groupUsers = SimStructs.baseStruct{iBase,1}.mcGroup{iGroup,1};
                for iUser = 1:length(groupUsers)
                    cUser = groupUsers(iUser,1);
                    for iBand = 1:nBands
                        Hsdp = cH{iBase,iBand}(:,:,cUser);
                        rX(cUser,iBand) = real(Hsdp * double(X{iBase,iBand}(:,iGroup)));
                        iX(cUser,iBand) = imag(Hsdp * double(X{iBase,iBand}(:,iGroup)));
                    end
                end
            end
            nEnabledAntenna = sum(double(binVar{iBase,1}));
            display(nEnabledAntenna);
        end
        gX = full(double(Gamma));bX = full(double(Beta));
    else
        display(solverOut);
        break;
    end
    
    if iIterateSCA < iterateSCAMax
        iIterateSCA = iIterateSCA + 1;
    else
        iterateSCA = 0;
    end
    
    objective = double(objective);
    if abs(objective - minPower) < epsilonT
        iterateSCA = 0;
    else
        minPower = objective;
    end
    
    display(objective);
    
end

SimParams.Debug.tempResource{2,1}{1,1} = rX;
SimParams.Debug.tempResource{3,1}{1,1} = iX;
SimParams.Debug.tempResource{4,1}{1,1} = bX;
SimParams.Debug.tempResource{4,2}{1,1} = gX;
SimParams.Debug.MultiCastSDPExchange = cell(nBases,nBands);
for iBase = 1:nBases
    for iBand = 1:nBands
        SimParams.Debug.MultiCastSDPExchange{iBase,iBand} = logical(int8(double(binVar{iBase,1})));
    end
end

end