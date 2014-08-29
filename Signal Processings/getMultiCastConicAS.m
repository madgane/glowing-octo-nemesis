
function [SimParams,SimStructs] = getMultiCastConicAS(SimParams,SimStructs)

initMultiCastVariables;
rX = SimParams.Debug.tempResource{2,1}{1,1};
iX = SimParams.Debug.tempResource{3,1}{1,1};

iterateSCA = 1;
iIterateSCA = 0;
minPower = 1e20;
iterateSCAMax = 50;

while iterateSCA
    
    gConstraints = [];
    X = cell(nBases,nBands);
    binVar = cell(nBases,nBands);
    aPowVar = cell(nBases,nBands);
    for iBand = 1:nBands
        for iBase = 1:nBases
            binVar{iBase,iBand} = binvar(SimParams.nTxAntenna,1,'full');
            aPowVar{iBase,iBand} = sdpvar(SimParams.nTxAntenna,1,'full');
            X{iBase,iBand} = sdpvar(SimParams.nTxAntenna,nGroupsPerCell(iBase,1),'full','complex');
        end
    end
    
    for iBand = 1:nBands
        for iBase = 1:nBases
            for iGroup = 1:nGroupsPerCell(iBase,1)
                groupUsers = SimStructs.baseStruct{iBase,1}.mcGroup{iGroup,1};
                for iUser = 1:length(groupUsers)
                    cUser = groupUsers(iUser,1);
                    Hsdp = cH{iBase,iBand}(:,:,cUser);
                    
                    tempVec = [sqrt(SimParams.N)];
                    riX = real(Hsdp * X{iBase,iBand}(:,iGroup));
                    iiX = imag(Hsdp * X{iBase,iBand}(:,iGroup));
                    tempExpression = rX(cUser,iBand)^2 + iX(cUser,iBand)^2 + 2 * (rX(cUser,iBand) * (riX - rX(cUser,iBand)) + iX(cUser,iBand) * (iiX - iX(cUser,iBand)));
                    
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
                    
                    gConstraints = [gConstraints , reqSINRPerUser(cUser,1) * (tempVec' * tempVec) - tempExpression <= 0];
                end
            end
        end
    end
    
    for iBase = 1:nBases
        for iBand = 1:nBands
            for iAntenna = 1:SimParams.nTxAntenna
                tempVector = [2 * X{iBase,iBand}(iAntenna,:), (aPowVar{iBase,iBand}(iAntenna,1) - binVar{iBase,iBand}(iAntenna,1))];
                gConstraints = [gConstraints, cone(tempVector,(aPowVar{iBase,iBand}(iAntenna,1) + binVar{iBase,iBand}(iAntenna,1)))];
            end
        end
    end
    
    objective = 0;
    for iBand = 1:nBands
        for iBase = 1:nBases
            objective = objective + aPowVar{iBase,iBand}' * aPowVar{iBase,iBand} + max((sum(binVar{iBase,iBand}) - SimParams.nTxAntennaEnabled),0) * 1e3;
        end
    end
    
    options = sdpsettings('verbose',0,'solver','Mosek');
    solverOut = solvesdp(gConstraints,objective,options);
    SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray,SimParams.iGroupArray) = solverOut.solvertime + SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray,SimParams.iGroupArray);
    
    if solverOut.problem == 0
        for iBand = 1:nBands
            for iBase = 1:nBases
                for iGroup = 1:nGroupsPerCell(iBase,1)
                    groupUsers = SimStructs.baseStruct{iBase,1}.mcGroup{iGroup,1};
                    for iUser = 1:length(groupUsers)
                        cUser = groupUsers(iUser,1);
                        Hsdp = cH{iBase,iBand}(:,:,cUser);
                        rX(cUser,iBand) = real(Hsdp * double(X{iBase,iBand}(:,iGroup)));
                        iX(cUser,iBand) = imag(Hsdp * double(X{iBase,iBand}(:,iGroup)));
                    end
                end
                nEnabledAntenna = sum(double(binVar{iBase,iBand}));
                display(nEnabledAntenna);
            end
        end        
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
    if abs(objective - minPower) < 1e-2
        iterateSCA = 0;
    else
        minPower = objective;
    end
    
    display(objective);
    
end

SimParams.Debug.tempResource{2,1}{1,1} = rX;
SimParams.Debug.tempResource{3,1}{1,1} = iX;
SimParams.Debug.MultiCastSDPExchange = cell(nBases,nBands);
for iBase = 1:nBases
    for iBand = 1:nBands
        SimParams.Debug.MultiCastSDPExchange{iBase,iBand} = logical(double(binVar{iBase,iBand}));
    end
end

end