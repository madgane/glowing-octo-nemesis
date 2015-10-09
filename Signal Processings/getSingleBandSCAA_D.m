
function [SimParams,SimStructs] = getSingleBandSCAA_D(SimParams,SimStructs)

initMultiCastVariables;
rX = SimParams.Debug.tempResource{2,1}{1,1};
iX = SimParams.Debug.tempResource{3,1}{1,1};

iterateSCA = 1;
iIterateSCA = 0;
minPower = 1e20;

binVar_P = cell(nBases,1);
invVariable = cell(nBases,1);
for iBase = 1:nBases
    binVar_P{iBase,1} = ones(SimParams.nTxAntenna,1);    
    invVariable{iBase,1} = ones(SimParams.nTxAntenna,1);
end

if SimParams.nTxAntennaEnabled == SimParams.nTxAntenna
    return;
end

while iterateSCA
    
    iBand = 1;
    gConstraints = [];
    X = cell(nBases,1);
    binVar = cell(nBases,1);
    aPowVar = cell(nBases,1);
    
    for iBase = 1:nBases
        binVar{iBase,1} = binvar(SimParams.nTxAntenna,1,'full');
        aPowVar{iBase,1} = sdpvar(SimParams.nTxAntenna,1,'full');
        X{iBase,iBand} = sdpvar(SimParams.nTxAntenna,nGroupsPerCell(iBase,1),'full','complex');
    end
    
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
    
    for iBase = 1:nBases
        for iAntenna = 1:SimParams.nTxAntenna
            tVec = [];
            for iBand = 1:nBands
                tVec = [tVec , X{iBase,iBand}(iAntenna,:)];
            end
            
            tempVector = [2 * tVec, (aPowVar{iBase,1}(iAntenna,1) - binVar{iBase,1}(iAntenna,1))];
            gConstraints = [gConstraints, cone(tempVector,(aPowVar{iBase,1}(iAntenna,1) + binVar{iBase,1}(iAntenna,1) * binVar_P{iBase,1}(iAntenna,1)))];
        end
        
        gConstraints = [gConstraints, 0 <= binVar{iBase,1} <= 1];
        gConstraints = [gConstraints, sum(binVar{iBase,1}) == SimParams.nTxAntennaEnabled];
    end
    
    objective = 0;
    for iBase = 1:nBases
        objective = objective + sum(aPowVar{iBase,1});
    end
    
    objective = objective * objWeight + sum(invVariable{iBase,1} .* (binVar{iBase,1} - binVar_P{iBase,1}));
    
    options = sdpsettings('verbose',0,'solver','Mosek');
    solverOut = optimize(gConstraints,objective,options);
    SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray) = solverOut.solvertime + SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray);
    
    if ((solverOut.problem == 0) || (solverOut.problem == 3) || (solverOut.problem == 4))
        
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
            
            binVar_P{iBase,1} = abs(value(binVar{iBase,1})) + sqrt(lowEpsilon);
            invVariable{iBase,1} = 1./binVar_P{iBase,1};            
            nEnabledAntenna = sum(double(binVar{iBase,1}));
            
        end
    else
        display(solverOut);
        break;
    end
        
    objective = double(objective);
    if (abs(objective - minPower) / abs(minPower)) < epsilonT
        if enableBreak
            iterateSCA = 0;
        end
    else
        minPower = objective;
    end
    
    if iIterateSCA < iterateSCAMax
        iIterateSCA = iIterateSCA + 1;
    else
        iterateSCA = 0;
    end
    
    fprintf('Enabled Antennas - \t');
    fprintf('%2.3f \t',value(binVar{iBase,1}));
    fprintf('\nUsing [%2.2f] Active Transmit Elements, Objective is - %f \n',nEnabledAntenna,objective);    
    
    if sum(abs(value(binVar{iBase,1})) < epsilonT) == (SimParams.nTxAntenna - SimParams.nTxAntennaEnabled)
        enableBreak = 1;
    end
    
end

SimParams.Debug.tempResource{2,1}{1,1} = rX;
SimParams.Debug.tempResource{3,1}{1,1} = iX;

SimParams.Debug.MultiCastSDPExchange = cell(nBases,nBands);
for iBase = 1:nBases
    for iBand = 1:nBands
        SimParams.Debug.MultiCastSDPExchange{iBase,iBand} = logical(int8(double(binVar{iBase,iBand})));
    end
end

end