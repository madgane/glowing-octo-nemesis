
function [SimParams,SimStructs] = getMultiCastConicBS_B(SimParams,SimStructs)

initMultiCastVariables;
rX = SimParams.Debug.tempResource{2,1}{1,1};
iX = SimParams.Debug.tempResource{3,1}{1,1};
bX = SimParams.Debug.tempResource{4,1}{1,1};

qExponent = 2;
iterateSCA = 1;
iIterateSCA = 0;
minPower = 1e20;

binVar_P = cell(nBases,1);
for iBase = 1:nBases
    binVar_P{iBase,1} = rand(SimParams.nTxAntenna,1);
end

if SimParams.nTxAntennaEnabled == SimParams.nAntennaArray
    return;
end

while iterateSCA
    
    gConstraints = [];
    X = cell(nBases,nBands);
    binVar = cell(nBases,1);
    aPowVar = cell(nBases,1);
    
    for iBase = 1:nBases
        binVar{iBase,1} = sdpvar(SimParams.nTxAntenna,1,'full');
        aPowVar{iBase,1} = sdpvar(SimParams.nTxAntenna,1,'full');
        for iBand = 1:nBands            
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
                for iBand = 1:nBands
                    Hsdp = cH{iBase,iBand}(:,:,cUser);
                    
                    tempVec = [sqrt(SimParams.N)];
                    riX = real(Hsdp * X{iBase,iBand}(:,iGroup));
                    iiX = imag(Hsdp * X{iBase,iBand}(:,iGroup));
                    fixedPoint = (rX(cUser,iBand)^2 + iX(cUser,iBand)^2) / bX(cUser,iBand);
                    tempExpression = fixedPoint + (2 / bX(cUser,iBand)) * (rX(cUser,iBand) * (riX - rX(cUser,iBand)) + iX(cUser,iBand) * (iiX - iX(cUser,iBand))) ...
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
                    
                    gConstraints = [gConstraints, Gamma(cUser,iBand) - tempExpression <= 0];
                    gConstraints = [gConstraints, (tempVec' * tempVec) - Beta(cUser,iBand) <= 0];
                end
                
                cGamma = Gamma(cUser,:) + 1;
                gConstraints = [gConstraints, gReqSINRPerUser(cUser,1) - geomean(cGamma(:)) <= 0];
            end
        end
    end
    
    tempSum = 0;
    for iBase = 1:nBases
        for iAntenna = 1:SimParams.nTxAntenna
            tVec = [];
            tempRHS = binVar_P{iBase,1}(iAntenna,1)^qExponent + qExponent * (binVar{iBase,1}(iAntenna,1) - binVar_P{iBase,1}(iAntenna,1)) * binVar_P{iBase,1}(iAntenna,1)^(qExponent - 1);
            for iBand = 1:nBands                
                tVec = [tVec , X{iBase,iBand}(iAntenna,:)];
            end
            tempVector = [2 * tVec, (aPowVar{iBase,1}(iAntenna,1) - tempRHS)];
            gConstraints = [gConstraints, cone(tempVector,(aPowVar{iBase,1}(iAntenna,1) + tempRHS))];
        end
        gConstraints = [gConstraints, 0 <= binVar{iBase,1} <= 1];
        gConstraints = [gConstraints, sum(binVar{iBase,1}) == SimParams.nTxAntennaEnabled];
%         gConstraints = [gConstraints, tempSum <= SimParams.nTxAntennaEnabled];
    end
    
    objective = 0;
    for iBase = 1:nBases
            objective = objective + aPowVar{iBase,1}' * aPowVar{iBase,1};
    end
    
    options = sdpsettings('verbose',0,'solver','Mosek');
    solverOut = optimize(gConstraints,objective,options);
    SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray) = solverOut.solvertime + SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray);
    
    if ((solverOut.problem == 0) || (solverOut.problem == 3) || (solverOut.problem == 4))
        
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
            binVar_P{iBase,1} = double(binVar{iBase,1});          
            nEnabledAntenna = sum(double(binVar{iBase,1}));
        end
        bX = full(double(Beta));
    else
        display(solverOut);
        break;
    end
    
    objective = double(objective);
    if (abs(objective - minPower) / abs(minPower)) < epsilonT
        iterateSCA = 0;
    else
        minPower = objective;
    end
    
    if iIterateSCA < iterateSCAMax
        iIterateSCA = iIterateSCA + 1;
    else
        iterateSCA = 0;
    end
    
    display(double(binVar{iBase,1})');
    fprintf('Using [%d] Active Transmit Elements, Total power required is - %f \n',nEnabledAntenna,objective);    
    
end

SimParams.Debug.tempResource{2,1}{1,1} = rX;
SimParams.Debug.tempResource{3,1}{1,1} = iX;
SimParams.Debug.tempResource{4,1}{1,1} = bX;

SimParams.Debug.MultiCastSDPExchange = cell(nBases,nBands);
for iBase = 1:nBases
    for iBand = 1:nBands
        SimParams.Debug.MultiCastSDPExchange{iBase,iBand} = logical(int8(double(binVar{iBase,1})));
    end
    display(double(binVar{iBase,1})');
end

end