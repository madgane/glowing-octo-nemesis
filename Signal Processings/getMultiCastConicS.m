
function [SimParams,SimStructs] = getMultiCastConicS(SimParams,SimStructs)

initMultiCastVariables;
rX = SimParams.Debug.tempResource{2,1}{1,1};
iX = SimParams.Debug.tempResource{3,1}{1,1};
bX = SimParams.Debug.tempResource{4,1}{1,1};

iterateSCA = 1;
iIterateSCA = 0;
sparseIterate = 1;

minPower = 1e20;
prevObjective = 1e20;

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

uVar = cell(nBases,1);
for iBase = 1:nBases
    uVar{iBase,1} = ones(nTxAntenna,1);
end

while sparseIterate
    
    while iterateSCA
    
        gConstraints = [];
        X = cell(nBases,nBands);
        aPowVar = cell(nBases,1);
        binVar = sdpvar(SimParams.nTxAntenna,nBases,'full');
        
        for iBase = 1:nBases
            for iBand = 1:nBands
                aPowVar{iBase,iBand} = sdpvar(SimParams.nTxAntenna,nBands,'full');
                X{iBase,iBand} = sdpvar(nTxAntenna,nGroupsPerCell(iBase,1),'full','complex');
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
                        
                        tempVec = [sqrt(SimParams.N)];
                        Hsdp = cH{iBase,iBand}(:,enabledAntenna{iBase,iBand},cUser);
                        
                        riX = real(Hsdp * X{iBase,iBand}(:,iGroup));
                        iiX = imag(Hsdp * X{iBase,iBand}(:,iGroup));
                        fixedPoint = (rX(cUser,iBand)^2 + iX(cUser,iBand)^2) / bX(cUser,iBand);
                        tempExpression = fixedPoint + (2 / bX(cUser,iBand)) * (rX(cUser,iBand) * (riX - rX(cUser,iBand)) + iX(cUser,iBand) * (iiX - iX(cUser,iBand))) ...
                            - (fixedPoint / bX(cUser,iBand)) * (Beta(cUser,iBand) - bX(cUser,iBand));
                        
                        for jBase = 1:nBases
                            for jGroup = 1:nGroupsPerCell(jBase,1)
                                Hsdp = cH{jBase,iBand}(:,enabledAntenna{jBase,iBand},cUser);
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
                    gConstraints = [gConstraints , gReqSINRPerUser(cUser,1) - geomean(cGamma(:)) <= 0];
                    
                end
            end
            
            for iBase = 1:nBases
                for iBand = 1:nBands
                    for iAntenna = 1:SimParams.nTxAntenna
                        tempVector = [X{iBase,iBand}(iAntenna,:)];
                        gConstraints = [gConstraints, cone(tempVector,aPowVar{iBase,1}(iAntenna,iBand))];
                    end
                end
                for iAntenna = 1:SimParams.nTxAntenna
                    gConstraints = [gConstraints, cone(aPowVar{iBase,1}(iAntenna,:),binVar(iAntenna,iBase))];
                end
            end
            
            
        end
        
        objective = 0;
        for iBase = 1:nBases
            objective = objective + sum(uVar{iBase,1} .* binVar(:,iBase));
        end
        
        options = sdpsettings('verbose',0,'solver','fmincon');
        solverOut = optimize(gConstraints,objective,options);
        SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray) = solverOut.solvertime + SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray);
        
        if ((solverOut.problem == 0) || (solverOut.problem == 3) || (solverOut.problem == 4))
            for iBand = 1:nBands
                for iBase = 1:nBases
                    for iGroup = 1:nGroupsPerCell(iBase,1)
                        groupUsers = SimStructs.baseStruct{iBase,1}.mcGroup{iGroup,1};
                        for iUser = 1:length(groupUsers)
                            cUser = groupUsers(iUser,1);
                            Hsdp = cH{iBase,iBand}(:,enabledAntenna{iBase,iBand},cUser);
                            rX(cUser,iBand) = real(Hsdp * double(X{iBase,iBand}(:,iGroup)));
                            iX(cUser,iBand) = imag(Hsdp * double(X{iBase,iBand}(:,iGroup)));
                        end
                    end
                end
            end
            bX = full(double(Beta));
        else
            display(solverOut);
            break;
        end
        
        objective = double(objective);
        iterateSCA = 0;
        
    end
    
    iterateSCA = 1;    
    for iBase = 1:nBases
        uVar{iBase,1} = 1./(double(binVar(:,iBase)) + epsilonT);
        display(uVar{iBase,1});
        display(double(binVar(:,iBase)));
    end
    
    if abs(double(objective) - prevObjective) < epsilonT
        sparseIterate = 0;
    else
        prevObjective = double(objective);
    end
    
end

SimParams.Debug.maxActiveAntenna = SimParams.nTxAntenna;

for iBase = 1:nBases
    SimParams.Debug.maxActiveAntenna = min(sum(double(binVar(:,iBase)) >= epsilonT),SimParams.Debug.maxActiveAntenna);
end

end