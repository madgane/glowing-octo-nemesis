
function [SimParams,SimStructs] = getMultiCastConicA(SimParams,SimStructs,ObjType)

initMultiCastVariables;
rX = SimParams.Debug.tempResource{2,1}{1,1};
iX = SimParams.Debug.tempResource{3,1}{1,1};

maxObj = 1e2;
iterateSCA = 1;
iIterateSCA = 0;
minPower = 1e20;

reqSINRPerUser = 2.^(QueuedPkts / nBands) - 1;
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

while iterateSCA
    
    gConstraints = [];
    X = cell(nBases,nBands);
    for iBand = 1:nBands
        for iBase = 1:nBases
            X{iBase,iBand} = sdpvar(nTxAntenna,nGroupsPerCell(iBase,1),'full','complex');
        end
    end
    
    switch ObjType
        case 'Dual'
            feasVariable = sdpvar(1,1);
    end
    
    for iBand = 1:nBands
        for iBase = 1:nBases
            for iGroup = 1:nGroupsPerCell(iBase,1)
                groupUsers = SimStructs.baseStruct{iBase,1}.mcGroup{iGroup,1};
                for iUser = 1:length(groupUsers)
                    cUser = groupUsers(iUser,1);
                    Hsdp = cH{iBase,iBand}(:,enabledAntenna{iBase,iBand},cUser);
                    
                    tempVec = [sqrt(SimParams.N)];
                    riX = real(Hsdp * X{iBase,iBand}(:,iGroup));
                    iiX = imag(Hsdp * X{iBase,iBand}(:,iGroup));
                    tempExpression = rX(cUser,iBand)^2 + iX(cUser,iBand)^2 + 2 * (rX(cUser,iBand) * (riX - rX(cUser,iBand)) + iX(cUser,iBand) * (iiX - iX(cUser,iBand)));
                    
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
                    
                    switch ObjType
                        case 'MP'
                            gConstraints = [gConstraints , reqSINRPerUser(cUser,1) * (tempVec' * tempVec) - tempExpression <= 0];
                        case 'Dual'
                            gConstraints = [gConstraints , reqSINRPerUser(cUser,1) * (tempVec' * tempVec) - tempExpression <= feasVariable];
                        case 'FC'
                            gConstraints = [gConstraints , reqSINRPerUser(cUser,1) * (tempVec' * tempVec) - tempExpression <= 0];
                    end
                    
                end
            end
        end
    end
    
    switch ObjType
        case 'FC'
            objective = [];
        case 'MP'            
            objective = 0;
            for iBand = 1:nBands
                for iBase = 1:nBases
                    objective = objective + (X{iBase,iBand}(:)' * X{iBase,iBand}(:));
                end
            end
        case 'Dual'
            objective = 0;
            for iBand = 1:nBands
                for iBase = 1:nBases
                    objective = objective + (X{iBase,iBand}(:)' * X{iBase,iBand}(:));
                end
            end
            objective = max(feasVariable) * maxObj + objective;
    end
    
    options = sdpsettings('verbose',0,'solver','Fmincon');
    solverOut = optimize(gConstraints,objective,options);
    SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray) = solverOut.solvertime + SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray);
    
    if solverOut.problem == 0
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
    else
        display(solverOut);
        if sum(strcmpi({'FC','Dual'},ObjType))
            rX = zeros(size(rX));iX = zeros(size(iX));
        else
            display(solverOut);
            break;
        end
    end
    
    switch ObjType
        case 'Dual'
            if (double(feasVariable) < 0)
                break;
            end
            
            fprintf('Setting maxObj - %4.2f, Feasible Variable - %f \n',maxObj,double(feasVariable));
            
        case 'MP'
            objective = double(objective);
            if (abs(objective - minPower) / abs(minPower)) < epsilonT
                iterateSCA = 0;
            else
                minPower = objective;
            end
            display(double(objective));
        case 'FC'
            yalmiperror(solverOut.problem);
            break;
    end
    
    if iIterateSCA < iterateSCAMax
        iIterateSCA = iIterateSCA + 1;
    else
        iterateSCA = 0;
    end
end

for iBase = 1:nBases
    for iBand = 1:nBands
        tempPrecoder = double(X{iBase,iBand});
        SimStructs.baseStruct{iBase,1}.PG{iBand,1} = zeros(SimParams.nTxAntenna,nGroupsPerCell(iBase,1));
        SimStructs.baseStruct{iBase,1}.PG{iBand,1}(enabledAntenna{iBase,iBand},:) = tempPrecoder;
    end
end

end