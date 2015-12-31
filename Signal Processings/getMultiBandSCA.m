
function [SimParams,SimStructs] = getMultiBandSCA(SimParams,SimStructs,ObjType)

initMultiCastVariables;
rX = SimParams.Debug.tempResource{2,1}{1,1};
iX = SimParams.Debug.tempResource{3,1}{1,1};
bX = SimParams.Debug.tempResource{4,1}{1,1};

fcCount = 0;
iterateSCA = 1;
iIterateSCA = 0;
minPower = 1e20;

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
    
    Beta = sdpvar(nUsers,nBands,'full');
    Gamma = sdpvar(nUsers,nBands,'full');

    switch ObjType
        case 'Dual'
            feasVariable = sdpvar(1,1);
        case 'MaxMin'
            fairnessVariable = sdpvar(1,1);
    end
    
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

                switch ObjType
                    case 'MP'
                        gConstraints = [gConstraints , gReqSINRPerUser(cUser,1) - geomean(cGamma(:)) <= 0];
                    case 'Dual'
                        gConstraints = [gConstraints , gReqSINRPerUser(cUser,1) - geomean(cGamma(:)) <= feasVariable];
                    case 'FC'
                        gConstraints = [gConstraints , gReqSINRPerUser(cUser,1) - geomean(cGamma(:)) <= 0];
                    case 'MaxMin'
                        gConstraints = [gConstraints , (fairnessVariable + 1) - geomean(cGamma(:)) <= 0];
                end
            end
        end
    end
    
    if (logical(sum(strcmpi(ObjType,{'FC','Dual'}))) && (SimParams.nTxAntennaEnabled <= singleBandCheck))
        for iGroup = 1:nGroupsPerCell(iBase,1)
            for iBand = 1:nBands
                if iBand ~= mod(iGroup - 1,nBands) + 1
                    gConstraints = [gConstraints, X{iBase,iBand}(:,iGroup)' * X{iBase,iBand}(:,iGroup) <= 0];
                end
            end
        end
    end
    
    if logical(sum(strcmpi(ObjType,'MaxMin')))
        for iBase = 1:nBases
            xVector = [];
            for iBand = 1:nBands
                xVector = [xVector ; X{iBase,iBand}(:)];
            end
            gConstraints = [gConstraints, xVector' * xVector <= sum(SimStructs.baseStruct{iBase,1}.sPower(1,:))];
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
            objective = feasVariable + objective * objWeight;
        case 'MaxMin'
            objective = -fairnessVariable;
    end
    
    options = sdpsettings('verbose',0,'solver','Mosek');
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
        rX = randn(size(rX));iX = randn(size(iX));
    end
    
    switch ObjType
        case 'Dual'
            if (double(feasVariable) < 0)
                fcCount = fcCount + 1;
                if fcCount > 0
                    break;
                end
            end  
            
            fprintf('Feasible Variable - %f \n',double(feasVariable));
            
        case {'MP','MaxMin'}
            objective = value(objective);
            if (abs(objective - minPower) / abs(minPower)) < epsilonT
                iterateSCA = 0;
            else
                minPower = objective;
            end
            
            if strcmpi(ObjType,'MaxMin')
                SimParams.Debug.minUserRate = value(fairnessVariable);
            end
            
            fprintf('Total Power Required - %4.4f \n',value(objective));
            
        case 'FC'
            if solverOut.problem == 0
                break;
            else
                yalmiperror(solverOut.problem);
            end
            
    end
    
    if iIterateSCA < iterateSCAMax
        iIterateSCA = iIterateSCA + 1;
    else
        iterateSCA = 0;
    end
   
end

if strcmpi(ObjType,'Dual')
    if value(feasVariable) > 0
        SimParams.Debug.SCA_initFailureFlag = 1;
    end
end

for iBase = 1:nBases
    for iBand = 1:nBands
        tempPrecoder = double(X{iBase,iBand});
        SimStructs.baseStruct{iBase,1}.PG{iBand,1} = zeros(SimParams.nTxAntenna,nGroupsPerCell(iBase,1));
        SimStructs.baseStruct{iBase,1}.PG{iBand,1}(enabledAntenna{iBase,iBand},:) = tempPrecoder;
    end
end

SimParams.Debug.tempResource{4,1}{1,1} = bX;
SimParams.Debug.tempResource{5,1}{1,1} = value(Gamma);

end