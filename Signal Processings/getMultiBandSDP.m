
function [SimParams,SimStructs] = getMultiBandSDP(SimParams,SimStructs)

initMultiCastVariables;

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

SimParams.Debug.groupRank = [];

bX = SimParams.Debug.tempResource{4,1}{1,1};
gX = SimParams.Debug.tempResource{5,1}{1,1};

while iterateSCA
    
    gConstraints = [];
    X = cell(nBases,nBands);
    for iBand = 1:nBands
        for iBase = 1:nBases
            X{iBase,iBand} = sdpvar(nTxAntenna,nTxAntenna,nGroupsPerCell(iBase,1),'hermitian','complex');
        end
    end
    
    b = sdpvar(nUsers,nBands,'full');
    Gamma = sdpvar(nUsers,nBands,'full');
        
    for iBand = 1:nBands
        for iBase = 1:nBases
            for iGroup = 1:nGroupsPerCell(iBase,1)
                groupUsers = SimStructs.baseStruct{iBase,1}.mcGroup{iGroup,1};
                for iUser = 1:length(groupUsers)
                    cUser = groupUsers(iUser,1);
                    tempB = SimParams.N;
                    for jBase = 1:nBases
                        for jGroup = 1:nGroupsPerCell(jBase,1)
                            Hsdp = cH{jBase,iBand}(:,enabledAntenna{jBase,iBand},cUser)' * cH{jBase,iBand}(:,enabledAntenna{jBase,iBand},cUser);
                            if and((iBase == jBase),(iGroup == jGroup))
                                tempX = trace(Hsdp * X{jBase,iBand}(:,:,jGroup));
                            else
                                tempB = tempB + trace(Hsdp * X{jBase,iBand}(:,:,jGroup));
                            end
                        end
                    end
                    gConstraints = [gConstraints, tempB <= b(cUser,iBand)];
                    taylorApprox = (gX(cUser,iBand) - bX(cUser,iBand))^2 + 2 * (gX(cUser,iBand) - bX(cUser,iBand)) * ((Gamma(cUser,iBand) - gX(cUser,iBand)) - (b(cUser,iBand) - bX(cUser,iBand)));
                    gConstraints = [gConstraints, 4 * tempX + taylorApprox >= (Gamma(cUser,iBand) + b(cUser,iBand))^2];
                end
                gConstraints = [gConstraints, X{iBase,iBand}(:,:,iGroup) >= 0];
            end
        end
    end
    
    cGamma = Gamma + 1;
    for iUser = 1:nUsers
        gConstraints = [gConstraints , gReqSINRPerUser(iUser,1) - geomean(cGamma(iUser,:)) <= 0];
%         gConstraints = [gConstraints , reqSINRPerUser(iUser,1) - Gamma(iUser,2) <= 0];
    end
    
    gConstraints = [gConstraints, Gamma >= 0];
    
    objective = 0;
    for iBand = 1:nBands
        for iBase = 1:nBases
            for iGroup = 1:nGroupsPerCell(iBase,1)
                objective = objective + real(trace(X{iBase,iBand}(:,:,iGroup)));
            end
        end
    end
    
    options = sdpsettings('verbose',0,'solver','Sedumi');
    solverOut = optimize(gConstraints,objective,options);
    SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray) = solverOut.solvertime + SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray);
    
    if ~((solverOut.problem == 0) || (solverOut.problem == 3) || (solverOut.problem == 4))
        display(yalmiperror(solverOut.problem));
        display('SDP Failed !');
        for iBase = 1:nBases
            for iBand = 1:nBands
                SimStructs.baseStruct{iBase,1}.P_SDP{iBand,1} = double(X{iBase,iBand});
                SimStructs.baseStruct{iBase,1}.PG{iBand,1} = zeros(nTxAntenna,nGroupsPerCell(iBase,1));
            end
        end
        gX = rand(nUsers,nBands);
        bX = ones(nUsers,nBands) + rand(nUsers,nBands);
    end
    
    if iIterateSCA < iterateSCAMax
        iIterateSCA = iIterateSCA + 1;
    else
        iterateSCA = 0;
    end
    
    bX = double(b);
    gX = double(Gamma);    

    txPower = 0;
    for iBand = 1:nBands
        for iGroup = 1:nGroupsPerCell(iBase,1)
            txPower = txPower + trace(double(X{iBase,iBand}(:,:,iGroup)));
        end
    end    
    
    if (abs(txPower - minPower) / abs(txPower)) < epsilonT
        iterateSCA = 0;
    else
        minPower = txPower;
    end

    fprintf('Transmit Power Required for Tx - %2.2f \n',txPower);
    
end

SimParams.Debug.SDP_ReqSINR_Band = value(Gamma);
SimParams.Debug.MultiCastSDPExchange = linspace(1,SimParams.nTxAntenna,SimParams.nTxAntenna);

end
