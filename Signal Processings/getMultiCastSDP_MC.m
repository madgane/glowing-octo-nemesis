
function [SimParams,SimStructs] = getMultiCastSDP_MC(SimParams,SimStructs)

initMultiCastVariables;

iterateSCA = 1;
iIterateSCA = 0;
minPower = 1e20;
nIterations = 50;

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

gX = zeros(nUsers,nBands);
bX = ones(nUsers,nBands);

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
    
    if solverOut.problem ~= 0
        display(yalmiperror(solverOut.problem));
        display('SDP Failed !');
        for iBase = 1:nBases
            for iBand = 1:nBands
                SimStructs.baseStruct{iBase,1}.P_SDP{iBand,1} = double(X{iBase,iBand});
                SimStructs.baseStruct{iBase,1}.PG{iBand,1} = zeros(nTxAntenna,nGroupsPerCell(iBase,1));
            end
        end
        return;
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

maxPower = 1e20;
Y = cell(nBases,nBands);

for iBase = 1:nBases
    for iBand = 1:nBands
        SimStructs.baseStruct{iBase,1}.P_SDP{iBand,1} = double(X{iBase,iBand});
    end
end

tY = cell(nBases,nBands,SimParams.nGroupArray);
SimParams.Debug.SDP_vars.groupRank = cell(nBases,1);
randSelection = cell(nBases,nBands,SimParams.nGroupArray);

for iBase = 1:nBases
    SimParams.Debug.SDP_vars.groupRank{iBase,1} = zeros(1,nGroupsPerCell(iBase,1));
    for iBand = 1:nBands
        dX = double(X{iBase,iBand});
        Y{iBase,iBand} = zeros(nTxAntenna,nGroupsPerCell(iBase,1));
        for iGroup = 1:nGroupsPerCell(iBase,1)
            [P, D] = eig(dX(:,:,iGroup));
            P = P(:,(diag(D) >= epsilonT));
            D = diag(D); D = D(D >= epsilonT);
            tY{iBase,iBand,iGroup} = P * sqrt(diag(D));
            SimParams.Debug.SDP_vars.groupRank{iBase,1}(1,iGroup) = length(D);
            SimParams.Debug.groupRank = [SimParams.Debug.groupRank, length(D)];
            randSelection{iBase,iBand,iGroup} = complex(rand(length(D),nIterations),rand(length(D),nIterations)) / sqrt(2);
        end
    end
end

xPrecoder = cell(nBases,nBands);
for iBase = 1:nBases
    for iBand = 1:nBands
        xPrecoder{iBase,iBand} = zeros(SimParams.nTxAntenna,nGroupsPerCell(iBase,1));
    end
end

for iIterate = 1:nIterations
    
    for iBase = 1:nBases
        for iBand = 1:nBands
            for iGroup = 1:nGroupsPerCell(iBase,1)
                Y{iBase,iBand}(:,iGroup) = tY{iBase,iBand,iGroup} * randSelection{iBase,iBand,iGroup}(:,iIterate);
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
                    tempSum = -SimParams.N * gX(cUser,iBand);
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
    
    options = sdpsettings('verbose',0,'solver','fmincon');
    solverOut = optimize(gConstraints,objective,options);
    SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray) = solverOut.solvertime + SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray);
    
    if solverOut.problem
        display(yalmiperror(solverOut.problem));
    end
    
    for iBase = 1:nBases
        grpP = double(grpPwr{iBase,1});
        for iBand = 1:nBands
            tempPrecoder = Y{iBase,iBand} * diag(sqrt(grpP(:,iBand)));
            xPrecoder{iBase,iBand}(enabledAntenna{iBase,iBand},:) = tempPrecoder;
        end
    end
    
    tX = cell2mat(xPrecoder);
    totalPower = norm(tX(:))^2;
    
    if totalPower < maxPower
        for iBase = 1:nBases
            for iBand = 1:nBands
                SimStructs.baseStruct{iBase,1}.PG{iBand,1} = xPrecoder{iBase,iBand};
            end
        end
        
        maxPower = totalPower;
    end
    
end


end
