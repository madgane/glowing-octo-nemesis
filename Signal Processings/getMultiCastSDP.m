
function [SimParams,SimStructs] = getMultiCastSDP(SimParams,SimStructs,nIterations)

cH = SimStructs.linkChan;
nBases = SimParams.nBases;
nBands = SimParams.nBands;

% Debug Buffers initialization

nUsers = SimParams.nUsers;
QueuedPkts = zeros(nUsers,1);
reqSINRPerUser = zeros(nUsers,1);

for iBase = 1:nBases
    linkedUsers = SimStructs.baseStruct{iBase,1}.linkedUsers;
    for iUser = 1:length(linkedUsers)
        cUser = linkedUsers(iUser,1);
        QueuedPkts(cUser,1) = SimStructs.userStruct{cUser,1}.trafficStats.backLogPkt;
        reqSINRPerUser(cUser,1) = 2^(QueuedPkts(cUser,1)) - 1;
    end
end

nGroupsPerCell = zeros(SimParams.nBases,1);
for iBase = 1:nBases
    nGroupsPerCell(iBase,1) = length(SimParams.mcGroups{iBase,1});
end

gConstraints = [];
X = cell(nBases,nBands);
for iBand = 1:nBands
    for iBase = 1:nBases
        X{iBase,iBand} = sdpvar(SimParams.nTxAntenna,SimParams.nTxAntenna,nGroupsPerCell(iBase,1),'hermitian','complex');
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
                        Hsdp = cH{jBase,iBand}(:,:,cUser)' * cH{jBase,iBand}(:,:,cUser);
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

options = sdpsettings('verbose',1,'solver','SDPT3');
solverOut = solvesdp(gConstraints,objective,options);
SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray,SimParams.iGroupArray) = solverOut.solvertime + SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray,SimParams.iGroupArray);

if solverOut.problem
    display(solverOut);
    display('SDP Failed !');
    for iBase = 1:nBases
        for iBand = 1:nBands
            SimStructs.baseStruct{iBase,1}.P_SDP{iBand,1} = zeros(SimParams.nTxAntenna,SimParams.nTxAntenna,nGroupsPerCell(iBase,1))
            SimStructs.baseStruct{iBase,1}.PG{iBand,1} = zeros(SimParams.nTxAntenna,nGroupsPerCell(iBase,1));
        end
    end
else
    
    maxPower = 1e20;
    Y = cell(nBases,nBands);
    beamPower = cell(iBase,1);
    
    for iBase = 1:nBases
        for iBand = 1:nBands
            SimStructs.baseStruct{iBase,1}.P_SDP{iBand,1} = double(X{iBase,iBand});
        end
    end
    
    for iIterate = 1:nIterations
        
        reIterate = 1;
        while reIterate
            for iBase = 1:nBases
                beamPower{iBase,1} = zeros(nGroupsPerCell(iBase,1),nBands);
                for iBand = 1:nBands
                    dX = double(X{iBase,iBand});
                    Y{iBase,iBand} = zeros(SimParams.nTxAntenna,nGroupsPerCell(iBase,1));
                    for iGroup = 1:nGroupsPerCell(iBase,1)
                        [P, D] = eig(dX(:,:,iGroup));
                        randSelection = complex(randn(SimParams.nTxAntenna,1),randn(SimParams.nTxAntenna,1));
                        Y{iBase,iBand}(:,iGroup) = P * sqrt(D) * randSelection / norm(randSelection);
                        beamPower{iBase,1}(iGroup,iBand) = norm(Y{iBase,iBand}(:,iGroup),2)^2;
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
                                    Hsdp = cH{jBase,iBand}(:,:,cUser);
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
        
        totalPower = 0;
        for iBase = 1:nBases
            beamPower{iBase,1} = beamPower{iBase,1} .* double(grpPwr{iBase,1});
            totalPower = sum(beamPower{iBase,1}) + totalPower;
        end
        
        if totalPower < maxPower
            for iBase = 1:nBases
                grpP = double(grpPwr{iBase,1});
                for iBand = 1:nBands
                    SimStructs.baseStruct{iBase,1}.PG{iBand,1} = Y{iBase,iBand} * diag(sqrt(grpP(:,iBand)));
                end
            end
            maxPower = totalPower;
        end
        
    end
    
end

end
