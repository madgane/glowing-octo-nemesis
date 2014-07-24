
function [SimParams,SimStructs] = getMultiCastPrecoders(SimParams,SimStructs)

epsilonT = 1e-5;
cH = SimStructs.linkChan;
nBases = SimParams.nBases;
nBands = SimParams.nBands;

% Debug Buffers initialization

SimParams.Debug.tempResource{2,SimParams.iDrop} = cell(SimParams.nUsers,1);
SimParams.Debug.tempResource{3,SimParams.iDrop} = cell(SimParams.nUsers,1);
SimParams.Debug.tempResource{4,SimParams.iDrop} = cell(SimParams.nUsers,SimParams.nBands);

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

underscore_location = strfind(SimParams.weightedSumRateMethod,'_');
if isempty(underscore_location)
    selectionMethod = SimParams.weightedSumRateMethod;
else
    selectionMethod = SimParams.weightedSumRateMethod(1:underscore_location-1);
end

nGroupsPerCell = zeros(SimParams.nBases,1);
for iBase = 1:nBases
    nGroupsPerCell(iBase,1) = length(SimParams.mcGroups{iBase,1});
end

switch selectionMethod
    
    case 'SDPMethod'
        
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
                                    tempSum = tempSum - trace(Hsdp * X{jBase,iBand}(:,:,jGroup));
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
    
        options = sdpsettings('verbose',1,'solver','sdpt3');
        solverOut = solvesdp(gConstraints,objective,options);
        
        if solverOut.problem
            display('Unable to Solve');
        else
            for iBand = 1:nBands
                SimStructs.baseStruct{iBase,1}.PG = cell(nBands,1);
                for iBase = 1:nBases                    
                    dX = double(X{iBase,iBand});
                    for iGroup = nGroupsPerCell(iBase,1)
                        [P D] = eig(dX(:,:,iGroup));
                        
                        
                        
                        SimStructs.baseStruct{iBase,1}.PG{iBand,1} = 
                    end
                end
            end
        end
        
        
        
        
    otherwise
        display('Unknown Algorithm !');    
        
end

for iUser = 1:nUsers
    SimParams.Debug.tempResource{2,SimParams.iDrop}{iUser,1} = SimParams.Debug.tempResource{2,SimParams.iDrop}{iUser,1};
    for iBand = 1:nBands
        SimParams.Debug.tempResource{4,SimParams.iDrop}{iUser,iBand} = SimParams.Debug.tempResource{4,SimParams.iDrop}{iUser,iBand};
        SimStructs.userStruct{iUser,1}.W{iBand,1} = W{iUser,iBand};
    end
end
