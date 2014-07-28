
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
        
        options = sdpsettings('verbose',1,'solver','sdpt3');
        solverOut = solvesdp(gConstraints,objective,options);
        
        if solverOut.problem
            display(solverOut);
            display('SDP Failed !');
            for iBase = 1:nBases
                SimStructs.baseStruct{iBase,1}.PG = cell(nBands,1);
                for iBand = 1:nBands
                    SimStructs.baseStruct{iBase,1}.PG{iBand,1} = zeros(SimParams.nTxAntenna,nGroupsPerCell(iBase,1));
                end
            end
        else
            reIterate = 1;
            while reIterate
                Y = cell(nBases,nBands);
                for iBand = 1:nBands
                    for iBase = 1:nBases
                        dX = double(X{iBase,iBand});
                        Y{iBase,iBand} = zeros(SimParams.nTxAntenna,nGroupsPerCell(iBase,1));
                        for iGroup = 1:nGroupsPerCell(iBase,1)
                            [P, D] = eig(dX(:,:,iGroup));
                            randSelection = ones(SimParams.nTxAntenna,1);%complex(randn(SimParams.nTxAntenna,1),randn(SimParams.nTxAntenna,1));
                            Y{iBase,iBand}(:,iGroup) = P * sqrt(D) * randSelection / norm(randSelection);
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
                
                options = sdpsettings('verbose',1,'solver','linprog');
                solverOut = solvesdp(gConstraints,objective,options);
                
                if ~solverOut.problem
                    reIterate = 0;
                    display('Linear Programming Failed !');
                end
            end
            
            for iBase = 1:nBases
                grpP = double(grpPwr{iBase,1});
                SimStructs.baseStruct{iBase,1}.PG = cell(nBands,1);
                for iBand = 1:nBands
                    SimStructs.baseStruct{iBase,1}.PG{iBand,1} = Y{iBase,iBand} * diag(sqrt(grpP(:,iBand)));
                end
            end
            
        end
        
    case 'ConicMethod'
        
        % Feasible point generation !
        
        gConstraints = [];
        X = cell(nBases,nBands);
        for iBand = 1:nBands
            for iBase = 1:nBases
                X{iBase,iBand} = sdpvar(SimParams.nTxAntenna,nGroupsPerCell(iBase,1),'full','complex');
            end
        end
        
        bX = sdpvar(nUsers,nBands,'full');
        
        for iBand = 1:nBands
            for iBase = 1:nBases
                for iGroup = 1:nGroupsPerCell(iBase,1)
                    groupUsers = SimStructs.baseStruct{iBase,1}.mcGroup{iGroup,1};
                    for iUser = 1:length(groupUsers)
                        cUser = groupUsers(iUser,1);
                        Hsdp = cH{iBase,iBand}(:,:,cUser);
                        
                        rX = real(Hsdp * X{iBase,iBand}(:,iGroup));
                        iX = imag(Hsdp * X{iBase,iBand}(:,iGroup));
                        tempExpression = rX^2 + iX^2;
                        gConstraints = [gConstraints , tempExpression - reqSINRPerUser(cUser,1) * bX(cUser,iBand) >= 0];
                        
                        tempSum = SimParams.N;
                        for jBase = 1:nBases
                            for jGroup = 1:nGroupsPerCell(jBase,1)
                                Hsdp = cH{jBase,iBand}(:,:,cUser);
                                if ~and((iBase == jBase),(iGroup == jGroup))
                                    rX = real(Hsdp * X{jBase,iBand}(:,jGroup));
                                    iX = imag(Hsdp * X{jBase,iBand}(:,jGroup));
                                    tempSum = tempSum + rX^2 + iX^2;
                                end
                            end
                        end
                        gConstraints = [gConstraints, cone(tempSum,bX(cUser,iBand)) <= 0];
                    end
                end
            end
        end
        
        objective = 0;
        for iBand = 1:nBands
            for iBase = 1:nBases
                for iGroup = 1:nGroupsPerCell(iBase,1)
                    objective = objective + real(trace(X{iBase,iBand}(:,iGroup) * X{iBase,iBand}(:,iGroup)'));
                end
            end
        end
        
        options = sdpsettings('verbose',1,'solver','fmincon');
        solverOut = solvesdp(gConstraints,objective,options);
        
        % Optimal Precoder Search !
        
        if solverOut.problem == 0
            rX = zeros(nUsers,nBands);
            iX = zeros(nUsers,nBands);
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
                    
                end
            end
        end
        
        iterateSCA = 1;
        while iterateSCA            
            
            gConstraints = [];
            X = cell(nBases,nBands);
            for iBand = 1:nBands
                for iBase = 1:nBases
                    X{iBase,iBand} = sdpvar(SimParams.nTxAntenna,nGroupsPerCell(iBase,1),'full','complex');
                end
            end
            
            bX = sdpvar(nUsers,nBands);
            
            for iBand = 1:nBands
                for iBase = 1:nBases
                    for iGroup = 1:nGroupsPerCell(iBase,1)
                        groupUsers = SimStructs.baseStruct{iBase,1}.mcGroup{iGroup,1};
                        for iUser = 1:length(groupUsers)
                            cUser = groupUsers(iUser,1);
                            Hsdp = cH{iBase,iBand}(:,:,cUser);
                            
                            riX = real(Hsdp * X{iBase,iBand}(:,iGroup));
                            iiX = imag(Hsdp * X{iBase,iBand}(:,iGroup));
                            tempExpression = rX(cUser,iBand) * (riX - rX(cUser,iBand)) + iX(cUser,iBand) * (iiX - iX(cUser,iBand));
                            gConstraints = [gConstraints , tempExpression - reqSINRPerUser(cUser,1) * bX(cUser,iBand)^2 >= 0];
                            
                            tempSum = sqrt(SimParams.N);
                            for jBase = 1:nBases
                                for jGroup = 1:nGroupsPerCell(jBase,1)
                                    Hsdp = cH{jBase,iBand}(:,:,cUser);
                                    if ~and((iBase == jBase),(iGroup == jGroup))
                                        riX = real(Hsdp * X{jBase,iBand}(:,jGroup));
                                        iiX = imag(Hsdp * X{jBase,iBand}(:,jGroup));
                                        tempSum = [tempSum , riX , iiX];
                                    end
                                end
                            end
                            gConstraints = [gConstraints, cone(tempSum,bX(cUser,iBand))];
                        end
                    end
                end
            end
            
            objective = 0;
            for iBand = 1:nBands
                for iBase = 1:nBases
                    for iGroup = 1:nGroupsPerCell(iBase,1)
                        objective = objective + real(trace(X{iBase,iBand}(:,iGroup) * X{iBase,iBand}(:,iGroup)'));
                    end
                end
            end
            
            options = sdpsettings('verbose',1,'solver','quadprog');
            solverOut = solvesdp(gConstraints,objective,options);
            
            if solverOut.problem == 0
                iterateSCA = 0;
            end
            
        end
        
        
    otherwise
        display('Unknown Algorithm !');
        
end

end
