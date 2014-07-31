
function [SimParams,SimStructs] = getMultiCastPrecoders(SimParams,SimStructs)

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

underscore_location = strfind(SimParams.DesignType,'_');
if isempty(underscore_location)
    selectionMethod = SimParams.DesignType;
else
    selectionMethod = SimParams.DesignType(1:underscore_location-1);
end

nGroupsPerCell = zeros(SimParams.nBases,1);
for iBase = 1:nBases
    nGroupsPerCell(iBase,1) = length(SimParams.mcGroups{iBase,1});
end

switch selectionMethod
    
    case 'SDPMethod'
        
        [SimParams,SimStructs] = getMultiCastSDP(SimParams,SimStructs,50);
        
    case 'ConicMethod'
        
        iterateSCAMax = 50;
        txPower = 2 * max(reqSINRPerUser) * SimParams.N;
        
        rX = randn(nUsers,nBands) * txPower;
        iX = randn(nUsers,nBands) * txPower;
        
        searchType = 'Dual';
        
        % Feasible point generation !
        
        switch searchType
            
            case 'Dual'
                
                objective = -1;
                iterateSCA = 1;
                iIterateSCA = 0;
                while iterateSCA
                    
                    gConstraints = [];
                    X = cell(nBases,nBands);
                    for iBand = 1:nBands
                        for iBase = 1:nBases
                            X{iBase,iBand} = sdpvar(SimParams.nTxAntenna,nGroupsPerCell(iBase,1),'full','complex');
                        end
                    end
                    oX = sdpvar(nUsers,nBands,'full');
                    for iBand = 1:nBands
                        for iBase = 1:nBases
                            for iGroup = 1:nGroupsPerCell(iBase,1)
                                groupUsers = SimStructs.baseStruct{iBase,1}.mcGroup{iGroup,1};
                                for iUser = 1:length(groupUsers)
                                    cUser = groupUsers(iUser,1);
                                    Hsdp = cH{iBase,iBand}(:,:,cUser);
                                    
                                    riX = real(Hsdp * X{iBase,iBand}(:,iGroup));
                                    iiX = imag(Hsdp * X{iBase,iBand}(:,iGroup));
                                    tempExpression = rX(cUser,iBand)^2 + iX(cUser,iBand)^2 + 2 * (rX(cUser,iBand) * (riX - rX(cUser,iBand)) + iX(cUser,iBand) * (iiX - iX(cUser,iBand)));
                                    tempSum = SimParams.N;
                                    for jBase = 1:nBases
                                        for jGroup = 1:nGroupsPerCell(jBase,1)
                                            Hsdp = cH{jBase,iBand}(:,:,cUser);
                                            if ~and((iBase == jBase),(iGroup == jGroup))
                                                riX = real(Hsdp * X{jBase,iBand}(:,jGroup));
                                                iiX = imag(Hsdp * X{jBase,iBand}(:,jGroup));
                                                tempSum = tempSum + riX^2 + iiX^2;
                                            end
                                        end
                                    end
                                    gConstraints = [gConstraints , reqSINRPerUser(cUser,1) * tempSum - tempExpression <= oX(cUser,iBand)];
                                end
                            end
                        end
                    end
                    
                    gConstraints = [gConstraints, oX <= 0];
                    options = sdpsettings('verbose',1,'solver','Gurobi');
                    solverOut = solvesdp(gConstraints,[],options);
                    SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray,SimParams.iGroupArray) = solverOut.solvertime + SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray,SimParams.iGroupArray);
                    
                    if solverOut.problem == 0
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
                    else
                        display(solverOut);
                        rX = rX / 2;iX = iX / 2;
                        continue;
                    end
                    
                    if iIterateSCA < iterateSCAMax
                        iIterateSCA = iIterateSCA + 1;
                    else
                        iterateSCA = 0;
                    end
                    
                    if double(objective) < 0
                        iterateSCA = 0;
                    end
                end
                
            case 'Nonlinear'
                
                gConstraints = [];
                X = cell(nBases,nBands);
                for iBand = 1:nBands
                    for iBase = 1:nBases
                        X{iBase,iBand} = sdpvar(SimParams.nTxAntenna,nGroupsPerCell(iBase,1),'full','complex');
                    end
                end
                
                for iBand = 1:nBands
                    for iBase = 1:nBases
                        for iGroup = 1:nGroupsPerCell(iBase,1)
                            groupUsers = SimStructs.baseStruct{iBase,1}.mcGroup{iGroup,1};
                            for iUser = 1:length(groupUsers)
                                cUser = groupUsers(iUser,1);
                                Hsdp = cH{iBase,iBand}(:,:,cUser);
                                rcX = real(Hsdp * X{iBase,iBand}(:,iGroup));
                                icX = imag(Hsdp * X{iBase,iBand}(:,iGroup));
                                tempExpression = rcX^2 + icX^2;
                                
                                tempSum = SimParams.N;
                                for jBase = 1:nBases
                                    for jGroup = 1:nGroupsPerCell(jBase,1)
                                        Hsdp = cH{jBase,iBand}(:,:,cUser);
                                        if ~and((iBase == jBase),(iGroup == jGroup))
                                            rcX = real(Hsdp * X{jBase,iBand}(:,jGroup));
                                            icX = imag(Hsdp * X{jBase,iBand}(:,jGroup));
                                            tempSum = tempSum + rcX^2 + icX^2;
                                        end
                                    end
                                end
                                gConstraints = [gConstraints , reqSINRPerUser(cUser,1) * tempSum - tempExpression <= 0];
                            end
                        end
                    end
                end
                
                options = sdpsettings('verbose',1,'solver','nomad');
                solverOut = solvesdp(gConstraints,[],options);
                SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray,SimParams.iGroupArray) = solverOut.solvertime + SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray,SimParams.iGroupArray);
                
                if solverOut.problem == 0
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
                else
                    display(solverOut);
                end
                
            case 'SDP'
                
                [SimParams,SimStructs] = getMultiCastSDP(SimParams,SimStructs,1);
                
                for iBase = 1:nBases
                    for iBand = 1:nBands
                        for iGroup = 1:nGroupsPerCell(iBase,1)
                            groupUsers = SimStructs.baseStruct{iBase,1}.mcGroup{iGroup,1};
                            for iUser = 1:length(groupUsers)
                                cUser = groupUsers(iUser,1);
                                Hsdp = cH{iBase,iBand}(:,:,cUser);
                                rX(cUser,iBand) = real(Hsdp * SimStructs.baseStruct{iBase,1}.PG{iBand,1}(:,iGroup));
                                iX(cUser,iBand) = imag(Hsdp * SimStructs.baseStruct{iBase,1}.PG{iBand,1}(:,iGroup));
                            end
                        end
                    end
                end
                
            otherwise
                
                rX = randn(nUsers,nBands) * txPower;
                iX = randn(nUsers,nBands) * txPower;
                display('Using Randomized data !');
                
        end
        
        iterateSCA = 1;
        iIterateSCA = 0;
        minPower = 1e20;
        while iterateSCA
            
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
                            
                            riX = real(Hsdp * X{iBase,iBand}(:,iGroup));
                            iiX = imag(Hsdp * X{iBase,iBand}(:,iGroup));
                            tempExpression = rX(cUser,iBand)^2 + iX(cUser,iBand)^2 + 2 * (rX(cUser,iBand) * (riX - rX(cUser,iBand)) + iX(cUser,iBand) * (iiX - iX(cUser,iBand)));
                            gConstraints = [gConstraints , reqSINRPerUser(cUser,1) * bX(cUser,iBand) * bX(cUser,iBand) - tempExpression <= 0];
                            tempVec = [sqrt(SimParams.N)];
                            for jBase = 1:nBases
                                for jGroup = 1:nGroupsPerCell(jBase,1)
                                    Hsdp = cH{jBase,iBand}(:,:,cUser);
                                    if ~and((iBase == jBase),(iGroup == jGroup))
                                        riX = real(Hsdp * X{jBase,iBand}(:,jGroup));
                                        iiX = imag(Hsdp * X{jBase,iBand}(:,jGroup));
                                        tempVec = [tempVec, riX, iiX];
                                    end
                                end
                            end
                            gConstraints = [gConstraints , tempVec * tempVec' - bX(cUser,iBand) <= 0];
                        end
                    end
                end
            end
            
            objective = 0;
            for iBand = 1:nBands
                for iBase = 1:nBases
                    objective = objective + (X{iBase,iBand}(:)' * X{iBase,iBand}(:));
                end
            end
            
            options = sdpsettings('verbose',1,'solver','Gurobi');
            solverOut = solvesdp(gConstraints,objective,options);
            SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray,SimParams.iGroupArray) = solverOut.solvertime + SimParams.solverTiming(SimParams.iPkt,SimParams.iAntennaArray,SimParams.iGroupArray);
            
            if solverOut.problem == 0
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
            else
                display(solverOut);
                txPower = 2 * max(reqSINRPerUser) * SimParams.N;
                rX = randn(nUsers,nBands) * txPower;
                iX = randn(nUsers,nBands) * txPower;
                continue;
            end
            
            if iIterateSCA < iterateSCAMax
                iIterateSCA = iIterateSCA + 1;
            else
                iterateSCA = 0;
            end
            
            if abs(double(objective) - minPower) < 1e-3
                iterateSCA = 0;
            else
                minPower = double(objective);
            end
        end
        
        for iBase = 1:nBases
            for iBand = 1:nBands
                SimStructs.baseStruct{iBase,1}.PG{iBand,1} = double(X{iBase,iBand});
            end
        end
        
    case 'KKTMethod'
        
        minPower = -10;
        iterateSCA = 1;
        iIterateSCA = 0;
        iterateSCAMax = 1000;
        dualGroup = zeros(nUsers,nBands);
                
        X = cell(nBases,nBands);
        Xt = cell(nBases,nBands);
        for iBand = 1:nBands
            for iBase = 1:nBases
                X{iBase,iBand} = complex(zeros(SimParams.nTxAntenna,nGroupsPerCell(iBase,1)),zeros(SimParams.nTxAntenna,nGroupsPerCell(iBase,1)));
                Xt{iBase,iBand} = complex(randn(SimParams.nTxAntenna,nGroupsPerCell(iBase,1)),randn(SimParams.nTxAntenna,nGroupsPerCell(iBase,1)));
            end
        end
        
        while iterateSCA            
            for dualIterate = 1:2               
                for iBand = 1:nBands
                    for iBase = 1:nBases
                        for iGroup = 1:nGroupsPerCell(iBase,1)
                            tempExpression = 0;
                            groupUsers = SimStructs.baseStruct{iBase,1}.mcGroup{iGroup,1};
                            for iUser = 1:length(groupUsers)
                                cUser = groupUsers(iUser,1);
                                Hsdp = cH{iBase,iBand}(:,:,cUser);
                                tempExpression = tempExpression + (Hsdp' * Hsdp) * Xt{iBase,iBand}(:,iGroup) * dualGroup(cUser,iBand);
                            end
                            tempSum = eye(SimParams.nTxAntenna);
                            for jBase = 1:nBases
                                for jGroup = 1:nGroupsPerCell(jBase,1)
                                    if ~and((iBase == jBase),(iGroup == jGroup))
                                        ifGroupUsers = SimStructs.baseStruct{jBase,1}.mcGroup{jGroup,1};
                                        for jUser = 1:length(ifGroupUsers)
                                            ifjUser = ifGroupUsers(jUser,1);
                                            Hsdp = cH{iBase,iBand}(:,:,ifjUser);
                                            tempSum = tempSum + reqSINRPerUser(ifjUser,1) * dualGroup(ifjUser,iBand) * (Hsdp' * Hsdp);
                                        end
                                    end
                                end
                            end
                            X{iBase,iBand}(:,iGroup) = pinv(tempSum) * tempExpression;
                        end
                    end
                end
                
                for iBand = 1:nBands
                    for iBase = 1:nBases
                        for iGroup = 1:nGroupsPerCell(iBase,1)
                            groupUsers = SimStructs.baseStruct{iBase,1}.mcGroup{iGroup,1};
                            for iUser = 1:length(groupUsers)
                                cUser = groupUsers(iUser,1);
                                Hsdp = cH{iBase,iBand}(:,:,cUser);
                                tempExpression = abs(Hsdp * X{iBase,iBand}(:,iGroup))^2;
                                tempSum = SimParams.N;
                                for jBase = 1:nBases
                                    for jGroup = 1:nGroupsPerCell(jBase,1)
                                        if ~and((iBase == jBase),(iGroup == jGroup))
                                            Hsdp = cH{jBase,iBand}(:,:,cUser);
                                            tempSum = tempSum + abs(Hsdp * X{jBase,iBand}(:,jGroup))^2;
                                        end
                                    end
                                end
                                dualGroup(cUser,iBand) = dualGroup(cUser,iBand) + 1e-3 * (tempSum * reqSINRPerUser(cUser,1) - tempExpression);
                                dualGroup(cUser,iBand) = max(0,dualGroup(cUser,iBand));
                            end
                        end
                    end
                end
                
                dX = cell2mat(X);
                disp(norm(dX(:),2));
                
            end
            
            Xt = X;
            if iIterateSCA < iterateSCAMax
                iIterateSCA = iIterateSCA + 1;
            else
                iterateSCA = 0;
            end
            
            dX = cell2mat(X);
            objPower = norm(dX(:),2);
            if abs(objPower - minPower) < 1e-6
                iterateSCA = 0;
            else
                minPower = objPower;
            end
        end
        
        for iBase = 1:nBases
            for iBand = 1:nBands
                SimStructs.baseStruct{iBase,1}.PG{iBand,1} = X{iBase,iBand};
            end
        end
        
        
    otherwise
        display('Unknown Precoding Method !');
        
end

end
