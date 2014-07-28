function [SimParams,SimStructs] = performGroupReception(SimParams,SimStructs)

singleNode = 1;
nUsers = SimParams.nUsers;
uscoreIndex = find(SimParams.pathLossModel == '_');

if isempty(uscoreIndex)
    uscoreIndex = length(SimParams.pathLossModel) + 1;
end

userGroups = zeros(nUsers,2);
groupUsers = cell(SimParams.nBases,1);
for iBase = 1:SimParams.nBases
    groupUsers{iBase,1} = SimStructs.baseStruct{iBase,1}.mcGroup;
end

for iUser = 1:nUsers
    userGroups(iUser,1) = SimStructs.userStruct{iUser,1}.baseNode;
    for iGroup = 1:length(groupUsers{userGroups(iUser,1)})
        if sum(iUser == groupUsers{userGroups(iUser,1),1}{iGroup,1})
            userGroups(iUser,2) = iGroup;
        end
    end
end

for iBand = 1:SimParams.nBands
    
    if strcmp(SimParams.DebugMode,'true')
        activeUsersDebug = cell(SimParams.nBases,1);
        for dBase = 1:SimParams.nBases
            activeUsersDebug{dBase,1} = sort(SimStructs.baseStruct{dBase,1}.assignedUsers{iBand,1})';
        end
        celldisp(activeUsersDebug);
    end
    
    linkChannel = SimStructs.actualChannel;
    for iUser = 1:SimParams.nUsers
        
        userActive = 0;
        performCooperation = 0;
        cUser = SimStructs.userStruct{iUser,1};
        
        baseNode = cUser.baseNode;
        neighNode = cUser.neighNode;
        
        if length(baseNode) ~= singleNode
            performCooperation = 1;
        end
        
        for kBase = 1:length(baseNode)
            if ~isempty(find(iUser == SimStructs.baseStruct{baseNode(1,kBase),1}.assignedUsers{iBand,1}))
                userActive = 1;
                break;
            end
        end
        
        Wmmse = cUser.W{iBand,1};
        if isempty(Wmmse)
            if strcmp(SimParams.DebugMode,'true')
                display('No RX Beamformer !');
            end
            Wmmse = eye(SimParams.nRxAntenna);
        end
        
        if userActive
            if strcmp(SimParams.pathLossModel(1:uscoreIndex(1,1) - 1),'3GPP')
                RoI = 10^(SimStructs.userStruct{iUser,1}.phyParams.restOfIF / 10);
                I = (SimParams.N + RoI) * eye(SimParams.nRxAntenna) * (Wmmse' * Wmmse);
            else
                I = SimParams.N * eye(SimParams.nRxAntenna) * (Wmmse' * Wmmse);
            end
        else
            I = 1;
        end
        
        S = 0;
        Nacc = 0;
        if userActive
            for iBase = 1:length(baseNode)
                
                baseIndex = baseNode(1,iBase);
                cBase = SimStructs.baseStruct{baseIndex,1};
                
                gP = cBase.PG{iBand,1};
                P = gP(:,userGroups(iUser,2));
                H = linkChannel{baseIndex,iBand}(:,:,iUser);
                S = Wmmse' * H * P + S;
                
                % Intra Group Calculation
                
                pIndices = userGroups(iUser,2) ~= (1:length(groupUsers{iBase,1}));
                P = gP(:,pIndices);
                if ~isempty(P)
                    N = Wmmse' * H * P;
                    if performCooperation
                        Nacc = Nacc + N;
                    else
                        I = I + N * N';
                    end
                end
                
            end
            
            if SimParams.nBases > 1
                if isempty(neighNode)
                    display('No Neighboring Groups');
                end
            end
            
            % Inter Group Interference Calculation
            
            for iBase = 1:length(neighNode)
                
                baseIndex = neighNode(1,iBase);
                cBase = SimStructs.baseStruct{baseIndex,1};
                
                P = cBase.PG{iBand,1};
                H = linkChannel{baseIndex,iBand}(:,:,iUser);
                
                if ~isempty(P)
                    N = Wmmse' * H * P;
                    I = I + N * N';
                end
                
            end
            
            if performCooperation
                I = I + Nacc * Nacc';
            end
            
            SimParams.Debug.receivedRSSI(:,:,iUser,iBand) = I;
            SimParams.Debug.activeStatus(iUser,iBand) = 1;
            
        else
            
            for iBase = 1:length(baseNode)
                
                baseIndex = baseNode(1,iBase);
                cBase = SimStructs.baseStruct{baseIndex,1};
                
                gP = cBase.PG{iBand,1};
                pIndices = iUser == cBase.assignedUsers{iBand,1};
                P = gP(:,pIndices);
                H = linkChannel{baseIndex,iBand}(:,:,iUser);
                S = Wmmse' * H * P + S;
                
                % Intra Stream Calculation
                
                pIndices = iUser ~= cBase.assignedUsers{iBand,1};
                P = gP(:,pIndices);
                
                if ~isempty(P)
                    N = Wmmse' * H * P;
                    if performCooperation
                        Nacc = Nacc + N;
                    else
                        I = I + N * N';
                    end
                end
                
            end
            
            if SimParams.nBases > 1
                if isempty(neighNode)
                    display('No Neighbors');
                end
            end
            
            % Inter Stream Calculation
            
            for iBase = 1:length(neighNode)
                
                baseIndex = neighNode(1,iBase);
                cBase = SimStructs.baseStruct{baseIndex,1};
                
                gP = cBase.PG{iBand,1};
                pIndices = iUser ~= cBase.assignedUsers{iBand,1};
                P = gP(:,pIndices);
                
                H = linkChannel{baseIndex,iBand}(:,:,iUser);
                
                if ~isempty(P)
                    N = Wmmse' * H * P;
                    I = I + N * N';
                end
                
            end
            
            if performCooperation
                I = I + Nacc * Nacc';
            end
            
            SimParams.Debug.receivedRSSI(:,:,iUser,iBand) = I;
            SimParams.Debug.activeStatus(iUser,iBand) = 0;
            
        end
        
        L = eye(size(I)) + pinv(I) * (S * S');
        xThrpt = real(log2(det(L)));
        
        if ~isnan(xThrpt)
            SimStructs.userStruct{iUser,1}.crThrpt = SimStructs.userStruct{iUser,1}.crThrpt + real(xThrpt);
            SimStructs.userStruct{iUser,1}.lastThrpt = xThrpt + SimStructs.userStruct{iUser,1}.lastThrpt;
            SimStructs.userStruct{iUser,1}.dropThrpt(SimParams.iDrop,1) = xThrpt;
            SimParams.Debug.resAllocation(SimParams.iDrop,iBand,iUser,SimParams.iSNR) = xThrpt;
        else
            SimStructs.userStruct{iUser,1}.dropThrpt(SimParams.iDrop,1) = 0;
            SimParams.Debug.resAllocation(SimParams.iDrop,iBand,iUser,SimParams.iSNR) = 0;
            
        end
        
        if userActive
            if ~isnan(xThrpt)
                if sign(xThrpt)
                    SimStructs.userStruct{iUser,1}.tAllocation = SimStructs.userStruct{iUser,1}.tAllocation + 1;
                end
            end
        end
    end
end

for iBand = 1:SimParams.nBands
    for iBase = 1:SimParams.nBases
        totPower = trace(SimStructs.baseStruct{iBase,1}.PG{iBand,1} * SimStructs.baseStruct{iBase,1}.PG{iBand,1}');
        SimParams.txPower(SimParams.iPkt,SimParams.iSNR,iBase,iBand) = totPower;
    end
end

