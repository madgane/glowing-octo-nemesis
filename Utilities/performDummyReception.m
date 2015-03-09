function [SimParams,SimStructs] = performDummyReception(varargin)

SimParams = varargin{1,1};
SimStructs = varargin{1,2};

switch nargin
    case 2
        xBandIndices = 1:SimParams.nBands;
    case 3
        xBandIndices = varargin{1,3};
end

singleNode = 1;
linkChannel = SimStructs.linkChan;

uscoreIndex = find(SimParams.pathLossModel == '_');
if isempty(uscoreIndex)
    uscoreIndex = length(SimParams.pathLossModel) + 1;
end

for xBand = 1:length(xBandIndices)
    
    iBand = xBandIndices(1,xBand);
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
            display('No RX Beamformer !');
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
                
                gP = cBase.P{iBand,1};
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
                
                gP = cBase.P{iBand,1};
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
            
        else
            
            for iBase = 1:length(baseNode)
                
                baseIndex = baseNode(1,iBase);
                cBase = SimStructs.baseStruct{baseIndex,1};
                
                gP = cBase.P{iBand,1};
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
                
                gP = cBase.P{iBand,1};
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
            
        end
        
        L = eye(size(I)) + (pinv(I) * (S * S'));
        xThrpt = real(log2(det(L))) *  SimParams.BITFactor;

        if ~isnan(xThrpt)
            SimParams.Debug.privateExchanges.resAllocation(iBand,iUser) = real(xThrpt);
        else
            SimParams.Debug.privateExchanges.resAllocation(iBand,iUser) = 0;
        end
        
    end
    
end
