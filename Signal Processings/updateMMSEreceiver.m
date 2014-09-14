
function [SimParams,SimStructs] = updateMMSEreceiver(SimParams,SimStructs,withPrecoder)

groupUsers = cell(SimParams.nBases,1);
groupStreams = cell(SimParams.nBases,1);

if ~withPrecoder
    [SimParams,SimStructs] = updateZFprecoders(SimParams,SimStructs,0);
end


for iBand = 1:SimParams.nBands
    
    for iBase = 1:SimParams.nBases
        groupUsers{iBase,1} = SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1}';
        groupStreams{iBase,1} = SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1}';
    end
    
    W = eye(SimParams.nRxAntenna);
    
    for iBase = 1:SimParams.nBases
        for iUser = 1:length(groupUsers{iBase,1})
            cUser = groupUsers{iBase,1}(1,iUser);
            cStream = groupStreams{iBase,1}(1,iUser);
            H = SimStructs.linkChan{iBase,iBand}(:,:,cUser);
            P = SimStructs.baseStruct{iBase,1}.P{iBand,1}(:,iUser);
            W = W + (H * P)' * (H * P);
            
            for jUser = 1:length(groupUsers{iBase,1})
                xUser = groupUsers{iBase,1}(1,jUser);
                if xUser ~= cUser
                    Px = SimStructs.baseStruct{iBase,1}.P{iBand,1}(:,jUser);
                    W = W + (H * Px)' * (H * Px);
                end
            end
            
            for jBase = 1:length(SimStructs.userStruct{cUser,1}.neighNode)
                ifBase = SimStructs.userStruct{cUser,1}.neighNode(jBase,1);
                Hk = SimStructs.linkChan{ifBase,iBand}(:,:,cUser);
                for kUser = 1:length(groupUsers{ifBase,1})
                    Pi = SimStructs.baseStruct{ifBase,1}.P{iBand,1}(:,kUser);
                    W = W + (Hk * Pi)' * (Hk * Pi);
                end
            end
            
            SimStructs.userStruct{cUser,1}.W{iBand,1}(:,cStream) = (H * P)' / W;
            
        end
    end
    
end

