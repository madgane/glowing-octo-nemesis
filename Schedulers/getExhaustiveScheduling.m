function [SimParams,SimStructs] = getExhaustiveScheduling(SimParams,SimStructs)

chScheduler = char(SimParams.SchedType);
uScoreIndex = find(chScheduler == '_');
if isempty(uScoreIndex)
    scheduleMethod = SimParams.SchedType;
else
    scheduleMethod = chScheduler(uScoreIndex(1,1) + 1:end);
end

switch scheduleMethod
    
    case 'MaxCapacity'
        
        for iBase = 1:SimParams.nBases
            
            uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
            kUsers = length(uIndices);
            
            userArray = getExhaustiveArray(kUsers,SimParams.muxRank);
            [nGroups,~] = size(userArray);
            
            for iBand = 1:SimParams.nBands
                
                chINVpwr = zeros(nGroups,1);
                eH = SimStructs.linkChan{iBase,iBand}(:,:,uIndices);
                
                for iGroup = 1:nGroups
                    
                    augH = [];
                    for iStream = 1:SimParams.muxRank
                        augH = [augH ; eH(:,:,userArray(iGroup,iStream))];
                    end
                    eP = pinv(augH);
                    chINVpwr(iGroup,1) = trace(eP' * eP);
                    
                end
                
                [~,minI] = min(chINVpwr);
                SimStructs.baseStruct{iBase}.assignedUsers{iBand,1} = uIndices(userArray(minI,:));
                SimStructs.baseStruct{iBase}.assignedStreams{iBand,1} = ones(SimParams.muxRank,1);
                
            end
            
        end
        
    case 'MinQueue'
        
        for iBase = 1:SimParams.nBases
            
            uIndices = SimStructs.baseStruct{iBase,1}.linkedUsers;
            kUsers = length(uIndices);
            
            QueueBacklogs = zeros(kUsers,1);
            userArray = getExhaustiveArray(kUsers,SimParams.muxRank);
            [nGroups,~] = size(userArray);
            futureQueueSize = zeros(nGroups,1);
            
            for iBand = 1:SimParams.nBands
                
                for iUser = 1:kUsers
                    QueueBacklogs(iUser,1) = SimStructs.userStruct{uIndices(iUser,1),1}.weighingFactor;
                end
                
                eH = SimStructs.linkChan{iBase,iBand}(:,:,uIndices);
                
                for iGroup = 1:nGroups
                    
                    augH = [];
                    for iStream = 1:SimParams.muxRank
                        augH = [augH ; eH(:,:,userArray(iGroup,iStream))];
                        Q(iStream,1) = QueueBacklogs(userArray(iGroup,iStream));
                    end
                    
                    eP = pinv(augH);
                    if SimParams.queueWt == 2
                        eP = performQueuedWF(eP,SimStructs.baseStruct{iBase,1}.sPower(1,iBand),Q);
                    else
                        eP = performWFAlgorithm(eP,SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
                    end
                    
                    for iUser = 1:kUsers
                        uIndex = find(iUser == userArray(iGroup,:));
                        if ~isempty(uIndex)
                            SNR = eH(:,:,iUser) * eP(:,uIndex);
                            futureQueueSize(iGroup,1) = futureQueueSize(iGroup,1) + max(QueueBacklogs(iUser,1) - log2(1 + SNR' * SNR),0);
                        else
                            futureQueueSize(iGroup,1) = futureQueueSize(iGroup,1) + QueueBacklogs(iUser,1);
                        end
                    end
                end
                
                [~,minI] = min(futureQueueSize);
                SimStructs.baseStruct{iBase}.assignedUsers{iBand,1} = uIndices(userArray(minI,:));
                SimStructs.baseStruct{iBase}.assignedStreams{iBand,1} = ones(SimParams.muxRank,1);
                
            end
            
        end
        
end

end

