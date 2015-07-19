
function [SimParams,SimStructs] = getGScheduling(SimParams,SimStructs)

nBands = SimParams.nBands;
nBases = SimParams.nBases;
maxRank = SimParams.maxRank;
cellUserIndices = cell(nBases,1);
cellUserQueues = cell(nBases,1);

for iBase = 1:nBases
    cellUserIndices{iBase,1} = SimStructs.baseStruct{iBase,1}.linkedUsers;
    cellUserQueues{iBase,1} = zeros(length(cellUserIndices{iBase,1}),1);
    for iUser = 1:length(cellUserIndices{iBase,1})
        cUser = cellUserIndices{iBase,1}(iUser,1);
        cellUserQueues{iBase,1}(iUser,1) = SimStructs.userStruct{cUser,1}.trafficStats.backLogPkt;
    end
end

charScheduling = char(SimParams.SchedType);
uscore_index = find(charScheduling == '_');
caseStudy = charScheduling(uscore_index + 1:end);

switch (caseStudy)
    
    case 'I-CELL'
        
        for iBand = 1:nBands
            for iBase = 1:nBases
                for jBase = 1:nBases
                    if iBase ~= jBase
                        for iUser = 1:length(cellUserIndices{jBase,1})
                            SimStructs.linkChan{iBase,iBand}(:,:,cellUserIndices{jBase,1}(iUser,1)) = zeros(SimParams.nRxAntenna,SimParams.nTxAntenna);
                        end
                    end
                end
            end
        end
        
    case 'TDM'

        if strcmpi(SimParams.sysMode,'true')
            modBases = SimParams.nSectors;
        else
            modBases = nBases;
        end

        for iBand = 1:nBands
            for iBase = 1:nBases
                if (mod(iBase - 1,modBases) == mod(SimParams.iDrop - 1,modBases))
                    SimStructs.baseStruct{iBase,1}.sPower(1,iBand) = SimStructs.baseStruct{iBase,1}.sPower(1,iBand) * modBases;
                else
                    SimStructs.baseStruct{iBase,1}.sPower(1,iBand) = SimStructs.baseStruct{iBase,1}.sPower(1,iBand) * 0;
                end
            end
        end
        
    case 'R-GROUP'
        
        xFrequency = SimParams.groupArrivalFreq / 2;
        if (mod((SimParams.iDrop - 1),xFrequency) == 0)
            for iBand = 1:nBands
                for iBase = 1:nBases
                    xUsers = floor(length(cellUserQueues{iBase,1}) / 2);
                    [~,qSortI] = sort(cellUserQueues{iBase,1},'descend'); 
                
                    if (mod(SimParams.iDrop - 1,nBases) == (iBase - 1))
                        selUsers = qSortI(1:xUsers);
                        uIndices = cellUserIndices{iBase,1}(selUsers,1);
                    else
                        selUsers = qSortI(xUsers + 1:end);
                        uIndices = cellUserIndices{iBase,1}(selUsers,1);
                    end
                
                    uIndices = repmat(uIndices',maxRank,1);uIndices = sort(uIndices(:));
                    SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1} = uIndices;
                    SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1} = repmat(linspace(1,maxRank,maxRank)',length(uIndices)/maxRank,1);
                end
            end
            
            SimParams.Debug.SchedBackup = cell(nBases,1);
            
            for iBase = 1:nBases
                SimParams.Debug.SchedBackup{iBase,1}.assignedUsers = SimStructs.baseStruct{iBase,1}.assignedUsers;
                SimParams.Debug.SchedBackup{iBase,1}.assignedStreams = SimStructs.baseStruct{iBase,1}.assignedStreams;
            end            
        else
            for iBase = 1:nBases
                SimStructs.baseStruct{iBase,1}.assignedUsers = SimParams.Debug.SchedBackup{iBase,1}.assignedUsers;
                SimStructs.baseStruct{iBase,1}.assignedStreams = SimParams.Debug.SchedBackup{iBase,1}.assignedStreams;
            end
        end
        
    case 'B-TDM'

        if (mod(SimParams.iDrop - 1,SimParams.groupArrivalFreq) == 0)
            randShfl = randi(SimParams.nBases,1,1);
            SimParams.Debug.randShfl = randShfl;
        end
        
        for iBand = 1:nBands
            for iBase = 1:nBases
                if (mod((iBand - 1 + SimParams.Debug.randShfl),nBases) ~= (iBase - 1))
                    SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1} = [];
                    SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1} = [];
                    SimStructs.baseStruct{iBase,1}.sPower(1,iBand) = SimStructs.baseStruct{iBase,1}.sPower(1,iBand) * 0;
                else
                    SimStructs.baseStruct{iBase,1}.sPower(1,iBand) = SimStructs.baseStruct{iBase,1}.sPower(1,iBand) * nBases;
                end
            end
        end
        
    otherwise
        
        
end
    
end
