function [SimParams,SimStructs] = systemLinking(SimParams,SimStructs)

uscoreIndex = find(SimParams.pathLossModel == '_');
if isempty(uscoreIndex)
    uscoreIndex = length(SimParams.pathLossModel) + 1;
end

if ~strcmp(SimParams.pathLossModel(1:uscoreIndex(1,1) - 1),'3GPP')
    
    for iBase = 1:SimParams.nBases
        SimStructs.baseStruct{iBase,1}.linkedUsers = [];
    end
    
    xBases = 1:SimParams.nBases;
    [~,maxI] = max(SimParams.PL_Profile,[],1);
    
    for iUser = 1:SimParams.nUsers
        cNode = maxI(1,iUser);
        SimStructs.userStruct{iUser,1}.baseNode = cNode;
        SimStructs.userStruct{iUser,1}.neighNode = find(xBases ~= cNode);
        SimStructs.baseStruct{cNode,1}.linkedUsers = [SimStructs.baseStruct{cNode,1}.linkedUsers ; iUser];
    end
    
else
    
    for iBase = 1:SimParams.nBases
        SimStructs.baseStruct{iBase,1}.linkedUsers = [];
    end
    
    for iUser = 1:SimParams.nUsers
        
        xCites = SimStructs.userStruct{iUser,1}.phyParams.listedCites;
        
        SimStructs.userStruct{iUser,1}.baseNode = xCites(1,1);
        SimStructs.userStruct{iUser,1}.neighNode = xCites(2:end,1)';
        SimStructs.baseStruct{xCites(1,1)}.linkedUsers = [SimStructs.baseStruct{xCites(1,1)}.linkedUsers ; iUser];
        
    end
    
end

if SimParams.multiCasting
    for iBase = 1:SimParams.nBases
        sIndex = 1;
        gPhaseShift = [0 90] * pi / 180;
        linkedUsers = SimStructs.baseStruct{iBase,1}.linkedUsers;
        cUsers = linkedUsers(randperm(length(linkedUsers)));
        SimStructs.baseStruct{iBase,1}.mcGroup = cell(length(SimParams.mcGroups{iBase,1}),1);
        for iGroup = 1:length(SimParams.mcGroups{iBase,1})
            eIndex = sum(SimParams.mcGroups{iBase,1}(1,1:iGroup));
            if iGroup ~= 1
                sIndex = sum(SimParams.mcGroups{iBase,1}(1,1:iGroup - 1)) + 1;
            end
            groupUserIndices = cUsers(sIndex:eIndex);
            SimStructs.baseStruct{iBase,1}.mcGroup{iGroup,1} = groupUserIndices;
            for iUser = 1:length(groupUserIndices)
                SimStructs.userStruct{groupUserIndices(iUser,1)}.phaseShift = -pi * sin(gPhaseShift(1,iGroup) + (iUser - 1) * pi / 90);
                SimStructs.userStruct{groupUserIndices(iUser,1)}.groupIndex = iGroup;
            end
        end
    end
end

end



