
function [SimParams,SimStructs] = getCoordinateScheduling(SimParams,SimStructs)

for iUser = 1:SimParams.nUsers
    SimStructs.userStruct{iUser,1}.baseNode = [];
    SimStructs.userStruct{iUser,1}.neighNode = [];
end

for iBand = 1:SimParams.nBands
    
    Haug = [];iIndex = 0;
    Ucell = cell(SimParams.nUsers,SimParams.nBases);
    uLocs = zeros(SimParams.nUsers * SimParams.maxRank * SimParams.nBases,3);
    
    for iUser = 1:SimParams.nUsers
        for iBase = 1:SimParams.nBases
            
            H = SimStructs.linkChan{iBase,iBand}(:,:,iUser);
            [U,~,~] = svd(H);M = U' * H;
            Ucell{iUser,iBase} = U;
            
            if SimParams.queueWt
                M = M.' * (SimStructs.userStruct{iUser,1}.weighingFactor);
            else
                M = M.' * sign(SimStructs.userStruct{iUser,1}.weighingFactor);
            end
            Haug = [Haug M(:,1:SimParams.maxRank)];
            
            for iStream = 1:SimParams.maxRank
                iIndex = iIndex + 1;
                uLocs(iIndex,:) = [iUser iStream iBase];
            end
            
        end
    end
    
    Hm = Haug;
    nIterations = 5;
    xIndex = SimParams.nUsers * SimParams.maxRank * SimParams.nBases;
    
    ppVolume = zeros(xIndex,1);
    Nspace = cell(SimParams.nBases,1);
    assIndex = cell(SimParams.nBases,SimParams.nBands);
    assUsers = cell(SimParams.nBases,SimParams.nBands);
    assStreams = cell(SimParams.nBases,SimParams.nBands);
    
    for iIteration = 1:nIterations
    
        Aspace = cell(SimParams.nBases,1);
        userIterIndex = zeros(SimParams.nBases,1);

        for iRank = 1:SimParams.muxRank
            for iUser = 1:xIndex
                cBase = uLocs(iUser,3);
                M = [Aspace{cBase,1} Nspace{cBase,1} Hm(:,iUser)];
                ppVolume(iUser,1) = real(det(M' * M));
            end
            
            rIndex = 1;continueAgain = 1;
            [~,sortI] = sort(ppVolume,'descend');maxI = sortI(1,1);
            
            while continueAgain
                maxI = sortI(rIndex,1);
                cUser = uLocs(maxI,1);cStream = uLocs(maxI,2);cBase = uLocs(maxI,3);
                [~,totalLength] = size(Aspace{cBase,1});totalLength = totalLength + 1;
                if totalLength > (SimParams.muxRank / SimParams.nBases)
                    continueAgain = 1;
                    rIndex = rIndex + 1;
                else
                    continueAgain = 0;
                    userIterIndex(cBase,1) = userIterIndex(cBase,1) + 1;
                end
            end
            
            Aspace{cBase,1} = [Aspace{cBase,1} Hm(:,maxI)];
            assIndex{cBase,iBand}(userIterIndex(cBase,1),1) = maxI;
            assUsers{cBase,iBand}(userIterIndex(cBase,1),1) = cUser;
            assStreams{cBase,iBand}(userIterIndex(cBase,1),1) = cStream;
                        
            Nspace = cell(SimParams.nBases,1);
            for jBase = 1:SimParams.nBases
                for iBase = 1:SimParams.nBases
                    if iBase ~= jBase
                        for iUser = 1:length(assUsers{jBase,1})
                            kUser = assUsers{jBase,1}(iUser,1);kStream = assStreams{jBase,1}(iUser,1);
                            M = Ucell{kUser,jBase}(:,kStream)' * SimStructs.linkChan{iBase,iBand}(:,:,kUser);
                            Nspace{iBase,1} = [Nspace{iBase,1} M.'];
                        end
                    end
                end
            end
        end  
    end
    
end

for iBase = 1:SimParams.nBases
    
    usersLinked = [];
    for iBand = 1:SimParams.nBands
        usersLinked = [usersLinked ; assUsers{iBase,iBand}];
        SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1} = assUsers{iBase,iBand};
        SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1} = assStreams{iBase,iBand};
    end
    
    SimStructs.baseStruct{iBase,1}.linkedUsers = unique(usersLinked);
end

iBand = 1;
ovBaseNodes = 1:SimParams.nBases;
for iUser = 1:SimParams.nUsers
    assNode = [];
    for iBase = 1:SimParams.nBases
        if ~isempty(find(iUser == assUsers{iBase,iBand}))
            assNode = [assNode iBase];
        end
    end
    
    if ~isempty(assNode)
        SimStructs.userStruct{iUser,1}.baseNode = assNode;
        SimStructs.userStruct{iUser,1}.neighNode = ovBaseNodes(1,(assNode ~= ovBaseNodes));
    else
        SimStructs.userStruct{iUser,1}.baseNode = [];
        SimStructs.userStruct{iUser,1}.neighNode = [];
    end
end
