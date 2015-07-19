function [SimParams,SimStructs] = userLayoutGeneration(SimParams,SimStructs)

losArray = zeros(SimParams.nUsers,1);
positionArray = cell(SimParams.nBases,1);

debugRSSI = zeros(SimParams.nUsers,1);
userLocIndices = cell(SimParams.nBases,1);
nCites = getCellsOverLayout(SimParams.nTiers,1);

xUser = 0;
hexSide = SimParams.sysConfig.ISD / sqrt(3);
eastRotRad = SimParams.sysConfig.layoutFeatures.layoutAngleFromEast * pi / 180;

for iCite = 1:nCites
    citeLocation = SimParams.wrapCellLocArray(iCite,1);
    for iSector = 1:SimParams.nSectors
        cCite = (iCite - 1) * SimParams.nSectors + iSector;
        while (length(userLocIndices{cCite,1}) < SimParams.perCiteUsers)
            reIterate = 1;
            userPosition = 0;
            while reIterate
                userPosition = getPointInRhombus(hexSide,eastRotRad);
                if abs(userPosition) > SimParams.sysConfig.layoutFeatures.minDistance
                    reIterate = 0;
                end
            end

            userPosition = citeLocation + userPosition;
            [rssiMeasures, losMeasures] = getRSSIMeasures(SimParams,userPosition);

            [maxRSSI, maxRSSIindex] = max(rssiMeasures,[],2);
            [sortV,sortI] = sort(maxRSSI,'descend');

            if maxRSSIindex(sortI(1,1),1) ~= 1
                continue;
            end

            if length(userLocIndices{sortI(1,1),1}) < SimParams.perCiteUsers
                xUser = xUser + 1;
                userLocIndices{sortI(1,1),1} = [userLocIndices{sortI(1,1),1} userPosition];
                SimStructs.userStruct{xUser,1}.phyParams.location = userPosition;

                SimStructs.userStruct{xUser,1}.phyParams.listedCites = sortI(1:(SimParams.nNeighbors + 1),1);
                SimStructs.userStruct{xUser,1}.phyParams.listedRSSI = sortV(1:(SimParams.nNeighbors + 1),1);

                for iNeighbor = 1:(SimParams.nNeighbors + 1)
                    SimStructs.userStruct{xUser,1}.losFading{sortI(iNeighbor,1),1} = losMeasures{sortI(iNeighbor,1),maxRSSIindex(sortI(iNeighbor,1))};
                end

                SimStructs.userStruct{xUser,1}.phyParams.restOfIF_Linear = sum(10.^(sortV((SimParams.nNeighbors + 2):end,1)./10));
                SimStructs.userStruct{xUser,1}.phyParams.restOfIF = 10 * log10(sum(10.^(sortV((SimParams.nNeighbors + 2):end,1)./10)));

                linRSSI = 10.^(sortV./10);
                debugRSSI(xUser,1) = 10 * log10(linRSSI(1,1) / (10^(SimParams.systemNoise / 10) + sum(linRSSI(2:end,1))));
                positionArray{SimStructs.userStruct{xUser,1}.phyParams.listedCites(1,1)} = [positionArray{SimStructs.userStruct{xUser,1}.phyParams.listedCites(1,1)} ; userPosition];
                losArray(xUser,1) = strcmp(SimStructs.userStruct{xUser,1}.losFading{sortI(1,1)},'true');

            end
        end
    end
end

SimParams.debugLayout.losArray = losArray;
SimParams.debugLayout.userRSSI = debugRSSI;
SimParams.debugLayout.positionArray = positionArray;

% figure(1);
% plotDebugLayout(SimParams);
% figure(2);
% cdfplot(debugRSSI);

end


% function [SimParams,SimStructs] = userLayoutGeneration(SimParams,SimStructs)
% 
% ROT_90 = (pi / 2);
% debugRSSI = zeros(SimParams.nUsers,1);
% positionArray = zeros(SimParams.nUsers,1);
% userLocIndices = cell(SimParams.nBases,1);
% nUsersOverCell = floor(SimParams.nUsers / SimParams.nBases);
% 
% xUser = 0;fillUsers = 1;
% hexSide = SimParams.sysConfig.ISD / sqrt(3);
% eastRotRad = SimParams.sysConfig.layoutFeatures.layoutAngleFromEast * pi / 180;
% dropUserSize = SimParams.sysConfig.ISD * (SimParams.nTiers - 1) + hexSide;
% 
% while fillUsers
%     
%     reIterate = 1;
%     while reIterate
%         userPosition = getPointInRhombus(dropUserSize,eastRotRad);
%         userPosition = userPosition * exp(sqrt(-1) * ROT_90);
%         separationM = min(min(abs(userPosition - SimParams.wrapCellLocArray(:,1))));
%         if separationM > SimParams.sysConfig.layoutFeatures.minDistance
%             reIterate = 0;            
%         end
%     end
%     
%     [rssiMeasures, losMeasures] = getRSSIMeasures(SimParams,userPosition);
%     
%     [maxRSSI, maxRSSIindex] = max(rssiMeasures,[],2);
%     [sortV,sortI] = sort(maxRSSI,'descend');
%     
%     if maxRSSIindex(sortI(1,1),1) ~= 1
%         continue;
%     end
%     
%     if length(userLocIndices{sortI(1,1),1}) < nUsersOverCell
%         xUser = xUser + 1;
%         userLocIndices{sortI(1,1),1} = [userLocIndices{sortI(1,1),1} userPosition];
%         SimStructs.userStruct{xUser,1}.phyParams.location = userPosition;
%         
%         SimStructs.userStruct{xUser,1}.phyParams.listedCites = sortI(1:(SimParams.nNeighbors + 1),1);
%         SimStructs.userStruct{xUser,1}.phyParams.listedRSSI = sortV(1:(SimParams.nNeighbors + 1),1);
%         
%         for iNeighbor = 1:(SimParams.nNeighbors + 1)
%             SimStructs.userStruct{xUser,1}.losFading{sortI(iNeighbor,1),1} = losMeasures{sortI(iNeighbor,1),maxRSSIindex(sortI(iNeighbor,1))};
%         end
%         
%         SimStructs.userStruct{xUser,1}.phyParams.restOfIF = 10 * log10(sum(10.^(sortV((SimParams.nNeighbors + 2):end,1)./10)));
%         
%         linRSSI = 10.^(sortV./10);
%         positionArray(xUser,1) = userPosition;
%         debugRSSI(xUser,1) = 10 * log10(linRSSI(1,1) / (10^(SimParams.systemNoise / 10) * SimParams.N + sum(linRSSI(2:end,1))));
%     end
%     
%     iCount = 0;
%     for iCell = 1:SimParams.nBases
%         if length(userLocIndices{iCell,1}) == nUsersOverCell
%             iCount = iCount + 1;
%         end            
%     end
%     
%     if iCount == SimParams.nBases
%         fillUsers = 0;
%     end
%     
% end
% 
% end


% function [SimParams,SimStructs] = userLayoutGeneration(SimParams,SimStructs)
%
% debugRSSI = zeros(SimParams.nUsers,1);
% userLocIndices = cell(SimParams.nBases,1);
% nUsersOverCell = floor(SimParams.nUsers / SimParams.nBases);
%
% xUser = 0;
% if SimParams.nSectors == 3
%     hexSide = SimParams.sysConfig.ISD / 3;
% else
%     hexSide = SimParams.sysConfig.ISD / 2;
% end
%
% dropDistance = (SimParams.nTiers - 1) * SimParams.sysConfig.ISD + hexSide;
% eastRotRad = SimParams.sysConfig.layoutFeatures.layoutAngleFromEast * pi / 180;
% systemNoise = SimParams.sysConfig.NoisePwr_dBm - 10 * log10(SimParams.sysConfig.usableTones);
%
% for iUser = 1:SimParams.nUsers
%
%     reIterate = 1;
%     while reIterate
%         userPosition = getPointInRhombus(dropDistance,0,eastRotRad,0);
%         minDistance = min(min(abs(SimParams.wrapCellLocArray - userPosition)));
%
%         if minDistance > SimParams.sysConfig.layoutFeatures.minDistance
%             reIterate = 0;
%         end
%     end
%
%     [rssiMeasures, losMeasures] = getRSSIMeasures(SimParams,userPosition);
%     [maxRSSI, maxRSSIindex] = max(rssiMeasures,[],2);
%     [sortV,sortI] = sort(maxRSSI,'descend');
%
%     if length(userLocIndices{sortI(1,1),1}) < nUsersOverCell
%         xUser = xUser + 1;
%         userLocIndices{sortI(1,1),1} = [userLocIndices{sortI(1,1),1} userPosition];
%         SimStructs.userStruct{xUser,1}.phyParams.location = userPosition;
%
%         SimStructs.userStruct{xUser,1}.phyParams.listedCites = sortI(1:(SimParams.nNeighbors + 1),1);
%         SimStructs.userStruct{xUser,1}.phyParams.listedRSSI = sortV(1:(SimParams.nNeighbors + 1),1);
%
%         for iNeighbor = 1:(SimParams.nNeighbors + 1)
%             SimStructs.userStruct{xUser,1}.losFading{sortI(iNeighbor,1),1} = losMeasures{sortI(iNeighbor,1),maxRSSIindex(sortI(iNeighbor,1))};
%         end
%
%         SimStructs.userStruct{xUser,1}.phyParams.restOfIF = 10 * log10(sum(10.^(sortV((SimParams.nNeighbors + 2):end,1)./10)));
%
%         linRSSI = 10.^(sortV./10);
%         debugRSSI(xUser,1) = 10 * log10(linRSSI(1,1) / (10^(systemNoise / 10) * SimParams.N + sum(linRSSI(2:end,1))));
%
%     end
% end
%
% end
