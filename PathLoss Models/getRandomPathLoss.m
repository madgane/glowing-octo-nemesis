function [SimParams] = getRandomPathLoss(SimParams)

uscoreIndex = find(SimParams.pathLossModel == '_');
if isempty(uscoreIndex)
    minDistPerc = 0.48;
else
    minDistPerc = str2double(SimParams.pathLossModel(uscoreIndex + 1:end));
end

SimParams.nSectors = 1;
currentPLModel = SimParams.pathLossModel;
SimParams.PL_Profile = zeros(SimParams.nBases,SimParams.nUsers);

SimParams.pathLossModel = '3GPP_UMi';
[SimParams] = configureLTEParams(SimParams);
eastRotRad = SimParams.sysConfig.layoutFeatures.layoutAngleFromEast * pi / 180;
hexSide = SimParams.sysConfig.ISD / sqrt(3);

hexAngleD = 60;
baseLocArray = zeros(SimParams.nBases,1);
userLocArray = zeros(SimParams.nUsers,1);

switch SimParams.nBases
    
    case 1
        SimParams.pathLossModel = currentPLModel;
        SimParams.PL_Profile = -rand(SimParams.nBases,SimParams.nUsers) * minDistPerc;
        return;
%         baseLocArray(1,1) = 0;
%         for iUser = 1:SimParams.nUsers
%             reIterate = 1;
%             while reIterate
%                 userPosition = getPointInRhombus(hexSide,eastRotRad);
%                 if abs(userPosition) > SimParams.sysConfig.layoutFeatures.minDistance
%                     reIterate = 0;
%                 end
%             end
%             userLocArray(iUser,1) = userPosition;
%         end
        
    case 2
        tierOffset = pi / 6;
        tierAngle = 2 * pi / 6;
        baseLocArray(1,1) = 0 + 0j;
        for iBase = 2:SimParams.nBases
            baseAngle = -tierOffset + (iBase - 1) * tierAngle;
            baseAngleD = round((baseAngle + tierOffset) * 180 / pi);
            isdDistance = abs(cosd(mod(baseAngleD,hexAngleD)));
            baseLocArray(iBase,1) = exp(sqrt(-1) * baseAngle) * SimParams.sysConfig.ISD * isdDistance;
        end
        
        xUser = 0;
        awayFromBS = SimParams.sysConfig.ISD * minDistPerc;
        
        for iCell = 1:SimParams.nBases            
            xCells = 1:SimParams.nBases;            
            xVector = (baseLocArray(iCell,1) - sum(baseLocArray(xCells ~= iCell,1)));
            xVector = xVector / norm(xVector);
            for iUser = 1:(SimParams.nUsers / SimParams.nBases)
                reIterate = 1;
                while reIterate
                    userPosition =rand * SimParams.sysConfig.ISD * 0.5;
                    if abs(userPosition) > awayFromBS
                        reIterate = 0;
                    end
                end
                xUser = xUser + 1;
                userLocArray(xUser,1) = -xVector * userPosition + baseLocArray(iCell,1);
            end
        end        
end

for iBase = 1:SimParams.nBases
    cBaseLoc = baseLocArray(iBase,1);
    for iUser = 1:SimParams.nUsers
        cUserLoc = userLocArray(iUser,1);
        separationM = abs(cBaseLoc - cUserLoc);
        xRSSI = evaluateLTE_PL(SimParams,separationM,'true','true');
        antennaGain = getAntennaPatterGain(cBaseLoc,cUserLoc,SimParams,0);
        
        SimParams.PL_Profile(iBase,iUser) = xRSSI + antennaGain;
    end
end

SimParams.PL_Profile = SimParams.PL_Profile - SimParams.systemNoise;
SimParams.pathLossModel = currentPLModel;

% figure;
% cdfplot(max(SimParams.PL_Profile) - 10 * log10(10.^(min(SimParams.PL_Profile)./10) + 1));

% figure;
% plot(baseLocArray,'^r');
% hold all;plot(userLocArray,'ob');

end