function [SimParams, SimStructs] = cellLayoutGeneration(SimParams,SimStructs)

RAD_30 = (pi / 6);
RAD_60 = (pi / 3);
RAD_90 = (pi / 2);

xCell = 0;
hexAngleD = 60;
tierOffset = RAD_30;

baseLocArray = zeros(getCellsOverLayout(SimParams.nTiers,1),1);

for iTier = 1:SimParams.nTiers
    if iTier == 1
        nCells = 1;
    else
        nCells = 2^(iTier - 2) * 6;
    end
    
    tierAngle = 2 * pi / nCells;
    
    for iCell = 1:nCells
        xCell = xCell + 1;
        baseAngle = -tierOffset + (iCell - 1) * tierAngle;
        
        baseAngleD = round((baseAngle + tierOffset) * 180 / pi);
        isdDistance = abs(cosd(mod(baseAngleD,hexAngleD)));
        
        baseLocArray(xCell,1) = exp(sqrt(-1) * baseAngle) * SimParams.sysConfig.ISD * (iTier - 1) * isdDistance;
    end
end

% Wrap cell model ...

nCells = 6;
wrapCellArray = zeros(length(baseLocArray),(nCells + 1));
wrapCellArray(:,1) = baseLocArray;

for iCell = 1:nCells
    yShift = (SimParams.nTiers - 1) * exp(-sqrt(-1) * RAD_30);xShift = SimParams.nTiers * exp(-sqrt(-1) * RAD_90);
    cmplxDisplacement = SimParams.sysConfig.ISD * (xShift + yShift);
    wrapCellArray(:,iCell + 1) = baseLocArray + cmplxDisplacement * exp(sqrt(-1) * RAD_60 * (iCell - 1));
end

SimParams.wrapCellLocArray = wrapCellArray;

% function [SimParams, SimStructs] = cellLayoutGeneration(SimParams,SimStructs)
%
% RAD_30 = (pi / 6);
% RAD_60 = (pi / 3);
% RAD_90 = (pi / 2);
% RAD_15 = (pi / 12);
%
% xCell = 0;
% hexAngleD = 60;
% tierOffset = RAD_30;
%
% baseLocArray = zeros(getCellsOverLayout(SimParams.nTiers,1),1);
%
% if SimParams.nSectors == 3
%
%     for iTier = 1:SimParams.nTiers
%         if iTier == 1
%             nCells = 1;
%         else
%             nCells = 2^(iTier - 2) * 6;
%         end
%
%         tierAngle = 2 * pi / nCells;
%
%         for iCell = 1:nCells
%             xCell = xCell + 1;
%             baseAngle = -tierOffset + (iCell - 1) * tierAngle;
%
%             baseAngleD = round((baseAngle + tierOffset) * 180 / pi);
%             isdDistance = abs(cosd(mod(baseAngleD,hexAngleD)));
%
%             baseLocArray(xCell,1) = exp(sqrt(-1) * baseAngle) * SimParams.sysConfig.ISD * (iTier - 1) * isdDistance;
%         end
%     end
%
%     % Wrap cell model ...
%
%     nCells = 6;
%     wrapCellArray = zeros(length(baseLocArray),(nCells + 1));
%     wrapCellArray(:,1) = baseLocArray;
%
%     for iCell = 1:nCells
%         yShift = (SimParams.nTiers - 1) * exp(-sqrt(-1) * RAD_30);xShift = SimParams.nTiers * exp(-sqrt(-1) * RAD_90);
%         cmplxDisplacement = SimParams.sysConfig.ISD * (xShift + yShift);
%         wrapCellArray(:,iCell + 1) = baseLocArray + cmplxDisplacement * exp(sqrt(-1) * RAD_60 * (iCell - 1));
%     end
%
% else
%
%     tierOffset = 0;
%     for iTier = 1:SimParams.nTiers
%         if iTier == 1
%             nCells = 1;
%         else
%             nCells = 2^(iTier - 2) * 6;
%         end
%
%         tierAngle = 2 * pi / nCells;
%
%         for iCell = 1:nCells
%             xCell = xCell + 1;
%             baseAngle = -tierOffset + (iCell - 1) * tierAngle;
%
%             baseAngleD = round((baseAngle + tierOffset) * 180 / pi);
%             isdDistance = abs(cosd(mod(baseAngleD,hexAngleD)));
%
%             baseLocArray(xCell,1) = exp(sqrt(-1) * baseAngle) * SimParams.sysConfig.ISD * (iTier - 1) * isdDistance;
%         end
%     end
%
%     % Wrap cell model ...
%
%     nCells = 6;
%     wrapCellArray = zeros(length(baseLocArray),(nCells + 1));
%     wrapCellArray(:,1) = baseLocArray;
%
%     for iCell = 1:nCells
%         yShift = ((SimParams.nTiers) / sqrt(3)) * exp(-sqrt(-1) * RAD_30);xShift = ((SimParams.nTiers) / sqrt(3)) * exp(-sqrt(-1) * RAD_30);
%         cmplxDisplacement = 1.25 * SimParams.sysConfig.ISD * (xShift + yShift);
%         wrapCellArray(:,iCell + 1) = baseLocArray + cmplxDisplacement * exp(sqrt(-1) * RAD_60 * (iCell - 1)) * exp(sqrt(-1) * (RAD_15)/2);
%     end
%
% end
%
% SimParams.wrapCellLocArray = wrapCellArray;
%
% % hold all;
% % plot(SimParams.wrapCellLocArray(:,1),'k.');
% % plot(SimParams.wrapCellLocArray(:,2),'b.');
% % plot(SimParams.wrapCellLocArray(:,3),'g.');
% %
% % plot(SimParams.wrapCellLocArray(:,4),'m.');
% % plot(SimParams.wrapCellLocArray(:,5),'y.');
% % plot(SimParams.wrapCellLocArray(:,6),'c.');
% % plot(SimParams.wrapCellLocArray(:,7),'r.');
%
