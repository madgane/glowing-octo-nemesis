
clc;
clear all;

nIters = 1000;
X = zeros(nIters,3);

SimStructs = 1;
SimParams.nTiers = 2;
SimParams.nSectors = 1;
SimParams.sysConfig.ISD = 300;
SimParams.nBases = getCellsOverLayout(SimParams.nTiers,SimParams.nSectors);
SimParams.nUsers = SimParams.nBases * 10;
SimParams.sysConfig.layoutFeatures.layoutAngleFromEast = 0;

[SimParams,~] = cellLayoutGeneration(SimParams,SimStructs);

% cellLayoutGeneration;
% mUserLocs = zeros(nIters,length(baseLocArray),3);

for iIter = 1:nIters
    X(iIter,1) = getPointInRhombus(SimParams.sysConfig.ISD/sqrt(3),0);
    X(iIter,2) = getPointInRhombus(SimParams.sysConfig.ISD/sqrt(3),0);
    X(iIter,3) = getPointInRhombus(SimParams.sysConfig.ISD/sqrt(3),0);
end

hold all;
plot(X(:,1),'r.');
plot(X(:,2),'b.');
plot(X(:,3),'m.');

plot(SimParams.wrapCellLocArray(:,1),'ro');

for iCell = 1:length(SimParams.wrapCellLocArray(:,1))
    plot(X(:,1) + SimParams.wrapCellLocArray(iCell,1),'b.');
end

for iCell = 1:length(SimParams.wrapCellLocArray(:,2))
    plot(X(:,1) + SimParams.wrapCellLocArray(iCell,2),'g.');
end

for iCell = 1:length(SimParams.wrapCellLocArray(:,3))
    plot(X(:,1) + SimParams.wrapCellLocArray(iCell,3),'r.');
end

for iCell = 1:length(SimParams.wrapCellLocArray(:,3))
    plot(X(:,1) + SimParams.wrapCellLocArray(iCell,4),'g.');
end

for iCell = 1:length(SimParams.wrapCellLocArray(:,3))
    plot(X(:,1) + SimParams.wrapCellLocArray(iCell,5),'r.');
end

for iCell = 1:length(SimParams.wrapCellLocArray(:,3))
    plot(X(:,1) + SimParams.wrapCellLocArray(iCell,6),'g.');
end

for iCell = 1:length(SimParams.wrapCellLocArray(:,3))
    plot(X(:,1) + SimParams.wrapCellLocArray(iCell,7),'r.');
end

break;

for wrapIndex = 1:7
    for iCell = 1:length(wrapCellArray(:,1))
        %     plot(X(:,1) + baseLocArray(iCell,1),'b.');
        %     plot(X(:,2) + baseLocArray(iCell,1),'g.');
        %     plot(X(:,3) + baseLocArray(iCell,1),'r.');
        
        mUserLocs(:,iCell,1) = X(:,1) + wrapCellArray(iCell,wrapIndex);
        mUserLocs(:,iCell,2) = X(:,2) + wrapCellArray(iCell,wrapIndex);
        mUserLocs(:,iCell,3) = X(:,3) + wrapCellArray(iCell,wrapIndex);
    end
    
    if wrapIndex == 1
        plot(mUserLocs(:,:,1),'b.');
        plot(mUserLocs(:,:,2),'g.');
        plot(mUserLocs(:,:,3),'r.');
    else
        plot(mUserLocs(:,:,1),'m.');
        plot(mUserLocs(:,:,2),'c.');
        plot(mUserLocs(:,:,3),'y.');
    end
end


% nCells = 6;
% tierOffset = pi / 6;
% tierAngle = 2 * pi / nCells;
%
% for iCell = 1:nCells
%     baseAngle = -tierOffset + iCell * tierAngle;
%     cmplxRotation = exp(sqrt(-1) * baseAngle) * SimParams.sysConfig.ISD * 4;
%     plot(mUserLocs(:,:,1) + cmplxRotation,'m.');
%     plot(mUserLocs(:,:,2) + cmplxRotation,'c.');
%     plot(mUserLocs(:,:,3) + cmplxRotation,'y.');
% end
%
% layoutParams.hBS = 25;
% layoutParams.hUT = 1.5;
% layoutParams.antTilt = 15;
% layoutParams.layoutAngleFromEast = 60;
%
% Y = zeros(nIters,1);
% for iIter = 1:nIters
%     Y(iIter,1) = getAntennaPatterGain(0+sqrt(-1)*0,X(iIter,1),layoutParams);
% end
%
% hold all;
% plot(angle(X) * 180 / pi,Y,'o')
%
%
% % layoutParams.layoutAngleFromEast = 0;
% % layoutParams.hUT = 0;
% % layoutParams.hBS = 0;
% %
% % nIters = 1000;
% % xUsers = zeros(nIters,3);
% % userLoc = zeros(nIters,1);
% %
% % userAngle = linspace(-pi,pi,nIters);
% %
% % baseLoc = 0 + 0j;
% % for iIter = 1:nIters
% %     userLoc(iIter,1) = exp(sqrt(-1) * userAngle(1,iIter)) * 100;
% %     for xSector = 1:3
% %         xUsers(iIter,xSector) = getAntennaPatterGain(baseLoc,userLoc(iIter,1),layoutParams,xSector);
% %     end
% % end
% %
% % hold all;
% % plot(userAngle * 180 / pi,xUsers(:,1));
% % plot(userAngle * 180 / pi,xUsers(:,2));
% % plot(userAngle * 180 / pi,xUsers(:,3));
% %
% %
