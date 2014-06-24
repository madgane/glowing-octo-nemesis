
function plotDebugLayout(SimParams)

markerStrArray = {'d','o','p','>','<','+','x'};

hold all;
plot(SimParams.wrapCellLocArray(:,1),'Marker','^','MarkerSize',8,'MarkerFaceColor','m','Color','m','LineStyle','none');

xCell = 0;
for iCell = 1:(SimParams.nBases / SimParams.nSectors)
    xCell = xCell + 1;    
    plot(SimParams.debugLayout.positionArray{xCell,1},'Marker',markerStrArray{1,mod(iCell,length(markerStrArray))+1},'MarkerSize',4,'MarkerFaceColor','b','Color','b','LineStyle','none');
    xCell = xCell + 1;
    plot(SimParams.debugLayout.positionArray{xCell,1},'Marker',markerStrArray{1,mod(iCell,length(markerStrArray))+1},'MarkerSize',4,'MarkerFaceColor','g','Color','g','LineStyle','none');
    xCell = xCell + 1;
    plot(SimParams.debugLayout.positionArray{xCell,1},'Marker',markerStrArray{1,mod(iCell,length(markerStrArray))+1},'MarkerSize',4,'MarkerFaceColor','r','Color','r','LineStyle','none'); 
end

