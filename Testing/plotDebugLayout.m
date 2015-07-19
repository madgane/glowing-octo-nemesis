
function plotDebugLayout(SimParams)

markerStrArray = {'d','o','p','>','<','v','^','s'};

hold all;
plot(SimParams.wrapCellLocArray(:,1),'Marker','^','MarkerSize',8,'MarkerFaceColor','m','Color','m','LineStyle','none');

colorRGB = rand(SimParams.nBases,3);

for xCell = 1:SimParams.nBases
    cColor = colorRGB(xCell,:);
    plot(SimParams.debugLayout.positionArray{xCell,1},'Marker',markerStrArray{1,mod(xCell,length(markerStrArray))+1},'MarkerSize',4,'MarkerFaceColor',cColor,'Color',cColor,'LineStyle','none');    
end
