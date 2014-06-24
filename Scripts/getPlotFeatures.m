
function [plotFeature] = getPlotFeatures()

markerType = cell(1,1);lineType = cell(1,1);lineColor = cell(1,1);

xIndex = 1;
lineType{xIndex,1} = '-';
lineColor{xIndex,1} = 'b';
markerType{xIndex,1} = '.';

xIndex = xIndex + 1;
lineType{xIndex,1} = ':';
lineColor{xIndex,1} = 'g';
markerType{xIndex,1} = 'o';

xIndex = xIndex + 1;
lineType{xIndex,1} = '-.';
lineColor{xIndex,1} = 'r';
markerType{xIndex,1} = 'x';

xIndex = xIndex + 1;
lineType{xIndex,1} = '--';
lineColor{xIndex,1} = 'c';
markerType{xIndex,1} = '+';

xIndex = xIndex + 1;
lineColor{xIndex,1} = 'm';
markerType{xIndex,1} = '*';

xIndex = xIndex + 1;
lineColor{xIndex,1} = 'k';
markerType{xIndex,1} = 's';

xIndex = xIndex + 1;
markerType{xIndex,1} = 'd';

xIndex = xIndex + 1;
markerType{xIndex,1} = 'v';

xIndex = xIndex + 1;
markerType{xIndex,1} = '>';

xIndex = xIndex + 1;
markerType{xIndex,1} = 'p';

xIndex = xIndex + 1;
markerType{xIndex,1} = 'h';


nLT = length(lineType);
nLC = length(lineColor);
nMarkers = length(markerType);

currIndex = randi([1 100],1,1);
cColor = lineColor{mod(currIndex - 1,nLC) + 1,1};
cLineType = lineType{mod(currIndex - 1,nLT) + 1,1};
cMarker = markerType{mod(currIndex - 1,nMarkers) + 1,1};
plotFeature = sprintf('%s%s%s',cLineType,cMarker,cColor);

end
