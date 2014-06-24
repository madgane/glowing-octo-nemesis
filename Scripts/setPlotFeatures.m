
function setPlotFeatures(plotHandle,randIndex)

colorCell = cell(1,1);
markerCell = cell(1,1);
lineTypeCell = cell(1,1);
lineWidthCell = cell(1,1);

rIndex = 1;
markerCell{rIndex,1} = '.';
lineTypeCell{rIndex,1} = '-';
lineWidthCell{rIndex,1} = 1;
colorCell{rIndex,1} = '[0,0,0.7]';

rIndex = rIndex + 1;
markerCell{rIndex,1} = 'o';
lineTypeCell{rIndex,1} = ':';
lineWidthCell{rIndex,1} = 2;
colorCell{rIndex,1} = '[0.7,0,0]';

rIndex = rIndex + 1;
markerCell{rIndex,1} = 'x';
lineTypeCell{rIndex,1} = '-.';
colorCell{rIndex,1} = '[0,0.8,0.8]';

rIndex = rIndex + 1;
markerCell{rIndex,1} = '+';
lineTypeCell{rIndex,1} = '--';
colorCell{rIndex,1} = '[0.8,0,0.8]';

rIndex = rIndex + 1;
markerCell{rIndex,1} = '*';
colorCell{rIndex,1} = '[0.75,0,0.75]';

rIndex = rIndex + 1;
markerCell{rIndex,1} = 's';
colorCell{rIndex,1} = '[0.5,0.5,0.5]';

rIndex = rIndex + 1;
markerCell{rIndex,1} = 'd';
colorCell{rIndex,1} = '[0,0.5,0]';

rIndex = rIndex + 1;
markerCell{rIndex,1} = 'v';
colorCell{rIndex,1} = '[0,0.75,0.75]';

rIndex = rIndex + 1;
markerCell{rIndex,1} = '>';
colorCell{rIndex,1} = '[0.6,0.2,0]';

rIndex = rIndex + 1;
markerCell{rIndex,1} = 'h';
colorCell{rIndex,1} = '[0.08,0.17,0.55]';

rIndex = rIndex + 1;
markerCell{rIndex,1} = 'p';
colorCell{rIndex,1} = '[0.68,0.47,0]';

rIndex = rIndex + 1;
colorCell{rIndex,1} = '[0.87,0.49,0]';

rIndex = rIndex + 1;
colorCell{rIndex,1} = '[0.85,0.7,1]';

rIndex = rIndex + 1;
colorCell{rIndex,1} = '[1,0.69,0.39]';

nColorCells = length(colorCell);
nMarkerCells = length(markerCell);
nLineTypeCells = length(lineTypeCell);
nLineWidthCells = length(lineWidthCell);

randI = randIndex;
colIndex = mod((randI - 1),nColorCells) + 1;
markIndex = mod((randI - 1),nMarkerCells) + 1;
linTIndex = mod((randI - 1),nLineTypeCells) + 1;
linWIndex = mod((randI - 1),nLineWidthCells) + 1;

set(plotHandle,'Marker',markerCell{markIndex,1});
set(plotHandle,'Color',eval(colorCell{colIndex,1}));
set(plotHandle,'LineStyle',lineTypeCell{linTIndex,1});
set(plotHandle,'LineWidth',lineWidthCell{linWIndex,1});
set(plotHandle,'MarkerFaceColor',eval(colorCell{colIndex,1}));

