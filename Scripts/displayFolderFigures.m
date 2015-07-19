
function displayFolderFigures(varargin)

if nargin == 1
    switchCase = 1;
    folderName = varargin{1};
else
    switchCase = varargin{2};
    folderName = varargin{1};
end
    

close all;
listOfFiles = dir(folderName);
rmpath(sprintf('%s\\Debug',pwd));

figColor = {'b','r',[0,0.7,0],'m',[0.7,0.7,0],[0.7,0,0.7]};
figLineType = {'-','-.','--',':'};
configParams.legendString = {};

fLength = length(figColor);
fType = length(figLineType);

combTypeA = repmat((1:fLength)',fType,1);
combTypeB = repmat((1:fType),fLength,1);

cFile = 0;
combType = [combTypeA, combTypeB(:)];

for iFile = 1:length(listOfFiles)    
    
    cFile = cFile + 1;
    fltIndex = combType(cFile,1);
    flwIndex = combType(cFile,2);
    configParams.figColor = figColor(1,fltIndex);
    configParams.lineType = figLineType{1,flwIndex};

    if listOfFiles(iFile).isdir
        if (sum(strcmpi(listOfFiles(iFile).name,{'.','..','MSE-C','MSE-L'})) == 0)
            cFolder = sprintf('%s/%s',folderName,listOfFiles(iFile).name);
            [configParams] = displayQueueStatusGlobal(cFolder,configParams,switchCase);
        end
    end
    
end

figure(1);legend(configParams.legendString);
figure(2);legend(configParams.legendString);
figure(3);legend(configParams.legendString);

