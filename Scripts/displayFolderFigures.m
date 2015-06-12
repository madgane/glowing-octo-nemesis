
function displayFolderFigures(folderName)

close all;
listOfFiles = dir(folderName);
rmpath(sprintf('%s/Debug',pwd));

figLineWidth = [0.5:0.5:2];
figLineType = {'-','-.',':','--'};
configParams.legendString = {};

for iFile = 1:length(listOfFiles)    
    fltIndex = mod(iFile - 1,(length(figLineType))) + 1;
    flwIndex = mod(iFile - 1,(length(figLineWidth))) + 1;
    configParams.lineWidth = figLineWidth(1,fltIndex);
    configParams.lineType = figLineType{1,flwIndex};
    if listOfFiles(iFile).isdir
        if (sum(strcmpi(listOfFiles(iFile).name,{'.','..','MSE-NM20-20','MSE-WM20-20','MSE-NM20-0','MSE-WM20-0','ADMM-WOS1','ADMM-WS1','CVX','TDM','MSE-L'})) == 0)
            cFolder = sprintf('%s/%s',folderName,listOfFiles(iFile).name);
            [configParams] = displayQueueStatusGlobal(cFolder,configParams);
        end
    else
        [configParams] = displayQueueStatusGlobal(folderName,configParams);
    end
end

figure(1);legend(configParams.legendString);
figure(2);legend(configParams.legendString);
figure(3);legend(configParams.legendString);