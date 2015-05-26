
function displayQueueStatusGlobal(folderName)

listOfFiles = dir(folderName);

for iFile = 1:length(listOfFiles)   
    if ~(listOfFiles(iFile).isdir)
        stringT = sprintf('%s/%s',folderName,listOfFiles(iFile).name);
        load(stringT);
        if ~isempty(strfind(stringT,'-'))
            display(SimParams.LegendName);
            SimParams.Log.Clock
        else
            display('Overall Duration');
            xConfig.Clock
        end
        clearvars SimParams SimStructs;
    end        
end