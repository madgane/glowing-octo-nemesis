
clc;
clear all;

userVariable = 1;
maxRand = 123513;
whatToDo = 'plot';
nQueueStability = 0;
scriptFile = 'outFile_xilinx_1.mat';

if isunix
    scriptFile = sprintf('..//..//simResults//%s',scriptFile);
else
    scriptFile = sprintf('..\\..\\simResults\\%s',scriptFile);
end

load(scriptFile);
nOptions = length(gParams);
gLegend = cell(nOptions,1);
titleString = cell(1,1);

if strcmp(whatToDo,'analyze')
    for iIndex = 1:nOptions
        dispString = sprintf('Printing %d Data',iIndex);disp(dispString);
        display(gParams{iIndex,1}.pathLossModel);
        dispString = sprintf('Number of Users - %d',gParams{iIndex,1}.nUsers);disp(dispString);
        dispString = sprintf('Number of Bases - %d',gParams{iIndex,1}.nBases);disp(dispString);
        display(gParams{iIndex,1}.dispString);
        display('----------------------------------------');
    end    
else

    iIndex = 1;
    display(gParams{iIndex,1}.pathLossModel);
    dispString = sprintf('Number of Users - %d',gParams{iIndex,1}.nUsers);disp(dispString);
    dispString = sprintf('Number of Bases - %d',gParams{iIndex,1}.nBases);disp(dispString);
    
    if strcmp(gParams{1,1}.queueMode,'false')
        for iIndex = 1:nOptions
            randIndex = randi([1 maxRand],1,1);
            display(gParams{iIndex,1}.pathLossModel);
            if ~userVariable
                plotHandle = plot(gParams{iIndex,1}.snrIndex,gResults{iIndex,1}.sumThrpt);setPlotFeatures(plotHandle,randIndex);hold all;
            else
                plotHandle = plot(gParams{iIndex,1}.userCount,gResults{iIndex,1}.sumThrpt);setPlotFeatures(plotHandle,randIndex);hold all;
            end
            gLegend{iIndex,1} = stringFormatting(gParams{iIndex,1},'legend',nQueueStability,userVariable);
        end
    else
        for iIndex = 1:nOptions
            if nQueueStability
                plotHandle = cell(2,1);
                figure(1);plotHandle{1,1} = plot(gParams{iIndex,1}.maxArrival,std(squeeze(gResults{iIndex,1}.queueBackLogs),1));hold all;
                figure(2);plotHandle{2,1} = plot(gParams{iIndex,1}.maxArrival,mean(squeeze(gResults{iIndex,1}.queueBackLogs),1));hold all;
            else
                plotHandle = cell(1,1);
                plotHandle{1,1} = plot(mean(squeeze(gResults{iIndex,1}.queueBackLogsOverTime(1,:,1,:)),1));hold all;
            end
            
            randIndex = randi([1 maxRand],1,1);
            for xIndex = 1:length(plotHandle)
                setPlotFeatures(plotHandle{xIndex,1},randIndex);
            end
            gLegend{iIndex,1} = stringFormatting(gParams{iIndex,1},'legend',nQueueStability);
        end
    end
    
    for xIndex = 1:length(plotHandle)
        
        if nQueueStability
            figure(xIndex);
            nQueueStability = xIndex;            
        end
        
        plotString = stringFormatting(gParams{iIndex,1},'plot',nQueueStability,userVariable);
        
        legend(gLegend);
        xlabel(plotString.xLabel);
        ylabel(plotString.yLabel);
        title(plotString.plotTitle);
        grid on;
    end
    
end
