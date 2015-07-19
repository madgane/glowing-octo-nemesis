
function [configParams] = displayQueueStatusGlobal(folderName,configParams,switchCase)

iIndex = 1;
xSystemConfig = cell(0);
listOfFiles = dir(folderName);

for iFile = 1:length(listOfFiles)
    if ~(listOfFiles(iFile).isdir)
        stringT = sprintf('%s/%s',folderName,listOfFiles(iFile).name);
        load(stringT);
        if ~isempty(strfind(listOfFiles(iFile).name,'-'))
            if ~isempty(find(SimParams.nExchangesOTA == [3,5,10,50]))
                display(SimParams.LegendName);
                SimParams.Log.Clock
                xSystemConfig{iIndex,1} = SimParams;
                iIndex = iIndex + 1;
            end
        else
            fprintf('Running on System - %s',xConfig.HostName);
            xConfig.Clock
        end
        clearvars SimParams SimStructs;
    end
end

if ~isempty(xSystemConfig)
    switch switchCase
        case 1
            configParams = displayFigures(xSystemConfig,configParams);
        case 2
            configParams = displayCDF(xSystemConfig,configParams);
    end
end

end

function configParams = displayFigures(cParams,configParams)

xL = length(configParams.legendString);
figColor = configParams.figColor;
figLineType = {configParams.lineType};

figLineWidth = {1,1.5};
figMarker = {'o','+','h','x','*'};
legendString = configParams.legendString;

fLwidth = length(figLineWidth);
fMarker = length(figMarker);

combTypeA = repmat((1:fLwidth),fMarker,1);
combTypeB = repmat((1:fMarker)',fLwidth,1);
combType = [combTypeB, combTypeA(:)];

for iParam = 1:length(cParams)
    
    xConfigParams = cParams{iParam,1};
    lUsers = xConfigParams.nUsers;
    lDrops = xConfigParams.nDrops + 1;
    lPackets = length(xConfigParams.maxArrival);
    
    packetServiceOverTime = zeros(lPackets,lUsers,lDrops);
    packetArrivalsOverTime = zeros(lPackets,lUsers,(lDrops - 1));
    queueBacklogsOverTime = zeros(lPackets,lUsers,lDrops);
    queueResiduesOverTime = zeros(lPackets,lUsers,lDrops);
    
    for iPkt = 1:lPackets
        if ~isempty(xConfigParams.Log.Clock.E{iPkt})
            for iUser = 1:lUsers
                packetServiceOverTime(iPkt,iUser,:) = xConfigParams.QueueInfo.packetServiceOverTime(1,iUser,iPkt,:);
                packetArrivalsOverTime(iPkt,iUser,:) = xConfigParams.QueueInfo.packetArrivalsOverTime(1,iUser,iPkt,:);
                queueBacklogsOverTime(iPkt,iUser,:) = xConfigParams.QueueInfo.queueBacklogsOverTime(1,iUser,iPkt,:);
                queueResiduesOverTime(iPkt,iUser,:) = xConfigParams.QueueInfo.queueResiduesOverTime(1,iUser,iPkt,:);
            end
        end
    end
    
    fcIndex = 1;
    fltIndex = 1;
    fmIndex = combType(iParam,1);
    flwIndex = combType(iParam,2);
    
    figure(1);hold on;box on;grid on;
    xlabel('avg. arrival pkts per user');
    ylabel('avg. transmitted pkts for all users');
    yValues = mean(sum(packetServiceOverTime,2),3);
    semilogy(xConfigParams.maxArrival,yValues,'Color',figColor{1,fcIndex},'LineWidth',figLineWidth{1,flwIndex},...
        'LineStyle',figLineType{1,fltIndex},'MarkerFaceColor',figColor{1,fcIndex},'Marker',figMarker{1,fmIndex},'MarkerSize',4);
    
    figure(2);hold on;box on;grid on;
    xlabel('avg. arrival pkts per user');
    ylabel('avg. backlogged pkts for all users');
    yValues = mean(sum(queueBacklogsOverTime,2),3);
    semilogy(xConfigParams.maxArrival,yValues,'Color',figColor{1,fcIndex},'LineWidth',figLineWidth{1,flwIndex},...
        'LineStyle',figLineType{1,fltIndex},'MarkerFaceColor',figColor{1,fcIndex},'Marker',figMarker{1,fmIndex},'MarkerSize',4);

    figure(3);hold on;box on;grid on;
    xlabel('avg. arrival pkts per user');
    ylabel('avg. residual pkts for all users');
    yValues = mean(sum(queueResiduesOverTime,2),3);
    semilogy(xConfigParams.maxArrival,yValues,'Color',figColor{1,fcIndex},'LineWidth',figLineWidth{1,flwIndex},...
        'LineStyle',figLineType{1,fltIndex},'MarkerFaceColor',figColor{1,fcIndex},'Marker',figMarker{1,fmIndex},'MarkerSize',4);
    
    legendString{1,(iParam + xL)} = xConfigParams.LegendName;
end

configParams.legendString = legendString;

figure(1);
figure(2);
figure(3);

end

function configParams = displayCDF(cParams,configParams)

xL = length(configParams.legendString);
figColor = configParams.figColor;
figLineType = {configParams.lineType};

figLineWidth = {1,1.5};
figMarker = {'o','+','h','x','*'};
legendString = configParams.legendString;

fLwidth = length(figLineWidth);
fMarker = length(figMarker);

combTypeA = repmat((1:fLwidth),fMarker,1);
combTypeB = repmat((1:fMarker)',fLwidth,1);
combType = [combTypeB, combTypeA(:)];

for iParam = 1:length(cParams)
    
    xConfigParams = cParams{iParam,1};
    lUsers = xConfigParams.nUsers;
    lDrops = xConfigParams.nDrops + 1;
    lPackets = length(xConfigParams.maxArrival);
    
    avgThroughput = zeros(lPackets,lUsers);
    fResidualPackets = zeros(lPackets,lUsers);
    packetServiceOverTime = zeros(lPackets,lUsers,lDrops);
    packetArrivalsOverTime = zeros(lPackets,lUsers,(lDrops - 1));
    queueBacklogsOverTime = zeros(lPackets,lUsers,lDrops);
    queueResiduesOverTime = zeros(lPackets,lUsers,lDrops);
    
    for iPkt = 1:lPackets
        if ~isempty(xConfigParams.Log.Clock.E{iPkt})
            for iUser = 1:lUsers
                avgThroughput(iPkt,iUser) = xConfigParams.Thrpt(1,iUser,iPkt) * (xConfigParams.nDrops / xConfigParams.iDrop);
                fResidualPackets(iPkt,iUser) = xConfigParams.QueueInfo.residualPkts(1,iUser,iPkt);
                packetServiceOverTime(iPkt,iUser,:) = xConfigParams.QueueInfo.packetServiceOverTime(1,iUser,iPkt,:);
                packetArrivalsOverTime(iPkt,iUser,:) = xConfigParams.QueueInfo.packetArrivalsOverTime(1,iUser,iPkt,:);
                queueBacklogsOverTime(iPkt,iUser,:) = xConfigParams.QueueInfo.queueBacklogsOverTime(1,iUser,iPkt,:);
                queueResiduesOverTime(iPkt,iUser,:) = xConfigParams.QueueInfo.queueResiduesOverTime(1,iUser,iPkt,:);
            end
        end
    end
    
    fcIndex = 1;
    fltIndex = 1;
    fmIndex = combType(iParam,1);
    flwIndex = combType(iParam,2);
    
    figure(1);hold on;box on;grid on;
    xlabel('avg. tx rate per RE');
    ylabel('prob. of users with gviven avg. tx rate');
    yValues = mean(mean(packetServiceOverTime,3),1);[ySeq,xSeq] = cdfcalc(yValues);
    plot(xSeq,ySeq(1:length(xSeq)),'Color',figColor{1,fcIndex},'LineWidth',figLineWidth{1,flwIndex},...
        'LineStyle',figLineType{1,fltIndex},'MarkerFaceColor',figColor{1,fcIndex},'Marker',figMarker{1,fmIndex},'MarkerSize',4);
    
    figure(2);hold on;box on;grid on;
    xlabel('avg. backlogged packets per RE');
    ylabel('prob. of users with given avg. backlogged packets');
    yValues = mean(mean(queueResiduesOverTime,3),1);[ySeq,xSeq] = cdfcalc(yValues);
    plot(xSeq,ySeq(1:length(xSeq)),'Color',figColor{1,fcIndex},'LineWidth',figLineWidth{1,flwIndex},...
        'LineStyle',figLineType{1,fltIndex},'MarkerFaceColor',figColor{1,fcIndex},'Marker',figMarker{1,fmIndex},'MarkerSize',4);
    
    figure(3);hold on;box on;grid on;
    xlabel('residual packets per RE');
    ylabel('prob. of users with given residual packets');
    yValues = mean(fResidualPackets,1);[ySeq,xSeq] = cdfcalc(yValues);
    plot(xSeq,ySeq(1:length(xSeq)),'Color',figColor{1,fcIndex},'LineWidth',figLineWidth{1,flwIndex},...
        'LineStyle',figLineType{1,fltIndex},'MarkerFaceColor',figColor{1,fcIndex},'Marker',figMarker{1,fmIndex},'MarkerSize',4);

    legendString{1,(iParam + xL)} = xConfigParams.LegendName;
end

configParams.legendString = legendString;

figure(1);
figure(2);
figure(3);

end


