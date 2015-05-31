
function displayQueueStatusGlobal(folderName)

close all;
iIndex = 1;
xSystemConfig = cell(1);
listOfFiles = dir(folderName);
rmpath(sprintf('%s/Debug',pwd));

for iFile = 1:length(listOfFiles)
    if ~(listOfFiles(iFile).isdir)
        stringT = sprintf('%s/%s',folderName,listOfFiles(iFile).name);
        load(stringT);
        if ~isempty(strfind(listOfFiles(iFile).name,'-'))
            display(SimParams.LegendName);
            SimParams.Log.Clock
            xSystemConfig{iIndex,1} = SimParams;
            iIndex = iIndex + 1;
        else
            fprintf('Running on System - %s',xConfig.HostName);
            xConfig.Clock
        end
        clearvars SimParams SimStructs;
    end
end

displayFigures(xSystemConfig);

end

function displayFigures(cParams)

figLineWidth = {1};
figLineType = {'-.','-'};
figColor = {'b','r','m',[0,0.6,0],[0,0.75,0.5],[0.3,0.7,0.9]};
figMarker = {'o','d','s','h','+','*'};
legendString = cell(1,1);

for iParam = 1:length(cParams)
    
    xConfigParams = cParams{iParam,1};
    lUsers = xConfigParams.nUsers;
    lDrops = xConfigParams.nDrops + 1;
    lPackets = length(xConfigParams.maxArrival);
    
    randInt = randi([1,100],1,4);
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
    
    fcIndex = mod(randInt(1,1) - 1,(length(figColor))) + 1;
    fmIndex = mod(randInt(1,2) - 1,(length(figMarker))) + 1;
    fltIndex = mod(randInt(1,3) - 1,(length(figLineType))) + 1;
    flwIndex = mod(randInt(1,4) - 1,(length(figLineWidth))) + 1;
    
    
    figure(1);hold on;box on;grid on;
    xlabel('avg. arrival pkts per user');
    ylabel('avg. transmitted pkts for all users');
    yValues = mean(sum(packetServiceOverTime,2),3);
    plot(xConfigParams.maxArrival,yValues,'Color',figColor{1,fcIndex},'LineWidth',figLineWidth{1,flwIndex},...
        'LineStyle',figLineType{1,fltIndex},'MarkerFaceColor',figColor{1,fcIndex},'Marker',figMarker{1,fmIndex},'MarkerSize',4);
    
    figure(2);hold on;box on;grid on;
    xlabel('avg. arrival pkts per user');
    ylabel('avg. backlogged pkts for all users');
    yValues = mean(sum(queueBacklogsOverTime,2),3);
    plot(xConfigParams.maxArrival,yValues,'Color',figColor{1,fcIndex},'LineWidth',figLineWidth{1,flwIndex},...
        'LineStyle',figLineType{1,fltIndex},'MarkerFaceColor',figColor{1,fcIndex},'Marker',figMarker{1,fmIndex},'MarkerSize',4);

    figure(3);hold on;box on;grid on;
    xlabel('avg. arrival pkts per user');
    ylabel('avg. residual pkts for all users');
    yValues = mean(sum(queueResiduesOverTime,2),3);
    plot(xConfigParams.maxArrival,yValues,'Color',figColor{1,fcIndex},'LineWidth',figLineWidth{1,flwIndex},...
        'LineStyle',figLineType{1,fltIndex},'MarkerFaceColor',figColor{1,fcIndex},'Marker',figMarker{1,fmIndex},'MarkerSize',4);
    
    legendString{1,iParam} = xConfigParams.LegendName;
end

figure(1);legend(legendString);
figure(2);legend(legendString);
figure(3);legend(legendString);

end