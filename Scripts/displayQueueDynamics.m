
function displayQueueDynamics(varargin)

close all;

switch nargin
    case 1
        load(varargin{1});
        randI = 0;
    case 3
        globalCount = 1;
        SimParamsCell{1} = varargin{1};
        SimStructsCell{1} = varargin{2};
        randI = varargin{3};
end

figLineWidth = {1};
figLineType = {'-.'};
figColor = {'b',[0 0.5 0],'r','m',[0 0.75 0.75],[0.5 0.5 0]};
figMarker = {'o','v','<','d','s','p','v','^'};
legendString = cell(1,globalCount + 1);

for iScheme = 1:globalCount
    
    clc;
    SimParams = SimParamsCell{iScheme,1};
    SimStructs = SimStructsCell{iScheme,1};
    displaySystemDetails;displayQueues(SimParams,SimStructs);
    
    fcIndex = mod(randI + iScheme - 1,(length(figColor))) + 1;
    fmIndex = mod(randI + iScheme - 1,(length(figMarker))) + 1;
    fltIndex = mod(randI + iScheme - 1,(length(figLineType))) + 1;
    flwIndex = mod(randI + iScheme - 1,(length(figLineWidth))) + 1;

    figure(1);hold on;box on;grid on;
    yValues = mean(squeeze(sum(squeeze(SimParams.QueueInfo.queueResiduesOverTime),1)),2);
    plot(SimParams.maxArrival,yValues,'Color',figColor{1,fcIndex},'LineWidth',figLineWidth{1,flwIndex},...
        'LineStyle',figLineType{1,fltIndex},'MarkerFaceColor',figColor{1,fcIndex},'Marker',figMarker{1,fmIndex});
    legendString{1,iScheme} = sprintf('OTA iterations - %d',SimParams.nExchangesOTA);
    
    figure(2);hold on;box on;grid on;
    yValues = mean(squeeze(sum(squeeze(SimParams.QueueInfo.queueBacklogsOverTime),1)),2);
    plot(SimParams.maxArrival,yValues,'Color',figColor{1,fcIndex},'LineWidth',figLineWidth{1,flwIndex},...
        'LineStyle',figLineType{1,fltIndex},'MarkerFaceColor',figColor{1,fcIndex},'Marker',figMarker{1,fmIndex});
    legendString{1,iScheme} = sprintf('OTA iterations - %d',SimParams.nExchangesOTA);
end

legendString{1,globalCount + 1} = 'Average Arrival Rate';

figure(1);
yValues = mean(squeeze(sum(squeeze(SimParams.QueueInfo.packetArrivalsOverTime),1)),2);
plot(SimParams.maxArrival,yValues,'Color','b','LineWidth',2,...
    'LineStyle','-','MarkerFaceColor','b','Marker','o');
legend(legendString);

figure(2);
yValues = mean(squeeze(sum(squeeze(SimParams.QueueInfo.packetArrivalsOverTime),1)),2);
plot(SimParams.maxArrival,yValues,'Color','b','LineWidth',2,...
    'LineStyle','-','MarkerFaceColor','b','Marker','o');
legend(legendString);

end
