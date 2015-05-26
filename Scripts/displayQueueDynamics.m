
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

runIndex = 1;
pktIndex = 6;
figLineWidth = {2};
figLineType = {'-.','-'};
figColor = {'b','r','m',[0,0.6,0],[0,0.75,0.5],[0.3,0.7,0.9]};
figMarker = {'o','d','s','h','+','*'};
legendString = cell(1,1);

for iScheme = 1:globalCount
    
    randInt = randi([1,100],1,4);
    SimParams = SimParamsCell{iScheme,1};
    SimStructs = SimStructsCell{iScheme,1};
    displaySystemDetails;displayQueues(SimParams,SimStructs);
    
    fcIndex = mod(randInt(1,1) - 1,(length(figColor))) + 1;
    fmIndex = mod(randInt(1,2) - 1,(length(figMarker))) + 1;
    fltIndex = mod(randInt(1,3) - 1,(length(figLineType))) + 1;
    flwIndex = mod(randInt(1,4) - 1,(length(figLineWidth))) + 1;

%     figure(1);hold on;box on;grid on;
%     yValues = mean(squeeze(sum(squeeze(SimParams.QueueInfo.queueResiduesOverTime(:,:,:,2:end)),1)),2);
%     plot(SimParams.maxArrival,yValues,'Color',figColor{1,fcIndex},'LineWidth',figLineWidth{1,flwIndex},...
%         'LineStyle',figLineType{1,fltIndex},'MarkerFaceColor',figColor{1,fcIndex},'Marker',figMarker{1,fmIndex},'MarkerSize',6);
    
    figure(2);hold on;box on;grid on;
    yValues = mean(squeeze(sum(squeeze(SimParams.QueueInfo.queueBacklogsOverTime(:,:,:,2:end)),1)),2);
    plot(SimParams.maxArrival,yValues,'Color',figColor{1,fcIndex},'LineWidth',figLineWidth{1,flwIndex},...
        'LineStyle',figLineType{1,fltIndex},'MarkerFaceColor',figColor{1,fcIndex},'Marker',figMarker{1,fmIndex},'MarkerSize',6);

%     figure(3);hold on;box on;grid on;
%     yValues = sum(squeeze(SimParams.QueueInfo.residualPkts),1);
%     plot(SimParams.maxArrival,yValues,'Color',figColor{1,fcIndex},'LineWidth',figLineWidth{1,flwIndex},...
%         'LineStyle',figLineType{1,fltIndex},'MarkerFaceColor',figColor{1,fcIndex},'Marker',figMarker{1,fmIndex},'MarkerSize',6);

%     figure(4);hold on;box on;grid on;
%     yValues = sum(squeeze(SimParams.QueueInfo.queueResiduesOverTime(end,:,pktIndex,:)));
%     plot(yValues,'Color',figColor{1,fcIndex},'LineWidth',1,...
%         'LineStyle',figLineType{1,fltIndex},'MarkerFaceColor',figColor{1,fcIndex},'Marker',figMarker{1,fmIndex},'MarkerSize',2);

%     legendString{1,runIndex} = SimParams.weightedSumRateMethod;
    legendString{1,runIndex} = strcat(sprintf('%s-%s-%d-%d',SimParams.weightedSumRateMethod,SimParams.additionalParams,SimParams.queueWt,SimParams.nExchangesOTA));

    runIndex = runIndex + 1;

end


figure(1);
legend(legendString);
xlabel('Average packet arrivals (in bits)');
ylabel('Average backlogged packets (in bits)');

figure(2);
legend(legendString);
xlabel('Average packet arrivals (in bits)');
ylabel('Average backlogged packets (in bits)');

figure(3);
legend(legendString);
xlabel('Average packet arrivals (in bits)');
ylabel('Total residual packets left in the system (in bits)');

figure(4);
yValues = sum(squeeze(SimParams.QueueInfo.packetArrivalsOverTime(end,:,pktIndex,:)));
plot(yValues,'Color',[0 0 0],'LineWidth',1,'LineStyle','-','MarkerFaceColor',[0 0 0],'Marker','o','MarkerSize',2);
legend(legendString);
xlabel('Time Slots');
ylabel('Total number of residual packets in the system (in bits)');

% figure(1);
% yValues = mean(squeeze(sum(squeeze(SimParams.QueueInfo.packetArrivalsOverTime),1)),2);
% plot(SimParams.maxArrival,yValues,'Color','b','LineWidth',2,...
%     'LineStyle','-','MarkerFaceColor','b','Marker','o');
% legend(legendString);
% xlabel('Average packet arrivals (in bits)');
% ylabel('Average backlogged packets (in bits)');
% 
% figure(2);
% yValues = mean(squeeze(sum(squeeze(SimParams.QueueInfo.packetArrivalsOverTime),1)),2);
% plot(SimParams.maxArrival,yValues,'Color','b','LineWidth',2,...
%     'LineStyle','-','MarkerFaceColor','b','Marker','o');
% legend(legendString);
% xlabel('Average packet arrivals (in bits)');
% ylabel('Average backlogged packets (in bits)');
% 
% figure(3);
% xlabel('Average packet arrivals (in bits)');
% ylabel('Total residual packets left in the system (in bits)');
% legend(legendString);
%        
% figure(4);
% yValues = sum(squeeze(SimParams.QueueInfo.packetArrivalsOverTime(end,:,pktIndex,:)));
% plot(yValues,'Color',[0 0.6 0],'LineWidth',1,'LineStyle','-','MarkerFaceColor',[0 0.6 0],'Marker','o','MarkerSize',4);
% legendString{1,iScheme} = SimParams.weightedSumRateMethod;
% 
% xlabel('Time Slots');
% ylabel('Total number of residual packets in the system (in bits)');
% legend(legendString);

end
