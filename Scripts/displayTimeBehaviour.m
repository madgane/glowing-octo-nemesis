
function displayTimeBehaviour(varargin)

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
legendString = cell(1,globalCount);

for iScheme = 1:globalCount
    
    clc;
    SimParams = SimParamsCell{iScheme,1};
    SimStructs = SimStructsCell{iScheme,1};
    displaySystemDetails;displayQueues(SimParams,SimStructs);
    
    fcIndex = mod(randI + iScheme - 1,(length(figColor))) + 1;
    fmIndex = mod(randI + iScheme - 1,(length(figMarker))) + 1;
    fltIndex = mod(randI + iScheme - 1,(length(figLineType))) + 1;
    flwIndex = mod(randI + iScheme - 1,(length(figLineWidth))) + 1;
    
    figure(1);hold all;
    yValues = sum(squeeze(SimParams.QueueInfo.queueBacklogsOverTime(end,:,end,:)));
    plot(yValues,'Color',figColor{1,fcIndex},'LineWidth',figLineWidth{1,flwIndex},...
        'LineStyle',figLineType{1,fltIndex},'MarkerFaceColor',figColor{1,fcIndex},'Marker',figMarker{1,fmIndex});

    figure(3);hold all;
    yValues = sum(squeeze(SimParams.QueueInfo.packetServiceOverTime(end,:,end,:)));
    plot(yValues,'Color',figColor{1,fcIndex},'LineWidth',figLineWidth{1,flwIndex},...
        'LineStyle',figLineType{1,fltIndex},'MarkerFaceColor',figColor{1,fcIndex},'Marker',figMarker{1,fmIndex});

    legendString{1,iScheme} = SimParams.weightedSumRateMethod;
end

figure(1);
box on;
legend(legendString);
modifyFigure;

figure(3);
box on;
legend(legendString);
modifyFigure;

end
