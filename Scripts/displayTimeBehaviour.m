
function displayTimeBehaviour(varargin)

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
figLineType = {'-','--','-.',':'};
figColor = {'b','g','r','m','c','k'};
figMarker = {'.','x','+','v','p','s','d','o'};
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
    
    reply = input('Do you want to display in plot Y/N [Y]:','s');
    
    if strcmpi(reply,'Y')
        yValues = sum(squeeze(SimParams.QueueInfo.queueBacklogsOverTime(end,:,end,:)));
        plot(yValues,'Color',figColor{1,fcIndex},'LineWidth',figLineWidth{1,flwIndex},...
            'LineStyle',figLineType{1,fltIndex},'MarkerFaceColor',figColor{1,fcIndex},'Marker',figMarker{1,fmIndex});        
    end
    
    hold all;
    legendString{1,iScheme} = SimParams.weightedSumRateMethod;

end

box on;
legend(legendString);

end
