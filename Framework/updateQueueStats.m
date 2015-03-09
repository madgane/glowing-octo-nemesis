
function [SimParams,SimStructs] = updateQueueStats(varargin)

SimParams = varargin{1,1};
SimStructs = varargin{1,2};

if nargin == 2
    
    for iUser = 1:SimParams.nUsers
        currentArrival = SimStructs.userStruct{iUser,1}.trafficConfig.currentArrival;
        currentResidual =  SimStructs.userStruct{iUser,1}.trafficStats.backLogPkt - SimStructs.userStruct{iUser,1}.lastThrpt;
        
        SimStructs.userStruct{iUser,1}.trafficStats.residualPkt = max(0,currentResidual);
        SimStructs.userStruct{iUser,1}.trafficStats.backLogPkt = max(0,currentResidual) + currentArrival;
        SimStructs.userStruct{iUser,1}.trafficStats.residuesOverTime(1,SimParams.iDrop) = max(0,currentResidual);
        SimStructs.userStruct{iUser,1}.trafficStats.backlogsOverTime(1,SimParams.iDrop) = max(0,currentResidual) + currentArrival;
        SimStructs.userStruct{iUser,1}.trafficHistory.pktService(1,SimParams.iDrop) = SimStructs.userStruct{iUser,1}.lastThrpt;
        SimStructs.userStruct{iUser,1}.lastThrpt = 0;
    end
    
else
    
    for iUser = 1:SimParams.nUsers
        currentArrival = 0;
        currentResidual =  SimStructs.userStruct{iUser,1}.trafficStats.backLogPkt - SimStructs.userStruct{iUser,1}.lastThrpt;
        
        SimStructs.userStruct{iUser,1}.trafficStats.residualPkt = max(0,currentResidual);
        SimStructs.userStruct{iUser,1}.trafficStats.backLogPkt = max(0,currentResidual) + currentArrival;
        SimStructs.userStruct{iUser,1}.trafficStats.residuesOverTime(1,SimParams.iDrop + 1) = max(0,currentResidual);
        SimStructs.userStruct{iUser,1}.trafficStats.backlogsOverTime(1,SimParams.iDrop + 1) = max(0,currentResidual) + currentArrival;
        SimStructs.userStruct{iUser,1}.trafficHistory.pktService(1,SimParams.iDrop + 1) = SimStructs.userStruct{iUser,1}.lastThrpt;
        SimStructs.userStruct{iUser,1}.lastThrpt = 0;
    end
    
end

end