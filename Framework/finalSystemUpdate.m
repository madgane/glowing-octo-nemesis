
[SimParams,SimStructs] = updateQueueStats(SimParams,SimStructs,1);

for iUser = 1:SimParams.nUsers
    SimParams.PFmetric(iSNR,iUser,iPkt) = SimStructs.userStruct{iUser}.PFmetric;
    SimParams.fairness(iSNR,iUser,iPkt) = SimStructs.userStruct{iUser}.tAllocation / SimParams.utilityScale;
    SimParams.Thrpt(iSNR,iUser,iPkt) = (SimStructs.userStruct{iUser}.crThrpt - 1) / (SimParams.nDrops * SimParams.nBands);
    
    SimParams.QueueInfo.residualPkts(iSNR,iUser,iPkt) = SimStructs.userStruct{iUser,1}.trafficStats.residualPkt;
    SimParams.QueueInfo.queueBacklogs(iSNR,iUser,iPkt) = SimStructs.userStruct{iUser,1}.trafficStats.backLogPkt;
    SimParams.QueueInfo.packetServiceOverTime(iSNR,iUser,iPkt,:) = SimStructs.userStruct{iUser,1}.trafficHistory.pktService;
    SimParams.QueueInfo.packetArrivalsOverTime(iSNR,iUser,iPkt,:) = SimStructs.userStruct{iUser,1}.trafficHistory.pktArrival;
    SimParams.QueueInfo.queueBacklogsOverTime(iSNR,iUser,iPkt,:) = SimStructs.userStruct{iUser,1}.trafficStats.backlogsOverTime;
    SimParams.QueueInfo.queueResiduesOverTime(iSNR,iUser,iPkt,:) = SimStructs.userStruct{iUser,1}.trafficStats.residuesOverTime;        
end
