
function [SimParams,SimStructs] = initializeBuffers(SimParams)

nSINRSamples = length(SimParams.snrIndex);
nPacketSamples = length(SimParams.maxArrival);
SimStructs.userStruct = cell(SimParams.nUsers,1);
SimStructs.baseStruct = cell(SimParams.nBases,1);

SimParams.maxRank = min(SimParams.nRxAntenna,SimParams.nTxAntenna);
SimParams.muxRank = min(SimParams.nTxAntenna,(SimParams.nRxAntenna * SimParams.nUsers));

SimParams.Thrpt = zeros(nSINRSamples,SimParams.nUsers,nPacketSamples);
SimParams.fairness = zeros(nSINRSamples,SimParams.nUsers,nPacketSamples);
SimParams.utilityScale = SimParams.nDrops * SimParams.muxRank * SimParams.nBands;

SimParams.sumRateInstant = zeros(nSINRSamples,SimParams.nDrops,nPacketSamples);
SimParams.QueueInfo.queueBacklogs = zeros(nSINRSamples,SimParams.nUsers,nPacketSamples);
SimParams.QueueInfo.packetServiceOverTime = zeros(nSINRSamples,SimParams.nUsers,nPacketSamples,(SimParams.nDrops + 1));
SimParams.QueueInfo.queueBacklogsOverTime = zeros(nSINRSamples,SimParams.nUsers,nPacketSamples,(SimParams.nDrops + 1));
SimParams.txPower = zeros(length(SimParams.maxArrival),length(SimParams.snrIndex),SimParams.nBases,SimParams.nBands);

end

