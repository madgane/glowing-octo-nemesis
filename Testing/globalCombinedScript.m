
function [SimParams,SimStructs,SimResults] = globalCombinedScript(SimParams)

nSINRSamples = length(SimParams.snrIndex);
nPacketSamples = length(SimParams.maxArrival);
SimStructs.userStruct = cell(SimParams.nUsers,1);
SimStructs.baseStruct = cell(SimParams.nBases,1);

SimParams.maxRank = min(SimParams.nRxAntenna,SimParams.nTxAntenna);
SimParams.muxRank = min(SimParams.nTxAntenna,(SimParams.nRxAntenna * SimParams.nUsers / SimParams.nBases));

SimParams.Thrpt = zeros(nSINRSamples,SimParams.nUsers,nPacketSamples);
utilityScale = SimParams.nDrops * SimParams.muxRank * SimParams.nBands;
SimParams.fairness = zeros(nSINRSamples,SimParams.nUsers,nPacketSamples);

queueBacklogs = zeros(nSINRSamples,SimParams.nUsers,nPacketSamples);
queueBacklogsOverTime = zeros(nSINRSamples,SimParams.nUsers,nPacketSamples,SimParams.nDrops);

for iPkt = 1:length(SimParams.maxArrival)
    
    SimParams.iPkt = iPkt;
    [SimParams,SimStructs] = fwkInitialization(SimParams,SimStructs);
    
    for iSNR = 1:length(SimParams.snrIndex)

        SimParams.N = 1;
        SimParams.iSNR = iSNR;
        SimParams.sPower = 10.^(SimParams.snrIndex(iSNR)/10);
        [SimParams,SimStructs] = systemInitialize(SimParams,SimStructs);
        [SimParams,SimStructs] = systemLinking(SimParams,SimStructs);

        % Resetting for every SNRs
        stream = RandStream.getGlobalStream;reset(stream);

        for iDrop = 1:SimParams.nDrops
            SimParams.iDrop = iDrop;
            [SimParams,SimStructs] = dropInitialize(SimParams,SimStructs);
            [SimParams,SimStructs] = getScheduledUsers(SimParams,SimStructs);
            [SimParams,SimStructs] = getPMatrix(SimParams,SimStructs);
            [SimParams,SimStructs] = performReception(SimParams,SimStructs);
        end

        for iUser = 1:SimParams.nUsers
            SimParams.PFmetric(iSNR,iUser,iPkt) = SimStructs.userStruct{iUser}.PFmetric;
            SimParams.fairness(iSNR,iUser,iPkt) = SimStructs.userStruct{iUser}.tAllocation / utilityScale;
            SimParams.Thrpt(iSNR,iUser,iPkt) = (SimStructs.userStruct{iUser}.crThrpt - 1) / (SimParams.nDrops * SimParams.nBands);
            queueBacklogs(iSNR,iUser,iPkt) = SimStructs.userStruct{iUser,1}.trafficStats.backLogPkt;
            queueBacklogsOverTime(iSNR,iUser,iPkt,:) = SimStructs.userStruct{iUser,1}.trafficStats.backlogsOverTime;
        end

        cState = sprintf('SINR completed - %d',SimParams.snrIndex(iSNR));disp(cState);
    end
    
end

testPktIndex = 1;testSINRIndex = 1;
SimParams.profiler.schX = SimParams.profiler.schX / (iSNR * iDrop);
   
fairnessS = squeeze(SimParams.fairness(:,:,testPktIndex));
throughputSumS = squeeze(SimParams.Thrpt(:,:,testPktIndex));
queueBacklogsOverTimeS = squeeze(queueBacklogsOverTime(testSINRIndex,:,testPktIndex,:));
    
SimParams.sumThrpt = sum(throughputSumS,2);
SimResults.sumCapacity = SimParams.sumThrpt;


JainMean = mean(throughputSumS,2).^2;JainVar = var(throughputSumS,0,2);
JainIndex_capacity = JainMean ./ (JainMean + JainVar);
SimResults.jainCapacity = JainIndex_capacity;

JainMean = mean(fairnessS,2).^2;JainVar = var(fairnessS,0,2);
JainIndex_utility = JainMean ./ (JainMean + JainVar);
SimResults.jainUtility = JainIndex_utility;

SimResults.queueBklgs = sum(queueBacklogsOverTimeS,1);
SimResults.snrIndices = SimParams.snrIndex;

end
