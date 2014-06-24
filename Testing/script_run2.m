
cd '..';

clc;clear;
pathAddition;

script_cells = cell(6,1);
script_file = './Scripts/Results/script_run_2.mat';

xIndex = 1;
script_cells{xIndex,1}.queueWt = 0;
script_cells{xIndex,1}.precoding = 'Best_ZF_Method';
script_cells{xIndex,1}.precType = 'PerformScheduling';
script_cells{xIndex,1}.schType = 'BDScheduling_SP';

xIndex = xIndex + 1;
script_cells{xIndex,1}.queueWt = 1;
script_cells{xIndex,1}.precoding = 'Best_ZF_Method';
script_cells{xIndex,1}.precType = 'PerformScheduling';
script_cells{xIndex,1}.schType = 'BDScheduling_SP';

xIndex = xIndex + 1;
script_cells{xIndex,1}.queueWt = 2;
script_cells{xIndex,1}.precoding = 'Best_ZF_Method';
script_cells{xIndex,1}.precType = 'PerformScheduling';
script_cells{xIndex,1}.schType = 'BDScheduling_SP';

xIndex = xIndex + 1;
script_cells{xIndex,1}.queueWt = 0;
script_cells{xIndex,1}.precoding = 'Best_WMMSE_Method';
script_cells{xIndex,1}.precType = 'StreamScheduling';
script_cells{xIndex,1}.schType = 'BDScheduling_SP';

xIndex = xIndex + 1;
script_cells{xIndex,1}.queueWt = 1;
script_cells{xIndex,1}.precoding = 'Best_WMMSE_Method';
script_cells{xIndex,1}.precType = 'StreamScheduling';
script_cells{xIndex,1}.schType = 'BDScheduling_SP';

xIndex = xIndex + 1;
script_cells{xIndex,1}.queueWt = 1;
script_cells{xIndex,1}.precoding = 'Best_WMMSE_Method';
script_cells{xIndex,1}.precType = 'PerformScheduling';
script_cells{xIndex,1}.schType = 'RRScheduling';

for iScript = 1:length(script_cells)
    
    SimParams.version = version;
    SimParams.DebugMode = 'false';
    SimParams.queueMode = 'true';
    
    SimParams.ChannelModel = 'IID';
    SimParams.pathLossModel = 'Random_40';
    SimParams.DopplerType = 'Constant_100';
    
    SimParams.weighingEqual = 'false';
    SimParams.queueWt = script_cells{iScript,1}.queueWt;
    SimParams.SchedType = script_cells{iScript,1}.schType;
    SimParams.PrecodingMethod = script_cells{iScript,1}.precoding;
    SimParams.weightedSumRateMethod = script_cells{iScript,1}.precType;
    
    display(SimParams.queueWt);
    display(SimParams.SchedType);
    display(SimParams.PrecodingMethod);
    display(SimParams.weightedSumRateMethod);
    
    SimParams.nDrops = 10000;
    SimParams.snrIndex = [20];
    
    SimParams.PF_dur = 40;
    SimParams.sampTime = 1e-3;
    SimParams.estError = 0.00;
    SimParams.fbFraction = 0.0;
    
    SimParams.nBands = 1;
    SimParams.nBases = 1;
    SimParams.nUsers = 50;
    
    SimParams.nTxAntenna = 4;
    SimParams.nRxAntenna = 1;
    
    SimParams.gracePeriod = 0;
    SimParams.arrivalDist = 'Constant';
    
    SimParams.maxArrival = [0.1:0.1:1];
    SimParams.FixedPacketArrivals = [10,10,10,10,10,10,1,1,1,1];
    SimParams.PL_Profile = [5 -inf 5 -inf 5 -inf 1e-20 0; -inf 5 -inf 5 -inf 5 0 1e-20];
    
    nSINRSamples = length(SimParams.snrIndex);
    nPacketSamples = length(SimParams.maxArrival);
    SimStructs.userStruct = cell(SimParams.nUsers,1);
    SimStructs.baseStruct = cell(SimParams.nBases,1);
    
    SimParams.maxRank = min(SimParams.nRxAntenna,SimParams.nTxAntenna);
    SimParams.muxRank = min(SimParams.nTxAntenna,(SimParams.nRxAntenna * SimParams.nUsers));
    
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
            resetRandomness;
            
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
            
            if strcmp(SimParams.DebugMode,'true')
                display(squeeze(queueBacklogs(iSNR,:,iPkt)));
            end
            
            cState = sprintf('SINR,Arrival - (%d,%2.2d)',SimParams.snrIndex(iSNR),SimParams.maxArrival(1,iPkt));disp(cState);
        end
        
    end
    
    testPktIndex = 1;testSINRIndex = 1;
    fairnessS = squeeze(SimParams.fairness(:,:,testPktIndex));
    throughputSumS = squeeze(SimParams.Thrpt(:,:,testPktIndex));
    
    queueBackLogsS = squeeze(queueBacklogs(testSINRIndex,:,:));
    queueBacklogsOverTimeS = squeeze(queueBacklogsOverTime(testSINRIndex,:,testPktIndex,:));
    
    SimResults.snrIndex = SimParams.snrIndex;
    
    if strcmp(SimParams.queueMode,'false')
        
        SimParams.sumThrpt = sum(throughputSumS,2);
        SimResults.sumThrpt = SimParams.sumThrpt;
        
        JainMean = mean(throughputSumS,2).^2;JainVar = var(throughputSumS,0,2);
        JainIndex_capacity = JainMean ./ (JainMean + JainVar);
        SimResults.jainCapacity = JainIndex_capacity;
        
        JainMean = mean(fairnessS,2).^2;JainVar = var(fairnessS,0,2);
        JainIndex_utility = JainMean ./ (JainMean + JainVar);
        SimResults.jainUtility = JainIndex_utility;
        
    else
        
        avgBacklogs = zeros(length(SimParams.maxArrival),1);
        for cPacket = 1:length(SimParams.maxArrival)
            cQueue = queueBacklogs(:,:,cPacket);cQueue = cQueue(:);
            avgBacklogs(cPacket,1) = sum(cQueue);
        end        
        
        SimResults.queueStats.backLogs = queueBacklogs;
        SimResults.queueStats.backLogHistory = queueBacklogsOverTime;
        SimResults.avgQueues = avgBacklogs;
        SimResults.nDrops = SimParams.nDrops;
        SimResults.queueBacklogs = sum(queueBacklogsOverTimeS,1);
        
    end
    
    hold all;
    plot(SimResults.avgQueues);
        
    if ~exist(script_file,'file')
        gIndex = 1;
        gResults = cell(1,1);gParams = cell(1,1);gStructs = cell(1,1);
        gResults{gIndex,1} = SimResults;gParams{gIndex,1} = SimParams;gStructs{gIndex,1} = SimStructs;
        save(script_file,'gResults','gParams','gStructs','gIndex');
    else
        load(script_file);gIndex = gIndex + 1;
        gResults{gIndex,1} = SimResults;gParams{gIndex,1} = SimParams;gStructs{gIndex,1} = SimStructs;
        save(script_file,'gResults','gParams','gStructs','gIndex');
    end
    
end


