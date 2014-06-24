
clc;clear;
script_cells = cell(1,1);

% xIndex = 1;
% script_cells{xIndex,1}.queueWt = 1;
% script_cells{xIndex,1}.precoding = 'Best_ZF_Method';
% script_cells{xIndex,1}.precType = 'PreScheduling';
% script_cells{xIndex,1}.schType = 'BDScheduling_SP';

% xIndex = xIndex + 1;
% script_cells{xIndex,1}.queueWt = 1;
% script_cells{xIndex,1}.precoding = 'Best_ZF_Method';
% script_cells{xIndex,1}.precType = 'PreScheduling';
% script_cells{xIndex,1}.schType = 'BDScheduling_RNS';

% xIndex = xIndex + 1;
% script_cells{xIndex,1}.queueWt = 1;
% script_cells{xIndex,1}.precoding = 'Best_ZF_Method';
% script_cells{xIndex,1}.precType = 'PreScheduling';
% script_cells{xIndex,1}.schType = 'PFBDScheduling_M-AHP';

% xIndex = xIndex + 1;
% script_cells{xIndex,1}.queueWt = 1;
% script_cells{xIndex,1}.precoding = 'Best_ZF_Method';
% script_cells{xIndex,1}.precType = 'PreScheduling';
% script_cells{xIndex,1}.schType = 'GreedyScheduling';

% xIndex = xIndex + 1;
% script_cells{xIndex,1}.queueWt = 1;
% script_cells{xIndex,1}.precoding = 'Best_ZF_Method';
% script_cells{xIndex,1}.precType = 'PreScheduling';
% script_cells{xIndex,1}.schType = 'RRScheduling';

xIndex = 1;
script_cells{xIndex,1}.queueWt = 1;
script_cells{xIndex,1}.precoding = 'Best_ZF_Method';
script_cells{xIndex,1}.precType = 'PreScheduling';
script_cells{xIndex,1}.schType = 'PFScheduling_BF';


SimParams.version = version;
SimParams.outFile = 'outFile_xilinx_1.mat';

pathAddition;
display(SimParams.outFile);

for iScript = 1:length(script_cells)
    
    SimParams.DebugMode = 'false';
    SimParams.queueMode = 'false';
    SimParams.weighingEqual = 'true';
    
    SimParams.ChannelModel = 'IID';
    SimParams.pathLossModel = 'CellEdge';
    SimParams.DopplerType = 'Constant_100';
    
    SimParams.queueWt = script_cells{iScript,1}.queueWt;
    SimParams.SchedType = script_cells{iScript,1}.schType;
    SimParams.PrecodingMethod = script_cells{iScript,1}.precoding;
    SimParams.weightedSumRateMethod = script_cells{iScript,1}.precType;
    
    SimParams.dispString = sprintf('%s-%s-%s',SimParams.SchedType,...
        SimParams.PrecodingMethod,SimParams.weightedSumRateMethod);
    display(SimParams.dispString);
    
    SimParams.nDrops = 1000;
    SimParams.snrIndex = [10];
    
    SimParams.PF_dur = 40;
    SimParams.sampTime = 1e-3;
    SimParams.estError = 0.00;
    SimParams.fbFraction = 0.0;
    
    SimParams.nBands = 1;
    SimParams.nBases = 1;
    SimParams.nUsers = 20;
    SimParams.userCount = [10:10:100];
    
    SimParams.nTxAntenna = 4;
    SimParams.nRxAntenna = 1;
    
    SimParams.mdpFactor = 0;
    SimParams.gracePeriod = 0;
    SimParams.arrivalDist = 'Constant';
    
    SimParams.maxArrival = 100;
    SimParams.FixedPacketArrivals = [10,10,10,10,10,10,1,1,1,1];
    SimParams.PL_Profile = [5 -inf 5 -inf 5 -inf 1e-20 0; -inf 5 -inf 5 -inf 5 0 1e-20];
    
    SimParams.maxRank = min(SimParams.nRxAntenna,SimParams.nTxAntenna);
    SimParams.muxRank = min(SimParams.nTxAntenna,(SimParams.nRxAntenna * SimParams.nUsers));
    
    nUserCount = length(SimParams.userCount);
    maxUserCnt = SimParams.userCount(1,nUserCount);
    
    nPacketSamples = length(SimParams.maxArrival);
    SimStructs.userStruct = cell(maxUserCnt,1);
    SimStructs.baseStruct = cell(SimParams.nBases,1);
    
    SimParams.Thrpt = zeros(nUserCount,SimParams.nUsers,nPacketSamples);
    utilityScale = SimParams.nDrops * SimParams.muxRank * SimParams.nBands;
    SimParams.fairness = zeros(nUserCount,maxUserCnt,nPacketSamples);
    
    queueBacklogs = zeros(nUserCount,maxUserCnt,nPacketSamples);
    queueBacklogsOverTime = zeros(nUserCount,maxUserCnt,nPacketSamples,SimParams.nDrops);
    SimParams.txPower = zeros(length(SimParams.maxArrival),length(SimParams.snrIndex),SimParams.nBases);
    
    for iPkt = 1:length(SimParams.maxArrival)
        
        SimParams.iPkt = iPkt;
                
        for iUserCount = 1:length(SimParams.userCount)
            
            iSNR = 1;
            SimParams.N = 1;
            SimParams.iSNR = iSNR;
            SimParams.nUsers = SimParams.userCount(1,iUserCount);
            SimParams.sPower = 10.^(SimParams.snrIndex(iSNR)/10);
                                    
            [SimParams,SimStructs] = fwkInitialization(SimParams,SimStructs);
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
                SimParams.PFmetric(iUserCount,iUser,iPkt) = SimStructs.userStruct{iUser}.PFmetric;
                SimParams.fairness(iUserCount,iUser,iPkt) = SimStructs.userStruct{iUser}.tAllocation / utilityScale;
                SimParams.Thrpt(iUserCount,iUser,iPkt) = (SimStructs.userStruct{iUser}.crThrpt - 1) / (SimParams.nDrops * SimParams.nBands);
                queueBacklogs(iUserCount,iUser,iPkt) = SimStructs.userStruct{iUser,1}.trafficStats.backLogPkt;
                queueBacklogsOverTime(iUserCount,iUser,iPkt,:) = SimStructs.userStruct{iUser,1}.trafficStats.backlogsOverTime;
            end
            
            if strcmp(SimParams.DebugMode,'true')
                display(squeeze(queueBacklogs(iUserCount,:,iPkt)));
            end
            
            cState = sprintf('Total Users completed - %d',SimParams.userCount(iUserCount));disp(cState);
        end
        
    end
    
    SimResults.avgTxPower = SimParams.txPower / SimParams.nDrops;
    
    if strcmp(SimParams.queueMode,'false')
        
        SimResults.sumThrpt = sum(SimParams.Thrpt(:,:,end),2);
        SimResults.thrptFairness = sum(SimParams.fairness(:,:,end),2);
        SimParams.sumThrpt = SimResults.sumThrpt;
        
        %         figure(1);hold all;
        %         plot(SimParams.snrIndex,SimParams.sumThrpt,'o-');
        %         xlabel('SNR in dB');ylabel('Sum Capacity in Bits/Sec/Hz');grid on;
        
        JainMean = mean(SimResults.sumThrpt,2).^2;JainVar = var(SimResults.sumThrpt,0,2);
        JainIndex_capacity = JainMean ./ (JainMean + JainVar);
        
        %     figure(2);hold all;
        %     plot(SimParams.snrIndex,JainIndex_capacity,markerS);
        %     xlabel('SNR in dB');ylabel('Capacity Deviation across Users in Bits/Sec/Hz');grid on;
        
        JainMean = mean(SimResults.thrptFairness,2).^2;JainVar = var(SimResults.thrptFairness,0,2);
        JainIndex_utility = JainMean ./ (JainMean + JainVar);
        
        %     figure(3);hold all;
        %     plot(SimParams.snrIndex,JainIndex_utility,markerS);
        %     xlabel('SNR in dB');ylabel('Utility Deviation across Users');grid on;
        
    else
        
        SimResults.queueBackLogs = queueBacklogs;
        SimResults.queueBackLogsOverTime = queueBacklogsOverTime;
        
        %     figure(5);hold all;
        %     plot(1:SimParams.nDrops,sum(squeeze(SimResults.queueBackLogsOverTime(end,:,end,:)),1));
        %     xlabel('Slot Index');ylabel('Queue Backlogs (pkts) over Time');grid on;
        
        %     plot(SimParams.maxArrival,sum(squeeze(SimResults.queueBackLogs(end,:,:)),1));
        %     xlabel('Average Arrival Rate');ylabel('Average Queue Size (pkts)');grid on;
        %     hold all;
        
    end
    
    if ~exist(SimParams.outFile,'file')
        gIndex = 1;
        gResults = cell(1,1);gParams = cell(1,1);gStructs = cell(1,1);
        gResults{gIndex,1} = SimResults;gParams{gIndex,1} = SimParams;gStructs{gIndex,1} = SimStructs;
        save(SimParams.outFile,'gResults','gParams','gStructs','gIndex');
    else
        load(SimParams.outFile);gIndex = gIndex + 1;
        gResults{gIndex,1} = SimResults;gParams{gIndex,1} = SimParams;gStructs{gIndex,1} = SimStructs;
        save(SimParams.outFile,'gResults','gParams','gStructs','gIndex');
    end
    
end


