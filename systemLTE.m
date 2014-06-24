
% -------------------------------------------------------------------------
% SRA - sum-rate plot, QA - queue analysis, STA - system throughput
% analysis, NRA - network rate analysis, SEA - spectral efficiency analysis
% -------------------------------------------------------------------------

clc;clear all;

SimParams.maxDebugCells = 4;
SimParams.version = version;
SimParams.outFile = 'outFile_x1.mat';
SimParams.plotMode = 'SEA';

preConfiguration;
SimParams.DebugMode = 'false';
SimParams.precoderWithIdealChn = 'false';

SimParams.ChannelModel = 'Jakes';
SimParams.pathLossModel = '3GPP_RMa';

SimParams.queueWt = 1;
SimParams.mdpFactor = 2;
SimParams.robustNoise = 0;

SimParams.weighingEqual = 'true';
SimParams.SchedType = 'SkipScheduling';
SimParams.PrecodingMethod = 'Best_QwtWSRM_Method';
SimParams.weightedSumRateMethod = 'JointAlloc';

SimParams.nDrops = 1;
SimParams.snrIndex = [0];

SimParams.PF_dur = 40;
SimParams.SFSymbols = 14;
SimParams.sampTime = 1e-3;
SimParams.estError = 0.00;
SimParams.fbFraction = 0.00;

SimParams.nBands = 1;
SimParams.nTiers = 0;
SimParams.nSectors = 3;
SimParams.nNeighbors = 2; % Number of neighbors to realize
SimParams.perCiteUsers = 4;

SimParams.nTxAntenna = 4;
SimParams.nRxAntenna = 1;
SimParams.ffrProfile_dB = zeros(1,SimParams.nBands);

SimParams.nBases = getCellsOverLayout(SimParams.nTiers,SimParams.nSectors);
SimParams.nUsers = SimParams.nBases * SimParams.perCiteUsers;

SimParams.gracePeriod = 0;
SimParams.arrivalDist = 'Constant_10';

SimParams.maxArrival = 20;
SimParams.FixedPacketArrivals = [10,10,10,10,10,10,1,1,1,1];

bufferInitializations;
if strcmp(SimParams.DebugMode,'true')
    keyboard;
end

for iPkt = 1:length(SimParams.maxArrival)
    
    SimParams.iPkt = iPkt;
    [SimParams,SimStructs] = fwkInitialization(SimParams,SimStructs);
    SimParams.N = 10^(SimParams.systemNoise / 10);
    
    for iSNR = 1:length(SimParams.snrIndex)
        
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
            
            if strcmp(SimParams.precoderWithIdealChn,'true')
                SimStructs.linkChan = SimStructs.actualChannel;
            end
            
            [SimParams,SimStructs] = getPMatrix(SimParams,SimStructs);
            [SimParams,SimStructs] = performReception(SimParams,SimStructs);
            
            for iUser = 1:SimParams.nUsers
                sumRateInstant(iSNR,iDrop,iPkt) = sumRateInstant(iSNR,iDrop,iPkt) + SimStructs.userStruct{iUser,1}.dropThrpt(iDrop,1);
            end
            
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
        
        cState = sprintf('SINR completed - %d',SimParams.snrIndex(iSNR));disp(cState);
    end
    
end

SimResults.avgTxPower = SimParams.txPower / SimParams.nDrops;

switch SimParams.plotMode
    
    case 'SRA'
        
        SimResults.sumThrpt = sum(SimParams.Thrpt(:,:,end),2);
        SimResults.thrptFairness = sum(SimParams.fairness(:,:,end),2);
        SimParams.sumThrpt = SimResults.sumThrpt;
        
        plotFigure(SimParams.snrIndex,SimParams.sumThrpt,1,'plot');
        xlabel('SNR in dB');ylabel('sum rate in bits/sec/Hz');
        
        JainMean = mean(SimParams.Thrpt,2).^2;JainVar = var(SimParams.Thrpt,0,2);
        JainIndex_capacity = JainMean ./ (JainMean + JainVar);
        
        plotFigure(SimParams.snrIndex,JainIndex_capacity,2,'plot');
        xlabel('SNR in dB');ylabel('Rate Deviation across Users in bits/sec/Hz');
        
        JainMean = mean(SimParams.fairness,2).^2;JainVar = var(SimParams.fairness,0,2);
        JainIndex_utility = JainMean ./ (JainMean + JainVar);
        
        plotFigure(SimParams.snrIndex,JainIndex_utility,3,'plot');
        xlabel('SNR in dB');ylabel('Network Utility Deviation across Users');
        
        
    case 'QA'
                        
        SimResults.queueBackLogs = queueBacklogs;
        SimResults.queueBackLogsOverTime = queueBacklogsOverTime;
        
        plotFigure(1:SimParams.nDrops,sum(squeeze(SimResults.queueBackLogsOverTime(end,:,end,:)),1),4,'plot');
        xlabel('Slot Index');ylabel('Queue Backlogs (pkts) over Time');grid on;
        
        plotFigure(1:SimParams.nDrops,std(squeeze(SimResults.queueBackLogsOverTime(end,:,end,:)),1),5,'plot');
        xlabel('Slot Index');ylabel('{\sigma_Q} Queue Backlogs (pkts) over Time');grid on;
        
        plotFigure(SimParams.maxArrival,sum(squeeze(SimResults.queueBackLogs(end,:,:)),1),6,'plot');
        xlabel('Average Arrival Rate');ylabel('Average Queue Size (pkts)');grid on;

        
    case 'STA'
        
        nT = 1e3;nPRB = 50;nREinPRB = 120;nTot = nT * nPRB * nREinPRB * 1e-6;
        
        hold all;
        plotFigure(SimParams.Thrpt(1,:,1) * nTot,1,1,'cdfplot');
        xlabel('Throughput in Mbps');
        ylabel('CDF of Throughput in Mbps');
        
    case 'NRA'
        
        plotFigure(1:SimParams.nDrops,sumRateInstant,1,'plot');
        
    case 'SEA'
        
        hold all;
        usableFraction = (120 / 168);
        actThrpt = SimParams.Thrpt(1,:,1) * usableFraction;
        plotFigure(actThrpt,1,1,'cdfplot');
        xlabel('Sum Rate in bits / channel-use');
        ylabel('CDF of Sum Rate');
        
        SimParams.Reports.cellSpectralEfficiency = sum(SimParams.Thrpt) / (SimParams.nBases * SimParams.sysConfig.subChnlBWHz);
        SimParams.Reports.cellSpectralEfficiency = usableFraction * SimParams.Reports.cellSpectralEfficiency / (SimParams.sampTime / SimParams.SFSymbols);
        
    otherwise
        
        display('Unknown print options !');        

        
end
