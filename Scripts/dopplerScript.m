
clc;clear all;
script_cells = cell(1,1);

% Ideal ....

xIndex = 1;
script_cells{xIndex,1}.mdpFactor = 0;
script_cells{xIndex,1}.fbFraction = 0.00;
script_cells{xIndex,1}.schedType = 'PFScheduling_SP';
script_cells{xIndex,1}.weightedSumRateMethod = 'StreamScheduling';

xIndex = xIndex + 1;
script_cells{xIndex,1}.mdpFactor = 0;
script_cells{xIndex,1}.fbFraction = 0.00;
script_cells{xIndex,1}.schedType = 'PFScheduling_BF';
script_cells{xIndex,1}.weightedSumRateMethod = 'StreamScheduling';

xIndex = xIndex + 1;
script_cells{xIndex,1}.mdpFactor = 0;
script_cells{xIndex,1}.fbFraction = 0.00;
script_cells{xIndex,1}.schedType = 'GreedyScheduling';
script_cells{xIndex,1}.weightedSumRateMethod = 'StreamScheduling';

% Feedback of 0.25 ....

xIndex = xIndex + 1;
script_cells{xIndex,1}.mdpFactor = 0;
script_cells{xIndex,1}.fbFraction = 0.25;
script_cells{xIndex,1}.schedType = 'PFScheduling_SP';
script_cells{xIndex,1}.weightedSumRateMethod = 'StreamScheduling';

xIndex = xIndex + 1;
script_cells{xIndex,1}.mdpFactor = 0;
script_cells{xIndex,1}.fbFraction = 0.25;
script_cells{xIndex,1}.schedType = 'PFScheduling_BF';
script_cells{xIndex,1}.weightedSumRateMethod = 'StreamScheduling';

xIndex = xIndex + 1;
script_cells{xIndex,1}.mdpFactor = 0;
script_cells{xIndex,1}.fbFraction = 0.25;
script_cells{xIndex,1}.schedType = 'GreedyScheduling';
script_cells{xIndex,1}.weightedSumRateMethod = 'StreamScheduling';

% MDP ....

xIndex = xIndex + 1;
script_cells{xIndex,1}.mdpFactor = 2;
script_cells{xIndex,1}.fbFraction = 0.25;
script_cells{xIndex,1}.schedType = 'PFScheduling_SP';
script_cells{xIndex,1}.weightedSumRateMethod = 'StreamScheduling';

xIndex = xIndex + 1;
script_cells{xIndex,1}.mdpFactor = 2;
script_cells{xIndex,1}.fbFraction = 0.25;
script_cells{xIndex,1}.schedType = 'PFScheduling_BF';
script_cells{xIndex,1}.weightedSumRateMethod = 'StreamScheduling';

xIndex = xIndex + 1;
script_cells{xIndex,1}.mdpFactor = 2;
script_cells{xIndex,1}.fbFraction = 0.25;
script_cells{xIndex,1}.schedType = 'GreedyScheduling';
script_cells{xIndex,1}.weightedSumRateMethod = 'StreamScheduling';

% MDP + Precoder ....

xIndex = xIndex + 1;
script_cells{xIndex,1}.mdpFactor = 2;
script_cells{xIndex,1}.fbFraction = 0.25;
script_cells{xIndex,1}.schedType = 'PFScheduling_SP';
script_cells{xIndex,1}.weightedSumRateMethod = 'PreScheduling';

xIndex = xIndex + 1;
script_cells{xIndex,1}.mdpFactor = 2;
script_cells{xIndex,1}.fbFraction = 0.25;
script_cells{xIndex,1}.schedType = 'PFScheduling_BF';
script_cells{xIndex,1}.weightedSumRateMethod = 'PreScheduling';

xIndex = xIndex + 1;
script_cells{xIndex,1}.mdpFactor = 2;
script_cells{xIndex,1}.fbFraction = 0.25;
script_cells{xIndex,1}.schedType = 'GreedyScheduling';
script_cells{xIndex,1}.weightedSumRateMethod = 'PreScheduling';

SimParams.version = version;
SimParams.outFile = 'outFile_doppler_3.mat';
pathAddition;

for iScript = 1:length(script_cells)
       
    SimParams.sysMode = 'false';
    SimParams.DebugMode = 'false';
    SimParams.queueMode = 'false';
    SimParams.precoderWithIdealChn = 'false';
    
    SimParams.ChannelModel = 'Jakes';
    SimParams.pathLossModel = 'Random_30';
    SimParams.DopplerType = 'Uniform_100';
    
    SimParams.queueWt = 1;
    SimParams.robustNoise = 1;
    
    SimParams.weighingEqual = 'true';
    SimParams.PrecodingMethod = 'Best_WMMSE_Method';
    
    SimParams.nDrops = 1000;
    SimParams.snrIndex = [-5:5:15];
    
    SimParams.PF_dur = 40;
    SimParams.SFSymbols = 14;
    SimParams.sampTime = 1e-3;
    SimParams.estError = 0.00;
    
    SimParams.SchedType = script_cells{iScript,1}.schedType;
    SimParams.mdpFactor = script_cells{iScript,1}.mdpFactor;
    SimParams.fbFraction = script_cells{iScript,1}.fbFraction;
    SimParams.weightedSumRateMethod = script_cells{iScript,1}.weightedSumRateMethod;
    
    chString = sprintf('Scheduling Type - %s, mdpFactor - %d, fbFraction - %d',...
                SimParams.SchedType,SimParams.mdpFactor,SimParams.fbFraction);
    display(chString);
            
    SimParams.nBands = 1;
    SimParams.nBases = 1;
    SimParams.nUsers = 20;
    
    SimParams.nTxAntenna = 4;
    SimParams.nRxAntenna = 1;
    
    SimParams.gracePeriod = 0;
    SimParams.arrivalDist = 'Constant';
    
    SimParams.maxArrival = 20;
    SimParams.FixedPacketArrivals = [10,10,10,10,10,10,1,1,1,1];
    SimParams.PL_Profile = [5 -inf 5 -inf 5 -inf 1e-20 0; -inf 5 -inf 5 -inf 5 0 1e-20];
    
    if strcmp(SimParams.sysMode,'true')
        SimParams.snrIndex = [0];
        SimParams.nBands = 1;
        SimParams.nBases = 57;
        SimParams.nUsers = 570;
    end
    
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
    SimParams.txPower = zeros(length(SimParams.maxArrival),length(SimParams.snrIndex),SimParams.nBases);
    
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
                
                if strcmp(SimParams.precoderWithIdealChn,'true')
                    SimStructs.linkChan = SimStructs.actualChannel;
                end
                
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
            
            cState = sprintf('SINR completed - %d',SimParams.snrIndex(iSNR));disp(cState);
        end
        
    end
    
    SimResults.avgTxPower = SimParams.txPower / SimParams.nDrops;
    
    if ~exist(SimParams.outFile,'file')
        gIndex = 1;
        gResults = cell(1,1);gParams = cell(1,1);gStructs = cell(1,1);
        gResults{gIndex,1} = SimResults;gParams{gIndex,1} = SimParams;gStructs{gIndex,1} = SimStructs;
        save(SimParams.outFile,'gResults','gParams','gStructs','gIndex');
    else
        load(char(SimParams.outFile));gIndex = gIndex + 1;
        gResults{gIndex,1} = SimResults;gParams{gIndex,1} = SimParams;gStructs{gIndex,1} = SimStructs;
        save(SimParams.outFile,'gResults','gParams','gStructs','gIndex');
    end

    continue;    
    
    if strcmp(SimParams.sysMode,'false')
        
        if strcmp(SimParams.queueMode,'false')
            
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
            
        else
            
            SimResults.queueBackLogs = queueBacklogs;
            SimResults.queueBackLogsOverTime = queueBacklogsOverTime;
            
            figure(5);
            hold all;
            plot(1:SimParams.nDrops,sum(squeeze(SimResults.queueBackLogsOverTime(end,:,end,:)),1));
            xlabel('Slot Index');ylabel('Queue Backlogs (pkts) over Time');grid on;
            
            figure(2);
            hold all;
            plot(1:SimParams.nDrops,std(squeeze(SimResults.queueBackLogsOverTime(end,:,end,:)),1));
            xlabel('Slot Index');ylabel('{\sigma_Q} Queue Backlogs (pkts) over Time');grid on;
            
            
            %     plot(SimParams.maxArrival,sum(squeeze(SimResults.queueBackLogs(end,:,:)),1));
            %     xlabel('Average Arrival Rate');ylabel('Average Queue Size (pkts)');grid on;
            %     hold all;
            
        end
        
        % if ~exist(SimParams.outFile,'file')
        %     gIndex = 1;
        %     gResults = cell(1,1);gParams = cell(1,1);gStructs = cell(1,1);
        %     gResults{gIndex,1} = SimResults;gParams{gIndex,1} = SimParams;gStructs{gIndex,1} = SimStructs;
        %     save(SimParams.outFile,'gResults','gParams','gStructs','gIndex');
        % else
        %     load(SimParams.outputFile);gIndex = gIndex + 1;
        %     gResults{gIndex,1} = SimResults;gParams{gIndex,1} = SimParams;gStructs{gIndex,1} = SimStructs;
        %     save(SimParams.outFile,'gResults','gParams','gStructs','gIndex');
        % end
        
    else
        
        nT = 1e3;nPRB = 50;nREinPRB = 120;nTot = nT * nPRB * nREinPRB * 1e-6;
        
        hold all;
        plotFigure(SimParams.Thrpt(1,:,1) * nTot,1,1,'cdfplot');
        xlabel('Throughput in Mbps');
        ylabel('CDF of Throughput in Mbps');
        
    end

end