% -------------------------------------------------------------------------
% SRA - sum-rate plot, QA - queue analysis, STA - system throughput
% analysis, NRA - network rate analysis
% -------------------------------------------------------------------------

clc;clear all;

saveContents = 'false';
SimParams.outFile = 'defaultFile';
SimParams.saveChannelInfo = 'false';
SimParams.channelSaveFolder = 'Results';

SimParams.maxDebugCells = 4;
SimParams.version = version;
SimParams.plotMode = 'CPlot';

prelimCheck;
preConfiguration;
SimParams.sysMode = 'false';
SimParams.DebugMode = 'false';
SimParams.precoderWithIdealChn = 'false';
SimParams.totalPwrDistOverSC = 'true';

SimParams.ChannelModel = 'IID';
SimParams.pathLossModel = 'Perturbed_3';
SimParams.DopplerType = 'Uniform_5';

SimParams.queueWt = 1;
SimParams.mdpFactor = 0;
SimParams.robustNoise = 0;

SimParams.weighingEqual = 'false';
SimParams.SchedType = 'SkipScheduling';
SimParams.PrecodingMethod = 'Best_QwtWSRM_Method';
SimParams.weightedSumRateMethod = 'GenAlloc';
SimParams.additionalParams = 'MMSE';

SimParams.nDrops = 1;
SimParams.snrIndex = [10];

SimParams.PF_dur = 40;
SimParams.SFSymbols = 14;
SimParams.sampTime = 1e-3;
SimParams.estError = 0.00;
SimParams.fbFraction = 0.00;

SimParams.nBands = 3;
SimParams.nBases = 3;
SimParams.nUsers = 9;

SimParams.nTxAntenna = 4;
SimParams.nRxAntenna = 2;
SimParams.ffrProfile_dB = zeros(1,SimParams.nBands);

SimParams.gracePeriod = 0;
SimParams.arrivalDist = 'Constant';

SimParams.maxArrival = 16;
SimParams.FixedPacketArrivals = [6];
SimParams.PL_Profile = [5 -inf 5 -inf 5 -inf 1e-20 0; -inf 5 -inf 5 -inf 5 0 1e-20];

if strcmp(SimParams.sysMode,'true')
    SimParams.snrIndex = [0];
    SimParams.nBands = 1;
    SimParams.nBases = 57;
    SimParams.nUsers = 570;
end

[SimParams,SimStructs] = initializeBuffers(SimParams);

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

            if strcmp(SimParams.DebugMode,'true')
                display(SimParams.Debug.activeStatus(:,1)');
            end
            
            [SimParams,SimStructs] = dropInitialize(SimParams,SimStructs);
            [SimParams,SimStructs] = getScheduledUsers(SimParams,SimStructs);
            
            if iDrop > 1
                displayQueues(SimParams,SimStructs,iDrop - 1);
            end
            
            if strcmp(SimParams.precoderWithIdealChn,'true')
                SimStructs.linkChan = SimStructs.actualChannel;
            end
            
            [SimParams,SimStructs] = getPMatrix(SimParams,SimStructs);
            [SimParams,SimStructs] = performReception(SimParams,SimStructs);
            
            for iUser = 1:SimParams.nUsers
                SimParams.sumRateInstant(iSNR,iDrop,iPkt) = SimParams.sumRateInstant(iSNR,iDrop,iPkt) + SimStructs.userStruct{iUser,1}.dropThrpt(iDrop,1);
            end
            
        end
        
        finalSystemUpdate;        
        if strcmp(SimParams.DebugMode,'true')
            display(squeeze(SimParams.QueueInfo.queueBacklogs(iSNR,:,iPkt)));
        end
        
        cState = sprintf('SINR completed - %d',SimParams.snrIndex(iSNR));disp(cState);
    end
    
end

SimResults.avgTxPower = SimParams.txPower / SimParams.nDrops;    
displayOutputs(SimParams,SimStructs);

if strcmp(saveContents,'true')
    
    cd Results;
    if exist(sprintf('%s.mat',SimParams.outFile),'file')
        load(SimParams.outFile);
        globalCount = globalCount + 1;
    else
        globalCount = 1;
        SimParamsCell = cell(1,1);
        SimStructsCell = cell(1,1);
    end
    
    SimParamsCell{globalCount,1} = SimParams;
    SimStructsCell{globalCount,1} = SimStructs;
    save(SimParams.outFile,'globalCount','SimParamsCell','SimStructsCell');
    cd ../
    
end

