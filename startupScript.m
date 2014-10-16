% -------------------------------------------------------------------------
% SRA - sum-rate plot, QA - queue analysis, STA - system throughput
% analysis, NRA - network rate analysis
% -------------------------------------------------------------------------

clc;clear all;

saveContents = 'false';
if strfind(saveContents,'true')
    addpath(genpath(pwd));
end

SimParams.outFile = 'systemFile';
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
SimParams.pathLossModel = 'Isolated';
SimParams.DopplerType = 'Constant_25';

SimParams.queueWt = 1;
SimParams.mdpFactor = 0;
SimParams.robustNoise = 0;

SimParams.weighingEqual = 'false';
SimParams.SchedType = 'SkipScheduling';
if 1
SimParams.PrecodingMethod = 'Best_QwtWSRMD_Method';
SimParams.weightedSumRateMethod = 'MSEKKTMethod';
else
SimParams.PrecodingMethod = 'Best_RTQWSRM_Method';
SimParams.weightedSumRateMethod = 'distMSEAlloc';
end
SimParams.additionalParams = 'MMSE';

SimParams.nExchangesOTA = 20;
SimParams.exchangeResetInterval = 10;
SimParams.nExchangesOBH = 1;

SimParams.nDrops = 1;
SimParams.snrIndex = [10];

SimParams.PF_dur = 40;
SimParams.SFSymbols = 14;
SimParams.sampTime = 1e-3;
SimParams.estError = 0.00;
SimParams.fbFraction = 0.00;

SimParams.nBands = 3;
SimParams.nBases = 2;
SimParams.nUsers = 10;

SimParams.nTxAntenna = 4;
SimParams.nRxAntenna = 2;
SimParams.ffrProfile_dB = zeros(1,SimParams.nBands);

SimParams.gracePeriod = 0;
SimParams.arrivalDist = 'SteadyFlow';

SimParams.maxArrival = 5;
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
            SimParams.distIteration = 1;
            
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

    xVal = fix(clock);
    SimParams.Log.date = date;
    SimParams.Log.clock = sprintf('%d:%d:%d',xVal(1,4),xVal(1,5),xVal(1,6));
    
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

