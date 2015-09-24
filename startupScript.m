% -------------------------------------------------------------------------
% SRA - sum-rate plot, QA - queue analysis, STA - system throughput
% analysis, NRA - network rate analysis
% -------------------------------------------------------------------------

clc;
clearvars;
clearvars -global;
clear workspace;

saveContents = 'false';
if strfind(saveContents,'true')
    updatePath;
end

SimParams.outFile = 'JournalData-1';
SimParams.saveChannelInfo = 'false';
SimParams.channelSaveFolder = 'Results';

SimParams.maxDebugCells = 5;
SimParams.version = version;
SimParams.plotMode = 'DispMCInfo';

SimParams.sysMode = 'false';
SimParams.cvxDisabled = 'true';
SimParams.multiCasting = 'true';

prelimCheck;
preConfiguration;
SimParams.sysMode = 'false';
SimParams.DebugMode = 'false';
SimParams.precoderWithIdealChn = 'false';
SimParams.totalPwrDistOverSC = 'true';

SimParams.ChannelModel = 'IID';
SimParams.pathLossModel = 'CellEdge';
SimParams.DopplerType = 'Constant_140';

SimParams.queueWt = 1;
SimParams.BITFactor = 1;
SimParams.mdpFactor = 0;
SimParams.robustNoise = 0;

SimParams.weighingEqual = 'true';
SimParams.SchedType = 'SkipScheduling';
SimParams.PrecodingMethod = 'Best_MultiCastBF_Method';
SimParams.DesignType = 'ConicMethodB';

SimParams.nExchangesOTA = 50;
SimParams.exchangeResetInterval = 1;
SimParams.nExchangesOBH = 1;

SimParams.maxArrival = [4];
SimParams.arrivalDist = 'SteadyFlow';

SimParams.PF_dur = 40;
SimParams.SFSymbols = 14;
SimParams.sampTime = 1e-3;
SimParams.estError = 0.00;
SimParams.fbFraction = 0.00;
SimParams.nSymbolsBIT = 1e100;

SimParams.nBands = 3;
SimParams.nBases = 1;
SimParams.nDrops = 1;
SimParams.snrIndex = [10];

SimParams.nTxAntenna = 4;
SimParams.nRxAntenna = 1;
SimParams.ffrProfile_dB = zeros(1,SimParams.nBands);

SimParams.maxArrival = 5;
SimParams.groupArrivalFreq = 1;
SimParams.arrivalDist = 'Constant';
SimParams.FixedPacketArrivals = [6];
SimParams.PL_Profile = [5 -inf 5 -inf 5 -inf 1e-20 0; -inf 5 -inf 5 -inf 5 0 1e-20];

if strcmp(SimParams.multiCasting,'true')
    SimParams.nGroupArray = 2;
    SimParams.usersPerGroup = 10;
    SimParams.nAntennaArray = 8;
    SimParams.nTxAntennaEnabledArray = [8:SimParams.nAntennaArray];
    
    SimParams.mcGroups = cell(SimParams.nBases,1);
    SimParams.totalTXpower_G = zeros(length(SimParams.maxArrival),length(SimParams.nTxAntennaEnabledArray));
    SimParams.totalTXpower_SDP = zeros(length(SimParams.maxArrival),length(SimParams.nTxAntennaEnabledArray));
    SimParams.solverTiming = zeros(length(SimParams.maxArrival),length(SimParams.nTxAntennaEnabledArray));
end

gXParams = cell(length(SimParams.nTxAntennaEnabledArray),1);
gXStructs = cell(length(SimParams.nTxAntennaEnabledArray),1);

for iAntennaArray = 1:length(SimParams.nTxAntennaEnabledArray)
    
    SimParams.iAntennaArray = iAntennaArray;
    nGroupsPerCell = SimParams.nGroupArray;
    for iBase = 1:SimParams.nBases
        SimParams.mcGroups{iBase,1} = [];
        for iGroup = 1:nGroupsPerCell
            SimParams.mcGroups{iBase,1} = [SimParams.mcGroups{iBase,1}, SimParams.usersPerGroup];
        end
    end
    
    SimParams.nRxAntenna = 1;
    SimParams.nTxAntenna = SimParams.nAntennaArray;
    SimParams.nUsers = SimParams.usersPerGroup * SimParams.nBases * nGroupsPerCell;
    SimParams.nTxAntennaEnabled = SimParams.nTxAntennaEnabledArray(1,SimParams.iAntennaArray);
    
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
                SimParams.distIteration = iDrop;
                
                if strcmp(SimParams.DebugMode,'true')
                    SimParams.Debug.activeStatus(:,1);
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
                if strcmp(SimParams.multiCasting,'true')
                    [SimParams,SimStructs] = performGroupReception(SimParams,SimStructs);
                else
                    [SimParams,SimStructs] = performReception(SimParams,SimStructs);
                end
                
                for iUser = 1:SimParams.nUsers
                    SimParams.sumRateInstant(iSNR,iDrop,iPkt) = SimParams.sumRateInstant(iSNR,iDrop,iPkt) + SimStructs.userStruct{iUser,1}.dropThrpt(iDrop,1);
                end
                
                [SimParams,SimStructs] = updateTransmitPower(SimParams,SimStructs);
                
            end
            
            finalSystemUpdate;
            if strcmp(SimParams.DebugMode,'true')
                display(squeeze(SimParams.QueueInfo.queueBacklogs(iSNR,:,iPkt)));
            end
            
        end
    end
    
    cState = sprintf('Solved for - [%d] - Active Transmit Antenna Elements ',SimParams.nTxAntennaEnabled);disp(cState);
    
    gXParams{iAntennaArray,1} = SimParams;
    gXStructs{iAntennaArray,1} = SimStructs;
    
end

SimResults.avgTxPower = SimParams.txPower / SimParams.nDrops;
displayOutputs(gXParams,gXStructs);

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

