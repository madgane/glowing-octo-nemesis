% -------------------------------------------------------------------------
% SRA - sum-rate plot, QA - queue analysis, STA - system throughput
% analysis, NRA - network rate analysis
% -------------------------------------------------------------------------

clc;clear all;

saveContents = 'false';
if strfind(saveContents,'true')
    if isunix
        addpath(genpath(pwd));
    else
        updatePath;
    end
end

%addpath(genpath(pwd));

runUserIndex = 0;
gUserIndices = 5:5:50;
gSimParams = cell(length(gUserIndices),1);
gSimStructs = cell(length(gUserIndices),1);

for runUserIndex = 1:length(gUserIndices)
    
    runUserCount = gUserIndices(1,runUserIndex);
    
    SimParams.outFile = 'ResultFileB';
    SimParams.saveChannelInfo = 'false';
    SimParams.channelSaveFolder = 'Results';
    
    SimParams.maxDebugCells = 4;
    SimParams.version = version;
    SimParams.plotMode = 'SRA';
    
    %prelimCheck;
    %preConfiguration;
    SimParams.sysMode = 'false';
    SimParams.DebugMode = 'false';
    SimParams.precoderWithIdealChn = 'false';
    SimParams.totalPwrDistOverSC = 'true';
    
    SimParams.ChannelModel = 'IID';
    SimParams.pathLossModel = 'CellEdge';
    SimParams.DopplerType = 'Uniform_25';
    
    SimParams.queueWt = 1;
    SimParams.BITFactor = 1;
    SimParams.mdpFactor = 0;
    SimParams.robustNoise = 0;
    
    SimParams.weighingEqual = 'true';
    SimParams.SchedType = 'BDScheduling_RNS';
    SimParams.PrecodingMethod = 'Best_ZF_Method';
    SimParams.weightedSumRateMethod = 'GenAlloc_1';
    SimParams.additionalParams = 'MMSE';
    
    SimParams.nExchangesOTA = 1;
    SimParams.exchangeResetInterval = 1;
    SimParams.nExchangesOBH = 1;
    
    SimParams.nDrops = 1000;
    SimParams.snrIndex = [10];
    
    SimParams.PF_dur = 40;
    SimParams.SFSymbols = 14;
    SimParams.sampTime = 1e-3;
    SimParams.estError = 0.00;
    SimParams.fbFraction = 0.00;
    SimParams.nSymbolsBIT = 100;
    
    SimParams.nBands = 1;
    SimParams.nBases = 1;
    SimParams.nUsers = runUserCount;
    
    SimParams.nTxAntenna = 4;
    SimParams.nRxAntenna = 4;
    SimParams.ffrProfile_dB = zeros(1,SimParams.nBands);
    
    display(SimParams);
    
    SimParams.maxArrival = 20;
    SimParams.groupArrivalFreq = 1;
    SimParams.arrivalDist = 'Uniform';
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
                SimParams.distIteration = iDrop;
                
                if strcmp(SimParams.DebugMode,'true')
                    display(SimParams.Debug.activeStatus(:,1)');
                end
                
                [SimParams,SimStructs] = dropInitialize(SimParams,SimStructs);
                [SimParams,SimStructs] = getScheduledUsers(SimParams,SimStructs);
                
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
    
    gSimParams{runUserIndex,1} = SimParams;
    gSimStructs{runUserIndex,1} = SimStructs;
end

sumThrpt = zeros(runUserIndex,1);
for iUserCount = 1:runUserIndex
    sumThrpt(iUserCount,1) = sum(gSimParams{iUserCount,1}.Thrpt(:,:,end),2);
end

%figStruct.N = 1;figStruct.P = 'plot';
%figStruct.X = gUserIndices;figStruct.Y = [sumThrpt(1:3);sumThrpt(5:end)];

%plotFigure(figStruct);
%xlabel('total Number of Users');ylabel('sum rate in bits/sec/Hz');

display(SimParams);
display(sumThrpt);

