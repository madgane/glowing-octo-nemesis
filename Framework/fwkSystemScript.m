
function fwkSystemScript(SimParams)

[SimParams,SimStructs] = initializeBuffers(SimParams);

SimParams.processID = feature('GetPid');
SimParams.Log.Date = cell(1,length(SimParams.maxArrival));
SimParams.Log.Clock.S = cell(1,length(SimParams.maxArrival));
SimParams.Log.Clock.E = cell(1,length(SimParams.maxArrival));


for iPkt = 1:length(SimParams.maxArrival)
    
    xValA = fix(clock);
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
            SimParams.distIteration = iDrop;
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
            
            finalSystemUpdate;
            xValB = fix(clock);
            if strcmp(SimParams.saveContents,'true')
                SimParams.Log.Date{1,iPkt} = date;
                SimParams.Log.Clock.S{1,iPkt} = sprintf('%d:%d:%d',xValA(1,4),xValA(1,5),xValA(1,6));
                SimParams.Log.Clock.E{1,iPkt} = sprintf('%d:%d:%d',xValB(1,4),xValB(1,5),xValB(1,6));
                outFile = sprintf('%s/%s.mat',SimParams.channelSaveFolder,SimParams.outFile);
                save(outFile,'SimParams','SimStructs');
            end
        end
        
        finalSystemUpdate;
        if strcmp(SimParams.DebugMode,'true')
            display(squeeze(SimParams.QueueInfo.queueBacklogs(iSNR,:,iPkt)));
        end
        
        cState = sprintf('(PKT,SINR,JobID) - (%d,%d,%d)',SimParams.maxArrival(1,iPkt),SimParams.snrIndex(1,iSNR),SimParams.xCount);
        disp(cState);
        
    end
    
    xValB = fix(clock);
    if strcmp(SimParams.saveContents,'true')
        SimParams.Log.Date{1,iPkt} = date;
        SimParams.Log.Clock.S{1,iPkt} = sprintf('%d:%d:%d',xValA(1,4),xValA(1,5),xValA(1,6));
        SimParams.Log.Clock.E{1,iPkt} = sprintf('%d:%d:%d',xValB(1,4),xValB(1,5),xValB(1,6));
        outFile = sprintf('%s/%s.mat',SimParams.channelSaveFolder,SimParams.outFile);
        save(outFile,'SimParams','SimStructs');
    end
    
end

end