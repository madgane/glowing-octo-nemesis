

function [SimParams,SimStructs] = systemInitialize(SimParams,SimStructs)

for iUser = 1:SimParams.nUsers
    SimStructs.userStruct{iUser,1}.crThrpt = 1;
    SimStructs.userStruct{iUser,1}.PFmetric = 0;
    SimStructs.userStruct{iUser,1}.lastThrpt = 0;
    SimStructs.userStruct{iUser,1}.userID = iUser;
    SimStructs.userStruct{iUser,1}.tAllocation = 0;
    SimStructs.userStruct{iUser,1}.W = cell(SimParams.nBands,1);
    SimStructs.userStruct{iUser,1}.dropThrpt = zeros(SimParams.nDrops,1);
    
    SimStructs.userStruct{iUser,1}.trafficStats.backLogPkt = 0;
    SimStructs.userStruct{iUser,1}.trafficConfig.bufferLength = 'Inf';
    SimStructs.userStruct{iUser,1}.trafficHistory.pktService = zeros(1,(SimParams.nDrops + 1));
    SimStructs.userStruct{iUser,1}.trafficStats.backlogsOverTime = zeros(1,(SimParams.nDrops + 1));
    SimStructs.userStruct{iUser,1}.trafficStats.residuesOverTime = zeros(1,(SimParams.nDrops + 1));
    
end

for iBase = 1:SimParams.nBases
    SimStructs.baseStruct{iBase,1}.profile = 0;
    SimStructs.baseStruct{iBase,1}.baseID = iBase;
    SimStructs.baseStruct{iBase,1}.allocGains = cell(SimParams.nBands,1);
    SimStructs.baseStruct{iBase,1}.allocPattern = cell(SimParams.nBands,1);
    SimStructs.baseStruct{iBase,1}.P = cell(SimParams.nBands,1);
end

SimStructs.linkChan = cell(SimParams.nBases,SimParams.nBands);
SimStructs.actualChannel = cell(SimParams.nBases,SimParams.nBands);
SimStructs.prevChan = cell(SimParams.nBases,SimParams.nBands);

if strcmp(SimParams.ChannelModel,'Jakes')
    for iUser = 1:SimParams.nUsers
        for iBase = 1:SimParams.nBases
            for iBand = 1:SimParams.nBands
                reset(SimStructs.JakesChStruct{iUser,iBase,iBand});
            end
        end
    end
end

for iBase = 1:SimParams.nBases
    for iBand = 1:SimParams.nBands
        SimStructs.actualChannel{iBase,iBand} = zeros(SimParams.nRxAntenna,SimParams.nTxAntenna,SimParams.nUsers);
        SimStructs.linkChan{iBase,iBand} = zeros(SimParams.nRxAntenna,SimParams.nTxAntenna,SimParams.nUsers);
    end
end

% Traffic Modeling

[SimParams,SimStructs] = generateTraffic(SimParams,SimStructs);

SimParams.Debug.activeStatus = zeros(SimParams.nUsers,SimParams.nBands);
SimParams.Debug.receivedRSSI = zeros(SimParams.nRxAntenna,SimParams.nRxAntenna,SimParams.nUsers,SimParams.nBands);
SimParams.Debug.tempResource = cell(SimParams.maxDebugCells,SimParams.nDrops);

SimParams.Debug.tempResource{2,1} = cell(SimParams.nUsers,1);
SimParams.Debug.tempResource{3,1} = cell(SimParams.nUsers,1);
SimParams.Debug.tempResource{4,1} = cell(SimParams.nUsers,SimParams.nBands);

SimParams.currentQueue = 100;
SimParams.Debug.dataExchange = cell(6,1);

% 1 used for per-antenna-power-constraint work
% 2 used for queue-weighted precoding scheme
% 3 used for queue-deviation procedure
% 4 used for queue-weighted precoding scheme

end