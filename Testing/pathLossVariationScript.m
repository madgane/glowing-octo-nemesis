
clc;clear all;

SimParams.version = version;
SimParams.outFile = 'outFile_x1.mat';

pathAddition;
SimParams.sysMode = 'false';
SimParams.DebugMode = 'false';
SimParams.queueMode = 'false';

SimParams.ChannelModel = 'IID';
SimParams.pathLossModel = 'Fixed';
SimParams.DopplerType = 'Uniform_100';

SimParams.queueWt = 1;
SimParams.mdpFactor = 0;
SimParams.robustNoise = 0;

SimParams.weighingEqual = 'true';
SimParams.SchedType = 'XScheduling_StreamSearch';
SimParams.PrecodingMethod = 'Best_WMMSE_Method';
SimParams.weightedSumRateMethod = 'StreamScheduling';

SimParams.nDrops = 100;
SimParams.snrIndex = [0];

SimParams.PF_dur = 40;
SimParams.sampTime = 1e-3;
SimParams.estError = 0.00;
SimParams.fbFraction = 0.00;

SimParams.nBands = 1;
SimParams.nBases = 2;
SimParams.nUsers = 20;

SimParams.nTxAntenna = 4;
SimParams.nRxAntenna = 1;

SimParams.gracePeriod = 0;
SimParams.arrivalDist = 'Constant';

SimParams.pathLossProfile = [0:5:20];
SimParams.uIndices = 1:SimParams.nBases:SimParams.nUsers;

SimParams.maxArrival = 10;
SimParams.FixedPacketArrivals = [10,10,10,10,10,10,1,1,1,1];

if strcmp(SimParams.sysMode,'true')
    SimParams.snrIndex = [0];
    SimParams.nBands = 1;
    SimParams.nBases = 57;
    SimParams.nUsers = 570;
end

nSINRSamples = length(SimParams.snrIndex);
nPacketSamples = length(SimParams.maxArrival);
nPathLossSamples = length(SimParams.pathLossProfile);

SimStructs.userStruct = cell(SimParams.nUsers,1);
SimStructs.baseStruct = cell(SimParams.nBases,1);

SimParams.maxRank = min(SimParams.nRxAntenna,SimParams.nTxAntenna);
SimParams.muxRank = min(SimParams.nTxAntenna,(SimParams.nRxAntenna * SimParams.nUsers));

SimParams.Thrpt = zeros(nSINRSamples,SimParams.nUsers,nPacketSamples,nPathLossSamples);
utilityScale = SimParams.nDrops * SimParams.muxRank * SimParams.nBands;
SimParams.fairness = zeros(nSINRSamples,SimParams.nUsers,nPacketSamples);

queueBacklogs = zeros(nSINRSamples,SimParams.nUsers,nPacketSamples);
queueBacklogsOverTime = zeros(nSINRSamples,SimParams.nUsers,nPacketSamples,SimParams.nDrops);
SimParams.txPower = zeros(length(SimParams.maxArrival),length(SimParams.snrIndex),SimParams.nBases);

for iPathLoss = 1:length(SimParams.pathLossProfile)
    
    SimParams.PL_Profile = zeros(SimParams.nBases,SimParams.nUsers);
    for iBase = 1:SimParams.nBases
        iBaseUsers = iBase:SimParams.nBases:SimParams.nUsers;
        
        bsSINR_Profile = ones(1,SimParams.nUsers) * 0;
        bsSINR_Profile(1,iBaseUsers) = ones(1,length(SimParams.uIndices)) * SimParams.pathLossProfile(1,iPathLoss);
        bsSINR_Profile(1,iBaseUsers) = bsSINR_Profile(1,iBaseUsers) + 1e-10;

        SimParams.PL_Profile(iBase,:) = bsSINR_Profile;
    end
   
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
                SimParams.Thrpt(iSNR,iUser,iPkt,iPathLoss) = (SimStructs.userStruct{iUser}.crThrpt - 1) / (SimParams.nDrops * SimParams.nBands);
                queueBacklogs(iSNR,iUser,iPkt) = SimStructs.userStruct{iUser,1}.trafficStats.backLogPkt;
                queueBacklogsOverTime(iSNR,iUser,iPkt,:) = SimStructs.userStruct{iUser,1}.trafficStats.backlogsOverTime;
            end
            
            if strcmp(SimParams.DebugMode,'true')
                display(squeeze(queueBacklogs(iSNR,:,iPkt)));
            end
            
            cState = sprintf('SINR completed - %d, PL Exponent Completed - %d',SimParams.snrIndex(iSNR),SimParams.pathLossProfile(1,iPathLoss));
            disp(cState);
        end
        
    end
    
end

pathLossThrpt = squeeze(SimParams.Thrpt(end,:,end,:));
SimParams.sumThrpt = sum(pathLossThrpt);

hold all;
plot(SimParams.pathLossProfile,SimParams.sumThrpt);
xlabel('pathloss variation of the first cell users');
ylabel('sum rate in bits/sec/Hz')



