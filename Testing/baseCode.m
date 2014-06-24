
SimParams.ExtRun = 'false';

if ~isequal(SimParams.ExtRun,'true')
    clc;clear;
    SimParams.ExtRun = 'false';
end

SimParams.ChannelModel = 'IID';

if ~isequal(SimParams.ExtRun,'true')
    SimParams.SchedType = 'PFBDScheduling_AHP';
    SimParams.PrecodingMethod = 'Best_ZF_Method';
end

SimParams.nBands = 1;
SimParams.estError = 0.00;

if ~isequal(SimParams.ExtRun,'true')
    SimParams.nBases = 1;
    SimParams.nUsers = 20;
end

SimParams.PF_dur = 40;
SimParams.maxRank = 4;
SimParams.nTxAntenna = 4;
SimParams.nRxAntenna = 1;
SimParams.muxRank = min(SimParams.nTxAntenna,(SimParams.nRxAntenna * SimParams.nUsers));

if ~isequal(SimParams.ExtRun,'true')
    SimParams.nDrops = 100;
    SimParams.snrIndex = [-10:10:40];    
    SimParams.PL_Profile = -rand(SimParams.nBases,SimParams.nUsers) * 0;
end

SimParams.Thrpt = zeros(length(SimParams.snrIndex),SimParams.nUsers);
utilityScale = SimParams.nDrops * SimParams.muxRank * SimParams.nBands;
SimParams.fairness = zeros(length(SimParams.snrIndex),SimParams.nUsers);
SimStructs.userStruct = cell(SimParams.nUsers,1);SimStructs.baseStruct = cell(SimParams.nBases,1);

profile on

for iSNR = 1:length(SimParams.snrIndex)
    
    SimParams.iSNR = iSNR;
    SimParams.sPower = 10.^(SimParams.snrIndex(iSNR)/10);
    [SimParams,SimStructs] = systemInitialize(SimParams,SimStructs);
    [SimParams,SimStructs] = systemLinking(SimParams,SimStructs);
    
    for iDrop = 1:SimParams.nDrops
        SimParams.iDrop = iDrop;
        [SimParams,SimStructs] = dropInitialize(SimParams,SimStructs);
        [SimParams,SimStructs] = getScheduledUsers(SimParams,SimStructs);
        [SimParams,SimStructs] = getPMatrix(SimParams,SimStructs);
        [SimParams,SimStructs] = performReception(SimParams,SimStructs);
    end
    
    for iUser = 1:SimParams.nUsers
        SimParams.PFmetric(iSNR,iUser) = SimStructs.userStruct{iUser}.PFmetric;
        SimParams.fairness(iSNR,iUser) = SimStructs.userStruct{iUser}.tAllocation / utilityScale;
        SimParams.Thrpt(iSNR,iUser) = (SimStructs.userStruct{iUser}.crThrpt - 1) / SimParams.nDrops;
    end
    
    cState = sprintf('SINR completed - %d',SimParams.snrIndex(iSNR));disp(cState);
end

profile off

% profile viewer

SimParams.profiler.schX = SimParams.profiler.schX / (iSNR * iDrop);

if ~isequal(SimParams.ExtRun,'true')

    markerS = 'o-';

    figure(1);hold all;
    SimParams.sumThrpt = sum(SimParams.Thrpt,2);
    plot(SimParams.snrIndex,SimParams.sumThrpt,markerS);
    xlabel('SNR in dB');ylabel('Sum Capacity in Bits/Sec/Hz');grid on;
 
%     figure(2);hold all;
%     JainMean = mean(SimParams.Thrpt,2).^2;JainVar = var(SimParams.Thrpt,0,2);
%     JainIndex_capacity = JainMean ./ (JainMean + JainVar);

%     plot(SimParams.snrIndex,JainIndex_capacity,markerS);
%     xlabel('SNR in dB');ylabel('Capacity Deviation across Users in Bits/Sec/Hz');grid on;

%     figure(3);hold all;
%     JainMean = mean(SimParams.fairness,2).^2;JainVar = var(SimParams.fairness,0,2);
%     JainIndex_utility = JainMean ./ (JainMean + JainVar);

%     plot(SimParams.snrIndex,JainIndex_utility,markerS);
%     xlabel('SNR in dB');ylabel('Network Utility Deviation across Users');grid on;

end
