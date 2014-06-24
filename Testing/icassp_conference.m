
clc;clear;

SimParams.capRun = 'false';
SimParams.ChannelModel = 'IID';
SimParams.pathLossModel = 'Random';
SimParams.weighingEqual = 'false';

SimParams.nDrops = 100;
SimParams.snrIndex = [25];

SimParams.PF_dur = 40;
SimParams.estError = 0.00;

SimParams.nBands = 1;
SimParams.nBases = 1;
SimParams.nUsers = 20;

SimParams.nTxAntenna = 4;
SimParams.nRxAntenna = 1;

SimParams.gracePeriod = 10;
SimParams.maxArrival = 0.5;
SimParams.arrivalDist = 'Constant';

SchedType = {'GreedyScheduling','BDScheduling_SP','PFScheduling_BF','PFScheduling_BF_M','Best_WMMSE_Method'};
PrecType = {'Best_WMMSE_Method','Best_WMMSE_Method','Best_WMMSE_Method','Best_WMMSE_Method','Best_WMMSE_Method','Best_WMMSE_Method'};
WSRMethod = {'StreamScheduling','StreamScheduling','StreamScheduling','StreamScheduling','StreamScheduling','PerformScheduling'};

nSchedSchemes = length(SchedType);
rxSimParams = cell(nSchedSchemes,1);
rxSimStructs = cell(nSchedSchemes,1);
rxSimResults = cell(nSchedSchemes,1);

for iSchedulingScheme = 1:length(SchedType)
    
    SimParams.SchedType = SchedType{1,iSchedulingScheme};
    SimParams.PrecodingMethod = PrecType{1,iSchedulingScheme};
    SimParams.weightedSumRateMethod = WSRMethod{1,iSchedulingScheme};
        
    display('Currently Running :');
    display(SimParams.SchedType);
    display(SimParams.PrecodingMethod);
    
    [rxSimParams{iSchedulingScheme,1},rxSimStructs{iSchedulingScheme,1},rxSimResults{iSchedulingScheme,1}] = globalCombinedScript(SimParams);
    
    display('---------------------');
    
end

save ICASSP_Results.mat rxSimParams rxSimStructs rxSimResults;

icassp_figures;

