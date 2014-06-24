
clc;
close all;clear all;

SimParams.ExtRun = 'true';
SimParams.userRun = 'true';

SimParams.nDrops = 10;
snr_iterates = [5,35];

for iRun = 1:length(snr_iterates)
    
    SimParams.snrIndex = snr_iterates(1,iRun);
    
    %% RoundRobin Scheduling
    
    display('RoundRobin Scheduling');
    SimParams.SchedType = 'RRScheduling';
    SimParams.PrecodingMethod = 'Best_MZF_Method';
    
    baseCode_1;
    
    markerS = '-d';LW = 1;Col = [0.0,0.0,0.9];
    globalPlotHold(SimParams,markerS,LW,Col);
    
    %% QR Scheduling
    
    display('QR Scheduling');
    SimParams.SchedType = 'BDScheduling';
    SimParams.PrecodingMethod = 'Best_MZF_Method';
    SimParams.caseStudy = 5;
    
    baseCode_1;
    
    markerS = '-v';LW = 1;Col = [0.0,0.8,0.0];
    globalPlotHold(SimParams,markerS,LW,Col);
    
    %% User Reduction Determinant Scheduling
    
    display('User Reduction Determinant Scheduling');
    SimParams.SchedType = 'PFBDScheduling';
    SimParams.PrecodingMethod = 'Best_MZF_Method';
    SimParams.caseStudy = 10;
    
    baseCode_1;
    
    markerS = '-p';LW = 1;Col = [0.9,0.0,0.0];
    globalPlotHold(SimParams,markerS,LW,Col);
    
    %% Projection Product Scheduling
    
    display('Projection Product Scheduling');
    SimParams.SchedType = 'BDScheduling';
    SimParams.PrecodingMethod = 'Best_MZF_Method';
    SimParams.caseStudy = -1;
    
    baseCode_1;
    
    markerS = '-o';LW = 1;Col = [1.0,0.4,0.9];
    globalPlotHold(SimParams,markerS,LW,Col);
    
    %% Reduced PF User Selction

    display('PF User Scheduling');
    SimParams.SchedType = 'PFScheduling';
    SimParams.PrecodingMethod = 'Best_MZF_Method';
    SimParams.caseStudy = 2;

    baseCode_1;

    markerS = '-.*';LW = 2;Col = [0.0,0.9,0.9];
    globalPlotHold(SimParams,markerS,LW,Col);
    
    %% Pairwise User Scheduling
    
    display('Pairwise User Scheduling');
    SimParams.SchedType = 'PFBDScheduling';
    SimParams.PrecodingMethod = 'Best_MZF_Method';
    SimParams.caseStudy = 0;
    
    baseCode_1;
    
    markerS = '-.+';LW = 2;Col = [0.5,0.4,0.3];
    globalPlotHold(SimParams,markerS,LW,Col);
    
    %% Norm Based User Selction
    
    display('Norm-based User Scheduling');
    SimParams.SchedType = 'GreedyScheduling';
    SimParams.PrecodingMethod = 'Best_MZF_Method';
    
    baseCode_1;
    
    markerS = '-.s';LW = 2;Col = [0.6,0.2,0.7];
    globalPlotHold(SimParams,markerS,LW,Col);
    
    
    %% Legend
    
    legendString = cell(7,1);
    legendString{1,1} = 'RoundRobin';legendString{2,1} = 'BD / QR based Scheduling';
    legendString{3,1} = 'ReducedSelection';legendString{4,1} = 'ProdOfPerp';
    legendString{5,1} = 'ReducedPF';legendString{6,1} = 'Pair-Wise';
    legendString{7,1} = 'Norm-based';
    
    if iRun == 1
        figure(4);legend(legendString);
    end
    
end
