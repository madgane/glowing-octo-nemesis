
clc;
clear all;close all;

SimParams.ExtRun = 'true';
SimParams.userRun = 'false';

SimParams.nBases = 1;
SimParams.nUsers = 50;
SimParams.nDrops = 1000;
SimParams.snrIndex = [-10:10:40];
SimParams.PL_Profile = -rand(SimParams.nBases,SimParams.nUsers) * 30;

%% RoundRobin Scheduling

display('RoundRobin Scheduling');
SimParams.SchedType = 'RRScheduling';
SimParams.PrecodingMethod = 'Best_MZF_Method';

baseCode;

markerS = '-d';LW = 1;Col = [0.0,0.0,0.9];
globalPlotHold(SimParams,markerS,LW,Col);

%% QR Scheduling

display('QR Scheduling');
SimParams.SchedType = 'BDScheduling';
SimParams.PrecodingMethod = 'Best_MZF_Method';
SimParams.caseStudy = 5;

baseCode;

markerS = '-v';LW = 1;Col = [0.0,0.8,0.0];
globalPlotHold(SimParams,markerS,LW,Col);

%% User Reduction Determinant Scheduling

display('User Reduction Determinant Scheduling');
SimParams.SchedType = 'PFBDScheduling';
SimParams.PrecodingMethod = 'Best_MZF_Method';
SimParams.caseStudy = 10;

baseCode;

markerS = '-p';LW = 1;Col = [0.9,0.0,0.0];
globalPlotHold(SimParams,markerS,LW,Col);

%% Projection Product Scheduling

display('Projection Product Scheduling');
SimParams.SchedType = 'BDScheduling';
SimParams.PrecodingMethod = 'Best_MZF_Method';
SimParams.caseStudy = -1;

baseCode;

markerS = '-o';LW = 1;Col = [0.4,0.9,0.2];
globalPlotHold(SimParams,markerS,LW,Col);

%% PF User Selction

display('PF User Scheduling');
SimParams.SchedType = 'PFScheduling';
SimParams.PrecodingMethod = 'Best_MZF_Method';
SimParams.caseStudy = 1;

baseCode;

markerS = '-.v';LW = 2;Col = [1.0,0.0,1.0];
globalPlotHold(SimParams,markerS,LW,Col);

%% Norm Based User Selction

display('Norm-based User Scheduling');
SimParams.SchedType = 'GreedyScheduling';
SimParams.PrecodingMethod = 'Best_MZF_Method';

baseCode;

markerS = '-.s';LW = 2;Col = [0.2,0.7,0.7];
globalPlotHold(SimParams,markerS,LW,Col);

%% Normalized QR Scheduling

display('Normalized QR Scheduling');
SimParams.SchedType = 'BDScheduling';
SimParams.PrecodingMethod = 'Best_MZF_Method';
SimParams.caseStudy = -3;

baseCode;

markerS = '-.*';LW = 2;Col = [0.0,0.9,0.9];
globalPlotHold(SimParams,markerS,LW,Col);

%% Reduced PF User Selction

display('PF User Scheduling');
SimParams.SchedType = 'PFScheduling';
SimParams.PrecodingMethod = 'Best_MZF_Method';
SimParams.caseStudy = 2;

baseCode;

markerS = '-.x';LW = 2;Col = [0.3,0.2,1.0];
globalPlotHold(SimParams,markerS,LW,Col);


%% Legend

legendString = cell(8,1);
legendString{1,1} = 'RoundRobin';legendString{2,1} = 'BD / QR based Selection';
legendString{3,1} = 'ReducedSelection';legendString{4,1} = 'ProdOfPerp';
legendString{5,1} = 'PFScheduling';legendString{6,1} = 'Norm-Based';
legendString{7,1} = 'NormalizedQR';legendString{8,1} = 'Reduced PF';
figure(1);legend(legendString);figure(2);legend(legendString);figure(3);legend(legendString);

