
clc;clear all;close all;
load('Results\defaultOutFile_x42.mat');

maxIterations = 700;
allocatedBuffer = 2;

linePalette = {'-.','-',':','-'};
markerPalette = {'*','o','.','+','d','s','p'};
colorPalette = {'b','g','r','m','c','k','y'};

for iBase = 1:SimParamsCell{1,1}.nBases
    figure(iBase);hold on;
end

figLineWidth = 1;
figLineType = linePalette{1,3};
for iBase = 1:SimParamsCell{1,1}.nBases
    figure(iBase);
    linkedUsers = SimStructsCell{1,1}.baseStruct{iBase,1}.linkedUsers;
    for iUser = 1:length(linkedUsers)
        cUser = linkedUsers(iUser,1);
        figColor = colorPalette{1,iUser};    
        queuedPkts = ones(maxIterations,1) * SimStructsCell{1,1}.userStruct{cUser,1}.trafficHistory.pktArrival(1,1);
        plot(queuedPkts,'Color',figColor,'LineWidth',figLineWidth,...
            'LineStyle',figLineType);
    end
end

figure(1);box on;
xlabel('Iteration index');ylabel('Rate achieved by users independently (in bits)');
legend('User - {1} Queue','User - {2} Queue','User - {3} Queue','User - {4} Queue');

figure(2);box on;
xlabel('Iteration index');ylabel('Rate achieved by users independently (in bits)');
legend('User - {5} Queue','User - {6} Queue','User - {7} Queue','User - {8} Queue');
for iCount = 1:globalCount
    
    SimParams = SimParamsCell{iCount,1};
    SimStructs = SimStructsCell{iCount,1};
    displaySystemDetails;displayChannel;displayQueues;
    
    if ~strcmp(SimParams.weightedSumRateMethod,'BandAlloc')
        figLineWidth = 2;
        figLineType = linePalette{1,iCount};
        maxIterations = max(length(SimParams.Debug.tempResource{allocatedBuffer,1}{1,1}),maxIterations);
        for iBase = 1:SimParams.nBases
            figure(iBase);
            linkedUsers = SimStructs.baseStruct{iBase,1}.linkedUsers;
            for iUser = 1:length(linkedUsers)
                cUser = linkedUsers(iUser,1);
                figColor = colorPalette{1,iUser};
                achRate = SimParams.Debug.tempResource{allocatedBuffer,1}{cUser,1};
                plot(achRate,'Color',figColor,'LineWidth',figLineWidth,...
                    'LineStyle',figLineType);
            end
        end
    end
end

figLineWidth = 2;
figLineType = linePalette{1,4};
figure(SimParamsCell{1,1}.nBases + 1);hold on;

for iCount = 1:globalCount

    figColor = colorPalette{1,iCount};
    SimParams = SimParamsCell{iCount,1};
    SimStructs = SimStructsCell{iCount,1};
    
    QdevMatrix = cell2mat(SimParams.Debug.tempResource{3,1});
    plot(sum(QdevMatrix),'Color',figColor,'LineWidth',figLineWidth,...
        'LineStyle',figLineType);
    
end

box on;
xlabel('Iteration index');ylabel('Queue deviation (in bits)');
