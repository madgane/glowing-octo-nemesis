function displayOutputs(SimParams, SimStructs)

switch SimParams.plotMode
    
    case 'SRA'
        
        SimResults.sumThrpt = sum(SimParams.Thrpt(:,:,end),2);
        SimResults.thrptFairness = sum(SimParams.fairness(:,:,end),2);
        SimParams.sumThrpt = SimResults.sumThrpt;
        
        figStruct.N = 1;figStruct.P = 'plot';
        figStruct.X = SimParams.snrIndex;figStruct.Y = SimParams.sumThrpt;
        
        plotFigure(figStruct);
        xlabel('SNR in dB');ylabel('sum rate in bits/sec/Hz');
        
        JainMean = mean(SimParams.Thrpt,2).^2;JainVar = var(SimParams.Thrpt,0,2);
        JainIndex_capacity = JainMean ./ (JainMean + JainVar);
        
        figStruct.N = 2;figStruct.P = 'plot';
        figStruct.X = SimParams.snrIndex;figStruct.Y = JainIndex_capacity;

%         plotFigure(figStruct);
%         xlabel('SNR in dB');ylabel('Rate Deviation across Users in bits/sec/Hz');
        
        JainMean = mean(SimParams.fairness,2).^2;JainVar = var(SimParams.fairness,0,2);
        JainIndex_utility = JainMean ./ (JainMean + JainVar);
        
        figStruct.N = 2;figStruct.P = 'plot';
        figStruct.X = SimParams.snrIndex;figStruct.Y = JainIndex_utility;

%         plotFigure(figStruct);
%         xlabel('SNR in dB');ylabel('Network Utility Deviation across Users');
        
        
    case 'QA'

        figStruct.N = 1;figStruct.P = 'plot';
        figStruct.X = 1:SimParams.nDrops;figStruct.Y = sum(squeeze(SimParams.QueueInfo.queueBackLogsOverTime(end,:,end,:)),1);
        
        plotFigure(figStruct);
        xlabel('Slot Index');ylabel('Queue Backlogs (pkts) over Time');grid on;
        
        %         plotFigure(1:SimParams.nDrops,std(squeeze(SimParams.QueueInfo.queueBackLogsOverTime(end,:,end,:)),1),5,'plot');
        %         xlabel('Slot Index');ylabel('{\sigma_Q} Queue Backlogs (pkts) over Time');grid on;
        
        %         plotFigure(SimParams.maxArrival,sum(squeeze(SimParams.QueueInfo.queueBackLogs(end,:,:)),1),6,'plot');
        %         xlabel('Average Arrival Rate');ylabel('Average Queue Size (pkts)');grid on;
        
    case 'STA'
        
        nT = 1e3;nPRB = 50;nREinPRB = 120;nTot = nT * nPRB * nREinPRB * 1e-6;

        figStruct.N = 1;figStruct.P = 'cdfplot';
        figStruct.Y = SimParams.Thrpt(1,:,1) * nTot;
        
        plotFigure(figStruct);
        xlabel('Throughput in Mbps');
        ylabel('CDF of Throughput in Mbps');
        
    case 'NRA'
        
        plotFigure(struct('Y',SimParams.sumRateInstant));
        
    case 'QInfo'
        
        clc;
        
        displaySystemDetails;
        displayChannel(SimParams,SimStructs);
        
        for iDrop = SimParams.nDrops:SimParams.nDrops
            displayQueues(SimParams,SimStructs,iDrop);
        end
        
    case 'QTimePlot'
        
        plotFigure(struct('Y',sum(squeeze(SimParams.QueueInfo.queueResiduesOverTime(end,:,end,:)))));
        
    case 'CPlot'
        
        profile on;
        if ~isfield(SimParams,'distDecompSteps')
            SimParams.distDecompSteps = 1;
        end
        
        displaySystemDetails;
        displayChannel(SimParams,SimStructs);
        displayQueues(SimParams,SimStructs);
        
        totalDeviation = cell2mat(SimParams.Debug.tempResource{3,1});
        totalThroughput = cell2mat(SimParams.Debug.tempResource{2,1});
        
        totalDeviation = sum(totalDeviation);
        totalThroughput = sum(totalThroughput);
        
        plotFigure(struct('Y',totalDeviation(1:SimParams.distDecompSteps:end),'N',1));
        xlabel('SCA Update Points');
        ylabel('Queue deviation in bits / channel use');
        
        plotFigure(struct('Y',totalThroughput(1:SimParams.distDecompSteps:end),'N',2));
        xlabel('SCA Update Points');
        ylabel('Sum rate in bits / channel use');   
        
        profile off;
        
    case 'DispSchedUsers'
        
        fprintf(1,'\n');
        for iBase = 1:SimParams.nBases
            for iBand = 1:SimParams.nBands
                fprintf(1,'eNodeB - %d, Band - %d, Scheduled Users - %s \n',iBase,iBand,mat2str(SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1}));
            end
        end

    case 'DispMCInfo'
        
        displayQueues(SimParams,SimStructs);
        for iBand = 1:SimParams.nBands
            for iBase = 1:SimParams.nBases
                fprintf('txPower on SC [%d] from BS [%d] - %f \t',iBand,iBase,SimParams.txPower(end,end,iBase,iBand));
            end
            fprintf('\n');
        end
        
        
    otherwise
        
        display('Simulation Completed without any display !');
        
end

end