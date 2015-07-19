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
        
        figStruct.N = 1;figStruct.P = 'cdfplot';
        figStruct.Y = SimParams.Thrpt(1,:,1);
        
        plotFigure(figStruct);
        xlabel('Spectral Efficiency in bits/RE');
        ylabel('CDF of Spectral Efficiency');
        
        figStruct.N = 2;figStruct.P = 'cdfplot';
        figStruct.Y = squeeze(SimParams.QueueInfo.queueBacklogs(1,:,1));
        
        plotFigure(figStruct);
        xlabel('Number of Residual Packets in bits');
        ylabel('CDF of Residual Packets');
        
        
    case 'NRA'
        
        plotFigure(struct('Y',SimParams.sumRateInstant));
        
    case 'QInfo'
        
        displaySystemDetails;
        displayChannel(SimParams,SimStructs);
        
        for iDrop = SimParams.nDrops:SimParams.nDrops
            displayQueues(SimParams,SimStructs,iDrop);
        end
        
    case 'QTimePlot'
        
        figStruct.N = 1;
        figStruct.Y = sum(squeeze(SimParams.QueueInfo.queueResiduesOverTime(end,:,end,:)));
        plotFigure(figStruct);
        xlabel('Time instant (T)');
        ylabel('Total number of residual packets {\chi} in bits/channel use');
        
        figStruct.N = 2;
        figStruct.Y = sum(squeeze(SimParams.QueueInfo.packetServiceOverTime(end,:,end,:)));
        plotFigure(figStruct); 
        xlabel('Time instant (T)');
        ylabel('Total number of transmitted packets in bits/channel use');

        figStruct.N = 3;
        figStruct.Y = sum(squeeze(SimParams.QueueInfo.queueBacklogsOverTime(end,:,end,:)));
        plotFigure(figStruct);
        xlabel('Time instant (T)');
        ylabel('Total number of backlogged packets {\chi} in bits/channel use');
        
    case 'CPlot'
        
        profile on;
        if ~isfield(SimParams,'distDecompSteps')
            SimParams.distDecompSteps = 1;%SimParams.nExchangesOBH;
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
           
        %displayChannels;
        
    case 'QSurfacePlot'
        
        ylabel('Slots');
        xlabel('Average packet arrivals in bits');
        zlabel('Total Packets arrived');
        plotFigure(struct('Y',mean(squeeze(sum(squeeze(SimParams.QueueInfo.packetArrivalsOverTime),1)),2),'N',1));
        
        xlabel('Average packet arrivals in bits');
        ylabel('Slots');
        zlabel('Total Packets backlogged Packets');
        plotFigure(struct('Y',mean(squeeze(sum(squeeze(SimParams.QueueInfo.queueBacklogsOverTime(:,:,:,2:end)),1)),2),'N',2));

        xlabel('Average packet arrivals in bits');
        ylabel('Slots');
        zlabel('Total Packets residual Packets');
        plotFigure(struct('Y',mean(squeeze(sum(squeeze(SimParams.QueueInfo.queueResiduesOverTime(:,:,:,2:end)),1)),2),'N',3));
        
        figStruct.N = 4;
        figStruct.Y = sum(squeeze(SimParams.QueueInfo.residualPkts),1);        
        plotFigure(figStruct);
        xlabel('Average packet arrivals in bits');
        ylabel('Slots');
        zlabel('Final Residual Packets in the system');

                       
    otherwise
        
        display('Simulation Completed without any display !');
        
end

end