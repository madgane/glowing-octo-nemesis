function displayOutputs(xParams, xStructs)

SimParams = xParams{1,1};
SimStructs = xStructs{1,1};

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
        
        baseSDP_Power = zeros(1,length(xParams));
        baseBF_Power = zeros(1,length(xParams));        
        
        for iAntennaArray = 1:length(xParams)
                   
            SimParams = xParams{iAntennaArray,1};
            SimStructs = xStructs{iAntennaArray,1};
            
            totalBFPower = 0;
            totalSDPPower = 0;
            displayQueues(SimParams,SimStructs);
            for iBand = 1:SimParams.nBands
                for iBase = 1:SimParams.nBases
                    sdpPower = 0;
                    for iGroup = 1:length(SimStructs.baseStruct{iBase,1}.mcGroup)
                        if ~isempty(SimStructs.baseStruct{iBase,1}.P_SDP{iBand,1})
                            sdpPower = sdpPower + real(trace(SimStructs.baseStruct{iBase,1}.P_SDP{iBand,1}(:,:,iGroup)));
                        end
                    end
                    bfPower = real(trace(SimStructs.baseStruct{iBase,1}.PG{iBand,1} * SimStructs.baseStruct{iBase,1}.PG{iBand,1}'));
                    fprintf('Transmit Power on SC [%d] from BS [%d] : BF Pwr - %f \t SDP Pwr (LB) - %f \n',iBand,iBase,bfPower,sdpPower);
                    
                    totalBFPower = totalBFPower + bfPower;
                    totalSDPPower = sdpPower + totalSDPPower;
                end
                fprintf('\n');
            end
            
            fprintf('Total Power SDP (LB) - %f, BF - %f \n',totalSDPPower,totalBFPower);
            baseSDP_Power(1,iAntennaArray) = db(totalSDPPower,'power');
            baseBF_Power(1,iAntennaArray) = db(totalBFPower,'power');
            
            if ((strcmpi(SimParams.DesignType,'MB-SCAS')) && (iAntennaArray == 1))
                SimParams.nTxAntennaEnabledArray = SimParams.nTxAntennaEnabled:SimParams.nAntennaArray;
                baseBF_Power = db(totalBFPower,'power') * ones(1,length(SimParams.nTxAntennaEnabledArray));
                break;
            end
            
        end
        
        if ~isempty(strfind(SimParams.DesignType,'SDP'))
            plotFigure(struct('X',SimParams.nTxAntennaEnabledArray,'Y',baseSDP_Power,'N',1));
        end
        plotFigure(struct('X',SimParams.nTxAntennaEnabledArray,'Y',baseBF_Power,'N',1));
        
        xlabel('Total Number of Active Antenna Elements ({X_T})');
        ylabel('Transmit Power in dB');
        
    case 'MCPolar'
        
        baseBF_Power = zeros(1,length(xParams));                
        for iAntennaArray = 1:length(xParams)
            
            SimParams = xParams{iAntennaArray,1};
            SimStructs = xStructs{iAntennaArray,1};
            
            iBand = 1;
            xAngle = linspace(0,2*pi,1025);
            xAG = zeros(length(xAngle),2);
            displayQueues(SimParams,SimStructs);
            for iAngle = 1:length(xAngle)
                xFFT = 0;
                for iBase = 1:SimParams.nBases
                    for iGroup = 1:length(SimStructs.baseStruct{iBase,1}.mcGroup)
                        gTheta = exp(sqrt(-1) * pi * sin(xAngle(1,iAngle)) * (0:SimParams.nTxAntenna-1)) / sqrt(SimParams.nAntennaArray);
                        xFFT = xFFT + abs(gTheta * SimStructs.baseStruct{iBase,1}.PG{iBand,1}(:,iGroup))^2;
                    end
                end
                xAG(iAngle,:) = [xAngle(1,iAngle), xFFT];
            end            
            
            baseBF_Power(iAntennaArray,1) = db(real(trace(SimStructs.baseStruct{iBase,1}.PG{iBand,1} * SimStructs.baseStruct{iBase,1}.PG{iBand,1}')),'power');
            polar(xAG(:,1),xAG(:,2));hold all;
            
        end
        
        fprintf('\nTotal transmit power required for current transmission - \n');
        fprintf('%3.4f\t',baseBF_Power);
        fprintf('\n');
        
    case 'DispMCFInfo'
        
        minGrtRate = zeros(1,length(xParams));
        
        for iAntennaArray = 1:length(xParams)
                   
            SimParams = xParams{iAntennaArray,1};
            SimStructs = xStructs{iAntennaArray,1};
            displayQueues(SimParams,SimStructs);
            
            userTxPkts = zeros(SimParams.nUsers,1);
            for iUser = 1:SimParams.nUsers
                userTxPkts(iUser,1) = sum(squeeze(SimParams.Debug.resAllocation(SimParams.nDrops,:,iUser,end)),2);
            end
            
            minGrtRate(1,iAntennaArray) = min(userTxPkts);
        end
        
        plotFigure(struct('X',SimParams.nTxAntennaEnabledArray,'Y',minGrtRate,'N',1));
        
        xlabel('Total Number of Active Antenna Elements ({X_T})');
        ylabel('Minimum rate achieved in nats');

        
    otherwise
        
        display('Simulation Completed without any display !');
        
end

end