
function [displayString] = stringFormatting(SimParams,LegendOrPlotLables,nQueueStability,userVariable)

if strcmp(LegendOrPlotLables,'legend')
    
    switch (SimParams.PrecodingMethod)
        
        case 'Best_ZF_Method'
            displayString = sprintf('%s','ZF');
            
        case 'Best_BF_Method'
            displayString = sprintf('%s','BF');
            
        case 'Best_MZF_Method'
            displayString = sprintf('%s','MZF');
            
        case 'Best_Network_Method'
            displayString = sprintf('%s','NW-BF');
            
        case 'Best_DP_Method'
            displayString = sprintf('%s','DP');
            
        case 'Best_MaxMinSINR_Method'
            displayString = sprintf('%s','max-min-SINR');
            
        case 'Best_WSR_Method'
            displayString = sprintf('%s','WSR');
            
        case 'Best_WMMSE_Method'
            displayString = sprintf('%s','WMMSE');
            
        case 'Best_TDMZF_Method'
            displayString = sprintf('%s','TDMZF');
            
        case 'Best_CZF_Method'
            displayString = sprintf('%s','CZF');
            
        otherwise
            display('Unknown Precoder Design !');
            
    end
    
    if strcmp(SimParams.PrecodingMethod,'Best_WMMSE_Method')
        
        switch SimParams.weightedSumRateMethod
            
            case 'PreScheduling'
                displayString = sprintf('%s-%s',displayString,'RedUsers');
                
            case 'PerformScheduling'
                displayString = sprintf('%s-%s',displayString,'OverLoaded');
                
            case 'DistScheduling'
                displayString = sprintf('%s-%s',displayString,'OverLoaded');
                
            case 'CNetworkBF'
                displayString = sprintf('%s-%s',displayString,'OverLoaded');
                
            case 'DNetworkBF'
                displayString = sprintf('%s-%s',displayString,'OverLoaded');
                
            case 'StreamScheduling'
                displayString = sprintf('%s-%s',displayString,'RedUsers');
                
            otherwise
                display('Unknown Weighted Sum Rate Option !');
                
        end
        
    end
    
    chScheduler = char(SimParams.SchedType);
    uScoreIndex = find(chScheduler == '_');
    if isempty(uScoreIndex)
        scheduleMethod = SimParams.SchedType;
    else
        scheduleMethod = chScheduler(1:uScoreIndex(1,1) - 1);
    end
    
    switch scheduleMethod
        
        case 'RRScheduling'
            displayString = sprintf('%s-%s',displayString,'RRSched');
            
        case 'RouletteScheduling'
            displayString = sprintf('%s-%s',displayString,'RoulSched');
            
        case 'GreedyScheduling'
            displayString = sprintf('%s-%s',displayString,'GreedySched');
            
        case 'PFScheduling'
            displayString = sprintf('%s-%s',displayString,'PFSched');
            if length(uScoreIndex) > 1
                displayString = sprintf('%s-%s',displayString,'50%ile');
            end
            
        case 'BDScheduling'
            displayString = sprintf('%s-%s',displayString,'BDSched');
            
        case 'PFBDScheduling'
            displayString = sprintf('%s-%s',displayString,'PFBDSched');
            
        case 'GeneticScheduling'
            displayString = sprintf('%s-%s',displayString,'GenSched');
            
        case 'ExScheduling'
            displayString = sprintf('%s-%s',displayString,'ExSched');
            
        case 'XScheduling'
            displayString = sprintf('%s-%s',displayString,'XSched');
            
        case 'NetworkScheduling'
            displayString = sprintf('%s-%s',displayString,'NWSched');
            
        case 'CoordScheduling'
            displayString = sprintf('%s-%s',displayString,'CoordSched');
            
        otherwise
            display('Unknown Scheduling Type');
    end
    
    if ~strcmp(SimParams.PrecodingMethod,'Best_WMMSE_Method')
        
        switch SimParams.queueWt
            
            case 0
                displayString = sprintf('%s-%s',displayString,'Act-WF-PA');
                
            case 1
                displayString = sprintf('%s-%s',displayString,'Qwt-WF-PA');
                
            case 2
                displayString = sprintf('%s-%s',displayString,'Qwt-QWSR-PA');
                
            otherwise
                display('Unknown WF PA Type');
        end
        
    else
        
        switch SimParams.queueWt
            
            case 0
                displayString = sprintf('%s-%s',displayString,'Act');
                
            case 1
                displayString = sprintf('%s-%s',displayString,'QwtUsers');
                
        end
        
    end
    
    if strcmp(scheduleMethod,'XScheduling')
        xScheme = chScheduler(uScoreIndex(1,1) + 1:end);
        switch xScheme
            case 'EqualShare'
                displayString = sprintf('%s-%s',displayString,'I-QR');
            case 'StreamSearch'
                displayString = sprintf('%s-%s',displayString,'J-QR');
        end
    end
    
else
    
    if strcmp(SimParams.queueMode,'false')
        if (~userVariable)
            plotLabel.plotTitle = sprintf('%s','Sum Rate Plot');
            plotLabel.xLabel = sprintf('%s','SNR in dB');
            plotLabel.yLabel = sprintf('%s','sum rate in bits/sec/Hz');
        else
            plotLabel.plotTitle = sprintf('%s','Sum Rate Plot');
            plotLabel.xLabel = sprintf('%s','Number of Users');
            plotLabel.yLabel = sprintf('%s','sum rate in bits/sec/Hz');
        end
    else
        if nQueueStability
            if nQueueStability == 1
                plotLabel.xLabel = sprintf('%s','arrival rate in bits (or) pkts / user');
                plotLabel.yLabel = sprintf('%s','std (\sigma_Q) of queue backlogs in bits (or) pkts');
                plotLabel.plotTitle = sprintf('%s','Deviation over Queue upon different arrivals');
            else
                plotLabel.xLabel = sprintf('%s','arrival rate in bits (or) pkts / user');
                plotLabel.yLabel = sprintf('%s','Expected queue backlogs in bits (or) pkts');
                plotLabel.plotTitle = sprintf('%s','Expected Queue Size over arrivals');
            end
        else
            plotLabel.xLabel = sprintf('%s','slot indices');
            plotLabel.yLabel = sprintf('%s','total backlogged packets');
            plotLabel.plotTitle = sprintf('%s','Queue Stability Plot');
        end
    end
    
    plotLabel.plotTitle = sprintf('%s-PLType-%s-Users-%d-BS-%d',plotLabel.plotTitle,SimParams.pathLossModel,...
        SimParams.nUsers,SimParams.nBases);
    
    displayString = plotLabel;
    
end

