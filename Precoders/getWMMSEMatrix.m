function [SimParams,SimStructs] = getWMMSEMatrix(SimParams,SimStructs)

switch SimParams.weightedSumRateMethod
    
    case 'PreScheduling'
        
        [SimParams,SimStructs] = getPreWeightedMMSEDesign(SimParams,SimStructs);

    case 'PerformScheduling'
        
        [SimParams,SimStructs] = getWeightedMMSEDesign(SimParams,SimStructs);
        
    case 'DistributedScheduling'
        
        [SimParams,SimStructs] = getDWeightedMMSEDesign(SimParams,SimStructs);
        
    case 'CNetworkBF'
        
        [SimParams,SimStructs] = getCNetworkBFWMMSEDesign(SimParams,SimStructs);
        
    case 'DNetworkBF'
        
        [SimParams,SimStructs] = getDNetworkBFWMMSEDesign(SimParams,SimStructs);
        
        
    case 'StreamScheduling'
        
        [SimParams,SimStructs] = getStrWeightedMMSEDesign(SimParams,SimStructs);
        
    case 'DistScheduling'
        
        [SimParams,SimStructs] = getDistStrWeightedMMSEDesign(SimParams,SimStructs);
        
    case 'IndepScheduling'
        
        [SimParams,SimStructs] = getIndStrWeightedMMSEDesign(SimParams,SimStructs);
                
    otherwise
        
        display('Unknown Weighted Sum Rate Option !');
                
end

if strcmp(SimParams.DebugMode,'true')
    
    displayArray = zeros(SimParams.nBands,SimParams.nBases + 1);
    
    for iBand = 1:SimParams.nBands
        for iBase = 1:SimParams.nBases
            displayArray(iBand,iBase) = trace(SimStructs.baseStruct{iBase,1}.P{iBand,1}' * SimStructs.baseStruct{iBase,1}.P{iBand,1});
        end
       
        displayArray(iBand,iBase + 1) = SimParams.sPower(1,iBand);
    end
    
    display(displayArray);
    
end

end
