function [SimParams,SimStructs] = getRobustBeamformers(SimParams,SimStructs)

for iBase = 1:SimParams.nBases
    for iBand = 1:SimParams.nBands
    
        gammaPrev = 0;
        epsilon = 1e-3;
        continueAgain = 1;
        robustParams.cBase = iBase;
        robustParams.cBand = iBand;
        gammaMax = 100;gammaMin = 0;
        
        while continueAgain
            
            gammaMean = (gammaMax + gammaMin) * 0.5;robustParams.cGamma = gammaMean;
            [cvxFeasibility,X] = performRobustBeamforming(SimParams,SimStructs,robustParams);
            
            if strcmp(cvxFeasibility,'Solved')
                gammaMin = gammaMean;
            else
                gammaMax = gammaMean;
            end
            
            if strcmp(cvxFeasibility,'Solved')
                if abs(gammaMean - gammaPrev) < epsilon
                    continueAgain = 0;
                end
            else
                gammaPrev = gammaMean;
            end
            
        end
        
        SimStructs.baseStruct{iBase}.P{iBand,1} = X * sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand) / trace(X' * X));
        
    end
end

end
