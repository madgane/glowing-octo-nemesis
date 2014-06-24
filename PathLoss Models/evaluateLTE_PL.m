function [varargout] = evaluateLTE_PL(SimParams,separationM,disableLOS,disableShadowing)

varargout = cell(1,nargout);
spLight = 3e8;isLOS = 'false';
plModel = char(SimParams.pathLossModel);
uscore_index = find(plModel == '_');
pathLossModel = plModel(uscore_index(1,1) + 1:end);
separationM = max(separationM,SimParams.sysConfig.layoutFeatures.minDistance);

Fc = SimParams.sysConfig.carrierFreqGHz;

switch pathLossModel
    
    case 'InH'
        
        if separationM <= 18
            prob_LOS = 1;
        elseif separationM < 37
            prob_LOS = exp(-(separationM - 18)/27);
        else
            prob_LOS = 0.5;
        end
        
        if strcmp(disableLOS,'true')
            prob_LOS = 0.0;
        end
        
        if rand < prob_LOS
            isLOS = 'true';
            pathLoss_dB = 16.9 * log10(separationM) + 32.8 + 20 * log10(Fc);
            shadowFading = randn * SimParams.sysConfig.shadowing.LOS;
        else
            isLOS = 'false';
            pathLoss_dB = 43.3 * log10(separationM) + 11.5 + 20 * log10(Fc);
            shadowFading = randn * SimParams.sysConfig.shadowing.NLOS;
        end
        
    case 'UMi'
        
        outDoorUsers = (rand < 0.5);
        hUT = SimParams.sysConfig.layoutFeatures.hUT;
        hBS = SimParams.sysConfig.layoutFeatures.hBS;
        
        prob_LOS = min((18/separationM),1) * (1 - exp (-separationM / 36)) + exp (-separationM / 36);
        breakDistance = (4 * (hBS - 1) * (hUT - 1) * Fc * 1e9) / spLight;
        
        if strcmp(disableLOS,'true')
            prob_LOS = 0.0;
        end

        if rand < prob_LOS
            isLOS = 'true';
            if separationM < breakDistance
                pathLoss_dB = 22.0 * log10(separationM) + 28.0 + 20 * log10(Fc);
            else
                pathLoss_dB = 40 * log10(separationM) + 7.8 - 18.0 * log10(hBS - 1) - 18.0 * log10(hUT - 1) + 2.0 * log10(Fc);
            end
            
            shadowFading = randn * SimParams.sysConfig.shadowing.LOS;
        else
            isLOS = 'false';
            pathLoss_dB = 36.7 * log10(separationM) + 22.7 + 26.0 * log10(Fc);
            shadowFading = randn * SimParams.sysConfig.shadowing.NLOS;
        end
        
        vehiclePenetrationLoss = 20 * outDoorUsers;
        shadowFading = shadowFading + vehiclePenetrationLoss;
        
    case 'UMa'
        
        W = rand * 45 + 5;h = rand * 45 + 5;
        hUT = SimParams.sysConfig.layoutFeatures.hUT;
        hBS = SimParams.sysConfig.layoutFeatures.hBS;
        
        prob_LOS = min((18/separationM),1) * (1 - exp (-separationM / 63)) + exp (-separationM / 63);
        breakDistance = (4 * (hBS - 1) * (hUT - 1) * Fc * 1e9) / spLight;
        
        if strcmp(disableLOS,'true')
            prob_LOS = 0.0;
        end

        if rand < prob_LOS
            isLOS = 'true';
            if separationM < breakDistance
                pathLoss_dB = 22.0 * log10(separationM) + 28.0 + 20 * log10(Fc);
            else
                pathLoss_dB = 40 * log10(separationM) + 7.8 - 18.0 * log10(hBS - 1) - 18.0 * log10(hUT - 1) + 2.0 * log10(Fc);
            end
            
            shadowFading = randn * SimParams.sysConfig.shadowing.LOS;
        else
            isLOS = 'false';
            pathLoss_dB = 161.04 - 7.1 * log10(W) + 7.5 * log10(h) - (24.37 - 3.7 * (h/hBS)^2) * log10(hBS) ...
                + (43.42 - 3.1 * log10(hBS)) * (log10(separationM) - 3) + 20 * log10(Fc) - (3.2 * (log10(11.75 * hUT))^2 - 4.97);
            shadowFading = randn * SimParams.sysConfig.shadowing.NLOS;
        end
        
        vehiclePenetrationLoss = 9 + 5 * randn;
        shadowFading = shadowFading + vehiclePenetrationLoss;
        
    case 'SMa'
        
        W = rand * 45 + 5;h = rand * 45 + 5;
        hUT = SimParams.sysConfig.layoutFeatures.hUT;
        hBS = SimParams.sysConfig.layoutFeatures.hBS;
        breakDistance = (2 * pi * hBS * hUT * Fc * 1e9) / spLight;
        
        if separationM <= 10
            prob_LOS = 1;
        else
            prob_LOS = exp(-(separationM - 10) / 200);
        end
        
        if strcmp(disableLOS,'true')
            prob_LOS = 0.0;
        end

        if rand < prob_LOS
            isLOS = 'true';
            if separationM < breakDistance
                pathLoss_dB = 20 * log10(40 * pi * separationM * Fc / 3) + min(0.03 * h^(1.72),10) * log10(separationM) ...
                    - min(0.044 * h^(1.72),14.77) + 0.002 * log10(h) * separationM;
                shadowFading = randn * SimParams.sysConfig.shadowing.LOS;
            else
                pathLoss_dB = 20 * log10(40 * pi * breakDistance * Fc / 3) + min(0.03 * h^(1.72),10) * log10(breakDistance) ...
                    - min(0.044 * h^(1.72),14.77) + 0.002 * log10(h) * breakDistance ...
                    + 40 * log10(separationM / breakDistance);     
                shadowFading = randn * 6;
            end
            
        else
            isLOS = 'false';
            pathLoss_dB = 161.04 - 7.1 * log10(W) + 7.5 * log10(h) - (24.37 - 3.7 * (h/hBS)^2) * log10(hBS) ...
                + (43.42 - 3.1 * log10(hBS)) * (log10(separationM) - 3) + 20 * log10(Fc) - (3.2 * (log10(11.75 * hUT))^2 - 4.97);
            shadowFading = randn * SimParams.sysConfig.shadowing.NLOS;
        end

        vehiclePenetrationLoss = 9 + 5 * randn;
        shadowFading = shadowFading + vehiclePenetrationLoss;        
        
    case 'RMa'
        
        W = rand * 45 + 5;h = rand * 45 + 5;
        hUT = SimParams.sysConfig.layoutFeatures.hUT;
        hBS = SimParams.sysConfig.layoutFeatures.hBS;
        breakDistance = (2 * pi * hBS * hUT * Fc * 1e9) / spLight;
        
        if separationM <= 10
            prob_LOS = 1;
        else
            prob_LOS = exp(-(separationM - 10) / 1000);
        end
        
        if strcmp(disableLOS,'true')
            prob_LOS = 0.0;
        end

        if rand < prob_LOS
            isLOS = 'true';
            if separationM < breakDistance
                pathLoss_dB = 20 * log10(40 * pi * separationM * Fc / 3) + min(0.03 * h^(1.72),10) * log10(separationM) ...
                    - min(0.044 * h^(1.72),14.77) + 0.002 * log10(h) * separationM;
                shadowFading = randn * SimParams.sysConfig.shadowing.LOS;
            else
                pathLoss_dB = 20 * log10(40 * pi * breakDistance * Fc / 3) + min(0.03 * h^(1.72),10) * log10(breakDistance) ...
                    - min(0.044 * h^(1.72),14.77) + 0.002 * log10(h) * breakDistance ...
                    + 40 * log10(separationM / breakDistance);     
                shadowFading = randn * 6;
            end
            
        else
            isLOS = 'false';
            pathLoss_dB = 161.04 - 7.1 * log10(W) + 7.5 * log10(h) - (24.37 - 3.7 * (h/hBS)^2) * log10(hBS) ...
                + (43.42 - 3.1 * log10(hBS)) * (log10(separationM) - 3) + 20 * log10(Fc) - (3.2 * (log10(11.75 * hUT))^2 - 4.97);
            shadowFading = randn * SimParams.sysConfig.shadowing.NLOS;
        end
        
        vehiclePenetrationLoss = 9 + 5 * randn;
        shadowFading = shadowFading + vehiclePenetrationLoss;        
        
end

if separationM > SimParams.sysConfig.layoutFeatures.maxDistance
    pathLoss_dB = 5e6;
end

if strcmp(disableShadowing,'true')
    shadowFading = 0;
end

powerCompensation = SimParams.sysConfig.BStransmitPwr_dBm + SimParams.sysConfig.userTerminalBG + SimParams.sysConfig.baseTerminalBG;
otherNoise = SimParams.sysConfig.baseTerminalNF + SimParams.sysConfig.cableLoss;
varargout{1,1} = powerCompensation - pathLoss_dB - shadowFading - otherNoise - 10 * log10(SimParams.sysConfig.usableTones);

if nargout > 1
    varargout{1,2} = isLOS;
    varargout{1,3} = pathLoss_dB + shadowFading;
end
