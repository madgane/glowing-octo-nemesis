
function [SimParams,SimStructs] = fwkInitialization(SimParams,SimStructs)

resetRandomness;

% Path Loss Model related code

plModel = char(SimParams.pathLossModel);
uscore_index = find(plModel == '_');

if ~isempty(uscore_index)
    pathLossModel = plModel(1:uscore_index(1,1) - 1);
else
    pathLossModel = plModel;
end

for iUser = 1:SimParams.nUsers
    SimStructs.userStruct{iUser,1}.losFading = cell(SimParams.nBases,1);
    for iBase = 1:SimParams.nBases
        SimStructs.userStruct{iUser,1}.losFading{iBase,1} = 'false';
    end
end

switch pathLossModel
    
    case 'Isolated'
        
        SimParams.PL_Profile = zeros(SimParams.nBases,SimParams.nUsers);
        for iUser = 1:SimParams.nUsers
            modNode = mod(iUser - 1,SimParams.nBases) + 1;
            SimParams.PL_Profile(:,iUser) = -1e5;
            SimParams.PL_Profile(modNode,iUser) = 0;
        end
        
    case 'CellEdge'
        
        SimParams.PL_Profile = zeros(SimParams.nBases,SimParams.nUsers);
        for iUser = 1:SimParams.nUsers
            modNode = mod(iUser - 1,SimParams.nBases) + 1;
            SimParams.PL_Profile(:,iUser) = -1/1e5;
            SimParams.PL_Profile(modNode,iUser) = 0;
        end
        
    case 'Random'
        
        [SimParams] = getRandomPathLoss(SimParams);
        
    case 'Perturbed'
        
        if ~isempty(uscore_index)
            randGain = str2double(plModel(uscore_index(1,1) + 1:end));
        else
            randGain = 1;
        end
        
        userVector = 1:SimParams.nUsers;
        shareUserCnt = SimParams.nUsers / SimParams.nBases;
        SimParams.PL_Profile = zeros(SimParams.nBases,SimParams.nUsers);
        
        for iBase = 1:SimParams.nBases
            userIndices = (1 + (iBase - 1) * shareUserCnt):(iBase * shareUserCnt);oBaseIndices = setxor(userIndices,userVector);
            SimParams.PL_Profile(iBase,userIndices) = SimParams.PL_Profile(iBase,userIndices) + rand(1,length(userIndices)) * randGain;
            SimParams.PL_Profile(iBase,oBaseIndices) = SimParams.PL_Profile(iBase,oBaseIndices) - rand(1,length(oBaseIndices)) * randGain;
        end
        
    case 'Fixed'
        
        SimParams.PL_Profile = SimParams.PL_Profile;
        
    case 'LTE'
        
        switch plModel(uscore_index(1,1) + 1:end)
            case 'UMi'
                loadParams = importdata('Utilities\UserStats_uMicro10.dat');
            case 'UMa'
                loadParams = importdata('Utilities\UserStats_uMacro10.dat');
        end
        
        loadParams = loadParams(1:570,:);
        SimParams.PL_Profile = -ones(SimParams.nBases,SimParams.nUsers) * Inf;
        for iCell = 1:SimParams.nBases
            cCell = iCell - 1;
            uLocs = find(cCell == loadParams(:,2));
            SimParams.PL_Profile(iCell,uLocs) = loadParams(uLocs,end);
        end
        
    case '3GPP'
        
        SimParams.PL_Profile = zeros(SimParams.nBases,SimParams.nUsers);
        
        [SimParams] = configureLTEParams(SimParams);
        [SimParams, SimStructs] = cellLayoutGeneration(SimParams,SimStructs);
        [SimParams, SimStructs] = userLayoutGeneration(SimParams,SimStructs);
                
    otherwise
        
        xdB = char(SimParams.pathLossModel);xI = strfind(xdB,'_');
        SimParams.PL_Profile = zeros(SimParams.nBases,SimParams.nUsers);
        for iUser = 1:SimParams.nUsers
            modNode = mod(iUser - 1,SimParams.nBases) + 1;
            SimParams.PL_Profile(:,iUser) = str2double(xdB(xI+1:end));
            SimParams.PL_Profile(modNode,iUser) = 0;
        end
end

% Doppler / Small scale related code

if strcmp(pathLossModel,'3GPP')
    dopplerRealizations = ones(SimParams.nUsers,1) * SimParams.sysConfig.userDoppler;
else
    dopplerType = char(SimParams.DopplerType);
    uscore_index = find(dopplerType == '_');
    
    if ~isempty(uscore_index)
        dopType = dopplerType(1:uscore_index(1,1) - 1);
        currDoppler = str2double(dopplerType(uscore_index(1,1) + 1:end));
    else
        dopType = dopplerType;
    end
    
    switch dopType
        case 'Uniform'
            dopplerRealizations = rand(SimParams.nUsers,1) * currDoppler;
        case 'Constant'
            dopplerRealizations = ones(SimParams.nUsers,1) * currDoppler;
    end
end

SimParams.userDoppler = dopplerRealizations;
SimStructs.JakesChStruct = cell(SimParams.nUsers,SimParams.nBases,SimParams.nBands);

if strcmp(pathLossModel,'3GPP')
    if strcmp(SimParams.ChannelModel,'Jakes')
        for iUser = 1:SimParams.nUsers
            currentDoppler = SimParams.userDoppler(iUser,1);
            xCites = SimStructs.userStruct{iUser,1}.phyParams.listedCites;
            for iBase = 1:SimParams.nBases
                if isempty(find(iBase == xCites))
                    continue;
                end
                 for iBand = 1:SimParams.nBands
                    SimStructs.JakesChStruct{iUser,iBase,iBand} = comm.MIMOChannel;
                    if strfind(version,'R2012')
                        SimStructs.JakesChStruct{iUser,iBase,iBand}.ReceiveCorrelationMatrix = eye(SimParams.nRxAntenna);
                        SimStructs.JakesChStruct{iUser,iBase,iBand}.TransmitCorrelationMatrix = eye(SimParams.nTxAntenna);
                    else
                        SimStructs.JakesChStruct{iUser,iBase,iBand}.SpatialCorrelation = false;
                    end
                    SimStructs.JakesChStruct{iUser,iBase,iBand}.SampleRate = SimParams.SFSymbols / SimParams.sampTime;
                    SimStructs.JakesChStruct{iUser,iBase,iBand}.MaximumDopplerShift = currentDoppler;
                    SimStructs.JakesChStruct{iUser,iBase,iBand}.NumTransmitAntennas = SimParams.nTxAntenna;
                    SimStructs.JakesChStruct{iUser,iBase,iBand}.NumReceiveAntennas = SimParams.nRxAntenna;
                    SimStructs.JakesChStruct{iUser,iBase,iBand}.NormalizePathGains = 1;
                    SimStructs.JakesChStruct{iUser,iBase,iBand}.PathGainsOutputPort = 1;
                    SimStructs.JakesChStruct{iUser,iBase,iBand}.AveragePathGains = 0;
                    if strcmp(SimStructs.userStruct{iUser,1}.losFading{iBase,1},'true')
                        SimStructs.JakesChStruct{iUser,iBase,iBand}.FadingDistribution = 'Rician';
                        kFactor = SimParams.sysConfig.Kfactor.avg + SimParams.sysConfig.Kfactor.std * randn;
                        SimStructs.JakesChStruct{iUser,iBase,iBand}.KFactor = 10^(kFactor/10);
                    end
                end
            end
        end
    end
else
    if strcmp(SimParams.ChannelModel,'Jakes')
        for iUser = 1:SimParams.nUsers
            currentDoppler = SimParams.userDoppler(iUser,1);
            for iBase = 1:SimParams.nBases
                for iBand = 1:SimParams.nBands
                    SimStructs.JakesChStruct{iUser,iBase,iBand} = comm.MIMOChannel;
                    if strfind(version,'R2012')
                        SimStructs.JakesChStruct{iUser,iBase,iBand}.ReceiveCorrelationMatrix = eye(SimParams.nRxAntenna);
                        SimStructs.JakesChStruct{iUser,iBase,iBand}.TransmitCorrelationMatrix = eye(SimParams.nTxAntenna);
                    else
                        SimStructs.JakesChStruct{iUser,iBase,iBand}.SpatialCorrelation = false;
                    end
                    SimStructs.JakesChStruct{iUser,iBase,iBand}.SampleRate = SimParams.SFSymbols / SimParams.sampTime;
                    SimStructs.JakesChStruct{iUser,iBase,iBand}.MaximumDopplerShift = currentDoppler;
                    SimStructs.JakesChStruct{iUser,iBase,iBand}.NumTransmitAntennas = SimParams.nTxAntenna;
                    SimStructs.JakesChStruct{iUser,iBase,iBand}.NumReceiveAntennas = SimParams.nRxAntenna;
                    SimStructs.JakesChStruct{iUser,iBase,iBand}.NormalizePathGains = 1;
                    SimStructs.JakesChStruct{iUser,iBase,iBand}.PathGainsOutputPort = 1;
                    SimStructs.JakesChStruct{iUser,iBase,iBand}.AveragePathGains = 0;
                    if strcmp(SimStructs.userStruct{iUser,1}.losFading{iBase,1},'true')
                        SimStructs.JakesChStruct{iUser,iBase,iBand}.FadingDistribution = 'Rician';
                        kFactor = SimParams.sysConfig.Kfactor.avg + SimParams.sysConfig.Kfactor.std * randn;
                        SimStructs.JakesChStruct{iUser,iBase,iBand}.KFactor = 10^(kFactor/10);
                    end
                end
            end
        end
    end
end

SimParams.updateFeedback = zeros(SimParams.nUsers,1);

for iUser = 1:SimParams.nUsers
    feedbackCycle = (1 / SimParams.userDoppler(iUser,1)) * SimParams.fbFraction;
    SimParams.updateFeedback(iUser,1) = round(feedbackCycle / SimParams.sampTime) + 1;
end

SimParams.feedbackOffset = randi([0,min(SimParams.updateFeedback) - 1],SimParams.nUsers,1);
SimParams.Debug.resAllocation = zeros(SimParams.nDrops,SimParams.nBands,SimParams.nUsers,length(SimParams.snrIndex));
SimParams.Debug.privateExchanges.resAllocation = zeros(SimParams.nBands,SimParams.nUsers);

end
