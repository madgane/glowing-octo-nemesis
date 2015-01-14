
function [SimParams,SimStructs] = dropInitialize(SimParams,SimStructs)

% Channel Generation + Addition of Path Loss

uscoreIndex = find(SimParams.pathLossModel == '_');
SimParams.elapsedFBDuration = zeros(SimParams.nUsers,1);

if isempty(uscoreIndex)
    uscoreIndex = length(SimParams.pathLossModel) + 1;
end

for iBand = 1:SimParams.nBands
    for iBase = 1:SimParams.nBases
        for iUser = 1:SimParams.nUsers
            
            if strcmp(SimParams.pathLossModel(1:uscoreIndex(1,1) - 1),'3GPP')
                xRSSI = SimStructs.userStruct{iUser,1}.phyParams.listedRSSI;
                xCites = SimStructs.userStruct{iUser,1}.phyParams.listedCites;
                if isempty(find(iBase == xCites))
                    continue;
                else
                    cRSSI = xRSSI(find(iBase == xCites),1);
                    PL = 10^(cRSSI / 20);
                end
            else
                PL = 10^(SimParams.PL_Profile(iBase,iUser) / 20);
            end
            
            estError = complex(randn(SimParams.nRxAntenna,SimParams.nTxAntenna),randn(SimParams.nRxAntenna,SimParams.nTxAntenna)) * sqrt(0.5 * SimParams.estError) * PL;
           
            switch (SimParams.ChannelModel)
                
                case 'AWGN'
                    
                    SimStructs.actualChannel{iBase,iBand}(:,:,iUser) = complex(ones(SimParams.nRxAntenna,SimParams.nTxAntenna), ...
                        zeros(SimParams.nRxAntenna,SimParams.nTxAntenna)) * PL;
                    SimStructs.prevChan{iBase,iBand}(:,:,iUser) = complex(ones(SimParams.nRxAntenna,SimParams.nTxAntenna), ...
                        zeros(SimParams.nRxAntenna,SimParams.nTxAntenna)) * sqrt(1 - SimParams.estError) * PL + estError;
                    SimStructs.linkChan{iBase,iBand}(:,:,iUser) = sqrt(1 - SimParams.estError) * SimStructs.actualChannel{iBase,iBand}(:,:,iUser) + estError;

                                        
                case 'IID'
                    
                    SimParams.elapsedFBDuration(iUser,1) = 0;
                    SimStructs.actualChannel{iBase,iBand}(:,:,iUser) = complex(randn(SimParams.nRxAntenna,SimParams.nTxAntenna), ...
                        randn(SimParams.nRxAntenna,SimParams.nTxAntenna)) * sqrt(0.5) * PL;                    
                    SimStructs.prevChan{iBase,iBand}(:,:,iUser) = complex(randn(SimParams.nRxAntenna,SimParams.nTxAntenna), ...
                        randn(SimParams.nRxAntenna,SimParams.nTxAntenna)) * sqrt(0.5) * sqrt(1 - SimParams.estError) * PL + estError;    
                    SimStructs.linkChan{iBase,iBand}(:,:,iUser) = sqrt(1 - SimParams.estError) * SimStructs.actualChannel{iBase,iBand}(:,:,iUser) + estError;
                    
                case 'Jakes'
                    
                    SimParams.Debug.feedbackUpdate(iUser,1) = 0;
                    if SimParams.iDrop == 1
                        [~,PathGains] = step(SimStructs.JakesChStruct{iUser,iBase,iBand},ones(SimParams.SFSymbols,SimParams.nTxAntenna));
                        SimStructs.prevChan{iBase,iBand}(:,:,iUser) = reshape(squeeze(PathGains(SimParams.SFSymbols,end,:,:)).',SimParams.nRxAntenna,SimParams.nTxAntenna) * PL * sqrt(1 - SimParams.estError) + estError;
                    else
                        SimStructs.prevChan{iBase,iBand}(:,:,iUser) = SimStructs.linkChan{iBase,iBand}(:,:,iUser);
                    end
                    
                    [~,PathGains] = step(SimStructs.JakesChStruct{iUser,iBase,iBand},ones(SimParams.SFSymbols,SimParams.nTxAntenna));
                    SimStructs.actualChannel{iBase,iBand}(:,:,iUser) = reshape(squeeze(PathGains(SimParams.SFSymbols,end,:,:)).',SimParams.nRxAntenna,SimParams.nTxAntenna) * PL;
                    
                    if ~mod((SimParams.iDrop + SimParams.feedbackOffset(iUser,1) - 1),SimParams.updateFeedback(iUser,1))
                        SimStructs.linkChan{iBase,iBand}(:,:,iUser) = sqrt(1 - SimParams.estError) * SimStructs.actualChannel{iBase,iBand}(:,:,iUser) + estError;
                        if iBase * iBand == 1
                            SimParams.Debug.feedbackUpdate(iUser,1) = 1;
                        end
                    end
                    
                    if SimParams.iDrop == 1
                        SimStructs.linkChan{iBase,iBand}(:,:,iUser) = sqrt(1 - SimParams.estError) * SimStructs.actualChannel{iBase,iBand}(:,:,iUser) + estError;
                        if iBase * iBand == 1
                            SimParams.Debug.feedbackUpdate(iUser,1) = 1;
                        end
                    end
                    
                    if SimStructs.userStruct{iUser,1}.baseNode == iBase
                        if SimParams.Debug.feedbackUpdate(iUser,1) == 1
                            SimParams.elapsedFBDuration(iUser,1) = 0;
                        else
                            SimParams.elapsedFBDuration(iUser,1) = SimParams.elapsedFBDuration(iUser,1) + 1;
                        end
                    end                    
                    
            end
        end
    end
end

alphaPF = SimParams.PF_dur^-1;
for iUser = 1:SimParams.nUsers
    SimStructs.userStruct{iUser,1}.W = cell(SimParams.nBands,1);
    SimStructs.userStruct{iUser,1}.trafficConfig.currentArrival = SimStructs.userStruct{iUser,1}.trafficHistory.pktArrival(1,SimParams.iDrop);
    SimStructs.userStruct{iUser,1}.PFmetric = (1 - alphaPF) * SimStructs.userStruct{iUser,1}.crThrpt + alphaPF * SimStructs.userStruct{iUser,1}.lastThrpt;
end

for iBase = 1:SimParams.nBases
    SimStructs.baseStruct{iBase}.assignedUsers = cell(SimParams.nBands,1);
    SimStructs.baseStruct{iBase}.assignedStreams = cell(SimParams.nBands,1);
    SimStructs.baseStruct{iBase}.P = cell(SimParams.nBands,1);
end

% Traffic Related Calculations !

[SimParams,SimStructs] = updateQueueStats(SimParams,SimStructs);

activeUsers = zeros(SimParams.nUsers,1);
weighingFactor = zeros(SimParams.nUsers,1);

for iUser = 1:SimParams.nUsers
    weighingFactor(iUser,1) = SimStructs.userStruct{iUser,1}.trafficStats.backLogPkt;
    activeUsers(iUser,1) = sign(SimStructs.userStruct{iUser,1}.trafficStats.backLogPkt);
end

for iUser = 1:SimParams.nUsers
    if strcmp(SimParams.weighingEqual,'true')
        SimStructs.userStruct{iUser,1}.weighingFactor = activeUsers(iUser,1);
    else
        SimStructs.userStruct{iUser,1}.weighingFactor = weighingFactor(iUser,1);
    end
end

for iUser = 1:SimParams.nUsers
    switch SimParams.mdpFactor
        case 0
        case 2
            chnInstant = SimParams.sampTime;
            mdpFactor = abs(besselj(0,(2 * pi * SimParams.userDoppler(iUser,1) * chnInstant)));
            SimStructs.userStruct{iUser,1}.weighingFactor = SimStructs.userStruct{iUser,1}.weighingFactor * (mdpFactor)^(SimParams.elapsedFBDuration(iUser,1));
        otherwise
            SimStructs.userStruct{iUser,1}.weighingFactor = SimStructs.userStruct{iUser,1}.weighingFactor * (SimParams.mdpFactor)^(SimParams.elapsedFBDuration(iUser,1));
    end
end

saveChannelInformation;
[SimParams,SimStructs] = updateFFRProfile(SimParams,SimStructs);

end
