
function [SimParams,SimStructs] = initializeSCAOperatingPoint(SimParams,SimStructs)

initMultiCastVariables;

enableSearch = 1;
if (SimParams.iAntennaArray ~= 1)
    enableSearch = 0;
end

if (SimParams.nTxAntennaEnabled == (singleBandCheck + 1))
    enableSearch = 1;
end

if enableSearch
    
    if SimParams.nBands == 1
        
        switch searchType
            case 'SDP'
                [SimParams,SimStructs] = getSingleBandSDP(SimParams,SimStructs,xRandIterations);
            case 'FC'
                [SimParams,SimStructs] = getSingleBandSCA(SimParams,SimStructs,searchType);
            case 'Dual'
                [SimParams,SimStructs] = getSingleBandSCA(SimParams,SimStructs,searchType);
            case 'KKT'
                [SimParams,SimStructs] = getKKTMCPrecodersInitialization(SimParams,SimStructs);
            otherwise
                for iBase = 1:nBases
                    for iBand = 1:nBands
                        SimStructs.baseStruct{iBase,1}.PG{iBand,1} = complex(randn(SimParams.nTxAntenna,nGroupsPerCell(iBase,1)),randn(SimParams.nTxAntenna,nGroupsPerCell(iBase,1)));
                    end
                end
                display('Using Randomized data !');
        end
        
        for iBase = 1:nBases
            for iBand = 1:nBands
                for iGroup = 1:nGroupsPerCell(iBase,1)
                    groupUsers = SimStructs.baseStruct{iBase,1}.mcGroup{iGroup,1};
                    for iUser = 1:length(groupUsers)
                        cUser = groupUsers(iUser,1);
                        Hsdp = cH{iBase,iBand}(:,:,cUser);
                        SimParams.Debug.tempResource{2,1}{1,1}(cUser,iBand) = real(Hsdp * SimStructs.baseStruct{iBase,1}.PG{iBand,1}(:,iGroup));
                        SimParams.Debug.tempResource{3,1}{1,1}(cUser,iBand) = imag(Hsdp * SimStructs.baseStruct{iBase,1}.PG{iBand,1}(:,iGroup));
                    end
                end
            end
        end
        
        [SimParams,SimStructs] = getSingleBandSCA(SimParams,SimStructs,'MP');
        
        for iBase = 1:nBases
            for iBand = 1:nBands
                for iGroup = 1:nGroupsPerCell(iBase,1)
                    groupUsers = SimStructs.baseStruct{iBase,1}.mcGroup{iGroup,1};
                    for iUser = 1:length(groupUsers)
                        cUser = groupUsers(iUser,1);
                        Hsdp = cH{iBase,iBand}(:,:,cUser);
                        SimParams.Debug.tempResource{2,1}{1,1}(cUser,iBand) = real(Hsdp * SimStructs.baseStruct{iBase,1}.PG{iBand,1}(:,iGroup));
                        SimParams.Debug.tempResource{3,1}{1,1}(cUser,iBand) = imag(Hsdp * SimStructs.baseStruct{iBase,1}.PG{iBand,1}(:,iGroup));
                    end
                end
            end
        end
        
    else
        
        switch searchType
            case 'FC'
                [SimParams,SimStructs] = getMultiBandSCA(SimParams,SimStructs,searchType);
            case 'Dual'
                [SimParams,SimStructs] = getMultiBandSCA(SimParams,SimStructs,searchType);
            otherwise
                for iBase = 1:nBases
                    for iBand = 1:nBands
                        SimStructs.baseStruct{iBase,1}.PG{iBand,1} = complex(randn(SimParams.nTxAntenna,nGroupsPerCell(iBase,1)),randn(SimParams.nTxAntenna,nGroupsPerCell(iBase,1)));
                    end
                end
                display('Using Randomized data !');
        end
        
        for iBase = 1:nBases
            for iBand = 1:nBands
                for iGroup = 1:nGroupsPerCell(iBase,1)
                    groupUsers = SimStructs.baseStruct{iBase,1}.mcGroup{iGroup,1};
                    for iUser = 1:length(groupUsers)
                        cUser = groupUsers(iUser,1);
                        Hsdp = cH{iBase,iBand}(:,:,cUser);
                        SimParams.Debug.tempResource{2,1}{1,1}(cUser,iBand) = real(Hsdp * SimStructs.baseStruct{iBase,1}.PG{iBand,1}(:,iGroup));
                        SimParams.Debug.tempResource{3,1}{1,1}(cUser,iBand) = imag(Hsdp * SimStructs.baseStruct{iBase,1}.PG{iBand,1}(:,iGroup));
                    end
                end
            end
        end
        
        [SimParams,SimStructs] = getMultiBandSCA(SimParams,SimStructs,'MP');
        
        for iBase = 1:nBases
            for iBand = 1:nBands
                for iGroup = 1:nGroupsPerCell(iBase,1)
                    groupUsers = SimStructs.baseStruct{iBase,1}.mcGroup{iGroup,1};
                    for iUser = 1:length(groupUsers)
                        cUser = groupUsers(iUser,1);
                        Hsdp = cH{iBase,iBand}(:,:,cUser);
                        SimParams.Debug.tempResource{2,1}{1,1}(cUser,iBand) = real(Hsdp * SimStructs.baseStruct{iBase,1}.PG{iBand,1}(:,iGroup));
                        SimParams.Debug.tempResource{3,1}{1,1}(cUser,iBand) = imag(Hsdp * SimStructs.baseStruct{iBase,1}.PG{iBand,1}(:,iGroup));
                    end
                end
            end
        end
        
    end
    
    SimParams.Debug.Retain = SimParams.Debug.tempResource;
else
    SimParams.Debug.tempResource = SimParams.Debug.Retain;    
end

display('Initialization point found !');

end