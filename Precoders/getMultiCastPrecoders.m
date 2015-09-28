
function [SimParams,SimStructs] = getMultiCastPrecoders(SimParams,SimStructs)

initMultiCastVariables;

% Debug Buffers initialization

SimParams.Debug.tempResource{2,SimParams.iDrop} = cell(SimParams.nUsers,1);
SimParams.Debug.tempResource{3,SimParams.iDrop} = cell(SimParams.nUsers,1);
SimParams.Debug.tempResource{4,SimParams.iDrop} = cell(SimParams.nUsers,SimParams.nBands);

underscore_location = strfind(SimParams.DesignType,'_');
if isempty(underscore_location)
    selectionMethod = SimParams.DesignType;
else
    selectionMethod = SimParams.DesignType(1:underscore_location-1);
    designType = SimParams.DesignType(underscore_location+1:end);
end

SimParams.Debug.tempResource{2,1}{1,1} = randn(nUsers,nBands);
SimParams.Debug.tempResource{3,1}{1,1} = randn(nUsers,nBands);
SimParams.Debug.tempResource{4,1}{1,1} = rand(nUsers,nBands) + 1;

if isfield(SimParams.Debug,'MultiCastSDPExchange')
    SimParams.Debug = rmfield(SimParams.Debug,'MultiCastSDPExchange');
end

switch selectionMethod
    
    case 'SB-SDP'
        
        [SimParams,SimStructs] = getSingleBandSDP(SimParams,SimStructs,xRandIterations);
        
    case 'SB-SCA'
        
        switch searchType
            case 'SDP'
                [SimParams,SimStructs] = getSingleBandSCA(SimParams,SimStructs,1);
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
        
        display('Initialization point found !');
        [SimParams,SimStructs] = getSingleBandSCA(SimParams,SimStructs,'MP');
        
    case 'KKTMethod'
        
        [SimParams,SimStructs] = getKKTMultiCastPrecoders(SimParams,SimStructs);
        
    case 'SB-SDPA'
        
        [SimParams,SimStructs] = getSingleBandSDPA(SimParams,SimStructs,xRandIterations);
        
    case 'SB-SCAA'
        
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
        
        display('Initialization point found !');
        
        switch designType
            
            case 'A'
                [SimParams,SimStructs] = getSingleBandSCAA_A(SimParams,SimStructs);
            case 'B'
                [SimParams,SimStructs] = getSingleBandSCAA_B(SimParams,SimStructs);
            case 'C'
                [SimParams,SimStructs] = getSingleBandSCAA_C(SimParams,SimStructs);
        end
        
        display('Antenna subset selected !');
        [SimParams,SimStructs] = getSingleBandSCA(SimParams,SimStructs,'MP');
        
    case 'MB-SDP'
        
        [SimParams,SimStructs] = getMultiBandSDP(SimParams,SimStructs);
        [SimParams,SimStructs] = getMultiBandSDP_SB(SimParams,SimStructs,xRandIterations);
        
    case 'MB-SDPA'
        
        [SimParams,SimStructs] = getMultiBandSDPA(SimParams,SimStructs);
        [SimParams,SimStructs] = getMultiBandSDP_SB(SimParams,SimStructs,xRandIterations);

        
    case 'MB-SCA'
        
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
        
        display('Initialization point found !');
        [SimParams,SimStructs] = getMultiBandSCA(SimParams,SimStructs,'MP');
        
    case 'MB-SCAA'
        
        searchType = 'Dual';
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
        
        display('Initialization point found !');
        switch designType
            
            case 'A'
                [SimParams,SimStructs] = getMultiBandSCAA_A(SimParams,SimStructs);
            case 'B'
                [SimParams,SimStructs] = getMultiBandSCAA_B(SimParams,SimStructs);
            case 'C'
                [SimParams,SimStructs] = getMultiBandSCAA_C(SimParams,SimStructs);
            case 'D'
                [SimParams,SimStructs] = getMultiBandSCAA_D(SimParams,SimStructs);
        end

        display('Antenna subset selected !');
        [SimParams,SimStructs] = getMultiBandSCA(SimParams,SimStructs,'MP');

    case 'MB-SCAS'
        
        if (SimParams.iAntennaArray == 1)

            [SimParams,SimStructs] = getMultiBandSCA(SimParams,SimStructs,'Dual');
            
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

            display('Initialization point found !');
            [SimParams,SimStructs] = getMultiBandSCAS(SimParams,SimStructs);
            
            display('Antenna subset selected !');
            [SimParams,SimStructs] = getMultiBandSCA(SimParams,SimStructs,'MP');
            
        else
            
            for iBase = 1:nBases
                for iBand = 1:nBands
                    SimStructs.baseStruct{iBase,1}.PG{iBand,1} = complex(randn(SimParams.nTxAntenna,nGroupsPerCell(iBase,1)),randn(SimParams.nTxAntenna,nGroupsPerCell(iBase,1)));
                end
            end
            
        end
        
    case 'MB-SCAE'
        
        [SimParams,SimStructs] = getMultiBandSCAE(SimParams,SimStructs);
       
    otherwise
        display('Unknown Precoding Method !');
        
end

end
