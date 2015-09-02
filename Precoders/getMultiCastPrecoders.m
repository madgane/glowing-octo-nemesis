
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
end

SimParams.Debug.tempResource{2,1}{1,1} = randn(nUsers,nBands);
SimParams.Debug.tempResource{3,1}{1,1} = randn(nUsers,nBands);
SimParams.Debug.tempResource{4,1}{1,1} = rand(nUsers,nBands) + 1;
SimParams.Debug.tempResource{4,1}{2,1} = zeros(nUsers,nBands);

switch selectionMethod
    
    case 'SDPMethod'
        
        [SimParams,SimStructs] = getMultiCastSDP(SimParams,SimStructs,0,250);
        
    case 'ConicMethod'
        
        searchType = 'SDP';
        switch searchType
            case 'SDP'
                [SimParams,SimStructs] = getMultiCastSDP(SimParams,SimStructs,0,1);
            case 'FC'
                [SimParams,SimStructs] = getMultiCastConicA(SimParams,SimStructs,searchType);
            case 'Dual'
                [SimParams,SimStructs] = getMultiCastConicA(SimParams,SimStructs,searchType);
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
        [SimParams,SimStructs] = getMultiCastConicA(SimParams,SimStructs,'MP');
        
    case 'KKTMethod'
        
        [SimParams,SimStructs] = getKKTMultiCastPrecoders(SimParams,SimStructs);
        
    case 'SDPASMethod'
        
        [SimParams,SimStructs] = getMultiCastSDPAS(SimParams,SimStructs,50);
        
    case 'ConicASMethod'
        
        searchType = 'FC';
        switch searchType
            case 'SDP'
                [SimParams,SimStructs] = getMultiCastSDP(SimParams,SimStructs,0,1);
            case 'FC'
                [SimParams,SimStructs] = getMultiCastConicA(SimParams,SimStructs,searchType);
            case 'Dual'
                [SimParams,SimStructs] = getMultiCastConicA(SimParams,SimStructs,searchType);
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
        [SimParams,SimStructs] = getMultiCastConicAS(SimParams,SimStructs);
        display('Antenna subset selected !');
        [SimParams,SimStructs] = getMultiCastConicA(SimParams,SimStructs,'MP');
        
    case 'ConicMethodB'
        
        searchType = 'FC';
        switch searchType
            case 'FC'
                [SimParams,SimStructs] = getMultiCastConicB(SimParams,SimStructs,searchType);
            case 'Dual'
                [SimParams,SimStructs] = getMultiCastConicB(SimParams,SimStructs,searchType);
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
        [SimParams,SimStructs] = getMultiCastConicB(SimParams,SimStructs,'MP');
        
    case 'ConicBSMethod'
        
        searchType = 'Dual';
        switch searchType
            case 'FC'
                [SimParams,SimStructs] = getMultiCastConicB(SimParams,SimStructs,searchType);
            case 'Dual'
                [SimParams,SimStructs] = getMultiCastConicB(SimParams,SimStructs,searchType);
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
        [SimParams,SimStructs] = getMultiCastConicBS(SimParams,SimStructs);
        display('Antenna subset selected !');
        [SimParams,SimStructs] = getMultiCastConicB(SimParams,SimStructs,'MP');
        
    otherwise
        display('Unknown Precoding Method !');
        
end

end
