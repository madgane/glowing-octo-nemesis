
function [SimParams,SimStructs] = getMultiBandSCAE(SimParams,SimStructs,ObjType)

if nargin == 2
    ObjType = '';
    exMinPower = Inf;
else
    exMinPower = -Inf;
end
initMultiCastVariables;

rEX = SimParams.Debug.tempResource{2,1}{1,1};
iEX = SimParams.Debug.tempResource{3,1}{1,1};
bEX = SimParams.Debug.tempResource{4,1}{1,1};

enabledAntennaIndices = getExhaustiveArray(SimParams.nAntennaArray,SimParams.nTxAntennaEnabled);
SimParams.Debug.MultiCastSDPExchange = cell(nBases,nBands);
nCombinations = size(enabledAntennaIndices,1);

if nBands == 1
    
    for iCombination = 1:nCombinations
    
        SimParams.Debug.SCA_initFailureFlag = 0;
        SimParams.Debug.tempResource{2,1}{1,1} = rEX;
        SimParams.Debug.tempResource{3,1}{1,1} = iEX;
        
        xAntennasEnabled = zeros(1,SimParams.nAntennaArray);
        xAntennasEnabled(1,enabledAntennaIndices(iCombination,:)) = 1;
        
        for iBase = 1:nBases
            for iBand = 1:nBands
                SimParams.Debug.MultiCastSDPExchange{iBase,iBand} = logical(xAntennasEnabled);
            end
        end
        
        display('Antenna subset selected as - ');
        display(logical(xAntennasEnabled));
        
        [SimParams,SimStructs] = getSingleBandSCA(SimParams,SimStructs,'Dual');
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
        
        txPower = 0;
        if (SimParams.Debug.SCA_initFailureFlag == 0)            
            [xParams,xStructs] = getSingleBandSCA(SimParams,SimStructs,'MP');
            for iBase = 1:nBases
                M = cell2mat(xStructs.baseStruct{iBase,1}.PG);
                txPower = txPower + norm(M(:))^2;
            end
            
            if (txPower < exMinPower)
                gStructs = xStructs;
                gParams = xParams;
                exMinPower = txPower;
            end
        end
    end
        
else
    
    for iCombination = 1:nCombinations
        
        SimParams.Debug.SCA_initFailureFlag = 0;
        SimParams.Debug.tempResource{2,1}{1,1} = rEX;
        SimParams.Debug.tempResource{3,1}{1,1} = iEX;
        SimParams.Debug.tempResource{4,1}{1,1} = bEX;
        
        xAntennasEnabled = zeros(1,SimParams.nAntennaArray);
        xAntennasEnabled(1,enabledAntennaIndices(iCombination,:)) = 1;
        
        for iBase = 1:nBases
            for iBand = 1:nBands
                SimParams.Debug.MultiCastSDPExchange{iBase,iBand} = logical(xAntennasEnabled);
            end
        end
        
        display('Antenna subset selected as - ');
        display(logical(xAntennasEnabled));
        
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

        txPower = 0;
        if (SimParams.Debug.SCA_initFailureFlag == 0)         
            if strcmpi(ObjType,'MaxMin')
                
                [xParams,xStructs] = getMultiBandSCA(SimParams,SimStructs,'MaxMin');
                
                if (xParams.Debug.minUserRate > exMinPower)
                    gStructs = xStructs;
                    gParams = xParams;
                    exMinPower = xParams.Debug.minUserRate;
                end

            else
                
                [xParams,xStructs] = getMultiBandSCA(SimParams,SimStructs,'MP');
                
                for iBase = 1:nBases
                    M = cell2mat(xStructs.baseStruct{iBase,1}.PG);
                    txPower = txPower + norm(M(:))^2;
                end
                
                if (txPower < exMinPower)
                    gStructs = xStructs;
                    gParams = xParams;
                    exMinPower = txPower;
                end

            end
            
        end
        
    end
    
end

if exist('gParams')
    SimParams = gParams;
    SimStructs = gStructs;
else
    for iBase = 1:nBases
        for iBand = 1:nBands
            SimStructs.baseStruct{iBase,1}.PG{iBand,1} = zeros(SimParams.nTxAntenna,nGroupsPerCell(iBase,1));
        end
    end
end

end