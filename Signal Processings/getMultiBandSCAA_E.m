
function [SimParams,SimStructs] = getMultiBandSCAA_E(SimParams,SimStructs)

initMultiCastVariables;

X = cell(nBands,1);
A = zeros(SimParams.nTxAntenna,1);

for iBase = 1:nBases
    for iBand = 1:nBands
        X{iBand,1} = SimParams.Debug.SCA.RetainX{iBand,1};
    end
    for iAntenna = 1:SimParams.nTxAntenna
        A(iAntenna,1) = A(iAntenna,1) + norm(X{iBand,1}(iAntenna,:))^2;
    end
end

[~,sortI] = sort(A,'descend');

xAntennasEnabled = zeros(1,SimParams.nAntennaArray);
xAntennasEnabled(1,sortI(1:SimParams.nTxAntennaEnabled)) = 1;

for iBase = 1:nBases
    for iBand = 1:nBands
        SimParams.Debug.MultiCastSDPExchange{iBase,iBand} = logical(xAntennasEnabled);
    end
end

display('Antenna subset selected as - ');
display(logical(xAntennasEnabled));

if SimParams.nTxAntennaEnabled ~= SimParams.nTxAntenna
    
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
    
    display('Initialization point and subset found !');
    
end

end
