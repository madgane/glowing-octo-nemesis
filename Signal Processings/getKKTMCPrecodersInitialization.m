function [SimParams, SimStructs] = getKKTMCPrecodersInitialization(SimParams,SimStructs)

initMultiCastVariables;

minPower = -10;
iterateSCA = 1;
iIterateSCA = 0;
iterateSCAMax = 5;

dualLambda = 10;
X = cell(nBases,nBands);
Xt = cell(nBases,nBands);
for iBand = 1:nBands
    for iBase = 1:nBases
        X{iBase,iBand} = complex(zeros(SimParams.nTxAntenna,nGroupsPerCell(iBase,1)),zeros(SimParams.nTxAntenna,nGroupsPerCell(iBase,1)));
        Xt{iBase,iBand} = complex(randn(SimParams.nTxAntenna,nGroupsPerCell(iBase,1)),randn(SimParams.nTxAntenna,nGroupsPerCell(iBase,1)));
    end
end

while iterateSCA

    for iBand = 1:nBands
        for iBase = 1:nBases
            for iGroup = 1:nGroupsPerCell(iBase,1)
                tempExpression = 0;
                groupUsers = SimStructs.baseStruct{iBase,1}.mcGroup{iGroup,1};
                for iUser = 1:length(groupUsers)
                    cUser = groupUsers(iUser,1);
                    Hsdp = cH{iBase,iBand}(:,:,cUser);
                    tempExpression = tempExpression + (Hsdp' * Hsdp) * Xt{iBase,iBand}(:,iGroup) * dualLambda;
                end
                tempSum = eye(SimParams.nTxAntenna);
                for jBase = 1:nBases
                    for jGroup = 1:nGroupsPerCell(jBase,1)
                        if ~and((iBase == jBase),(iGroup == jGroup))
                            ifGroupUsers = SimStructs.baseStruct{jBase,1}.mcGroup{jGroup,1};
                            for jUser = 1:length(ifGroupUsers)
                                ifjUser = ifGroupUsers(jUser,1);
                                Hsdp = cH{iBase,iBand}(:,:,ifjUser);
                                tempSum = tempSum + reqSINRPerUser(ifjUser,1) * dualLambda * (Hsdp' * Hsdp);
                            end
                        end
                    end
                end
                X{iBase,iBand}(:,iGroup) = tempSum \ tempExpression;
            end
        end
    end
        
    Xt = X;
    if iIterateSCA < iterateSCAMax
        iIterateSCA = iIterateSCA + 1;
    else
        iterateSCA = 0;
    end
    
    dX = cell2mat(X);
    objPower = norm(dX(:),2);
    display(objPower);
    if abs(objPower - minPower) < 1e-8
        iterateSCA = 0;
    else
        minPower = objPower;
    end
end

for iBase = 1:nBases
    for iBand = 1:nBands
        SimStructs.baseStruct{iBase,1}.PG{iBand,1} = X{iBase,iBand};
    end
end

end