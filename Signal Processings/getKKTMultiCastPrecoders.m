function [SimParams, SimStructs] = getKKTMultiCastPrecoders(SimParams,SimStructs)

initMultiCastVariables;

minPower = -10;
iterateSCA = 1;
iIterateSCA = 0;
iterateSCAMax = 100;

X = cell(nBases,nBands);
Xt = cell(nBases,nBands);
for iBand = 1:nBands
    for iBase = 1:nBases
        X{iBase,iBand} = complex(zeros(SimParams.nTxAntenna,nGroupsPerCell(iBase,1)),zeros(SimParams.nTxAntenna,nGroupsPerCell(iBase,1)));
        Xt{iBase,iBand} = complex(randn(SimParams.nTxAntenna,nGroupsPerCell(iBase,1)),randn(SimParams.nTxAntenna,nGroupsPerCell(iBase,1)));
    end
end

dualGroup = zeros(nUsers,nBands);

while iterateSCA

    stepSize = 1e-2;
    for dualIterate = 1:100
        for iBand = 1:nBands
            for iBase = 1:nBases
                for iGroup = 1:nGroupsPerCell(iBase,1)
                    tempExpression = 0;
                    groupUsers = SimStructs.baseStruct{iBase,1}.mcGroup{iGroup,1};
                    for iUser = 1:length(groupUsers)
                        cUser = groupUsers(iUser,1);
                        Hsdp = cH{iBase,iBand}(:,:,cUser);
                        tempExpression = tempExpression + (Hsdp' * Hsdp) * Xt{iBase,iBand}(:,iGroup) * dualGroup(cUser,iBand);
                    end
                    tempSum = eye(SimParams.nTxAntenna);
                    for jBase = 1:nBases
                        for jGroup = 1:nGroupsPerCell(jBase,1)
                            if ~and((iBase == jBase),(iGroup == jGroup))
                                ifGroupUsers = SimStructs.baseStruct{jBase,1}.mcGroup{jGroup,1};
                                for jUser = 1:length(ifGroupUsers)
                                    ifjUser = ifGroupUsers(jUser,1);
                                    Hsdp = cH{iBase,iBand}(:,:,ifjUser);
                                    tempSum = tempSum + reqSINRPerUser(ifjUser,1) * dualGroup(ifjUser,iBand) * (Hsdp' * Hsdp);
                                end
                            end
                        end
                    end
                    X{iBase,iBand}(:,iGroup) = pinv(tempSum) * tempExpression;
                end
            end
        end
        
        for iBand = 1:nBands
            for iBase = 1:nBases
                for iGroup = 1:nGroupsPerCell(iBase,1)
                    groupUsers = SimStructs.baseStruct{iBase,1}.mcGroup{iGroup,1};
                    for iUser = 1:length(groupUsers)
                        cUser = groupUsers(iUser,1);
                        Hsdp = cH{iBase,iBand}(:,:,cUser);
                        tempExpression = abs(Hsdp * Xt{iBase,iBand}(:,iGroup))^2 + 2 * real(Xt{iBase,iBand}(:,iGroup)' * Hsdp' * Hsdp * (X{iBase,iBand}(:,iGroup) - Xt{iBase,iBand}(:,iGroup)));
                        tempSum = SimParams.N;
                        for jBase = 1:nBases
                            for jGroup = 1:nGroupsPerCell(jBase,1)
                                if ~and((iBase == jBase),(iGroup == jGroup))
                                    Hsdp = cH{jBase,iBand}(:,:,cUser);
                                    tempSum = tempSum + abs(Hsdp * X{jBase,iBand}(:,jGroup))^2;
                                end
                            end
                        end
                        dualGroup(cUser,iBand) = dualGroup(cUser,iBand) + stepSize * (tempSum * reqSINRPerUser(cUser,1) - tempExpression);
                        dualGroup(cUser,iBand) = max(0,dualGroup(cUser,iBand));
                    end
                end
            end
        end
        
        dX = cell2mat(X);
        totalPower = dX(:)'*dX(:);
        display(totalPower);
       
    end
    
    Xt = X;
    if iIterateSCA < iterateSCAMax
        iIterateSCA = iIterateSCA + 1;
    else
        iterateSCA = 0;
    end
    
    dX = cell2mat(X);
    objPower = norm(dX(:),2);
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