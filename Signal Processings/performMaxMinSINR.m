
function [SimParams,SimStructs] = performMaxMinSINR(SimParams,SimStructs)

epsilon = 1e-2;
groupUsers = cell(SimParams.nBases,1);
groupStreams = cell(SimParams.nBases,1);

for iBand = 1:SimParams.nBands
    
    xLength = 0;
    Qgamma = cell(SimParams.nBases);
    
    for iBase = 1:SimParams.nBases
        groupUsers{iBase,1} = SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1}';
        groupStreams{iBase,1} = SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1}';
        xLength = xLength + length(groupUsers{iBase,1});
    end
    
    gammaMax = 1000;gammaMin = 0;
    [SimParams,SimStructs] = updateMMSEreceiver(SimParams,SimStructs,0);
    
    while 1
        
        gamma = (gammaMax + gammaMin) * 0.5;
        
        for iBase = 1:SimParams.nBases
            Qgamma{iBase,1} = ones(1,length(groupUsers{iBase,1})) * sqrt(1 + 1 / gamma);
        end
        
        cvx_begin
        
        cvx_quiet('true');
        
        variable X(SimParams.nTxAntenna,length(groupUsers{1,1})) complex;
        if SimParams.nBases > 1
            variable Y(SimParams.nTxAntenna,length(groupUsers{2,1})) complex;
        end
        
        minimize (0)
        
        subject to
        
        for iBase = 1:SimParams.nBases
            for iUser = 1:length(groupUsers{iBase,1})
                cUser = groupUsers{iBase,1}(1,iUser);
                cStream = groupStreams{iBase,1}(1,iUser);
                H = SimStructs.linkChan{iBase,iBand}(:,:,cUser);
                W = SimStructs.userStruct{cUser,1}.W{iBand,1}(:,:,cStream);
                IFarray = [X' * H' * W];
                for jBase = 1:length(SimStructs.userStruct{cUser,1}.neighNode)
                    ifBase = SimStructs.userStruct{cUser,1}.neighNode(jBase,1);
                    Hk = SimStructs.linkChan{ifBase,iBand}(:,:,cUser);
                    IFarray = [IFarray ; Y' * Hk' * W];
                end
                
                IFarray = [IFarray ; 1];
                if iBase == 1
                    Darray = Qgamma{iBase,1}(1,iUser) * W * H * X(:,iUser);
                else
                    Darray = Qgamma{iBase,1}(1,iUser) * W * H * Y(:,iUser);
                end
                
                {IFarray , Darray} <In> complex_lorentz(xLength + 1);
                
            end
        end
        
        norm(X(:)) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
        
        if SimParams.nBases > 1
            norm(Y(:)) <= sqrt(SimStructs.baseStruct{iBase,1}.sPower(1,iBand));
        end
        
        cvx_end
        
        if strcmp(cvx_status,'Solved')
            if abs(gammaMin - gammaMax) < epsilon
                break;
            end
            gammaMin = gamma;
        else
            gammaMax = gamma;
        end
        
    end
    
    for iBase = 1:SimParams.nBases
        if iBase == 1
            SimStructs.baseStruct{iBase,1}.P{iBand,1} = X;
        else
            SimStructs.baseStruct{iBase,1}.P{iBand,1} = Y;
        end
    end
    
       
end
