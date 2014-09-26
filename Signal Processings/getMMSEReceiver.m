function [varargout] = getMMSEReceiver(SimParams,SimStructs,rxType,bsIndex)

linkedUsers = cell(SimParams.nBases,1);
for iBase = 1:SimParams.nBases
    linkedUsers{iBase,1} = SimStructs.baseStruct{iBase,1}.linkedUsers;
end

W0 = cell(SimParams.nUsers,SimParams.nBands);

for iBand = 1:SimParams.nBands
    for iBase = 1:SimParams.nBases
        if strcmpi('Self',rxType)
            if iBase ~= bsIndex
                continue;
            end
        end
        for iUser = 1:length(linkedUsers{iBase,1})
            cUser = linkedUsers{iBase,1}(iUser,1);
            
            if strcmpi(rxType,'Ones')
                W = ones(SimParams.nRxAntenna,SimParams.maxRank);
                SimStructs.userStruct{cUser,1}.W{iBand,1} = W;
                W0{cUser,iBand} = W;
                continue;
            end

            R = SimParams.N;
            for jBase = 1:SimParams.nBases
                for jUser = 1:length(linkedUsers{jBase,1})
                    switch rxType
                        case 'Self'
                            if (iBase == jBase)
                                H = SimStructs.linkChan{jBase,iBand}(:,:,cUser);
                                R = R + H * SimStructs.baseStruct{iBase,1}.P{iBand,1}(:,:,jUser) * SimStructs.baseStruct{iBase,1}.P{iBand,1}(:,:,jUser)' * H';
                            end
                        case 'Joint'
                            H = SimStructs.linkChan{jBase,iBand}(:,:,cUser);
                            R = R + H * SimStructs.baseStruct{iBase,1}.P{iBand,1}(:,:,jUser) * SimStructs.baseStruct{iBase,1}.P{iBand,1}(:,:,jUser)' * H';
                        otherwise
                            display('Unknown RX Type !');
                    end
                end
            end
            H = SimStructs.linkChan{iBase,iBand}(:,:,cUser);
            W = R \ (H * SimStructs.baseStruct{iBase,1}.P{iBand,1}(:,:,iUser));
            SimStructs.userStruct{cUser,1}.W{iBand,1} = W;
            W0{cUser,iBand} = W;
        end
    end
end

if nargout == 1
    varargout = cell(1,1);
    varargout{1,1} = W0;
else
    varargout = cell(1,2);
    varargout{1,1} = SimParams;
    varargout{1,2} = SimStructs;
end

end

