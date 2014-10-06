
function [SimParams,SimStructs] = initializeSCApoint(SimParams,SimStructs,bsIndex)

switch SimStructs.baseStruct{bsIndex,1}.selectionType
    case {'BF','Ones'}
        precInit = 'BF';
    case 'Last'
        precInit = 'Last';
end

cH = SimStructs.linkChan;
nBands = SimParams.nBands;
linkedUsers = SimStructs.baseStruct{bsIndex,1}.linkedUsers;
kUsers = length(SimStructs.baseStruct{bsIndex,1}.linkedUsers);

W0 = cell(SimParams.nUsers,nBands);
B0 = zeros(SimParams.maxRank,kUsers,nBands);
G0 = zeros(SimParams.maxRank,kUsers,nBands);
addNoise = zeros(SimParams.maxRank,SimParams.nUsers,nBands);
M0 = zeros(SimParams.nTxAntenna,SimParams.maxRank,kUsers,nBands);

for iBand = 1:nBands
    for iUser = 1:SimParams.nUsers
        W0{iUser,iBand} = SimStructs.userStruct{iUser,1}.pW{iBand,1};
    end
end

switch precInit
    
    case 'BF'
        for iBand = 1:nBands
            for iUser = 1:kUsers
                cUser = linkedUsers(iUser,1);
                [~,~,V] = svd(cH{bsIndex,iBand}(:,:,cUser));
                M0(:,:,iUser,iBand) = V(:,1:SimParams.maxRank);
            end
        end
        
    case 'Ones'
        for iBand = 1:nBands
            M0(:,:,:,iBand) = complex(ones(SimParams.nTxAntenna,SimParams.maxRank,kUsers),ones(SimParams.nTxAntenna,SimParams.maxRank,kUsers));
        end
        
    case 'Last'
        for iBand = 1:nBands
            M0(:,:,:,iBand) = SimParams.Debug.globalExchangeInfo.P{bsIndex,iBand};
            for iUser = 1:kUsers
                cUser = linkedUsers(iUser,1);
                for jBase = 1:SimParams.nBases
                    if jBase ~= bsIndex
                        addNoise(:,cUser,iBand) = SimParams.Debug.globalExchangeInfo.gI{jBase,1}(:,cUser,iBand).^2;
                    end
                end
            end
        end
end

for iBand = 1:nBands
    for iUser = 1:kUsers
        cUser = linkedUsers(iUser,1);
        for iLayer = 1:SimParams.maxRank
            B0(iLayer,iUser,iBand) = SimParams.N * norm(W0{cUser,iBand}(:,iLayer))^2 + addNoise(iLayer,cUser,iBand);
            for jUser = 1:kUsers
                if iUser ~= jUser
                    B0(iLayer,iUser,iBand) = B0(iLayer,iUser,iBand) + norm(W0{cUser,iBand}(:,iLayer)' * cH{bsIndex,iBand}(:,:,cUser) * M0(:,:,jUser,iBand),2)^2;
                else
                    B0(iLayer,iUser,iBand) = B0(iLayer,iUser,iBand) + norm(W0{cUser,iBand}(:,iLayer)' * cH{bsIndex,iBand}(:,:,cUser) * M0(:,iLayer~=1:SimParams.maxRank,iUser,iBand),2)^2;
                end
            end
            G0(iLayer,iUser,iBand) = norm(W0{cUser,iBand}(:,iLayer)' * cH{bsIndex,iBand}(:,:,cUser) * M0(:,iLayer,iUser,iBand),2)^2;
        end
    end
end

G0 = G0 ./ B0;
T0 = log2(1 + G0);

SimParams.Debug.globalExchangeInfo.funcOut{1,bsIndex} = M0;
SimParams.Debug.globalExchangeInfo.funcOut{2,bsIndex} = B0;
SimParams.Debug.globalExchangeInfo.funcOut{3,bsIndex} = G0;
SimParams.Debug.globalExchangeInfo.funcOut{4,bsIndex} = T0;
SimParams.Debug.globalExchangeInfo.funcOut{5,bsIndex} = W0;

end

