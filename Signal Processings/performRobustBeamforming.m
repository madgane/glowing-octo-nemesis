
function [optFeasibility,optPrecoders] = performRobustBeamforming(SimParams,SimStructs,robustParams)

pickUsers = SimStructs.baseStruct{robustParams.cBase,1}.assignedUsers{robustParams.cBand,1};
pickStreams = SimStructs.baseStruct{robustParams.cBase,1}.assignedStreams{robustParams.cBand,1};
kUsers = length(pickUsers);

H = cell(kUsers,1);
G = cell(kUsers,1);
Q = cell(kUsers,1);

for iUser = 1:kUsers
    cUser = pickUsers(iUser,1);cStream = pickStreams(iUser,1);
    Ch = SimStructs.linkChan{robustParams.cBase,robustParams.cBand}(:,:,cUser);
    [U,D,V] = svd(Ch); H{iUser,1} = U(:,cStream)' * D(cStream,cStream) * V(:,cStream);
    G{iUser,1} = robustParams.cGamma;% * (2^(SimStructs.userStruct{cUser,1}.trafficStats.backLogPkt) - 1);
    Q{iUser,1} = eye(SimParams.nTxAntenna) * 1e10;
end

cvx_quiet(false);

cvx_begin sdp

variable X(SimParams.nTxAntenna * kUsers,SimParams.nTxAntenna * kUsers) complex
variable lambda(kUsers,1)
expression U

minimize 0

subject to

for iUser = 1:kUsers
    
    for jUser = 1:kUsers
        if jUser ~= iUser
            sIndex = (jUser - 1) * SimParams.nTxAntenna + 1;
            eIndex = sIndex + SimParams.nTxAntenna - 1;
            U = X(sIndex:eIndex,sIndex:eIndex) + U;
        end
    end
    
    sIndex = (iUser - 1) * SimParams.nTxAntenna + 1;
    eIndex = sIndex + SimParams.nTxAntenna - 1;    
    U = X(sIndex:eIndex,sIndex:eIndex) - U * G{iUser,1};
    
    U11 = U + lambda(iUser,1) * Q{iUser,1};
    U12 = U * H{iUser,1};U21 = H{iUser,1}' * U';
    U22 = H{iUser,1}' * U * H{iUser,1} - lambda(iUser,1) - SimParams.N * G{iUser,1};
    
    [U11 U12 ; U21 U22] == hermitian_semidefinite(SimParams.nTxAntenna + 1);
    
end

lambda >= 0;
real(trace(X)) <= SimStructs.baseStruct{iBase,1}.sPower(1,iBand);

cvx_end

optFeasibility = cvx_status;
optPrecoders = zeros(SimParams.nTxAntenna,kUsers);

if strcmp(optFeasibility,'Solved')
    X = full(X);
    for iUser = 1:kUsers
        sIndex = (iUser - 1) * SimParams.nTxAntenna + 1;
        eIndex = sIndex + SimParams.nTxAntenna - 1;    
        [P,D] = eig(X(sIndex:eIndex,sIndex:eIndex));
        optPrecoders(:,iUser) = P * sqrt(D) * randn(SimParams.nTxAntenna,1);
    end
end

end




