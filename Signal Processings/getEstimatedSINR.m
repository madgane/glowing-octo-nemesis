
function SINR = getEstimatedSINR(H,interH,intraP,interP,iUser)

W = getEstimatedWMMSE(H,interH,intraP,interP,iUser);


[~,nUsers] = size(intraP);
[nNeighbors,~] = size(interH);

P = intraP(:,iUser);
sigPWR = norm(W' * H * P)^2;
noisePWR = 0;
    
for jUser = 1:nUsers
    if iUser ~= jUser
        noisePWR = noisePWR + norm(W' * H * intraP(:,jUser))^2;
    end
end

for jBase = 1:nNeighbors
    [~,nRx] = size(interP{jBase,1});
    for kUser = 1:nRx
        noisePWR = noisePWR + norm(W' * interH{jBase,1} * interP{jBase,1}(:,kUser))^2;
    end
end

noisePWR = noisePWR + trace(W' * W);
SINR = sigPWR / noisePWR;

if isnan(SINR)
    SINR = 0;
end
    

