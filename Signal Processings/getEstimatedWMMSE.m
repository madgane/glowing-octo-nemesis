
function W = getEstimatedWMMSE(H,interH,intraP,interP,iUser)

[~,nUsers] = size(intraP);
[nNeighbors,~] = size(interH);

P = intraP(:,iUser);
Heq = (H * P) * (H * P)';
    
for jUser = 1:nUsers
    if iUser ~= jUser
        Heq = Heq + (H * intraP(:,jUser)) * (H * intraP(:,jUser))';
    end
end

for jBase = 1:nNeighbors
    [~,jbUsers] = size(interP{jBase,1});
    for kUser = 1:jbUsers
        Heq = Heq + (interH{jBase,1} * interP{jBase,1}(:,kUser)) * (interH{jBase,1} * interP{jBase,1}(:,kUser))';
    end
end

[nRX,~] = size(H);
W = (H * P)' / (Heq + eye(nRX));W = W / norm(W);

end

