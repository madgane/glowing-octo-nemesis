
function W = calculateMMSEWvector(Harray,Pmatrix,iUser)

Heq = 0;
[~,nUsers] = size(Pmatrix);
[nReceive,~] = size(Harray{iUser,1});
W = (Harray{iUser,1} * Pmatrix(:,iUser))';
for kUser = 1:nUsers
    Xval = (Harray{kUser,1} * Pmatrix(:,kUser));
    Heq = Heq + Xval * Xval';
end

W = W / (Heq + eye(nReceive));
W = W' / norm(W);

end
