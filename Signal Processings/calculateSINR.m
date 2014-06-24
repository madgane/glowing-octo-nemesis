
function gammaG = calculateSINR(Harray,Pmatrix,Wvec,iUser,Pvector)

[~,nUsers] = size(Pmatrix);

if (nargin == 3)
    Pvector = ones(nUsers,1);
end

Npwr = 0;
H = Harray{iUser,1};
Spwr = (Wvec' * H * Pmatrix(:,iUser))' * (Wvec' * H * Pmatrix(:,iUser)) * Pvector(iUser,1);

for kUser = 1:nUsers
    if kUser ~= iUser
        Npwr = Npwr + (Wvec' * H * Pmatrix(:,kUser))' * (Wvec' * H * Pmatrix(:,kUser)) * Pvector(kUser,1);        
    end
end

Npwr = Npwr + trace(Wvec * Wvec');
gammaG = (Spwr / Npwr);

end

    