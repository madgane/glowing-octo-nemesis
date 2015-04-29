function [U1, D1, V1] = getIterateSVD(X,nIter)

%X = X';
iIter = 1;
[nR,nC] = size(X);
Q = cell(nIter,1);R = cell(nIter,1);

if 1

while iIter <= nIter
    
    if iIter == 1
        [Q{iIter,1}, R{iIter,1}] = qr(X);
    else
        [Q{iIter,1}, R{iIter,1}] = qr(R{iIter - 1,1}');
    end
    iIter = iIter + 1;
end

U = eye(nC);V = eye(nR);
for iIter = 1:nIter
    if mod(iIter - 1,2) == 0
        V = V * Q{iIter,1};
    else
        U = Q{iIter,1}' * U;
    end
end

D = R{nIter,1};
V1 = V;U1 = U';D1 = D;
    
end

U1 = eye(nR,nR);
V1 = eye(nC,nC);

A = X;
for i = 1:nIter
		
	[Q,R] = qr(A);
    Q = Q';R = R';
	U1 = Q * U1;
	A = R;
    
	[Q,R] = qr(A);
    Q = Q';R = R';
	V1 = Q * V1;
	A = R;

end

D1 = A;

end




