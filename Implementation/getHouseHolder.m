% function [Q,R] = getHouseHolder(A)
% % Compute the QR decomposition of an m-by-n matrix A using
% % Householder transformations.
% [m,n] = size(A);
% R = A; % Transformed matrix so far
% 
%     for j = 1:min(m,n)-1
%         
%         I = eye((m+1-j));
%         Q(:,:,j) = eye(m); 
%         
%         normx = norm(R(j:end,j));
%         
%         % u = x + ae1 
%         % a = -e^(i*argx(k)) * ||x||
%         u = R(j:end,j);
%         u(1) = u(1)- exp(i*angle(R(j,j))) * normx;
%         %v = u/||u||
%         v = u/norm(u);
%         
%         %Q = I -(1 + (x'v)/(v'x))*vv'   
%         Q(j:end,j:end,j) = I -(1+((u'*v)/(v'*u)))*v*v';
%         R(j:end,j:end) = Q(j:end,j:end,j)*R(j:end,j:end);
%     end
%         Q_tmp = Q(:,:,1);
%         for j = 2:min(m,n)-1
%             Q_tmp = Q_tmp * Q(:,:,j);
%         end
%         Q=Q_tmp;
%         R = (Q' * A);          
% end


function [Q,R] = getHouseHolder(A)

R = A;
[M,N] = size(A);

for iCol = 1:min(M,N)-1
    
    q = eye(min(M,N));
    X = R(iCol:end,iCol);
	U = X;
    
    normX = norm(X);
    U(1,1) = U(1,1) - normX;
    V = U;
    
    normU = norm(V)^2;    
    fproductA = (1 + ((X' * V)/(V' * X))) / normU;
    fproductA = (1 + (X'*V)*(X'*V)) / normU;
    
    q(iCol:end,iCol:end) = eye(length(V)) - fproductA * (V * V');
    R = q * R;
    
    if iCol == 1
        Q = q;
    else
        Q = q * Q;
    end
        
end

end


