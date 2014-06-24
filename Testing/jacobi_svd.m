function [U,S,V]=jacobi_svd(A)
% [U S V]=jacobi_svd(A)
% A is original square matrix
% Singular values come back in S (diag matrix)
% orig matrix = U*S*V'
%
% One-sided Jacobi algorithm for SVD
% lawn15: Demmel, Veselic, 1989, 
% Algorithm 4.1, p. 32

TOL=1.e-8;
MAX_STEPS=40;

n=size(A,1);
U=[real(A) -imag(A) ; imag(A) real(A)];
V=eye(n);
for steps=1:MAX_STEPS
  converge=0;
  for j=2: 2*n
    for k=1:j-1
      % compute [alpha gamma;gamma beta]=(k,j) submatrix of U'*U
      
      alpha = U(:,k)' *  U(:,k);  %might be more than  1 line
      beta = U(:,j)' *  U(:,j);   %might be more than 1 line
      gamma = U(:,j)' * U(:,k)  %might be more than 1 line
      converge=max(converge,abs(gamma)/sqrt(alpha*beta));
      
      % compute Jacobi rotation that diagonalizes 
      %    [alpha gamma;gamma beta]
      if gamma ~= 0
        zeta=(beta-alpha)/(2*gamma);
        t=sign(zeta)/(abs(zeta)+sqrt(1+zeta^2));
      else
        % if gamma=0, then zeta=infinity and t=0
        t=0;
      end
      c=???
      s=???
      
      % update columns k and j of U
      T=U(:,k);
      U(:,k)=c*T-s*U(:,j);
      U(:,j)=s*T+c*U(:,j);
      
      % update matrix V of right singular vectors

      ???

    end
  end
  if converge < TOL
    break;
  end
end
if steps >= MAX_STEPS
  error('jacobi_svd failed to converge!');
end

% the singular values are the norms of the columns of U
% the left singular vectors are the normalized columns of U
for j=1:n
  singvals(j)=norm(U(:,j));
  U(:,j)=U(:,j)/singvals(j);
end
S=diag(singvals);

