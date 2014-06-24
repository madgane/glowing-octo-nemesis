
clc;
clear all;

Nt = 4;Nk = 10;

augH = [];
H = complex(randn(Nk,Nt),randn(Nk,Nt));

for iK = 1:Nt
    
    if iK == 1
        projMatrix = eye(Nt);
    else
        projMatrix = eye(Nt) - augH' * inv(augH * augH') * augH;
    end
    
    
    for iN = 1:Nk
        projValue(iN,1) = norm(H(iN,:) * projMatrix);
    end
    
    [~,sortI] = sort(projValue,'descend');
    augH = [augH ; H(sortI(1,1),:)];
    
end



