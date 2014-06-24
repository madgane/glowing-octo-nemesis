
clc;
clear all;

usersOverSystem_actual = [10:10:100];
itUsers = length(usersOverSystem_actual);

sumCAdd = zeros(itUsers,4);
sumCMult = zeros(itUsers,4);

Nt = 4;
Nr = 1;
nPRB = 20;
nSec = 1e3;
nCmplxToRealMult = 4;
nCmplxToRealAdd = 2;

nCalcMult = nPRB * nSec * nCmplxToRealMult;
nCalcAdd = nPRB * nSec * nCmplxToRealAdd;

usersOverSystem = usersOverSystem_actual * Nr;

% SVD complexity

if Nr == 1
    svdComplexOperations = 0;
    svdOperations = repmat((svdComplexOperations * usersOverSystem_actual)',1,4);
else
    svdComplexOperations = 4 * (Nr * 2)^2 * (Nt * 2) + 22 * (Nt * 2)^3;
    svdOperations = repmat((svdComplexOperations * usersOverSystem_actual)',1,4);
end

% SUS algorithm using Nullspace

gIndex = 1;

for itUser = 1:itUsers
    for iAntenna = 1:Nt
        
        if iAntenna == 1
            for kUser = 1:usersOverSystem(itUser)
                
                [sumPerUser,prodPerUser] = est_flops([1,Nt],[Nt,1]);
                sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser;
                sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser;
                
            end
        else
            for kUser = 1:(usersOverSystem(itUser) - (iAntenna - 1))
                
                % For Matrix Multiplication
                
                [sumPerUser,prodPerUser] = est_flops([Nt,Nt],[Nt,1],1);
                sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser;
                sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser;
                
            end
        end
        
        if iAntenna < Nt
            
            % Inverse Calculation
            
            % Inner product Nt x iAntenna
            [sumPerUser,prodPerUser] = est_flops([iAntenna,Nt],[Nt,iAntenna],0);
            sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser;
            sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser;
            
            % Inverse of iAntenna x iAntenna
            [sumPerUser,prodPerUser] = est_flops([iAntenna,iAntenna],[iAntenna,iAntenna],0);
            sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser * iAntenna * iAntenna;
            sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser * iAntenna * iAntenna;
            
            % Null space calculation
            
            % First Product
            [sumPerUser,prodPerUser] = est_flops([Nt,iAntenna],[iAntenna,iAntenna],0);
            sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser;
            sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser;
            
            % Second Product
            [sumPerUser,prodPerUser] = est_flops([Nt,iAntenna],[iAntenna,Nt],0);
            sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser;
            sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser;
            
            sumCAdd(itUser,gIndex) = Nt + sumCAdd(itUser,gIndex);
            
        end
        
    end
end

gIndex = 2;
with_memory = 0;

% Reduced Nullspace Selection

for itUser = 1:itUsers
    for iAntenna = 1:Nt
        
        if iAntenna == 1
            for kUser = 1:usersOverSystem(itUser)
                
                [sumPerUser,prodPerUser] = est_flops([1,Nt],[Nt,1]);
                sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser;
                sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser;
                
            end
        else
            for kUser = 1:(usersOverSystem(itUser) - (iAntenna - 1))
                
                % For Matrix Multiplication
                
                if with_memory
                    
                    [sumPerUser,prodPerUser] = est_flops([1,Nt],[Nt,1],1);
                    sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser + 1;
                    sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser + 1;
                    
                else
                    
                    [sumPerUser,prodPerUser] = est_flops([(iAntenna - 1),Nt],[Nt,1],1);
                    sumPerUser = sumPerUser + (iAntenna - 1);
                    prodPerUser = prodPerUser + (iAntenna - 1);
                    sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser;
                    sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser;
                    
                end
                
            end
        end
        
        if iAntenna < Nt
            
            % Unit Norm Calculation
            
            if with_memory
                
                [sumPerUser,prodPerUser] = est_flops([1,Nt],[Nt,1],0);
                sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser + 1;
                sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser + 1;
                
            else
                
                [sumPerUser,prodPerUser] = est_flops([1,Nt],[Nt,1],0);
                sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + (sumPerUser + 1) * iAntenna;
                sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + (prodPerUser + 1) * iAntenna;
                
            end
            
        end
        
    end
end

gIndex = 3;
with_memory = 1;

% Reduced Nullspace Selection

for itUser = 1:itUsers
    for iAntenna = 1:Nt
        
        if iAntenna == 1
            for kUser = 1:usersOverSystem(itUser)
                
                [sumPerUser,prodPerUser] = est_flops([1,Nt],[Nt,1]);
                sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser;
                sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser;
                
            end
        else
            for kUser = 1:(usersOverSystem(itUser) - (iAntenna - 1))
                
                % For Matrix Multiplication
                
                if with_memory
                    
                    [sumPerUser,prodPerUser] = est_flops([1,Nt],[Nt,1],1);
                    sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser + 1;
                    sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser + 1;
                    
                else
                    
                    [sumPerUser,prodPerUser] = est_flops([(iAntenna - 1),Nt],[Nt,1],1);
                    sumPerUser = sumPerUser + (iAntenna - 1);
                    prodPerUser = prodPerUser + (iAntenna - 1);
                    sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser;
                    sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser;
                    
                end
                
            end
        end
        
        if iAntenna < Nt
            
            % Unit Norm Calculation
            
            if with_memory
                
                [sumPerUser,prodPerUser] = est_flops([1,Nt],[Nt,1],0);
                sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser + 1;
                sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser + 1;
                
            else
                
                [sumPerUser,prodPerUser] = est_flops([1,Nt],[Nt,1],0);
                sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + (sumPerUser + 1) * iAntenna;
                sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + (prodPerUser + 1) * iAntenna;
                
            end
            
        end
        
    end
end

gIndex = 4;
   
% Greedy Approach

for itUser = 1:itUsers
    
    totComb = 0;
    for iAntenna = 1:Nt
        
        totComb = usersOverSystem(itUser) - iAntenna + 1;
        [sumPerUser,prodPerUser] = est_flops([1,Nt],[Nt,1],0);
        
        sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser * totComb;
        sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser * totComb;
        
    end
    
end

searchComplexity = zeros(length(usersOverSystem),4);
for iUser = 1:length(usersOverSystem)
    kUsers = usersOverSystem(1,iUser);
    searchComplexity(iUser,:) = Nt * kUsers;
end

sumPlot = (sumCAdd * nCalcAdd + sumCMult * nCalcMult + svdOperations + searchComplexity) * 1e-9;
semilogy(usersOverSystem_actual,sumPlot);
box on;

legendString = cell(4,1);
legendString{1,1} = 'QR based search';legendString{2,1} = 'Reduced Null Space';
legendString{3,1} = 'Reduced Nullspace (with memory)';legendString{4,1} = 'Greedy';
legend(legendString);

xlabel('Number of Users in the System');
ylabel('Number of G-FLOPS');

for iUser = 1:4
    sumPlot(:,iUser)
end
