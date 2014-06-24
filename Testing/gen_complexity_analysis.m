
clc;

Nt = 4;
usersOverSystem = [10:10:100];
itUsers = length(usersOverSystem);

sumCAdd = zeros(itUsers,7);
sumCMult = zeros(itUsers,7);

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

% Pairwise User Selection

% for itUser = 1:itUsers
%     for iAntenna = 1:Nt
%         
%         sumUsers = usersOverSystem(itUser) - (iAntenna - 1);
%         psdUsers = sumUsers * (sumUsers + 1) * 0.5;
%         
%         % For Pairwise Matrix Multiplication
%         
%         [sumPerUser,prodPerUser] = est_flops([(1 + iAntenna),Nt],[Nt,(1 + iAntenna)]);
%         
%         prodPerUser = prodPerUser + ((2 + iAntenna - 1)^3 / 3);
%         sumPerUser = sumPerUser + ((1 + iAntenna)^2 / 3);
%         
%         sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser * psdUsers;
%         sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser * psdUsers;
%         
%         if iAntenna < Nt
%             
%             % Inverse Calculation
%             
%             % Inner product Nt x iAntenna
%             [sumPerUser,prodPerUser] = est_flops([iAntenna,Nt],[Nt,iAntenna],0);
%             sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser;
%             sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser;
%             
%             % Inverse of iAntenna x iAntenna
%             [sumPerUser,prodPerUser] = est_flops([iAntenna,iAntenna],[iAntenna,iAntenna],0);
%             sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser * iAntenna * iAntenna;
%             sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser * iAntenna * iAntenna;
%             
%             % Null space calculation
%             
%             % First Product
%             [sumPerUser,prodPerUser] = est_flops([Nt,iAntenna],[iAntenna,iAntenna],0);
%             sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser;
%             sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser;
%             
%             % Second Product
%             [sumPerUser,prodPerUser] = est_flops([Nt,iAntenna],[iAntenna,Nt],0);
%             sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser;
%             sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser;
%             
%             sumCAdd(itUser,gIndex) = Nt * Nt + sumCAdd(itUser,gIndex);
%             
%         end
%         
%     end
% end

gIndex = 4;

% M-AHP 

for itUser = 1:itUsers
    
    for iAntenna = 1:Nt
        
        sumUsers = usersOverSystem(itUser) - (iAntenna + 1);
        lowerTGLusers = sumUsers * (sumUsers + 1) * 0.5 - sumUsers;
        
        % Inner product
        
        [s0,p0] = est_flops([(iAntenna),Nt],[Nt,(iAntenna)],0);
        [s1,p1] = est_flops([(iAntenna + 1),Nt],[Nt,(iAntenna + 1)],0);
        
        s1 = s1 + (iAntenna + 1)^3;s0 = s0 + (iAntenna)^3;
        p1 = p1 + (iAntenna + 1)^3;p0 = p0 + (iAntenna)^3;
        
        if iAntenna == Nt
            sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + s0 * sumUsers;
            sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + p0 * sumUsers;
        else
            sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + s1 * lowerTGLusers + s0 * sumUsers;
            sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + p1 * lowerTGLusers + p0 * sumUsers;
            
            [s,p] = est_flops([sumUsers,sumUsers],[sumUsers,sumUsers]);
            sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + s * 2;
            sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + p * 2;
        end
        
    end
    
end         


gIndex = 5;

% M-AHP 

for itUser = 1:itUsers
    
    for iAntenna = 1:Nt
        
        sumUsers = usersOverSystem(itUser) - (iAntenna + 1);
        lowerTGLusers = sumUsers * (sumUsers + 1) * 0.5 - sumUsers;
        
        % Inner product
        
        [s0,p0] = est_flops([(iAntenna),Nt],[Nt,(iAntenna)],0);
        [s1,p1] = est_flops([(iAntenna + 1),Nt],[Nt,(iAntenna + 1)],0);
        
        s1 = s1 + (iAntenna + 1)^3;s0 = s0 + (iAntenna)^3;
        p1 = p1 + (iAntenna + 1)^3;p0 = p0 + (iAntenna)^3;
        
        if iAntenna == Nt
            sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + s0 * sumUsers;
            sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + p0 * sumUsers;
        else
            sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + s1 * lowerTGLusers + s0 * sumUsers;
            sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + p1 * lowerTGLusers + p0 * sumUsers;
        end
    
    end
    
end         
   
% Exhaustive using Nullspace
gIndex = 6;

for itUser = 1:itUsers
    
    totComb = factorial(usersOverSystem(itUser)) / ...
        (factorial(Nt) * factorial(usersOverSystem(itUser) - Nt));
    
    [sumPerUser,prodPerUser] = est_flops([Nt,Nt],[Nt,Nt],0);
    sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser * Nt^3;
    sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser * Nt^3;
    
    sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser * totComb;
    sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser * totComb;
    
end


gIndex = 7;

for itUser = 1:itUsers
    
    totComb = 0;
    for iAntenna = 1:Nt
        
        totComb = usersOverSystem(itUser) - iAntenna + 1;
        [sumPerUser,prodPerUser] = est_flops([1,Nt],[Nt,1],0);
        
        sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser * totComb;
        sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser * totComb;
        
    end
    
end

sumPlot = sumCAdd + sumCMult;
semilogy(usersOverSystem,sumPlot);
box on;

legendString = cell(7,1);
legendString{1,1} = 'BD Scheduling';legendString{2,1} = 'PIP Selection';
legendString{3,1} = 'PIP Selection (with memory)';legendString{4,1} = 'PWS Scheme-1';
legendString{5,1} = 'PWS Scheme';legendString{6,1} = 'Exhaustive';legendString{7,1} = 'Greedy';
legend(legendString);

xlabel('Number of Users in the System');
ylabel('Number of Complex Operations in log scale');


