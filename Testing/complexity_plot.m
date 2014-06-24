
clc;
clear all;

Nt = 4;
usersOverSystem = [Nt:100];
itUsers = length(usersOverSystem);

sumCAdd = zeros(itUsers,4);
sumCMult = zeros(itUsers,4);

% SUS algorithm using Nullspace

gIndex = 1;

for itUser = 1:itUsers
    for iAntenna = 1:Nt
        
        if iAntenna == 1
            for kUser = 1:usersOverSystem(itUser)
                
                prodPerUser = Nt;
                sumPerUser = Nt - 1;
                
                sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser;
                sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser;
                
            end
        else       
            for kUser = 1:(usersOverSystem(itUser) - (iAntenna - 1))
                
                % For Matrix Multiplication
                
                prodPerUser = Nt * Nt + Nt;
                sumPerUser = Nt * (Nt - 1) + Nt - 1;
                
                sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser;
                sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser;

            end
        end
        
        if iAntenna ~= Nt

            % Inverse Calculation

            prodPerUser = Nt * Nt * Nt * Nt;
            sumPerUser = (Nt - 1) * Nt * Nt;

            sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser;
            sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser;

            % Null space calculation

            prodPerUser = Nt * iAntenna * iAntenna + Nt * Nt * iAntenna;
            sumPerUser = (iAntenna - 1) * iAntenna + (iAntenna - 1) * Nt;

            sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser;
            sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser;
            
        end

    end
end

gIndex = 2;

% Reduced Nullspace Selection

for itUser = 1:itUsers
    for iAntenna = 1:Nt
        
        if iAntenna == 1
            for kUser = 1:usersOverSystem(itUser)
                
                prodPerUser = Nt;
                sumPerUser = Nt - 1;
                
                sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser;
                sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser;
                
            end
        else       
            for kUser = 1:(usersOverSystem(itUser) - (iAntenna - 1))
                
                % For Matrix Multiplication
                
                prodPerUser = Nt * (iAntenna - 1) + Nt;
                sumPerUser = (iAntenna - 1) * (Nt - 1) + Nt - 1;
                
                sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser;
                sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser;
                
                % Norm calculation for each user
                
                prodPerUser = (iAntenna - 1);sumPerUser = iAntenna;                
                sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser;
                sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser;

            end
        end
        
        if iAntenna ~= Nt

            % Unit Norm Calculation

            prodPerUser = Nt;
            sumPerUser = (Nt - 1);

            sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser;
            sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser;
            
        end

    end
end

gIndex = 3;

% Pairwise User Selection

for itUser = 1:itUsers
    for iAntenna = 1:Nt
        
        sumUsers = usersOverSystem(itUser) - (iAntenna - 1);
        psdUsers = sumUsers * (sumUsers + 1) * 0.5;
        
        % For Pairwise Matrix Multiplication
        
        prodPerUser = Nt * (2 + iAntenna - 1) * Nt + ((2 + iAntenna - 1)^3 / 3);
        sumPerUser = (Nt - 1) * (2 + iAntenna - 1);
        
        sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser * psdUsers;
        sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser * psdUsers;
        
        if iAntenna ~= Nt

            % Inverse Calculation

            prodPerUser = Nt * Nt * Nt * Nt;
            sumPerUser = (Nt - 1) * Nt * Nt;

            sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser;
            sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser;

            % Null space calculation

            prodPerUser = Nt * iAntenna * iAntenna + Nt * Nt * iAntenna;
            sumPerUser = (iAntenna - 1) * iAntenna + (iAntenna - 1) * Nt;

            sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser;
            sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser;

            % Dominant Eigen Mode
            
            prodPerUser = sumUsers * sumUsers * sumUsers;
            sumPerUser = (sumUsers - 1) * sumUsers * sumUsers;

            sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser;
            sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser;
            
            
        end

    end
end

% Exhaustive using Nullspace
gIndex = 4;

for itUser = 1:itUsers
    
    totComb = factorial(usersOverSystem(itUser)) / (factorial(Nt) * factorial(usersOverSystem(itUser) - Nt));
       
    prodPerUser = Nt * Nt * Nt + (Nt^3 / 3);
    sumPerUser = (Nt - 1) * Nt;
    
    sumCAdd(itUser,gIndex) = sumCAdd(itUser,gIndex) + sumPerUser * totComb;
    sumCMult(itUser,gIndex) = sumCMult(itUser,gIndex) + prodPerUser * totComb;
    
end

sumPlot = sumCAdd + sumCMult;
plot(usersOverSystem,sumPlot(:,3:4));



