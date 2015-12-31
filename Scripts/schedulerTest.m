% Testing Scheduling Algorithm

nUsers = 10;
nReceive = 4;
nTransmit = 8;
nDrops = 100;

G = dftmtx(nTransmit);

for iSNR = -5:5:40
    
    H = complex(randn(nReceive,nTransmit,nUsers,nDrops),randn(nReceive,nTransmit,nUsers,nDrops)) / sqrt(2);
    
    for iDrop = 1:nDrops
    
        % Scheduler Algorithm
        
        M = [];
        
        switch algoType
            
            case 'Random'
                
                chosenUsers = randi(1,nUsers,1,nMIMOUsers);
                for xUser = chosenUsers
                    [U,~,~] = svd(H(:,:,xUser,iDrop);
                    M = [M, (U(:,1)' * H(:,:,xUser,iDrop)).'];
                end
                        
        
        % Precoder Design
        
        Z = pinv(M * M') * M';
        
        % Rate calculation
        
        
    
    end
    
end
