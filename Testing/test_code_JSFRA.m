
clc;
clear all;

cvx_solver('Mosek');

nBases = 2;
nBands = 2;
nUsers = 8;

nTxAntenna = 4;
nRxAntenna = 1;

userServingBS = [1 1 1 1 2 2 2 2];
H = complex(randn(nRxAntenna,nTxAntenna,nUsers,nBands,nBases),randn(nRxAntenna,nTxAntenna,nUsers,nBands,nBases)) / sqrt(2);

SNR_dB = 10;
currentIterate = 1;
iteratePoints_SCA = 20;

qExponent = 1;
sPower = 10^(SNR_dB/10);
userWeights = ones(nUsers,1);
userQueues = randi([5,10],nUsers,1);

associatedUsers = cell(nBases,1);
for iBase = 1:nBases
    associatedUsers{iBase,1} = find(iBase == userServingBS)';
end

sumRate = zeros(iteratePoints_SCA,1);
queueDev = zeros(iteratePoints_SCA,1);

p_o = 2*rand(nUsers,nBands);
b_o = 10*rand(nUsers,nBands);
q_o = zeros(nUsers,nBands);

while 1
    
    cvx_begin
    
    expressions p(nUsers,nBands) q(nUsers,nBands)
    variable M(nTxAntenna,nUsers,nBands) complex
    variables t(nUsers,nBands) b(nUsers,nBands) g(nUsers,nBands)
    variables userObjective(nUsers,1) epiObjective
    
    minimize(epiObjective)
    
    subject to
    
    for iUser = 1:nUsers
        userWeights(iUser,1) * abs(userQueues(iUser,1) - sum(vec(t(iUser,:)))) <= userObjective(iUser,1);
    end
    
    epiObjective >= norm(userObjective,qExponent);
    
    for iBase = 1:nBases
        for iBand = 1:nBands
            for iUser = 1:length(associatedUsers{iBase,1})
                
                cUser = associatedUsers{iBase,1}(iUser,1);
                intVector = 1;
                
                for jBase = 1:nBases
                    currentH = H(:,:,cUser,iBand,jBase);
                    for jUser = 1:length(associatedUsers{jBase,1})
                        rUser = associatedUsers{jBase,1}(jUser,1);
                        if rUser ~= cUser
                            intVector = [intVector ; currentH * M(:,rUser,iBand)];
                        end
                    end
                end
                
                norm(intVector,2) <= sqrt(b(cUser,iBand));
                log(1 + g(cUser,iBand)) >= t(cUser,iBand) * log(2);
                
                currentH = H(:,:,cUser,iBand,iBase);
                p(cUser,iBand) = real(currentH * M(:,cUser,iBand));
                q(cUser,iBand) = imag(currentH * M(:,cUser,iBand));
                
                (p_o(cUser,iBand)^2 + q_o(cUser,iBand)^2) / (b_o(cUser,iBand)) + ...
                    (2 / b_o(cUser,iBand)) * (p_o(cUser,iBand) * (p(cUser,iBand) - p_o(cUser,iBand))) + ...
                    (2 / b_o(cUser,iBand)) * (q_o(cUser,iBand) * (q(cUser,iBand) - q_o(cUser,iBand))) - ...
                    (p_o(cUser,iBand)^2 + q_o(cUser,iBand)^2) / (b_o(cUser,iBand)^2) * ...
                    (b(cUser,iBand) - b_o(cUser,iBand)) >= g(cUser,iBand);
                
            end
        end
        norm(vec(M(:,associatedUsers{iBase,1},:)),2) <= sqrt(sPower);
    end
    
    cvx_end
    
    sumRate(currentIterate,1) = sum(t(:));
    queueDev(currentIterate,1) = norm(userObjective,1);
    
    if strfind(cvx_status,'Solved')
%         for iBand = 1:nBands
%             for iBase = 1:nBases
%                 for iUser = 1:length(associatedUsers{iBase,1})
%                     cUser = associatedUsers{iBase,1}(iUser,1);
%                     currentH = H(:,:,cUser,iBand,iBase);
%                     p_o(cUser,iBand) = real(currentH * M(:,cUser,iBand));
%                     q_o(cUser,iBand) = imag(currentH * M(:,cUser,iBand));
%                 end
%             end
%         end
%         
%         for iBase = 1:nBases
%             for iBand = 1:nBands
%                 for iUser = 1:length(associatedUsers{iBase,1})
%                     cUser = associatedUsers{iBase,1}(iUser,1);
%                     intVector = 1;
%                     for jBase = 1:nBases
%                         currentH = H(:,:,cUser,iBand,jBase);
%                         for jUser = 1:length(associatedUsers{jBase,1})
%                             rUser = associatedUsers{jBase,1}(jUser,1);
%                             if rUser ~= cUser
%                                 intVector = [intVector ; currentH * M(:,rUser,iBand)];
%                             end
%                         end
%                     end
%                     b_o(cUser,iBand) = norm(intVector)^2;
%                 end
%             end
%         end
        
        p_o = p;
        q_o = q;
        b_o = b;

    else
        p_o = p_o / 2;
        q_o = q_o / 2;
        b_o = b_o * 2;
    end
    
    if currentIterate <= iteratePoints_SCA
        currentIterate = currentIterate + 1;
    else
        break;
    end
    
end

figure(1);
plot(sumRate,'b-o','LineWidth',2);
xlabel('SCA Points');ylabel('Sum Rate in bits/sec/Hz');
title('Network Throughput');
figure(2);
plot(queueDev,'r-o','LineWidth',2);
xlabel('SCA Points');ylabel('Total Queued Packets Remaining in bits');
title('Backlogged Packets in Bits');
