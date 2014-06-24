
clc;
clear all;close all;

Nt = 4;
users_iterated = 10:10:500;
analysis_gmac = zeros(length(users_iterated),4);


for iIndex = 1:length(users_iterated)
    
    K = users_iterated(1,iIndex);
    analysis_gmac(iIndex,1) = (Nt^2 + K * Nt^2 + Nt^2 * (Nt - 1) + sum((1:(Nt-1)).^3)) * Nt - 1;

    analysis_gmac(iIndex,2) = sum(Nt * (1:(Nt))) * K + Nt * Nt;
    
    Kr = 0.5 * K;
    analysis_gmac(iIndex,3) = (Nt^2 + Kr * Nt^2 + Nt^2 * (Nt - 1) + sum((1:(Nt-1)).^3)) * Nt - 1;
    
    analysis_gmac(iIndex,4) = sum(Nt * (1:(Nt))) * Kr + Nt * Nt;
    
end

plot(users_iterated,analysis_gmac * 2)