
function globalPlotHold(SimParams,markerS,LW,Col)

if ~isequal(SimParams.userRun,'true')
    
    figure(1);hold on;
    SimParams.sumThrpt = sum(SimParams.Thrpt,2);
    plot(SimParams.snrIndex,SimParams.sumThrpt,markerS,'LineWidth',LW,'Color',Col);
    xlabel('SNR in dB');ylabel('Sum Capacity in Bits/Sec/Hz');grid on;

    figure(2);hold on;
    JainMean = mean(SimParams.Thrpt,2).^2;JainVar = var(SimParams.Thrpt,0,2);
    JainIndex_capacity = JainMean ./ (JainMean + JainVar);

    plot(SimParams.snrIndex,JainIndex_capacity,markerS,'LineWidth',LW,'Color',Col);
    xlabel('SNR in dB');ylabel('Capacity Deviation across Users in Bits/Sec/Hz');grid on;

    figure(3);hold on;
    JainMean = mean(SimParams.fairness,2).^2;JainVar = var(SimParams.fairness,0,2);
    JainIndex_utility = JainMean ./ (JainMean + JainVar);

    plot(SimParams.snrIndex,JainIndex_utility,markerS,'LineWidth',LW,'Color',Col);
    xlabel('SNR in dB');ylabel('Network Utility Deviation across Users');grid on;

else
    
    figure(4);hold on;
    plot(SimParams.userIdices,SimParams.Thrpt,markerS,'LineWidth',LW,'Color',Col);
    xlabel('Number of Users');ylabel('Sum Capacity in Bits/Sec/Hz');grid on;
    hold off;
    
end



end