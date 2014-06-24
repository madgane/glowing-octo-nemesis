
clc;
clear all;close all;

load '.\..\..\MATLAB\outFile_doppler_4.mat';
LineVariable = {'-','--','-.',':'};
MarkerVariable = {'o','+','p','d','s','^','>'};
ColorVariable = {'b','g','r','m','c','k'};

LineVariable = {'-','-','-',':',':',':','-.','-.','-.'};
MarkerVariable = {'o','d','p','o','d','p','o','d','p'};
ColorVariable = {'b','g','r','b','g','r','b','g','r'};

for iPlot = 1:gIndex
    
    figLineWidth = 1;
    figColor = ColorVariable{1,mod(iPlot - 1,length(ColorVariable)) + 1};
    figLineType = LineVariable{1,mod(iPlot - 1,length(LineVariable)) + 1};
    figMarker = MarkerVariable{1,mod(iPlot - 1,length(MarkerVariable)) + 1};
    
    xValues = gParams{iPlot,1}.snrIndex;
    
    yValues = sum(gParams{iPlot,1}.Thrpt(:,:,end),2);
    thrptFairness = sum(gParams{iPlot,1}.fairness(:,:,end),2);
    
    figure(1);hold all;
    plot(xValues,yValues,'Color',figColor,'LineWidth',figLineWidth,...
        'LineStyle',figLineType,'MarkerFaceColor',figColor,'Marker',figMarker);
    
    xlabel('SNR in dB');ylabel('sum rate in bits/sec/Hz');
    
    JainMean = mean(gParams{iPlot,1}.Thrpt,2).^2;JainVar = var(gParams{iPlot,1}.Thrpt,0,2);
    yValues = JainMean ./ (JainMean + JainVar);
    
    figure(2);hold all;
    plot(xValues,yValues,'Color',figColor,'LineWidth',figLineWidth,...
        'LineStyle',figLineType,'MarkerFaceColor',figColor,'Marker',figMarker);
    xlabel('SNR in dB');ylabel('Jain Index for user rates');
    
    JainMean = mean(gParams{iPlot,1}.fairness,2).^2;JainVar = var(gParams{iPlot,1}.fairness,0,2);
    yValues = JainMean ./ (JainMean + JainVar);
    
    figure(3);hold all;
    plot(xValues,yValues,'Color',figColor,'LineWidth',figLineWidth,...
        'LineStyle',figLineType,'MarkerFaceColor',figColor,'Marker',figMarker);
    xlabel('SNR in dB');ylabel('Network Utility Deviation across Users');
    
end

legendVariable = cell(gIndex,1);
for iIndex = 1:gIndex
    switch gParams{iIndex,1}.fbFraction
        case 0.00
            currLegend = 'Ideal';
        case 0.25
            currLegend = 'Fb 0.25';
    end
    
    [~,~,currLegend] = getScheduledUsers(gParams{iIndex,1},gStructs{iIndex,1},currLegend);
    
    switch gParams{iIndex,1}.mdpFactor
        case 0
        case 2
            currLegend = sprintf('%s - %s',currLegend,'with MDP');
    end
    
    legendVariable{iIndex,1} = currLegend;
    
end

for iPlot = 1:3
    figure(iPlot);
    legend(legendVariable);
end



