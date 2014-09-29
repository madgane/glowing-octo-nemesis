
clc;
clear all;close all;

uiopen('Y:\Latex\glowing-octo-meme\Figures\Linux\fig-9-2.fig',1);
H = get(gca);F = get(gcf);
figStruct = get(get(gca,'Children'));

lineCell = {'-','-.'};
figColorCell = {[0,0,0.75],[0,0.75,0],[0.75,0,0],[1,0,1],[0.25,0.25,0.5],[0.15,0.45,0.15],[0.75,0.5,0.35]};
markerCell = {'o','s','d','h','v','*','^','p','>'};

figure(2);
figMarkerSize = 4;
hold all;grid on;box on;

nMarkers = 12;
nCurves = length(figStruct);
legendString = cell(nCurves,1);
xCurveLegend = zeros(nCurves,1);
for iCurve = 1:nCurves
    nJump = nCurves + iCurve;
    jCurve = nCurves - iCurve + 1;
    randInterval = randperm(length(figStruct(iCurve).XData),nMarkers);
    figColor = figColorCell{mod(iCurve - 1,length(figColorCell)) + 1};
    lineType = lineCell{mod(iCurve - 1,length(lineCell)) + 1};
    marker = markerCell{mod(iCurve - 1,length(markerCell)) + 1};
    plot(figStruct(iCurve).XData,figStruct(iCurve).YData,'Color',figColor,'LineWidth',figStruct(iCurve).LineWidth,...
        'LineStyle',lineType,'MarkerFaceColor',figColor,'Marker','none','MarkerSize',figMarkerSize);
    plot(figStruct(iCurve).XData(randInterval),figStruct(iCurve).YData(randInterval),'Color',figColor,'LineWidth',figStruct(iCurve).LineWidth,...
        'LineStyle','none','MarkerFaceColor',figColor,'Marker',marker,'MarkerSize',figMarkerSize);
    xCurveLegend(jCurve,1) = plot(figStruct(iCurve).XData(randInterval(1,1)),figStruct(iCurve).YData(randInterval(1,1)),'Color',figColor,'LineWidth',figStruct(iCurve).LineWidth,...
        'LineStyle',lineType,'MarkerFaceColor',figColor,'Marker',marker,'MarkerSize',figMarkerSize);
    legendString{jCurve,1} = figStruct(iCurve).DisplayName;
end

legend(xCurveLegend,legendString);
xlabel('SCA Update Points');
ylabel('Queue deviation (\chi) in bits / channel use');

saveas(gcf,'Y:\Latex\glowing-octo-meme\Figures\Linux\fig-9-2-2.eps');
