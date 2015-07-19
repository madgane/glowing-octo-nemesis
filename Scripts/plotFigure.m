
function plotFigure(figStruct)

rng('shuffle');
figMarkerCell = {'s','o','*','v','d','+','p','h','>','^','<','.','x'};
xMarker = randi(1000,1);xMarker = mod(xMarker - 1,length(figMarkerCell)) + 1;


figColor = rand(1,3);
figMarker = figMarkerCell{1,xMarker};
figLineType = '-.';
figLineWidth = 1;
figMarkerSize = 2;
    
if ~isfield(figStruct,'N')
    figStruct.N = 1;
end

if ~isfield(figStruct,'X')
    figStruct.X = 1:length(figStruct.Y);
end

if ~isfield(figStruct,'P')
    figStruct.P = 'plot';
end

figure(figStruct.N);hold all;grid on;box on;

switch figStruct.P
    
    case 'plot'
        
        plot(figStruct.X,figStruct.Y,'Color',figColor,'LineWidth',figLineWidth,...
            'LineStyle',figLineType,'MarkerFaceColor',figColor,'Marker',figMarker,'MarkerSize',figMarkerSize);
        
    case 'semilogy'
        
        semilogy(figStruct.X,figStruct.Y,'Color',figColor,'LineWidth',figLineWidth,...
            'LineStyle',figLineType,'MarkerFaceColor',figColor,'Marker',figMarker,'MarkerSize',figMarkerSize);
        
    case 'cdfplot'
        
        [yValues, xValues] = cdfcalc(figStruct.Y);        
        stairs(xValues,yValues(2:end),'Color',figColor,'LineWidth',figLineWidth,...
            'LineStyle',figLineType,'MarkerFaceColor',figColor,'Marker',figMarker,'MarkerSize',figMarkerSize);

end
        