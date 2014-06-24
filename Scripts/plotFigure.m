
function plotFigure(figStruct)

figColor = 'k';
figMarker = '.';
figLineType = '-';
figLineWidth = 1;
    
if ~isfield(figStruct,'N')
    figStruct.N = 1;
end

if ~isfield(figStruct,'X')
    figStruct.X = 1:length(figStruct.Y);
end

if ~isfield(figStruct,'P')
    figStruct.P = 'plot';
end

figure(figStruct.N);hold all;grid on;

switch figStruct.P
    
    case 'plot'
        
        plot(figStruct.X,figStruct.Y,'Color',figColor,'LineWidth',figLineWidth,...
            'LineStyle',figLineType,'MarkerFaceColor',figColor,'Marker',figMarker);
        
    case 'semilogy'
        
        semilogy(figStruct.X,figStruct.Y,'Color',figColor,'LineWidth',figLineWidth,...
            'LineStyle',figLineType,'MarkerFaceColor',figColor,'Marker',figMarker);
        
    case 'cdfplot'
        
        [yValues, xValues] = cdfcalc(figStruct.Y);        
        stairs(xValues,yValues(2:end),'Color',figColor,'LineWidth',figLineWidth,...
            'LineStyle',figLineType,'MarkerFaceColor',figColor,'Marker',figMarker);

end
        