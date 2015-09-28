
function plotFigure(figStruct)

rng('shuffle');
ovMarker = {'o','v','d','+','*','p','s'};

figColor = rand(1,3);
figMarker = ovMarker{1,randi(length(ovMarker),1,1)};
figLineType = '-.';
figLineWidth = 1;
figMarkerSize = 6;
    
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
        