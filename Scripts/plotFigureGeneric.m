
function plotFigureGeneric(varargin)

figColor = 'm';
figMarker = 'o';
figLineType = '-.';
figLineWidth = 1;

switch nargin
    case 1
        yValues = varargin{1,1};
        xValues = 1:length(yValues);
        figPage = 1;figType = 'plot';
        
    case 2
        xValues = varargin{1,1};
        if ischar(varargin{1,2})
            figType = varargin{1,2};
            figPage = 1;
            yValues = 1:length(xValues);
        else
            yValues = varargin{1,2};
            figPage = 1;figType = 'plot';
        end
        
    case 3
        xValues = varargin{1,1};yValues = varargin{1,2};
        if ischar(varargin{1,3})
            figType = varargin{1,3};
            figPage = 1;
        else
            figPage = varargin{1,3};
            figType = 'plot';
        end
        
    case 4
        xValues = varargin{1,1};yValues = varargin{1,2};
        figPage = varargin{1,3};figType = varargin{1,4};
        
    otherwise
        display('invalid argument !');
        
end
    
figure(figPage);hold all;grid on;

switch figType
    
    case 'plot'
        
        plot(xValues,yValues,'Color',figColor,'LineWidth',figLineWidth,...
            'LineStyle',figLineType,'MarkerFaceColor',figColor,'Marker',figMarker);
        
    case 'semilogy'
        
        semilogy(xValues,yValues,'Color',figColor,'LineWidth',figLineWidth,...
            'LineStyle',figLineType,'MarkerFaceColor',figColor,'Marker',figMarker);
        
    case 'cdfplot'
        
        [yValues, xValues] = cdfcalc(xValues);        
        stairs(xValues,yValues(2:end),'Color',figColor,'LineWidth',figLineWidth,...
            'LineStyle',figLineType,'MarkerFaceColor',figColor,'Marker',figMarker);

end
        