function [xChromozome , yChromozome] = xInversionFunction(GenStruct,xChromozome,yChromozome)

xLocs = 2;

if GenStruct.pInversion < rand(1,1)
    return;
end

indInversionX = sort(randi([1,GenStruct.xLen],1,xLocs));
indInversionY = sort(randi([1,GenStruct.xLen],1,xLocs));

uChromozome = xChromozome(:,indInversionX(1,1):indInversionX(1,2));
vChromozome = yChromozome(:,indInversionY(1,1):indInversionY(1,2));
xChromozome(:,indInversionX(1,1):indInversionX(1,2)) = fliplr(uChromozome);
yChromozome(:,indInversionY(1,1):indInversionY(1,2)) = fliplr(vChromozome);

end
