function [xChromozome , yChromozome] = xOverFunction(GenStruct,newGenChromozomes)

randProb = rand(1,GenStruct.xPop) .* GenStruct.xFitness;
[~,sortI] = sort(randProb,'descend');

aChromozome = newGenChromozomes{sortI(1,1),1};
bChromozome = newGenChromozomes{sortI(1,2),1};
xChromozome = aChromozome;yChromozome = bChromozome;

if GenStruct.pXover < rand(1,1)
    return;
end    

xChromozome(:,GenStruct.XoverPoints) = bChromozome(:,GenStruct.XoverPoints);
yChromozome(:,GenStruct.XoverPoints) = aChromozome(:,GenStruct.XoverPoints);
