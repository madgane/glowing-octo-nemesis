function [xChromozome , yChromozome] = xMutationFunction(GenStruct,xChromozome,yChromozome)

betaOne = 1.2;betaTwo = 10;
pMutation = (betaOne + betaTwo * (GenStruct.xGeneFitnessMean / GenStruct.xGeneFitnessDeviation))^(-1);

prbMutationX = rand(GenStruct.nGenes,GenStruct.xLen);prbMutationY = rand(GenStruct.nGenes,GenStruct.xLen);
prbMutationX = prbMutationX > pMutation;prbMutationY = prbMutationY > pMutation;

if GenStruct.pMutation < rand(1,1)
    return;
end

xChromozome = xor(xChromozome,prbMutationX);
yChromozome = xor(yChromozome,prbMutationY);

end
