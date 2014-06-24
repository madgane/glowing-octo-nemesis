function fitnessGain = fitnessFunction_2(GenStruct,newGenChromozomes)

fitnessGain = zeros(1,GenStruct.xPop);

for iPop = 1:GenStruct.xPop
    for iGene = 1:GenStruct.nGenes
        Z = GenStruct.CXprecoders{iGene,1};
        fitnessM = Z(:,newGenChromozomes{iPop,1}(iGene,:) == 1);
        fitnessGain(1,iPop) = abs(det(fitnessM));
    end
end
