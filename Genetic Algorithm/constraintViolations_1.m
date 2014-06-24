
function geneArray = constraintViolations_1(GenStruct,geneArray)

minConst = GenStruct.constraint_arguments(1,1);
maxConst = GenStruct.constraint_arguments(1,2);

for iPhenon = 1:GenStruct.xPop
    for iGene = 1:GenStruct.nGenes
        
        while 1
            if sum(geneArray{iPhenon,1}(iGene,:)) < minConst
                rIndices = randi([1,GenStruct.xLen],1,1);
                geneArray{iPhenon,1}(iGene,rIndices) = 1;
            else
                break;
            end
        end

        while 1
            if sum(geneArray{iPhenon,1}(iGene,:)) > maxConst
                oneLocs = find(geneArray{iPhenon,1}(iGene,:) == 1);
                rIndices = randi([1,length(oneLocs)],1,1);
                geneArray{iPhenon,1}(iGene,oneLocs(1,rIndices)) = 0;
            else
                break;
            end
        end
        
    end
end 
    
end
