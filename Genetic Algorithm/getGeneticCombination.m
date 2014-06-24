function bestPhenotype = getGeneticCombination(GenStruct)

% Population Initialization and Constraint Evaluations
geneArray = GenStruct.globalXInit(GenStruct);
newGenChromozomes = geneArray;

% Fitness Calculation
GenStruct.xFitness = GenStruct.geneFitnessFunc(GenStruct,newGenChromozomes);

for iGen = 1:GenStruct.xGen    
    
    GenStruct.xGeneFitnessMean = mean(GenStruct.xFitness);
    GenStruct.xGeneFitnessDeviation = std(GenStruct.xFitness);
    
    for iParent = 1:2:(GenStruct.xPop - 2)
        [xChromozome , yChromozome] = xOverFunction(GenStruct,newGenChromozomes);
        [xChromozome , yChromozome] = xMutationFunction(GenStruct,xChromozome,yChromozome);
        [xChromozome , yChromozome] = xInversionFunction(GenStruct,xChromozome,yChromozome);
        geneArray{iParent,1} = xChromozome;
        geneArray{iParent + 1,1} = yChromozome;
    end
    
    [~,sortI] = sort(GenStruct.xFitness,'descend');
    geneArray{iParent + 2,1} = newGenChromozomes{sortI(1,1),1};
    geneArray{iParent + 3,1} = newGenChromozomes{sortI(1,2),1};
    
    newGenChromozomes = geneArray;
    newGenChromozomes = GenStruct.geneConstFunc(GenStruct,newGenChromozomes);
    GenStruct.xFitness = GenStruct.geneFitnessFunc(GenStruct,newGenChromozomes);
    
end

[~,maxI] = max(GenStruct.xFitness);
bestPhenotype = newGenChromozomes{maxI,1};
