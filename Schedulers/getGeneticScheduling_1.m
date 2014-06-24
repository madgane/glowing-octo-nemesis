function [SimParams,SimStructs] = getGeneticScheduling_1(SimParams,SimStructs)

GenStruct.xGen = 5;
GenStruct.xPop = 10;
GenStruct.xLen = SimParams.nUsers;
GenStruct.nGenes = SimParams.nBases;

GenStruct.pXover = 1;
GenStruct.pMutation = 1;
GenStruct.pInversion = 0.25;
GenStruct.XoverPoints = 1:2:10;

GenStruct.globalXInit = @initializeFunction_1;
GenStruct.geneConstFunc = @constraintViolations_1;
GenStruct.geneFitnessFunc = @fitnessFunction_1;

GenStruct.SimParams = SimParams;
GenStruct.SimStructs = SimStructs;
GenStruct.constraint_arguments = [1,SimParams.nTxAntenna];

for iBand = 1:SimParams.nBands
    
    GenStruct.iBand = iBand;
    bestPhenoType = getGeneticCombination(GenStruct);

    for iBase = 1:SimParams.nBases
        xValue = find(bestPhenoType == 1);
        userIdx = SimStructs.baseStruct{iBase,1}.linkedUsers;
        SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1} = userIdx(xValue);
        SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1} = ones(length(xValue),1);
    end
   
end


