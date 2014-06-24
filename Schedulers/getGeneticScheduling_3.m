
function [SimParams,SimStructs] = getGeneticScheduling_3(SimParams,SimStructs)

N = SimParams.nBands;
E = SimParams.nBases;
K = SimParams.nUsers;
Nt = SimParams.nTxAntenna;

GPStruct.Ng = 5;
GPStruct.Np = 10;
childPair = cell(2,1);
GPStruct.maxUsers = Nt;
GPStruct.xOverPattern = 1:2:K;
GPStruct.xGeneFitness = cell(N,1);
bestGeneGroup = cell(SimParams.nBands,1);
GeneGroupHistory = cell((GPStruct.Ng + 1),SimParams.nBands);

for iBand = 1:N
    
    % Entering into GP
    GPStruct.cBand = iBand;
    GPStruct.geneGroup = getRandChromosomes(E,K,Nt,GPStruct.Np);
    GeneGroupHistory{1,iBand} = GPStruct.geneGroup;
        
    for iGeneration = 1:GPStruct.Ng
    
        % Metrics calculation
    
        GPStruct.xGeneFitness = getChromosomeWeights(SimParams,SimStructs,GPStruct);        
        GPStruct.xGeneFitnessMean = mean(GPStruct.xGeneFitness);
        GPStruct.xGeneFitnessDeviation = std(GPStruct.xGeneFitness);
        
        for iParents = 1:((GPStruct.Np / 2) - 1)
                  
            % Crossover
            [xChromozome yChromozome] = getChildChromosomes(GPStruct);
    
            % Mutation
            [childPair{1,1} childPair{2,1}] = performMutation(xChromozome,yChromozome,GPStruct);
            
            % Inversion
            [childPair{1,1} childPair{2,1}] = performInversion(childPair{1,1},childPair{2,1});

            % Constraint Qualification
            [xyChromozomes] = performConstraintQualification(childPair,GPStruct);    
            GeneGroupHistory{iGeneration + 1,iBand}(:,:,(iParents * 2 - 1)) = xyChromozomes{1,1};
            GeneGroupHistory{iGeneration + 1,iBand}(:,:,(iParents * 2)) = xyChromozomes{2,1};
            
        end
    
        % Elitism
        [~,srtI] = sort(GPStruct.xGeneFitness,'descend');
        GeneGroupHistory{iGeneration + 1,iBand}(:,:,(iParents * 2 + 1)) = GPStruct.geneGroup(:,:,srtI(1,1));
        GeneGroupHistory{iGeneration + 1,iBand}(:,:,(iParents * 2 + 2)) = GPStruct.geneGroup(:,:,srtI(1,2));
        GPStruct.geneGroup = GeneGroupHistory{iGeneration + 1,iBand};
    
    end
    
    GPStruct.xGeneFitness = getChromosomeWeights(SimParams,SimStructs,GPStruct);
    [~,maxI] = max(GPStruct.xGeneFitness);    
    bestGeneGroup{iBand,1} = GPStruct.geneGroup(:,:,maxI);    
    
end

for iBand = 1:SimParams.nBands
    presentIndex = bestGeneGroup{iBand,1};
    for iBase = 1:SimParams.nBases
        xValue = find(presentIndex(iBase,:) == 1);
        userIdx = SimStructs.baseStruct{iBase,1}.linkedUsers;
        SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1} = userIdx(xValue);
        SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1} = ones(length(xValue),1);
    end
end



end
