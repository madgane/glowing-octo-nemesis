function [SimParams,SimStructs] = getGeneticScheduling_4(SimParams,SimStructs)

GenStruct.xGen = 10;
GenStruct.xPop = 10;
GenStruct.xLen = SimParams.nUsers / SimParams.nBases;
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
GenStruct.constraint_arguments = [SimParams.nTxAntenna,SimParams.nTxAntenna];

GenStruct.CXprecoders = cell(GenStruct.nGenes,1);

for iBand = 1:SimParams.nBands
    
    GenStruct.iBand = iBand;
    for iBase = 1:SimParams.nBases
        Z = [];
        linkedUsers = SimStructs.baseStruct{iBase,1}.linkedUsers;
        Dmatrix = zeros(length(linkedUsers));
        
        for iUser = 1:length(linkedUsers)
            cUser = linkedUsers(iUser,1);
            [~, D, V] = svd(SimStructs.linkChan{iBase,iBand}(:,:,cUser));
            Z = [Z V(:,1)];
            Dmatrix(:,iUser) = D(1,1);
        end                
        GenStruct.CXprecoders{iBase,1} = abs(Z' * Z);
    end
    
    caseType = 1;
    
    switch (caseType)
        
        case 1
            Z = squeeze(SimStructs.linkChan{iBase,iBand}(:,:,linkedUsers));
            
            GenStruct.CXprecoders{iBase,1} = Z;
            bestPhenoType = getGeneticCombination(GenStruct);
             
        case 2
            Z = GenStruct.CXprecoders{1,1};X = 1./Z;
            Z = sum(X^4,2);Z = Z / sum(Z);[~,sortI] = sort(Z,'descend');
            bestPhenoType = zeros(1,GenStruct.xLen);bestPhenoType(1,sortI(1:4)) = 1;    
            
        case 3
            Z = GenStruct.CXprecoders{1,1} .* Dmatrix;
            X = diag(diag(Z)) + 1./Z - diag(diag(1./Z));
            Z = sum(X^4,2);Z = Z / sum(Z);[~,sortI] = sort(Z,'descend');
            bestPhenoType = zeros(1,GenStruct.xLen);bestPhenoType(1,sortI(1:4)) = 1;    
            
        case 4
            Xmatrix = Dmatrix' ./ Dmatrix;
            X = 1./GenStruct.CXprecoders{1,1} .* Xmatrix;
            Z = sum(X^4,2);Z = Z / sum(Z);[~,sortI] = sort(Z,'descend');
            bestPhenoType = zeros(1,GenStruct.xLen);bestPhenoType(1,sortI(1:4)) = 1;
            
    end            
    
    for iBase = 1:SimParams.nBases
        xValue = find(bestPhenoType(iBase,:) == 1);
        userIdx = SimStructs.baseStruct{iBase,1}.linkedUsers;
        SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1} = userIdx(xValue);
        SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1} = ones(length(xValue),1);
    end
   
end


