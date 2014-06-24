function fitnessWeights = fitnessFunction_1(GenStruct,newGenChromozomes)

iBand = GenStruct.iBand;
userLinkage = cell(GenStruct.nGenes,1);
fitnessWeights = zeros(1,GenStruct.xLen);

for iBase = 1:GenStruct.nGenes
    userLinkage{iBase,1} = GenStruct.SimStructs.baseStruct{iBase,1}.linkedUsers';
end

for iPop = 1:GenStruct.xPop
    
    xUsers = cell(GenStruct.nGenes,1);
    groupChannel = cell(GenStruct.nGenes,GenStruct.constraint_arguments(1,2));
    
    for iGene = 1:GenStruct.nGenes
        xUsers{iGene,1} = userLinkage{iGene,1}(1,find(newGenChromozomes{iPop,1}(iGene,:) == 1));
        for iUser = 1:length(xUsers{iGene,1});
            groupChannel{iGene,iUser} = GenStruct.SimStructs.linkChan{iGene,iBand}(:,:,xUsers{iGene,1}(1,iUser));
        end
    end
    
    basePrecoders = getZFprecoders(groupChannel,max(GenStruct.SimParams.sPower));
    
    interH = cell(GenStruct.nGenes - 1,1);
    interP = cell(GenStruct.nGenes - 1,1);
    
    for iBase = 1:GenStruct.nGenes
        intraP = basePrecoders{iBase,1};
        for iUser = 1:length(xUsers{iBase,1})
            kRunBase = 0;cUser = xUsers{iBase,1}(1,iUser);
            for kBase = 1:GenStruct.nGenes
                if kBase ~= iBase
                    kRunBase = kRunBase + 1;
                    interH{kRunBase,1} = GenStruct.SimStructs.linkChan{kBase,iBand}(:,:,cUser);
                    interP{kRunBase,1} = basePrecoders{kBase,1};                
                end
            end
            
            compoundSINR = getEstimatedSINR(...
                GenStruct.SimStructs.linkChan{iBase,iBand}(:,:,cUser),interH,intraP,interP,iUser);
            
            fitnessWeights(1,iPop) = fitnessWeights(1,iPop) + log2(1 + compoundSINR);
            
        end
    end
end

if (sum(fitnessWeights < 0))
    display('Error !! Negative fitness weights ! ');
end



