
iBase = 1;

for iBand = 1:SimParams.nBands
    for iUser = 1:SimParams.nUsers
        fprintf(1,'UserID - %d, SBlock - %d, Norm - %f \n',iUser,iBand,norm(SimStructs.linkChan{iBase,iBand}(:,:,iUser)));
        display(SimStructs.linkChan{iBase,iBand}(:,:,iUser));
    end
end
