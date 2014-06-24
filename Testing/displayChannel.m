function displayChannel(SimParams,SimStructs)

fprintf('\n');
display('Displaying System Channel (last)');
display('--------------------------------');

ChannelMatrix = zeros(SimParams.nUsers,SimParams.nBands);

for iUser = 1:SimParams.nUsers    
    cCell = SimStructs.userStruct{iUser,1}.baseNode;    
    for iBand = 1:SimParams.nBands        
        ChannelMatrix(iUser,iBand) = norm(SimStructs.actualChannel{cCell,iBand}(:,:,iUser));        
    end
end

display(ChannelMatrix);