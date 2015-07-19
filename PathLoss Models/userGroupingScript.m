
ifCellsOfEachUser = cell(nUsers,1);
ifUsersOfEachCell = cell(nBases,3);

for iBase = 1:nBases
    cBase = SimStructs.baseStruct{iBase,1};
    xUsers = cBase.linkedUsers;
    
    for iUser = 1:length(xUsers)
        xUser = xUsers(iUser,1);
        cUser = SimStructs.userStruct{xUser,1};
        userNeighbors = cUser.neighNode;
        for jBase = 1:length(userNeighbors)
            ifUsersOfEachCell{userNeighbors(1,jBase),1} = [ifUsersOfEachCell{userNeighbors(1,jBase),1}, xUser];
            ifUsersOfEachCell{userNeighbors(1,jBase),2} = [ifUsersOfEachCell{userNeighbors(1,jBase),2}, iUser];
            ifUsersOfEachCell{userNeighbors(1,jBase),3} = [ifUsersOfEachCell{userNeighbors(1,jBase),3}, iBase];
        end
        
        ifCellsOfEachUser{xUser,1} = userNeighbors;        
    end       
end