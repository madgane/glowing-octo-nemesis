function basePrecoders = getZFprecoders(groupChannel,totalP)

[nBases,nUsers] = size(groupChannel);
basePrecoders = cell(nBases,1);

for iBase = 1:nBases
    augH = [];
    for iUser = 1:nUsers
        if ~isempty(groupChannel{iBase,iUser})
            augH = [augH ; groupChannel{iBase,iUser}];
        end
    end
    
    X = pinv(augH);
    basePrecoders{iBase,1} = performWFAlgorithm(X,totalP);
end






