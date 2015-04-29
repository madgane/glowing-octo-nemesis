
function crossCheck(inCell,caseA,caseB)

equivChannel = [];
nMatrix = length(inCell);

switch caseA
    case 0
        [M,N] = size(inCell{1});
    case 1
        [N,M] = size(inCell{1});
end

for iMatrix = 1:nMatrix
    switch caseA
        case 0
            [U,~,~] = svd(inCell{iMatrix,1});
            tempChannel = U' * inCell{iMatrix,1};
            equivChannel = [equivChannel, transpose(tempChannel)];
        case 1
            equivChannel = [equivChannel, inCell{iMatrix,1}];
    end
end

Xsum = sqrt(sum(abs(equivChannel).^2,1));
display([uint16(Xsum * 2^15);uint32(Xsum * 2^15)]);

switch caseB
    case 0
        [~,sortA] = sort(Xsum,2,'descend');
        X = floor((sortA - 1) / M);
        Y = floor(mod(sortA - 1,min(M,N)));
    case 1
        [~,~,sortA] = qr(equivChannel,'vector');
        X = floor((sortA - 1) / M);
        Y = floor(mod(sortA - 1,min(M,N)));
end

display([X(1:N)' Y(1:N)']);

