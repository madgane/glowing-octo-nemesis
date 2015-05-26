function [xConfig] = parseXLFile(xConfig)

xStruct = cell(1,1);
[M,N] = size(xConfig.txt);

for iField = 1:M

    if and(any(xConfig.raw{iField,1}),~sum(isnan(xConfig.raw{iField,1})))
        for iModel = 2:N
            if ischar(xConfig.raw{iField,iModel})
                stringT = sprintf('%s{%d,1}.%s = ''%s'';','xStruct',(iModel - 1),xConfig.raw{iField,1},xConfig.raw{iField,iModel});
            else
                stringT = sprintf('%s{%d,1}.%s = %d;','xStruct',(iModel - 1),xConfig.raw{iField,1},xConfig.raw{iField,iModel});
            end
            eval(stringT);
        end
    end
    
end

xConfig.xStruct = xStruct;

end