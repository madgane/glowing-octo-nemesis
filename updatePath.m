
currentDirectories = genpath(pwd);
if isunix
    colonLocs = strfind(currentDirectories,':');
else
    colonLocs = strfind(currentDirectories,';');
end
colonLocs = [0 colonLocs];
xPath = [];

for iLen = 2:length(colonLocs)
    if isempty(strfind(currentDirectories(colonLocs(1,iLen - 1) + 1:colonLocs(1,iLen) - 1),'git'))
        display(currentDirectories(colonLocs(1,iLen - 1) + 1:colonLocs(1,iLen) - 1));
        xPath = [xPath,currentDirectories(colonLocs(1,iLen - 1) + 1:colonLocs(1,iLen) - 1),';'];
    end
end

addpath(char(xPath));