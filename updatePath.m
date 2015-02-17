
currentDirectories = genpath(pwd);
colonLocs = strfind(currentDirectories,';');
colonLocs = [0 colonLocs];
xPath = [];

for iLen = 2:length(colonLocs) - 2
    if isempty(strfind(currentDirectories(colonLocs(1,iLen - 1) + 1:colonLocs(1,iLen) - 1),'git'))
        xPath = [xPath,currentDirectories(colonLocs(1,iLen - 1) + 1:colonLocs(1,iLen) - 1),';'];
    end
end

addpath(char(xPath));