
epsilonT = 1e-5;
cH = SimStructs.linkChan;
nBases = SimParams.nBases;
nBands = SimParams.nBands;
maxRank = SimParams.maxRank;
globalMode = SimParams.totalPwrDistOverSC;

vector = @(x)(x(:));
usersPerCell = zeros(nBases,1);
cellUserIndices = cell(nBases,1);
cellNeighbourIndices = cell(nBases,1);

% Debug Buffers initialization

rankArray = linspace(1,maxRank,maxRank);
if isempty(SimParams.Debug.tempResource{2,SimParams.iDrop})
    SimParams.Debug.tempResource{2,SimParams.iDrop} = cell(SimParams.nUsers,1);
    SimParams.Debug.tempResource{3,SimParams.iDrop} = cell(SimParams.nUsers,1);
    SimParams.Debug.tempResource{4,SimParams.iDrop} = cell(SimParams.nUsers,SimParams.nBands);
end

for iBase = 1:nBases
    cellUserIndices{iBase,1} = unique(SimStructs.baseStruct{iBase,1}.linkedUsers);
    usersPerCell(iBase,1) = length(cellUserIndices{iBase,1});
end

nUsers = sum(usersPerCell);
QueuedPkts = zeros(nUsers,1);

for iBase = 1:nBases
    for iUser = 1:usersPerCell(iBase,1)
        cUser = cellUserIndices{iBase,1}(iUser,1);
        QueuedPkts(cUser,1) = SimStructs.userStruct{cUser,1}.trafficStats.backLogPkt;
    end
end

for iBase = 1:nBases
    for jBase = 1:nBases
        if jBase ~= iBase
            cellNeighbourIndices{iBase,1} = [cellNeighbourIndices{iBase,1} ; cellUserIndices{jBase,1}];
        end
    end
end

userWts = ones(SimParams.nUsers,1);
underscore_location = strfind(SimParams.weightedSumRateMethod,'_');
if isempty(underscore_location)
    qExponent = 1;
    selectionMethod = SimParams.weightedSumRateMethod;
else
    qExponent = str2double(SimParams.weightedSumRateMethod(underscore_location + 1:end));
    selectionMethod = SimParams.weightedSumRateMethod(1:underscore_location-1);
end

if ~strcmpi(SimParams.SchedType,'SkipScheduling')
    if strcmpi(SimParams.Debug.reSchedule,'true')
        SimParams.Debug.schTable = cell(nBases,nBands);
        for iBase = 1:nBases
            for iBand = 1:nBands
                SimParams.Debug.schTable{iBase,iBand}.assignedUsers = SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1};
                SimParams.Debug.schTable{iBase,iBand}.assignedStreams = SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1};
            end
        end
        [SimParams,SimStructs] = getSkipScheduling(SimParams,SimStructs);
        SimParams.Debug.reSchedule = 'false';
    end
end

