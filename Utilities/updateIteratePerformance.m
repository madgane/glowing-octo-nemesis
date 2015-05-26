
function [SimParams,SimStructs] = updateIteratePerformance(varargin)

switch nargin
    case 2
        SimParams = varargin{1,1};
        SimStructs = varargin{1,2};
    otherwise
        SimParams = varargin{1,1};
        SimStructs = varargin{1,2};
        cellM = varargin{1,3};
        W = varargin{1,4};
end

nBases = SimParams.nBases;
nBands = SimParams.nBands;
nUsers = SimParams.nUsers;

QueuedPkts = zeros(nUsers,1);
usersPerCell = zeros(nBases,1);
cellUserIndices = cell(nBases,1);

for iBase = 1:nBases
    for iBand = 1:nBands
        cellUserIndices{iBase,1} = [cellUserIndices{iBase,1} ; SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1}];
    end
    cellUserIndices{iBase,1} = unique(cellUserIndices{iBase,1});
    usersPerCell(iBase,1) = length(cellUserIndices{iBase,1});
end

for iBase = 1:nBases
    for iUser = 1:usersPerCell(iBase,1)
        cUser = cellUserIndices{iBase,1}(iUser,1);
        QueuedPkts(cUser,1) = SimStructs.userStruct{cUser,1}.trafficStats.backLogPkt;
    end
end

underscore_location = strfind(SimParams.weightedSumRateMethod,'_');
if isempty(underscore_location)
    qExponent = 1;
    selectionMethod = SimParams.weightedSumRateMethod;
else
    qExponent = str2double(SimParams.weightedSumRateMethod(underscore_location + 1:end));
    selectionMethod = SimParams.weightedSumRateMethod(1:underscore_location-1);
end

if nargin ~= 2    
    
    if iscell(cellM)
        
        [~,xBands] = size(cellM);
        
        if xBands ~= nBands
            for iBase = 1:nBases
                for iBand = 1:nBands
                    SimStructs.baseStruct{iBase,1}.P{iBand,1} = cellM{iBase,1}(:,:,:,iBand);
                end
            end
        else
            for iBase = 1:nBases
                for iBand = 1:nBands
                    SimStructs.baseStruct{iBase,1}.P{iBand,1} = cellM{iBase,1}(:,:,:,iBand);
                end
            end
        end
        
    else
        
        for iBase = 1:nBases
            for iBand = 1:nBands
                P = [];
                for iUser = 1:usersPerCell(iBase,1)
                    cUser = cellUserIndices{iBase,1}(iUser,1);
                    P = [P cellM(:,:,cUser,iBand)];
                end
                SimStructs.baseStruct{iBase,1}.P{iBand,1} = P;
            end
        end
        
    end
    
    for iUser = 1:nUsers
        for iBand = 1:nBands
            SimStructs.userStruct{iUser,1}.W{iBand,1} = W{iUser,iBand};
        end
    end

end

SimParams.Debug.privateExchanges.resAllocation = zeros(nBands,nUsers);
[SimParams,SimStructs] = performDummyReception(SimParams,SimStructs);
tBandUser = SimParams.Debug.privateExchanges.resAllocation;

qVector = zeros(nUsers,1);
for iUser = 1:nUsers
    qDeviation = max(QueuedPkts(iUser,1) - sum(tBandUser(:,iUser)),0);
    qVector(iUser,1) = QueuedPkts(iUser,1) - sum(tBandUser(:,iUser));
    SimParams.Debug.tempResource{2,SimParams.iDrop}{iUser,1} = [SimParams.Debug.tempResource{2,SimParams.iDrop}{iUser,1} sum(tBandUser(:,iUser))];
    SimParams.Debug.tempResource{3,SimParams.iDrop}{iUser,1} = [SimParams.Debug.tempResource{3,SimParams.iDrop}{iUser,1} qDeviation];
    for iBand = 1:nBands
        SimParams.Debug.tempResource{4,SimParams.iDrop}{iUser,iBand} = [SimParams.Debug.tempResource{4,SimParams.iDrop}{iUser,iBand} tBandUser(iBand,iUser)];
    end
end

tempSum = sum(cell2mat(SimParams.Debug.tempResource{2,SimParams.iDrop}));
tempQueue = sum(cell2mat(SimParams.Debug.tempResource{3,SimParams.iDrop}));
fprintf('(R,Q,|Q - R|) - (%.2f, %.2f, %.2f) \n',tempSum(end),tempQueue(end),norm(qVector,qExponent));

% saveGlobals(SimParams,SimStructs,'runPrintFile');

SimParams.currentQueue = tempQueue(end);

end