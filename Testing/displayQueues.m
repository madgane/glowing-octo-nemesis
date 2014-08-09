
function displayQueues(SimParams,SimStructs,iDrop)

if nargin == 2
    iDrop = SimParams.nDrops;
end

underscore_location = strfind(SimParams.weightedSumRateMethod,'_');
if isempty(underscore_location)
    qExponent = 1;
else
    qExponent = str2double(SimParams.weightedSumRateMethod(underscore_location + 1:end));
end

fprintf('\n');
fprintf('Displaying User Queue status for drop - %2d \n',iDrop);
display('------------------------------------------');

printLatexScript = 'false';
Queues = zeros(SimParams.nUsers,1);
presentQueues = zeros(SimParams.nUsers,1);
txPkts = zeros(SimParams.nUsers,SimParams.nBands);

for iUser = 1:SimParams.nUsers
    Queues(iUser,1) = SimStructs.userStruct{iUser,1}.trafficStats.backlogsOverTime(1,iDrop);
    presentQueues(iUser,1) = SimStructs.userStruct{iUser,1}.trafficStats.backlogsOverTime(1,(iDrop + 1));
    txPkts(iUser,:) = squeeze(SimParams.Debug.resAllocation(iDrop,:,iUser,end));
end

servedPkts = sum(txPkts,2);
Qdeviation = sum(max((Queues - servedPkts),0));
QueueMatrix = [txPkts Queues servedPkts presentQueues];

display(QueueMatrix);
fprintf('Queue Deviation - %f bits \n',Qdeviation);
fprintf('Total Tx Rate - %f bits \n',sum(txPkts(:)));
fprintf('Objective Function - %f bits \n',norm((Queues - servedPkts),qExponent));

% Displaying in Latex import format

if strcmp(printLatexScript,'true')
    
    [nRows,nCols] = size(QueueMatrix);

    for iRow = 1:nRows
        for iCol = 1:nCols
            fprintf('& %3.2f \t',QueueMatrix(iRow,iCol));
        end
        fprintf('\n');
    end
    
end

