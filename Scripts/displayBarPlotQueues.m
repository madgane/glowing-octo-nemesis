function displayBarPlotQueues(SimParams,SimStructs,iDrop)

if nargin == 2
    iDrop = SimParams.nDrops;
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

bar(QueueMatrix);
