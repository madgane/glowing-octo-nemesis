
fprintf('\n');
display('Displaying System Configuration');
display('-------------------------------');

fprintf('\n');
fprintf('Number of BS - %d \n',SimParams.nBases);
fprintf('Number of Users - %d \n',SimParams.nUsers);
fprintf('Number of Bands - %d \n',SimParams.nBands);
fprintf('Max Pkt Arrivals - %d \n',SimParams.maxArrival);

fprintf('\n');
fprintf('Total Receive Antennas - %d \n',SimParams.nRxAntenna);
fprintf('Total Transmit Antennas - %d \n',SimParams.nTxAntenna);

fprintf('\n');
fprintf('Traffic Model - %s \n',SimParams.arrivalDist);
fprintf('Pathloss Model - %s \n',SimParams.pathLossModel);

fprintf('\n');
fprintf('Scheduling Type - %s \n',SimParams.SchedType);
fprintf('Precoding Method - %s \n',SimParams.PrecodingMethod);
fprintf('Precoding Sub Model - %s \n \n',SimParams.weightedSumRateMethod);

fprintf('\n');
fprintf('PL Profile - \n');
disp(SimParams.PL_Profile);
