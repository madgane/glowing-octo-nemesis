
function [SimParams,SimStructs] = generateTraffic(SimParams,SimStructs)

enStatToolBox = 'true';

queueModel = char(SimParams.arrivalDist);
uscore_index_q = find(queueModel == '_');

if ~isempty(uscore_index_q)
    queueType = queueModel(1:uscore_index_q(1,1) - 1);
    maxPktArrival = str2double(queueModel(uscore_index_q(1,1) + 1:end));
else
    queueType = queueModel;
    maxPktArrival = SimParams.maxArrival(1,SimParams.iPkt);
end

switch queueType
    case 'Uniform'
        dRange = floor(maxPktArrival * 0.5);
        randArrival = randi([(maxPktArrival - dRange),(maxPktArrival + dRange)],1,SimParams.nUsers);
    case 'Constant'
        randArrival = ones(1,SimParams.nUsers) * maxPktArrival;
    case 'Fixed'
        randArrival = SimParams.FixedPacketArrivals;
    case 'ConstFixed'
        randArrival = SimParams.FixedPacketArrivals;
    case 'SteadyFlow'
        randArrival = ones(1,SimParams.nUsers) * maxPktArrival;
end

SimParams.avgPktValues = randArrival;
surplusPkts = mod(SimParams.nDrops,SimParams.groupArrivalFreq);

if surplusPkts ~= 0
    nSamples = SimParams.groupArrivalFreq - surplusPkts;
else
    nSamples = 0;
end

for iUser = 1:SimParams.nUsers
   
    cLambda = SimParams.avgPktValues(1,iUser);
    SimStructs.userStruct{iUser,1}.trafficConfig.avgArrRate = cLambda;

    if strcmp(enStatToolBox,'true')
        poissonArrivals = random('Poisson',cLambda,1,SimParams.nDrops);
    else
        poissonArrivals = getPoisson(cLambda,1,SimParams.nDrops);
    end
    
    if strcmp(queueType,'ConstFixed')
    poissonArrivals = SimParams.avgPktValues(1,iUser) * ones(1,length(poissonArrivals));
    end
    
    if strcmp(queueType,'SteadyFlow')
    poissonArrivals = SimParams.avgPktValues(1,iUser) * ones(1,length(poissonArrivals));
    end
    
    poissonArrivals = [poissonArrivals, zeros(1,nSamples)];    
    xArrivals = reshape(poissonArrivals,SimParams.groupArrivalFreq,length(poissonArrivals) /SimParams.groupArrivalFreq);
    xArrivals = upsample(sum(xArrivals,1),SimParams.groupArrivalFreq);
    SimStructs.userStruct{iUser,1}.trafficHistory.pktArrival = reshape(xArrivals(1:SimParams.nDrops),1,SimParams.nDrops);
    
end


end