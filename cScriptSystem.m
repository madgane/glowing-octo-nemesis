
clc;
clear;

tic;
updatePath;
groupConfig;
preConfiguration;

xParams = cell(1,1);
xStructs = cell(1,1);

xConfig.workBook = 'SModel-A';
xConfig.fileName = 'Utilities/systemConfig.xlsx';
[xConfig.num, xConfig.txt, xConfig.raw] = xlsread(xConfig.fileName,xConfig.workBook);
xConfig = parseXLFile(xConfig);

xConfig.defaultFolderName = sprintf('Results/10July2015/%s',xConfig.workBook);
xConfig.defaultFileName = 'GlobalData';

if ~exist(xConfig.defaultFolderName,'dir')
    mkdir(xConfig.defaultFolderName);
end

xVal = fix(clock);
xConfig.Clock.S = sprintf('%d:%d:%d',xVal(1,4),xVal(1,5),xVal(1,6));
outFile = sprintf('%s/%s.mat',xConfig.defaultFolderName,xConfig.defaultFileName);

for xCount = 1:length(xConfig.xStruct)
    
    xParams{xCount,1}.xCount = xCount;
    xParams{xCount,1}.saveContents = 'true';
    xParams{xCount,1}.outFile = sprintf('%s-%d',xConfig.defaultFileName,xCount);
    xParams{xCount,1}.saveChannelInfo = 'false';
    xParams{xCount,1}.channelSaveFolder = xConfig.defaultFolderName;
    xParams{xCount,1}.LegendName = xConfig.xStruct{xCount,1}.LegendName;
    
    xParams{xCount,1}.maxDebugCells = 5;
    xParams{xCount,1}.version = version;
    xParams{xCount,1}.plotMode = 'NoDisplay';
    
    xParams{xCount,1}.sysMode = 'true';
    xParams{xCount,1}.DebugMode = 'false';
    xParams{xCount,1}.precoderWithIdealChn = 'false';
    xParams{xCount,1}.totalPwrDistOverSC = 'true';
    
    xParams{xCount,1}.ChannelModel = 'Jakes';
    xParams{xCount,1}.pathLossModel = '3GPP_UMi';
        
    xParams{xCount,1}.queueWt = 0;
    xParams{xCount,1}.mdpFactor = 0;
    xParams{xCount,1}.robustNoise = 0;
    xParams{xCount,1}.weighingEqual = 'false';

    xParams{xCount,1}.PF_dur = 40;
    xParams{xCount,1}.SFSymbols = 14;
    xParams{xCount,1}.sampTime = 1e-3;
    xParams{xCount,1}.estError = 0.00;
    xParams{xCount,1}.fbFraction = 0.00;

    xParams{xCount,1}.SchedType = xConfig.xStruct{xCount,1}.SchedType;
    xParams{xCount,1}.PrecodingMethod = xConfig.xStruct{xCount,1}.PrecodingMethod;
    xParams{xCount,1}.weightedSumRateMethod = xConfig.xStruct{xCount,1}.weightedSumRateMethod;
    xParams{xCount,1}.additionalParams = xConfig.xStruct{xCount,1}.additionalParams;
    
    xParams{xCount,1}.nExchangesOTA = xConfig.xStruct{xCount,1}.nExchangesOTA;
    xParams{xCount,1}.exchangeResetInterval = xConfig.xStruct{xCount,1}.exchangeResetInterval;
    xParams{xCount,1}.nExchangesOBH = xConfig.xStruct{xCount,1}.nExchangesOBH;
    
    xParams{xCount,1}.nDrops = xConfig.xStruct{xCount,1}.nDrops;
    
    xParams{xCount,1}.BITFactor = 1;
    xParams{xCount,1}.nSymbolsBIT = xConfig.xStruct{xCount,1}.nSymbolsBIT;
    
    xParams{xCount,1}.nTxAntenna = xConfig.xStruct{xCount,1}.nTxAntenna;
    xParams{xCount,1}.nRxAntenna = xConfig.xStruct{xCount,1}.nRxAntenna;
    
    xParams{xCount,1}.maxArrival = 6;
    xParams{xCount,1}.groupArrivalFreq = 5;
    xParams{xCount,1}.arrivalDist = 'Constant';
    xParams{xCount,1}.FixedPacketArrivals = [6];
    
    xParams{xCount,1}.nSectors = 3;
    xParams{xCount,1}.nTiers = xConfig.xStruct{xCount,1}.nTiers;
    xParams{xCount,1}.perCiteUsers = xConfig.xStruct{xCount,1}.nUsersPerCite;
    
    xParams{xCount,1}.nBands = xConfig.xStruct{xCount,1}.nBands;
    xParams{xCount,1}.nBases = getCellsOverLayout(xParams{xCount,1}.nTiers,xParams{xCount,1}.nSectors);
    xParams{xCount,1}.nUsers = xParams{xCount,1}.perCiteUsers * xParams{xCount,1}.nBases;
    
    xParams{xCount,1}.nNeighbors = xParams{xCount,1}.nBases - 1;
    xParams{xCount,1}.snrIndex = 0 - 10 * log10(xParams{xCount,1}.nBands);
    xParams{xCount,1}.ffrProfile_dB = zeros(1,xParams{xCount,1}.nBands);    
end

pCluster = parcluster('local');
[~,xConfig.HostName] = system('hostname');
xConfig.HostName = strtrim(xConfig.HostName);
save(outFile,'xParams','xStructs','xConfig');

parfor xCount = 1:length(xParams)
    fwkSystemScript(xParams{xCount,1});
end
    
xVal = fix(clock);
xConfig.Date = date;
xConfig.Clock.E = sprintf('%d:%d:%d',xVal(1,4),xVal(1,5),xVal(1,6));
save(outFile,'xParams','xStructs','xConfig');

toc;
display('Completed Running the Simulation !');
