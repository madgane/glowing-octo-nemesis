function [rssiMeasures, losMeasures] = getRSSIMeasures(SimParams,userPosition)

wrapModelLayout = 7;
nCites = getCellsOverLayout(SimParams.nTiers,1);
rssiMeasures = zeros(SimParams.nBases,wrapModelLayout);
losMeasures = cell(SimParams.nBases,wrapModelLayout);

for iWrapMode = 1:wrapModelLayout
    for iCite = 1:nCites
        basePosition = SimParams.wrapCellLocArray(iCite,iWrapMode);
        separationM = abs(userPosition - basePosition);
        
        [xRSSI, tempLOS] = evaluateLTE_PL(SimParams,separationM,'true','false');
        
        for iSector = 1:SimParams.nSectors
            currentSite = (iCite - 1) * SimParams.nSectors + iSector;
            losMeasures{currentSite,iWrapMode} = tempLOS;
            antennaGain = getAntennaPatterGain(basePosition,userPosition,SimParams,iSector);
            rssiMeasures(currentSite,iWrapMode) = xRSSI + antennaGain;
        end
        
    end
end

end