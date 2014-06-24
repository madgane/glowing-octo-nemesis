function totCells = getCellsOverLayout(nTiers,nSectors)

totCells = 0;
for iTier = 1:nTiers
    if iTier == 1
        nCells = 1;
    else
        nCells = 2^(iTier - 2) * 6;
    end
    
    totCells = totCells + nCells * nSectors;    
end

end
