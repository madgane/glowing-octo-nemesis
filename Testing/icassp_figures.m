
load ICASSP_Results.mat;

nSchemes = length(rxSimParams);

if strcmp(rxSimParams{1,1}.capRun,'true')
    hold all;
    for iScheme = 1:nSchemes
        plot(rxSimParams{iScheme,1}.snrIndex,rxSimResults{iScheme,1}.sumCapacity);
    end
else
    hold all;
    for iScheme = 1:nSchemes
        plot(rxSimResults{iScheme,1}.queueBklgs);
    end
end

