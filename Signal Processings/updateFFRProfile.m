function [SimParams, SimStructs] = updateFFRProfile(SimParams,SimStructs)

if isfield(SimParams,'ffrProfile_dB')
    ffrProfile = 10.^(0.1 * SimParams.ffrProfile_dB);
    ffrProfile = SimParams.nBands * ffrProfile ./ sum(ffrProfile);
    tempProfile = SimParams.sPower * ffrProfile;
else
    tempProfile = ones(1,SimParams.nBands) * SimParams.sPower;
end

for iBase = 1:SimParams.nBases
    SimStructs.baseStruct{iBase,1}.sPower = circshift(tempProfile',((iBase - 1) * SimParams.nBases + SimParams.iDrop - 1))';
end

end
