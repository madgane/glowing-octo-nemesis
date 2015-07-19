function [SimParams,SimStructs,varargout] = getScheduledUsers(SimParams,SimStructs,varargin)

chScheduler = char(SimParams.SchedType);
uScoreIndex = find(chScheduler == '_');
if isempty(uScoreIndex)
    scheduleMethod = SimParams.SchedType;
else
    scheduleMethod = chScheduler(1:uScoreIndex(1,1) - 1);
end

if nargin == 2
    
    switch scheduleMethod
        
        case 'RRScheduling'
            [SimParams,SimStructs] = getRoundRobinScheduling(SimParams,SimStructs);
            
        case 'RouletteScheduling'
            [SimParams,SimStructs] = getRouletteWheelScheduling(SimParams,SimStructs);
            
        case 'GreedyScheduling'
            [SimParams,SimStructs] = getGreedyScheduling(SimParams,SimStructs);
            
        case 'PFScheduling'
            [SimParams,SimStructs] = getFairnessScheduling(SimParams,SimStructs);
            
        case 'BDScheduling'
            [SimParams,SimStructs] = getBDScheduling(SimParams,SimStructs);
            
        case 'PFBDScheduling'
            [SimParams,SimStructs] = getPFBDScheduling(SimParams,SimStructs);
            
        case 'ExScheduling'
            [SimParams,SimStructs] = getExhaustiveScheduling(SimParams,SimStructs);
            
        case 'XScheduling'
            [SimParams,SimStructs] = getXScheduling(SimParams,SimStructs);
            
        case 'NetworkScheduling'
            [SimParams,SimStructs] = getNetworkScheduling(SimParams,SimStructs);
            
        case 'CoordScheduling'
            [SimParams,SimStructs] = getCoordinateScheduling(SimParams,SimStructs);
            
        case 'SkipScheduling'
            [SimParams,SimStructs] = getSkipScheduling(SimParams,SimStructs);
            
        case 'GScheduling'
            [SimParams,SimStructs] = getGScheduling(SimParams,SimStructs);
            
        otherwise
            display('Unknown Scheduling Type');
    end
    
    for iBand = 1:SimParams.nBands
        for iBase = 1:SimParams.nBases
            SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1} = SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1}(SimStructs.baseStruct{iBase,1}.assignedUsers{iBand,1} ~= 0);
            SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1} = SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1}(SimStructs.baseStruct{iBase,1}.assignedStreams{iBand,1} ~= 0);
        end
    end

else
    
    switch scheduleMethod
        
        case 'RRScheduling'
            varargout{1,1} = sprintf('%s - %s',varargin{1,1},'RR sched.');
            
        case 'RouletteScheduling'
            varargout{1,1} = sprintf('%s - %s',varargin{1,1},'Random sched.');
            
        case 'GreedyScheduling'
            varargout{1,1} = sprintf('%s - %s',varargin{1,1},'Greedy sched.');
            
        case 'PFScheduling'
            
            charScheduling = char(SimParams.SchedType);
            uscore_index = find(charScheduling == '_');
            
            if ~isempty(uscore_index)
                fairnessType = charScheduling(uscore_index(1,1) + 1:end);
                switch fairnessType
                    case {'BF','AF'}
                        varargout{1,1} = sprintf('%s - %s',varargin{1,1},'PF sched.');
                    case 'PF'
                        varargout{1,1} = sprintf('%s - %s',varargin{1,1},'perc. PF sched.');
                    case 'SP'
                        varargout{1,1} = sprintf('%s - %s',varargin{1,1},'QR-PF sched.');
                end
            end
            
        case 'BDScheduling'
            varargout{1,1} = sprintf('%s - %s',varargin{1,1},'QR/SP sched.');
            
        case 'PFBDScheduling'
            varargout{1,1} = sprintf('%s - %s',varargin{1,1},'AHP sched.');
            
        case 'GeneticScheduling'
            [SimParams,SimStructs] = getGeneticScheduling_4(SimParams,SimStructs);
            
        case 'ExScheduling'
            varargout{1,1} = sprintf('%s - %s',varargin{1,1},'Exh. sched.');
            
        case 'XScheduling'
            varargout{1,1} = sprintf('%s - %s',varargin{1,1},'Joint sched.');
            
        case 'NetworkScheduling'
            [SimParams,SimStructs] = getNetworkScheduling(SimParams,SimStructs);
            
        case 'CoordScheduling'
            varargout{1,1} = sprintf('%s - %s',varargin{1,1},'Coord. sched.');
            
        case 'SkipScheduling'
            
        otherwise
            display('Unknown Scheduling Type');
    end
    
end
