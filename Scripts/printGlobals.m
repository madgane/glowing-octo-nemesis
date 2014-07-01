
function printGlobals(varargin)

if nargin == 0
	load('Results/runPrintFile.mat');
else
	load(sprintf('Results/%s.mat',varargin{1}));
end

displayOutputs(systemStructures{1},systemStructures{2});

end