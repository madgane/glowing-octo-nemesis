
function saveGlobals(varargin)

strFileName = sprintf('Results/%s.mat',varargin{1,nargin});

if exist(strFileName,'file')
	system(sprintf('rm %s',strFileName));
end

systemStructures = varargin;
save(strFileName,'systemStructures');

end