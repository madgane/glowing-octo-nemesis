
function saveGlobals(varargin)

if isunix
	strFileName = sprintf('Results/%s.mat',varargin{1,nargin});
	if exist(strFileName,'file')
		system(sprintf('rm %s',strFileName));
	end
else
	strFileName = sprintf('Results\\%s.mat',varargin{1,nargin});
	if exist(strFileName,'file')
		system(sprintf('del %s',strFileName));
	end
end

systemStructures = varargin;
save(strFileName,'systemStructures');

end