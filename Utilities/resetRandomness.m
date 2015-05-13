
function resetRandomness

versionNumber = version;
switch(versionNumber(end-6:end-2))
    
    case {'R2010'}
        stream = RandStream.getDefaultStream;reset(stream);
    case {'R2011'} 
        stream = RandStream.getGlobalStream;reset(stream);
    case {'R2012','R2013','R2014','R2015'}
        rng('default');
    otherwise
        rand('seed',0);randi('seed',0);randn('seed',0);
end
  
end