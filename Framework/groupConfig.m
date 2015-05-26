
currentPath = path;
currentDirectory = pwd;

if isempty(strfind(currentPath,'cvx'))
    if isunix
        cd ~/codes/solvers/cvx_linux;
        cvx_setup;
    else
        cd ..\solvers\cvx;
        cvx_setup;
    end
    cd(currentDirectory);
end

display('Added Path Variables !');
