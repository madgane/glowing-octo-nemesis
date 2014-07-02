
currentPath = path;

if isempty(strfind(currentPath,'glowing-octo-nemesis'))
    addpath(genpath(pwd));
end

if isempty(strfind(currentPath,'cvx'))
    if isunix
        cd ~/codes/solvers/cvx_linux;
        cvx_setup;
        cd ~/codes/finna-be-bugfixes;
    else
        cd ..\solvers\cvx;
        cvx_setup;
        cd ..\..\glowing-octo-nemesis;
    end
end

display('Added Path Variables !');

