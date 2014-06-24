
processID = feature('GetPid');
display(processID);

if isunix
    cvx_quiet('true');
    cvx_solver('Sedumi');
else
    cvx_quiet('true');
    cvx_solver('Mosek');
end