
processID = feature('GetPid');
display(processID);

if isunix
    cvx_quiet('true');
    cvx_solver('Mosek_2');
else
    cvx_quiet('true');
    cvx_solver('Mosek');
end