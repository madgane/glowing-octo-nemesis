
processID = feature('GetPid');
display(processID);

if isunix
    cvx_quiet('true');
    cvx_solver('Sedumi');
    cvx_expert('true');
else
    cvx_quiet('true');
    cvx_solver('Mosek_2');
    cvx_expert('true');
end