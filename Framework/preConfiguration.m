
processID = feature('GetPid');
display(processID);

if isunix
    cvx_quiet('true');
    cvx_solver('SEDUMI');
    cvx_expert('true');
else
    cvx_quiet('true');
    cvx_solver('Mosek');
    cvx_expert('true');
end