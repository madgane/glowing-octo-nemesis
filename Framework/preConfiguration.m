
processID = feature('GetPid');
display(processID);

if isfield(SimParams,'cvxDisabled')
    if strcmpi(SimParams.cvxDisabled,'false')
        cvx_quiet('true');
        cvx_expert('true');
        if isunix
            cvx_solver('SDPT3');
        else
            cvx_solver('Mosek');
        end
    end
end