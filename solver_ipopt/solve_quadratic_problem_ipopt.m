function [x, fval, status] = solve_quadratic_problem_ipopt(H, g, c, Aineq, bineq, Aeq, ...
                                           beq, lb, ub, x0)


    dim = size(g, 1);
    assert(size(H, 1) == dim && size(H, 2) == dim);

    
    f_ipopt.objective = @(x) quadratic(H, g, c, x);
    f_ipopt.gradient = @(x) H*x + g;
    f_ipopt.hessian = @(x, sigma, lambda) sparse(tril(H));
    f_ipopt.hessianstructure = @() sparse(tril(ones(size(H))));

    if isempty(lb) && ~isempty(ub)
        lb = -inf(dim, 1);
    elseif isempty(ub) && ~isempty(lb)
        ub = inf(dim, 1);
    end
    ipopt_options.lb = lb;
    ipopt_options.ub = ub;
    if ~isempty(Aineq) && ~isempty(Aeq)
        A = [Aineq;
             Aeq];
        con_ub = [bineq;
                  beq];
        con_lb = [-inf(size(bineq));
                  beq];
        ipopt_options.ipopt.jac_d_constant = 'yes';
        ipopt_options.ipopt.jac_c_constant = 'yes';
    elseif ~isempty(Aineq)
        A = Aineq;
        con_ub = bineq;
        con_lb = -inf(size(bineq));
        ipopt_options.ipopt.jac_d_constant = 'yes';
    elseif ~isempty(Aeq)
        A = Aeq;
        con_ub = beq;
        con_lb = beq;
        ipopt_options.ipopt.jac_c_constant = 'yes';
    else
        A = [];
    end
    if ~isempty(A)
        ipopt_options.cu = con_ub;
        ipopt_options.cl = con_lb;
        f_ipopt.constraints = @(x) A*x;
        f_ipopt.jacobian = @(x) sparse(A);
        f_ipopt.jacobianstructure = @(x) sparse(ones(size(A)));
    end

    ipopt_options.ipopt.print_level = 0;
    ipopt_options.ipopt.hessian_constant = 'yes';
        
    %ipopt_options.ipopt.nlp_scaling_method = 'equilibration-based';
    ipopt_options.ipopt.nlp_scaling_method = 'gradient-based';
    %ipopt_options.ipopt.acceptable_iter = 100;
    %ipopt_options.ipopt.tol = 1e-9;
    %ipopt_options.ipopt.compl_inf_tol = 1e-5;
    
    % ipopt_options.ipopt.dual_inf_tol = 0.5;
    % ipopt_options.ipopt.constr_viol_tol = solver_constr_tol;
    
    [x, info] = ipopt(x0, f_ipopt, ipopt_options);
    fval = f_ipopt.objective(x);
    %exitflag = info.status;
    if info.status <= -10
        'Debug';
    end
    if info.status == 2 || info.status > 3
        status = -200;
    else
        status = 0;
    end

end