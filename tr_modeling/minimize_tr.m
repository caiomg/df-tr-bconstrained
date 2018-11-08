function [x, fval, exitflag] = minimize_tr(polynomial, x_tr_center, ...
                                           radius, bl, bu)

    matlab_solver = true;
    dim = size(x_tr_center, 1);
    if isempty(bl)
        bl = -inf(dim, 1);
    end
    if isempty(bu)
        bu = inf(dim, 1);
    end
    tol_tr = 10*eps(max(1, norm(x_tr_center, inf)));

    % TR bounds
    bl_tr = x_tr_center - radius;
    bu_tr = x_tr_center + radius;
    % Joining TR bounds and decision variable bounds
    bl_mod = max(bl, bl_tr) + tol_tr;
    bu_mod = min(bu, bu_tr) - tol_tr;

    % Restoring feasibility at TR center
    bl_active = x_tr_center <= bl;
    bu_active = x_tr_center >= bu;
    bl_mod(bl_active) = bl(bl_active);
    bu_mod(bu_active) = bu(bu_active);
    
    
    [c, g, H] = get_matrices(polynomial);
    f = @(x) quadratic(H, g, c, x);
    
    % Getting away from a stationary point
    if norm(g + H*x_tr_center) > 1e-4
        x0 = x_tr_center;
    else
        x0 = (bl_mod + bu_mod)/2;
        if norm(g + H*x0) < 1e-4
            % A bit more effort
            [V, D] = eig(H); % eig is more accurate than eigs & fast enough
            [~, min_eig] = min(diag(D));
            v = V(:, min_eig);
            nz = abs(v) < 1e-5;
            
            alpha_l = (bl_mod - x0)./V(:, min_eig);
            alpha_u = (bu_mod - x0)./V(:, min_eig);
            x0 = x0 + min(min(alpha_l(~nz), alpha_u(~nz)))*V(:, min_eig);
        end
    end

    solver_x_tol = min(1e-8, 1e-2*min(min(bu_mod - bl_mod)));
    solver_constr_tol = min(1e-10, 1e-2*min(min(bu_mod - bl_mod)));
    if ~matlab_solver && exist('ipopt', 'file') ~= 3
        % Override option
        matlab_solver = true;
    end
    if matlab_solver
        if norm(H, inf) > 10*eps(norm(g))
            fmincon_options = optimoptions(@fmincon, 'Display', 'off', ...
                                           'Algorithm', 'interior-point', ...
                                           'SpecifyObjectiveGradient', true);%, ...
                                           %'StepTolerance', solver_x_tol, ...
                                           %'OptimalityTolerance', 1e-6);
                                           %'ConstraintTolerance', solver_constr_tol);
            [x, fval, exitflag] = fmincon(f, x0, [], [], [], [], ...
                                          bl_mod, bu_mod, [], ...
                                          fmincon_options);
            if exitflag < 0 || norm(x) < 0.001*radius
                fmincon_options = optimoptions(@fmincon, 'Display', 'off', ...
                                               'Algorithm', 'active-set', ...
                                               'SpecifyObjectiveGradient', true, ...
                                               'StepTolerance', solver_x_tol, ...
                                               'OptimalityTolerance', 1e-6);
                                               %'ConstraintTolerance', solver_constr_tol);
                [x2, fval2, exitflag2] = fmincon(f, x0, [], [], [], [], ...
                                              bl_mod, bu_mod, [], ...
                                              fmincon_options);
                 if (fval2 < fval && exitflag2 >= 0) || exitflag < 0
                     x = x2;
                     fval = fval2;
                     exitflag = exitflag2;
                 end
             end
        else
            linprog_problem.f = g;
            linprog_problem.solver = 'linprog';
            linprog_problem.lb = bl_mod;
            linprog_problem.ub = bu_mod;
            linprog_problem.options.Display = 'off';
            [x, ~, exitflag, output] = linprog(linprog_problem);
            fval = f(x);
        end
        if exitflag >= 0
             exitflag = 0;
        else
             warning('cmg:tr_minimization_failed', 'TR minimization failed');
        end
    else
        f_ipopt.objective = f;
        f_ipopt.gradient = @(x) H*x + g;
        f_ipopt.hessian = @(x, sigma, lambda) sparse(tril(H));
        f_ipopt.hessianstructure = @() sparse(tril(ones(size(H))));

        ipopt_options.lb = bl_mod;
        ipopt_options.ub = bu_mod;
        ipopt_options.ipopt.print_level = 0;
        ipopt_options.ipopt.hessian_constant = 'yes';
        ipopt_options.ipopt.acceptable_iter = 100;
        ipopt_options.ipopt.tol = 1e-9;
        ipopt_options.ipopt.compl_inf_tol = 1e-5;
        % ipopt_options.ipopt.dual_inf_tol = 0.5;
        % ipopt_options.ipopt.constr_viol_tol = solver_constr_tol;

        [x, info] = ipopt(x0, f_ipopt, ipopt_options);
        fval = f(x);
        exitflag = 0;
    end
end

    