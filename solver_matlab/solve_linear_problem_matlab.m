function [x, fval, status] = solve_linear_problem_matlab(f, c, Aineq, bineq, Aeq, beq, lb, ub, x0)
    
    
    linprog_problem.solver = 'linprog';
    linprog_problem.f = f;
    linprog_problem.lb = lb;
    linprog_problem.ub = ub;
    linprog_problem.Aineq = Aineq;
    linprog_problem.bineq = bineq;
    linprog_problem.Aeq = Aeq;
    linprog_problem.beq = beq;
    linprog_problem.options = optimoptions('linprog', ...
                                           'Display', 'off', ...
                                           'Algorithm', 'dual-simplex');
    [x, fval, exitflag, output] = linprog(linprog_problem);
    if exitflag >= 0
        status = 0;
    else
        status = -1;
    end
end
