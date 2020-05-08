function [x, fval, info] = solve_quadratic_problem(H, g, c, Aineq, bineq, Aeq, beq, lb, ub, x0)

    [x, fval, info] = solve_quadratic_problem_matlab(H, g, c, Aineq, ...
                                                    bineq, Aeq, beq, ...
                                                    lb, ub, x0);


end