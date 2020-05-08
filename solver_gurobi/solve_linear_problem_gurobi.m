function [x, fval, status] = solve_linear_problem_gurobi(g, c, Aineq, ...
                                                      bineq, Aeq, ...
                                                      beq, lb, ub, x0)


    n_ineqs = size(bineq, 1);
    n_eqs = size(beq, 1);
    dim = size(g, 1);

    model.obj = g;
    model.objcon = c;
    if n_ineqs == 0 && n_eqs == 0
        A = sparse(zeros(0, dim));
    else
        A = sparse([Aineq;
                     Aeq]);
    end
    model.A = A;
    model.rhs = [bineq;
               beq];
    model.sense = [repmat('<', n_ineqs, 1);
                   repmat('=', n_eqs, 1)];
    if isempty(lb)
        % We have to always set lb: Gurobi defaults to lb = 0
        lb = -inf(dim, 1);
    end
    model.lb = lb;
    if ~isempty(ub)
        model.ub = ub;
    end
    
    params.OutputFlag = 0;

    try
        result = gurobi(model, params);
    catch thiserror
        rethrow(thiserror);
    end
    if isfield(result, 'x')
        x = result.x;
        fval = result.objval;
        H = zeros(dim);
        fval_safeguard = quadratic(H, g, c, x);
        if abs(fval - fval_safeguard) - sqrt(eps(max(1, abs(fval)))) > 0
            warning('cmg:different_objective', ...
                    'Diference in objectives: % 8.3g', ...
                    abs(fval - fval_safeguard));
        end
    else
        x = [];
        fval = [];
    end
    if strcmp(result.status, 'OPTIMAL')
        status = 0;
    else
        status = -1;
    end

end