function [x, fval, exitflag] = minimize_tr(polynomial, x_tr_center, ...
                                           radius, bl, bu)

    dim = size(x_tr_center, 1);
    if isempty(bl)
        bl = -inf(dim, 1);
    end
    if isempty(bu)
        bu = inf(dim, 1);
    end

    % TR bounds
    bl_tr = x_tr_center - radius;
    bu_tr = x_tr_center + radius;
    % Joining TR bounds and decision variable bounds
    bl_mod = max(bl, bl_tr);
    bu_mod = min(bu, bu_tr);

    % Restoring feasibility at TR center
    bl_active = x_tr_center <= bl;
    bu_active = x_tr_center >= bu;
    x_tr_center(bl_active) = bl(bl_active);
    x_tr_center(bu_active) = bu(bu_active);
    
    [c, g, H] = get_matrices(polynomial);
    
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

    [x, fval, exitflag] = solve_quadratic_problem(H, g, c, [], [], ...
                                                  [], [], bl_mod, bu_mod, x0);
    x = project_to_bounds(x, bl, bu);
end
