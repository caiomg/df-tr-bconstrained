function [trial_point, trial_decrease] = solve_tr_subproblem(model, bl, bu, options)
% SOLVE_TR_SUBPROBLEM --- approximately solves trust-region
% subproblem.
%


    obj_pol = model.modeling_polynomials{1};
    x_tr_center = model.points_abs(:, model.tr_center);
    obj_pol = shift_polynomial(obj_pol, -x_tr_center); % Shift to origin
    radius = model.radius;
    [trial_point, trial_fval, exitflag] = minimize_tr(obj_pol, x_tr_center, radius, bl, bu);

    current_fval = model.fvalues(1, model.tr_center);

    trial_decrease = current_fval - trial_fval;
    
    if current_fval <= trial_fval
       1; 
    end
    
    tol_interp = max(1e-8, eps(max(1, max(model.fvalues(1, :))))*1e3);
    n_points = size(model.points_abs, 2);
    for k = 1:n_points
        val = evaluate_polynomial(obj_pol, model.points_abs(:, k));
        error_interp = abs(val - model.fvalues(1, k));
        if error_interp > tol_interp
            1;
        end
    end

end
