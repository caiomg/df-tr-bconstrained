function [trial_point, trial_decrease] = solve_tr_subproblem(model, bl, bu, options)
% SOLVE_TR_SUBPROBLEM --- approximately solves trust-region
% subproblem.
%


    obj_pol = model.modeling_polynomials{1};
    x_tr_center = model.center_point();
    obj_pol = shift_polynomial(obj_pol, -x_tr_center); % Shift to origin
    radius = model.radius;
    [trial_point, trial_fval, exitflag] = minimize_tr(obj_pol, x_tr_center, radius, bl, bu);

    current_fval = model.center_fvalues(1);

    trial_decrease = current_fval - trial_fval;
    
    if trial_decrease <= 0
       warning('cmg:trial_not_decrease', ...
               'Trial step does not provide decrease');
    end

end
