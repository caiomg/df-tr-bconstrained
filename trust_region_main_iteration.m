function [x, fval] = trust_region_main_iteration(funcs, model, bl, bu, options)
% TRUST_REGION_MAIN_ITERATION - 
%   


tol_radius = options.tol_radius;
tol_f = options.tol_f;
eps_c = options.eps_c;
eta_0 = options.eta_0;
eta_1 = options.eta_1;
gamma_0 = 0.0625;
gamma_1 = options.gamma_dec;
gamma_2 = options.gamma_inc;
radius_factor = options.radius_factor;

iter_max = options.iter_max;

radius_max = options.radius_max;
print_level = options.print_level;

debug_on = options.debug;
    

rho = 0;
iter = 1;
fval_current = model.center_fvalues();
x_current = model.center_point();
sum_rho = 0;
sum_rho_sqr = 0;
delay_reduction = 0;


for iter = 1:iter_max
    if (model.radius < tol_radius)
        break
    end
    if true || is_lambda_poised(model, options)
        % Move among points that are part of the model
        model = move_to_best_point(model, bl, bu);
        fval_current = model.center_fvalues();
        x_current = model.center_point();
    end
    model.modeling_polynomials = compute_polynomial_models(model);

    if debug_on
        err_model = check_interpolation(model);
    end
    % Criticality step -- if we are possibly close to the optimum
    if norm(measure_criticality(model, bl, bu)) <= eps_c
        [model, crit_measure] = criticality_step(model, funcs, bl, bu, options);
        criticality_step_performed = true;
        if crit_measure < tol_f
            break;
        end
    else
        criticality_step_performed = false;
    end
    iteration_model_fl = is_lambda_poised(model, options);

    % Print summary
    if print_level >= 1
        print_iteration(iter, fval_current, rho, model.radius, model.number_of_points());
    end
    
    % Compute step
    [trial_point, predicted_red] = solve_tr_subproblem(model, bl, bu, options);
    trial_step = trial_point - x_current;

    if ((predicted_red < tol_radius*1e-2) || ...
        (predicted_red < tol_radius*abs(fval_current) && ...
         norm(trial_step) < tol_radius) || ...
        (predicted_red < tol_f*abs(fval_current)*1e-3))
        rho = -inf;
        [model, mchange_flag] = ensure_improvement(model, funcs, bl, bu, options);
    else
        % Evaluate objective at trial point
        [fval_trial, f_succeeded] = evaluate_new_fvalues(funcs, trial_point);
        
        if f_succeeded
            % Actual reduction
            ared = fval_current - fval_trial;
            % Agreement factor
            rho = ared/(predicted_red);
        else
            % Step failed
            rho = -inf;
        end

        % Acceptance of the trial point
        if rho > eta_1
            % Successful iteration
            fval_current = fval_trial;
            x_current = trial_point;
            % Including this new point as the TR center
            [model, mchange_flag] = change_tr_center(model, trial_point, ...
                                                     fval_trial, ...
                                                     options);
            % this mchange_flag is not being used (rho > eta_1)
            if ~iteration_model_fl && mchange_flag == 4
                % Had to rebuild a model that wasn't even Fully
                % Linear
                % This shouldn't happen
                [model, mchange_flag] = ensure_improvement(model, ...
                                                           funcs, ...
                                                           bl, bu, options);
                % this mchange_flag is not being used (rho > eta_1)
            end
        else
             [model, mchange_flag] = try_to_add_point(model, ...
                                                      trial_point, ...
                                                      fval_trial, ...
                                                      funcs, bl, bu, ...
                                                      options);
             % if mchange_flag == 4, we had to rebuild the model
             % and the radius will be reduced
        end
        sum_rho = sum_rho + rho;
        sum_rho_sqr = sum_rho_sqr + rho^2;
    end


    % From time to time a step may end a bit outside the TR
    step_size = min(model.radius, norm(trial_step, inf));
    % Radius update
    if rho > eta_1
        radius_inc = max(1, gamma_2*(step_size/model.radius));
        model.radius = min(radius_inc*model.radius, radius_max);
    elseif (iteration_model_fl && (rho == -inf || mchange_flag == 4 || criticality_step_performed))
        % A good model should have provided a better point
        % We reduce the radius, since the error is related to the
        % radius
        % | f(x) - m(x) | < K*(radius)^2
        % rho == -inf -> too short step size
        % mchange_flag == 4 -> Couldn't add point, had to rebuild model
        if model.radius <= 2*tol_radius/gamma_1
            delay_reduction = delay_reduction + 1; 
        else
            delay_reduction = 0;
        end
        if delay_reduction >= 3 || delay_reduction == 0 || criticality_step_performed
            gamma_dec = gamma_1;
            model.radius = gamma_dec*model.radius;
            delay_reduction = 0;
        end
    end
end
x = model.center_point();
fval = model.center_fvalues();
