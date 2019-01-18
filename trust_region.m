function [x, fval] = trust_region(funcs, initial_points, initial_fvalues, ...
                                  bl, bu, options)
% TRUST_REGION - Derivative-free trust-region algorithm
%   


defaultoptions = struct('tol_radius', 1e-5, 'tol_f', 1e-6, ...
                        'eps_c', 1e-5, 'eta_0', 0, 'eta_1', 0.05, ...
                        'pivot_threshold', 1/16, 'add_threshold', 100, ...
                        'exchange_threshold', 1000,  ...
                        'initial_radius', 1, 'radius_max', 1e3, ...
                        'radius_factor', 6, ...
                        'gamma_inc', 2, 'gamma_dec', 0.5, ...
                        'criticality_mu', 100, 'criticality_beta', 10, ...
                        'criticality_omega', 0.5, 'basis', 'diagonal hessian', ...
                        'iter_max', 10000, 'print_level', 0);

if nargin < 6
    options = [];
end
options = parse_trust_region_inputs(options, defaultoptions);
if nargin < 5
    bu = [];
end
if nargin < 4
    bl = [];
end
if nargin < 3
    initial_fvalues = [];
end


tol_radius = options.tol_radius;
tol_f = options.tol_f;
eps_c = options.eps_c;
eta_0 = options.eta_0;
eta_1 = options.eta_1;
gamma_0 = 0.0625;
gamma_1 = options.gamma_dec;
gamma_2 = options.gamma_inc;
radius_factor = options.radius_factor;
rel_pivot_threshold = options.pivot_threshold;

iter_max = options.iter_max;

initial_radius = options.initial_radius;
radius_max = options.radius_max;
print_level = options.print_level;

dimension = size(initial_points, 1);

if (~isempty(bl) && ~isempty(find(initial_points(:, 1) < bl, 1))) || ...
        (~isempty(bu) && ~isempty(find(initial_points(:, 1) > bu, 1)))
     if isempty(initial_fvalues)
         % Replace
         initial_points(:, 1) = project_to_bounds(initial_points(:, 1), bl, bu);
     else
         % Add
         initial_points = [project_to_bounds(initial_points(:, 1), bl, bu), ...
                           initial_points];
         initial_fvalues(:, 1) = evaluate_new_fvalues(funcs, initial_points(:, 1));
     end
end

n_initial_points = size(initial_points, 2);
if n_initial_points == 1
    % Finding a random second point
    old_seed = rng('default');
    second_point = rand(size(initial_points));
    rng(old_seed);
    while norm(second_point, inf) < rel_pivot_threshold
        % Second point must not be too close
        second_point = 2*second_point;
    end
    second_point = (second_point - 0.5)*initial_radius;
    second_point = initial_points(:, 1) + second_point;
    if ~isempty(bu)
        second_point = min(bu, second_point);
    end
    if ~isempty(bl)
        second_point = max(bl, second_point);
    end
    initial_points(:, 2) = second_point;
    n_initial_points = 2;
end
% Calculating function values for other points of the set
if length(initial_fvalues) < n_initial_points
    for k = 1:n_initial_points
        [initial_fvalues(:, k), succeeded] = evaluate_new_fvalues(funcs, initial_points(:, k));
        if ~succeeded
           error('cmg:bad_starting_point', 'Bad starting point');
        end
    end
end


% Initializing model structure
model = tr_model(initial_points, initial_fvalues, initial_radius);
model = rebuild_model(model, options);
model = move_to_best_point(model, bl, bu);
%basis = band_prioritizing_basis(size(model.points_shifted, 1));
model.modeling_polynomials = compute_polynomial_models(model);
if size(model.points_abs, 2) < 2
    [model, exitflag] = ensure_improvement(model, funcs, bl, bu, options);
end


rho = 0;
iter = 1;
fval_current = model.fvalues(1);
x_current = model.points_abs(:, model.tr_center);
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
        fval_current = model.fvalues(1, model.tr_center);
        x_current = model.points_abs(:, model.tr_center);
    end
    model.modeling_polynomials = compute_polynomial_models(model);
    err_model = check_interpolation(model);
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
        print_iteration(iter, fval_current, rho, model.radius, size(model.points_abs, 2));
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
        fval_trial = evaluate_new_fvalues(funcs, trial_point);
        
        % Actual reduction
        ared = fval_current - fval_trial;
        % Agreement factor
        rho = ared/(predicted_red);

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
    iter = iter + 1;
end
x = model.points_abs(:, model.tr_center);
fval = model.fvalues(1, model.tr_center);

