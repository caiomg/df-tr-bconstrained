function [model, measure] = criticality_step(model, funcs, bl, bu, options)
% CRITICALITY_STEP -- ensures model is sufficiently poised and with
% a radius comparable to the gradient
%


mu = options.criticality_mu; % factor between radius and
                             % criticality measure
omega = options.criticality_omega; % factor used to reduce radius
beta = options.criticality_beta; % to ensure the final radius
                                 % reduction is not drastic
tol_radius = options.tol_radius; % tolerance of TR algorithm
tol_f = options.tol_f;

initial_radius = model.radius;
while ~is_lambda_poised(model, options) || is_old(model, options)
    model = ensure_improvement(model, funcs, bl, bu, options);
    model.modeling_polynomials = compute_polynomial_models(model);
end
measure = norm(measure_criticality(model, bl, bu));
while (model.radius > mu*measure)
    model.radius = omega*model.radius;
    while ~is_lambda_poised(model, options) || is_old(model, options)
        model = ensure_improvement(model, funcs, bl, bu, options);
        model.modeling_polynomials = compute_polynomial_models(model);
    end

    measure = norm(measure_criticality(model, bl, bu));
    if (model.radius < tol_radius || (beta*measure < tol_f) && ...
            model.radius < 100*tol_radius)
        % Better break.
        % Not the end of this algorithm, but satisfies stopping
        % condition for outer algorithm anyway...
        break;
    end
end

% The final radius is increased not to make the reduction drastic
model.radius = min(max(model.radius, beta*measure), initial_radius);
                       
end