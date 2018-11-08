tr_options = struct('tol_radius', 1e-5, 'tol_f', 1e-5, ...
                       'eps_c', 1e-4, 'eta_0', 0, 'eta_1', 0.05, ...
                       'gamma_inc', 2, 'gamma_dec', 0.5, ...
                        'initial_radius', 1, 'radius_max', 1e3, ...
                        'criticality_mu', 100, 'criticality_beta', 10, ...
                        'criticality_omega', 0.5, 'basis', 'full quadratic', ...
                        'pivot_threshold', 0.001, 'radius_factor', 6, ...
                        'xi_imp', 1.1);
                    

f = @(x) norm(x)^2;
bl = [-10; -10];
bu = [10; 10];

x0 = ones(2, 1);
fvalues = f(x0);
radius = 1;
model = tr_model(x0, fvalues, radius)

model = rebuild_model(model, tr_options)
%%
[model, success] = improve_model_nfp(model, {f}, bl, bu, tr_options)

%%
model.radius = 2;
model = add_point(model, x0 + 2, f(x0+2), tr_options)

%%
                    