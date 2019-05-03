unconstrained_problems = {'HS1', 'HS2', 'HS3', 'HS4', 'HS5', 'HS25', ...
                          'HS38', 'HS45'};

n_problems = length(unconstrained_problems);
fmincon_options = optimoptions('fmincon', 'Display', 'off', ...
                               'Algorithm', 'sqp', ...
                               'SpecifyObjectiveGradient', false);

terminate_cutest_problem()
clear global problem_path_cutest problem_name_cutest problem_data_cutest
global problem_data_cutest
clear results_unconstrained

warning('off', 'cmg:ill_conditioned_system')
warning('off', 'cmg:trial_not_decrease');

results_unconstrained(n_problems).fval_matlab = [];
results_unconstrained(n_problems).fval_trust = [];
results_unconstrained(n_problems).fcount_matlab = [];
results_unconstrained(n_problems).fcount_trust = [];
results_unconstrained(n_problems).name = '';

for n = 1:n_problems
    
    terminate_cutest_problem()

    problem_name = unconstrained_problems{n};
    results_unconstrained(n).name = problem_name;
    
    prob = setup_cutest_problem(problem_name, '../my_problems/');

    % Objective
    f_obj = @(x) get_cutest_objective(x);
    counter = evaluation_counter(f_obj);
    counter.set_max_count(10000);
    f = @(x) counter.evaluate(x);

    x0 = prob.x;
    bl = prob.bl;
    bu = prob.bu;

    [x_fmin, fvalue_fmin] = fmincon(f, x0, [], [], [], [], ...
                                          bl, bu, [], ...
                                          fmincon_options);
    f_count_fmincon = counter.get_count();
    counter.reset_count();
    
    results_unconstrained(n).ref.fx = fvalue_fmin;
    results_unconstrained(n).ref.count = f_count_fmincon;
    results_unconstrained(n).ref.viol = 0;

    options = [];

    try
        [x_trust, fvalue_trust] = trust_region({f}, x0, f(x0), bl, bu, options);
        f_count_trust = counter.get_count();
        results_unconstrained(n).test.fx = fvalue_trust;
        results_unconstrained(n).test.count = f_count_trust;
        results_unconstrained(n).test.viol = 0;
    catch this_exception
        results_unconstrained(n).test.count = counter.get_count();
        results_unconstrained(n).test.exception = this_exception;
    end

    counter.reset_count();
    fprintf(1, format_test_result(results_unconstrained(n)));

    terminate_cutest_problem()

end