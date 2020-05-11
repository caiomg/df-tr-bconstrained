unconstrained_problems = {...
    'AKIVA', 'BEALE', 'BOXBODLS', 'BRKMCC', 'BROWNBS', 'CLIFF', 'CUBE', ...
'DANWOODLS', 'DENSCHNA', 'DENSCHNB', 'DENSCHNC', 'DENSCHNF', 'EXPFIT', ...
'GBRAINLS', 'HAIRY', 'HIMMELBB', 'HIMMELBG', 'HIMMELBH', 'JENSMP', ...
'MARATOSB', 'MEXHAT', 'MISRA1ALS', 'MISRA1BLS', 'MISRA1CLS', ...
'MISRA1DLS', 'POWELLBSLS', 'ROSENBR', 'S308', 'SISSER', 'ZANGWIL2', ...
'BARD', 'BENNETT5LS', 'BOX3', 'CHWIRUT1LS', 'CHWIRUT2LS', 'DENSCHND', ...
'DENSCHNE', 'ECKERLE4LS', 'ENGVAL2', 'GAUSSIAN', 'GROWTHLS', 'GULF', ...
'HATFLDD', 'HATFLDE', 'HELIX', 'HIELOW', 'LSC2LS', 'MEYER3', 'MGH10LS', ...
'RAT42LS', 'ALLINITU', 'BROWNDEN', 'HIMMELBF', 'KOWOSB', 'MGH09LS', ...
'PALMER5D', 'RAT43LS', 'ROSZMAN1LS', 'KIRBY2LS', 'MGH17LS', 'OSBORNEA', ...
'BIGGS6', 'LANCZOS1LS', 'LANCZOS2LS', 'LANCZOS3LS', 'PALMER5C', ...
'HAHN1LS', 'PALMER1D', 'THURBERLS', 'GAUSS1LS', 'GAUSS2LS', 'GAUSS3LS', ...
'PALMER1C', 'PALMER2C', 'PALMER3C', 'PALMER4C', 'PALMER6C', 'PALMER7C', ...
'PALMER8C', 'VIBRBEAM', 'VESUVIALS', 'VESUVIOLS', 'VESUVIOULS', ...
'ENSOLS', 'OSBORNEB', ...
};

n_problems = length(unconstrained_problems);
fminunc_options = optimoptions('fminunc', 'Display', 'off', ...
                               'SpecifyObjectiveGradient', false, ...
                               'Algorithm', 'quasi-newton');

terminate_cutest_problem()
clear global problem_path_cutest problem_name_cutest problem_data_cutest
global problem_data_cutest
clear results_unconstrained

warning('off', 'cmg:ill_conditioned_system')
warning('off', 'cmg:trial_not_decrease');
warning('off', 'cmg:geometry_degenerating');

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

    [x_fmin, fvalue_fmin] = fminunc(f, x0, fminunc_options);
    f_count_fmincon = counter.get_count();
    counter.reset_count();
    
    results_unconstrained(n).ref.fx = fvalue_fmin;
    results_unconstrained(n).ref.count = f_count_fmincon;
    results_unconstrained(n).ref.viol = 0;

    options = [];

    try
        [x_trust, fvalue_trust] = trust_region({f}, x0, f(x0), [], [], options);
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