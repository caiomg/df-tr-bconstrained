function polynomials = recompute_polynomial_models(shifted_points, fvalues, basis)

    [dim, points_num] = size(shifted_points);
    n_interpolating_functions = size(fvalues, 1);
    linear_basis_size = dim + 1;
    basis_size = length(basis);
    
    if points_num <= linear_basis_size
        % Linear interpolation

        % Calculate coefficients of polynomial model:
        M = zeros(points_num, linear_basis_size);
        for m = 1:points_num
            for n = 1:linear_basis_size
                M(m, n) = evaluate_polynomial(basis(n), shifted_points(:, m));
            end
        end
        
        [Q, R] = qr(M');
        l_opts.LT = true;
        u = linsolve(R(1:points_num, :)', fvalues', l_opts);
        linear_coefficients = Q(:, 1:points_num)*u;
        quadratic_coefficients = zeros(basis_size - linear_basis_size, ...
                                       n_interpolating_functions);
        coefficients_basis = [linear_coefficients;
                              quadratic_coefficients];
        %% DEBUG
        a_test = M\fvalues';
        linsys_error = linear_coefficients - a_test;
        if norm(linsys_error, inf) > 1e-5
            warning('cmg:interp_error', 'Possible interpolation error');
        end
        %%
        
    else
        % Quadratic interpolation

        % Calculate coefficients of polynomial model:
        M = zeros(points_num, basis_size);
        for m = 1:points_num
            for n = 1:basis_size
                M(m, n) = evaluate_polynomial(basis(n), shifted_points(:, m));
            end
        end
        M0 = M(:, 1);
        Ml = M(2:end, 2:linear_basis_size);
        Mq = M(2:end, linear_basis_size+1:end);
        
        % Constant coefficient
        constant_coefficient = fvalues(:, 1)';
        
        Fv = fvalues(:, 2:end)' - M0(2:end)*constant_coefficient;
        % Minimum norm solution
        Mmn = [Mq*Mq', Ml;
               Ml', zeros(linear_basis_size-1:linear_basis_size-1)];
        fmn = [Fv;
               zeros(linear_basis_size-1, n_interpolating_functions)];
        sym_opts.SYM = true;
        [sol_mn, M_rcond] = linsolve(Mmn, fmn, sym_opts);
        
        % Linear coefficients
        linear_coefficients = sol_mn(points_num:end, :);
        
        % Quadratic coefficients
        quadratic_coefficients = Mq'*sol_mn(1:points_num-1, :);

        coefficients_basis = [constant_coefficient;
                              linear_coefficients;
                              quadratic_coefficients];
                         
        %% TESTING
        if M_rcond < eps()*1e4
            % This shouldn't happen
            warning('cmg:badly_conditioned_system', 'Badly conditioned system');
        end
        f_test = @(x) 0.5*(norm(x(linear_basis_size+1:end))^2);
        b1 = fvalues(1, :)';
        a0 = M\b1;
        fmincon_opts = optimoptions('fmincon', 'Display', 'off', ...
                                       'SpecifyObjectiveGradient', false);
%         a_test = fmincon(f_test, a0, [], [], M, b1, [], [], [], fmincon_opts);
%         coeff_mn_error = coefficients_basis(:, 1) - a_test;
%         coefficients_basis(:, 1) = a_test; %REMOVE!!!!!!!!
        N = null(Ml');
        try
            chol(N'*(Mq*Mq')*N);
        catch this_error
            rethrow(this_error);
        end
        %%

    end


    for nf = n_interpolating_functions:-1:1
        c_polynomial = polynomial_zero(dim);
        for m = 1:basis_size
           c_polynomial = add_p(c_polynomial, ...
                                    multiply_p(basis(m), coefficients_basis(m, nf)));
        end
        polynomials{nf} = c_polynomial;
    end