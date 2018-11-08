function [model, exitflag] = improve_model(model, ff, bl, bu, options, complete_only)
% IMPROVE_MODEL improves model to be lambda-poised (for some
% lambda)

    dim = size(model.points, 1);
    if isempty(bl)
        bl = -inf(dim, 1);
    end
    if isempty(bu)
        bu = inf(dim, 1);
    end
    if nargin < 6
        complete_only = false;
    end

    if strcmp(options.basis, 'dummy')
        return;
    end
    basis = model.basis;
    basis_size = length(basis);

    % pivot_threshold is constant through the algorithm. Ensures we can
    % obtain lambda-poised models when needed
    tol_pivot = options.pivot_threshold;
    radius_factor = options.poised_radius_factor;
    xi_imp = options.xi_imp;

    n_interpolating_functions = length(ff);

    radius = model.radius;
    points_abs = model.points;

    keep_points = discard_points_only(points_abs, radius*radius_factor);
    points_abs = points_abs(:, keep_points);
    fvalues = model.fvalues(:, keep_points);

    center = points_abs(:, 1);

    % Scaling of points and bounds
    [points_scaled, scale_factor_x]  = scale_interpolation_points(points_abs, ...
                                                      center, radius);
    shift_point = @(x) (x - center)/scale_factor_x;
    unshift_point = @(x) project_to_bounds(x*scale_factor_x + center, bl, bu);

    bl_scaled = shift_point(bl);
    bl_scaled = min(bl_scaled, 0); % Making sure to include center
    bu_scaled = shift_point(bu);
    bu_scaled = max(bu_scaled, 0); % Making sure to include center

    % A tolerance for the shifting, so new points indeed fall inside TR
    tol_shift = 10*eps(max(1, max(abs(center))))/min(1, scale_factor_x);
    if tol_shift > 1e-3
        tol_shift = 1e-3; % DANGER, still
    end
    larger_radius = 1 - tol_shift;

    %%
    linear_basis_size = dim + 1;
    basis_lin = basis(1:linear_basis_size);
    basis_quad = basis(linear_basis_size+1:end);

    quadratic_pivots_num = 0;
    interpolation_set_changed = false;

    %%
    while true
        [indices_pts, pivot_absvalues_lin, linear_pivot_polynomials, poly_perm, ...
              quadratic_pivot_polynomials] = check_geometry_linear_mn(points_scaled, ...
                                                          basis_lin, ...
                                                          basis_quad, tol_pivot);
        % Reordering information
        points_abs = points_abs(:, indices_pts);
        points_scaled = points_scaled(:, indices_pts);
        fvalues = fvalues(:, indices_pts);
        % Finding point to change
        pivot_index = find(pivot_absvalues_lin < tol_pivot, 1);
        linear_pivots_num = sum(pivot_absvalues_lin >= tol_pivot);


        % reordering...
        linear_pivot_polynomials = linear_pivot_polynomials(poly_perm);
        if ~isempty(pivot_index) && ~complete_only
            chosen_poly = linear_pivot_polynomials(pivot_index);
            tol_pivot_replacing = max(tol_pivot, ...
                                      pivot_absvalues_lin(pivot_index)*xi_imp);

            [pivot_found, new_point_abs, new_fvalues] = ...
                find_point_for_pivot(chosen_poly, bl_scaled, bu_scaled, ...
                                     larger_radius, ff, tol_pivot_replacing, ...
                                     unshift_point);
            if pivot_found
                new_point = shift_point(new_point_abs);
                points_scaled = [points_scaled(:, 1:pivot_index-1), ...
                                 new_point, ...
                                 points_scaled(:, pivot_index:end)];
                points_abs = [points_abs(:, 1:pivot_index-1), ...
                              new_point_abs, ...
                              points_abs(:, pivot_index:end)];
                fvalues = [fvalues(:, 1:pivot_index-1), ...
                           new_fvalues, ...
                           fvalues(:, pivot_index:end)];

                new_value = evaluate_polynomial(chosen_poly, new_point);
                pivot_absvalues_lin(pivot_index) = abs(new_value);
                % I should double check interpolation set here
                interpolation_set_changed = true;
                linear_pivots_num = linear_pivots_num + 1;
                %% DEBUGING
                [indices_pts2, pivot_absvalues2, linear_pivot_polynomials2, ...
                  quadratic_pivot_polynomials2] = check_geometry_linear_mn(points_scaled, ...
                                                              basis_lin(poly_perm), ...
                                                              basis_quad, tol_pivot);
                %%
                1;
            else
                warning('cmg:bad_fvalue', 'Bad f value'); %changeme
                status = false;
            end
        end
        break
        if linear_pivots_num == linear_basis_size
            break
        end
    end
%%

    % Check quadratic terms
    if linear_pivots_num < linear_basis_size
        pivot_absvalues = pivot_absvalues_lin;
    else
        while true
            [indices_pts_quad, pivot_absvalues_quad, ...
             quadratic_pivot_polynomials, polynomials_permutation] = ...
              check_geometry_quadratic_mn(points_scaled(:, ...
                                                        linear_basis_size+1:end), ...
                                          quadratic_pivot_polynomials, tol_pivot);

            points_abs(:, linear_basis_size + 1:end) = points_abs(:, ...
                                                              indices_pts_quad ...
                                                              + linear_basis_size);
            points_scaled(:, linear_basis_size + 1:end) = points_scaled(:, ...
                                                              indices_pts_quad ...
                                                              + linear_basis_size);
            fvalues(:, linear_basis_size + 1:end) = fvalues(:, ...
                                                            indices_pts_quad ...
                                                            + linear_basis_size);
            quadratic_pivots_num = sum(pivot_absvalues_quad > tol_pivot);
            quadratic_pivot_polynomials = ...
                quadratic_pivot_polynomials(polynomials_permutation);

            if ~interpolation_set_changed && ~complete_only
                % Geometry improvement
                pivot_index = find(pivot_absvalues_quad < tol_pivot, 1);
                point_index = pivot_index + linear_basis_size;
                if ~isempty(pivot_index)
                    chosen_poly = quadratic_pivot_polynomials(pivot_index);
                    tol_pivot_replacing = max(tol_pivot, ...
                                              pivot_absvalues_quad(pivot_index)*xi_imp);
                    %%%
                    [pivot_found, new_point_abs, new_fvalues] = ...
                        find_point_for_pivot(chosen_poly, bl_scaled, ...
                                             bu_scaled, larger_radius, ...
                                             ff, tol_pivot_replacing, ...
                                             unshift_point);
                    if pivot_found
                        new_point = shift_point(new_point_abs);
                        points_scaled = [points_scaled(:, 1:point_index-1), ...
                                         new_point, ...
                                         points_scaled(:, point_index:end)];
                        points_abs = [points_abs(:, 1:point_index-1), ...
                                      new_point_abs, ...
                                      points_abs(:, point_index:end)];
                        fvalues = [fvalues(:, 1:point_index-1), ...
                                   new_fvalues, ...
                                   fvalues(:, point_index:end)];

                        new_value = evaluate_polynomial(chosen_poly, new_point);
                        pivot_absvalues_quad(pivot_index) = abs(new_value);
                        % I should double check interpolation set here
                        interpolation_set_changed = true;
                        quadratic_pivots_num = quadratic_pivots_num + 1;
                        %% DEBUGING
                        [indices_pts_quad2, pivot_absvalues_quad2, ...
                         quadratic_pivot_polynomials2, ...
                         polynomials_permutation2] = ...
                            check_geometry_quadratic_mn(points_scaled(:, ...
                                                                      linear_basis_size+1:end), quadratic_pivot_polynomials, tol_pivot);
                        %%
                        1;
                    else
                        warning('cmg:bad_fvalue', 'Bad f value'); %changeme
                        status = false;
                    end                
                    %%%
                end
            end
            pivot_absvalues = [pivot_absvalues_lin, pivot_absvalues_quad];
            break
            if quadratic_pivots_num >= 1
                break
            end
        end
    end


%%%
    linear_interpolation = (quadratic_pivots_num == 0);
    if linear_interpolation
        % Discarding points not used for interpolation
        points_used_num = linear_pivots_num;
        points_abs = points_abs(:, 1:points_used_num);
        points_scaled = points_scaled(:, 1:points_used_num);
        fvalues = fvalues(:, 1:points_used_num);

        % Calculate coefficients of polynomial model:
        M = zeros(points_used_num);
        for m = 1:points_used_num
            for n = 1:linear_basis_size
                M(m, n) = evaluate_polynomial(basis(n), points_scaled(:, m));
            end
        end
    %         [L, U, P] = lu(M);
        l_opts.LT = true;
        u_opts.UT = true;
    %         uf = linsolve(L, P*fvalues', l_opts);

        [Q, R] = qr(M');
        u = linsolve(R(1:points_used_num, :)', fvalues', l_opts);
        linear_coefficients = Q(:, 1:points_used_num)*u;
        quadratic_coefficients = zeros(length(basis_quad), 1);
        %% DEBUG
        a_test = M\fvalues';
        linsys_error = linear_coefficients - a_test;
        if norm(linsys_error, inf) > 1e5
            warning('cmg:interp_error', 'Possible interpolation error');
        end
        %%
        
        coefficients_basis = [linear_coefficients;
                              quadratic_coefficients];
        pivot_polynomials = linear_pivot_polynomials;
    else
        points_used_num = linear_pivots_num + quadratic_pivots_num;
        points_abs = points_abs(:, 1:points_used_num);
        points_scaled = points_scaled(:, 1:points_used_num);
        fvalues = fvalues(:, 1:points_used_num);

        % Calculate coefficients of polynomial model:
        M = zeros(points_used_num);
        for m = 1:points_used_num
            for n = 1:basis_size
                M(m, n) = evaluate_polynomial(basis(n), points_scaled(:, m));
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
        linear_coefficients = sol_mn(points_used_num:end, :);
        
        % Quadratic coefficients
        quadratic_coefficients = Mq'*sol_mn(1:points_used_num-1, :);

        coefficients_basis = [constant_coefficient;
                              linear_coefficients;
                              quadratic_coefficients];
        pivot_polynomials = [linear_pivot_polynomials, ...
                             quadratic_pivot_polynomials];
                         
        %% TESTING
        if M_rcond < sqrt(eps)
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
        %%
        N = null(Ml');
        try
            chol(N'*(Mq*Mq')*N);
        catch this_error
            rethrow(this_error);
        end

    end

    model_polynomial = polynomial_zero(dim);
    for m = 1:basis_size
       model_polynomial = add_p(model_polynomial, ...
                                multiply_p(basis(m), coefficients_basis(m, 1)));
    end
    for nf = n_interpolating_functions:-1:2
        c_polynomial = polynomial_zero(dim);
        for m = 1:basis_size
           c_polynomial = add_p(c_polynomial, ...
                                    multiply_p(basis(m), coefficients_basis(m, nf)));
        end
        other_polynomials{nf-1} = c_polynomial;
    end

    model.model_polynomial = model_polynomial;
    if n_interpolating_functions > 1
        model.other_polynomials = other_polynomials;
    end
    model.points = points_abs;
    model.fvalues = fvalues;
    model.scale_factor_x = scale_factor_x;
    model.pivot_absvalues = pivot_absvalues(pivot_absvalues > 0);
    model.pivot_polynomials = pivot_polynomials;

    model.smallest_pivot = min(model.pivot_absvalues);
    if model.smallest_pivot >= tol_pivot
        model.poised_radius = radius;
        model.poised_center = model.points(:, 1);
    end

    exitflag = 0;






