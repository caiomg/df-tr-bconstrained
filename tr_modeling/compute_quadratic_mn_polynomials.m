function polynomials = compute_quadratic_mn_polynomials(points, center_i, fvalues)

    [dim, points_num] = size(points);
    functions_num = size(fvalues, 1);

    
    points_shifted = zeros(dim, points_num-1);
    fvalues_diff = zeros(functions_num, points_num-1);
    m2 = 1;
    for m = [1:center_i-1, center_i+1:points_num]
        points_shifted(:, m2) = points(:, m) - points(:, center_i);
        fvalues_diff(:, m2) = fvalues(:, m) - fvalues(:, center_i);
        m2 = m2 + 1;
    end
    
    M = points_shifted'*points_shifted;
    M = 0.5*(M.^2) + M;
    
    
    % Solve symmetric system
    sym_opts.SYM = true;
    [mult_mn, M_rcond] = linsolve(M, fvalues_diff', sym_opts);
    if M_rcond < 1e4*eps(1)
        % This shouldn't happen
        warning('cmg:badly_conditioned_system', 'Badly conditioned system');
    end
    
    for n = 1:functions_num
        g = zeros(dim, 1);
        H = zeros(dim);
        
        for m = 1:points_num - 1
            g = g + mult_mn(m, n)*points_shifted(:, m);
            H = H + mult_mn(m, n)*(points_shifted(:, m)*points_shifted(:, m)');
        end
        c = fvalues(n, center_i);

        polynomials{n} = matrices_to_polynomial(c, g, H);
    end

end
