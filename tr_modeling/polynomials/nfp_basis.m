function basis = nfp_basis(dimension)

    % Number of terms
    poly_num =  (dimension + 1)*(dimension + 2)/2;
    linear_size = dimension + 1;
    
    % Calculating basis of polynomials
    basis(poly_num) = polynomial_zero(dimension); % initializing
    for k = 1:linear_size
        coefficients_current = zeros(poly_num, 1);
        coefficients_current(k) = 1;
        basis(k).coefficients = coefficients_current;
        basis(k).dimension = dimension;
    end
    
    % Quadratic entries
    c0 = 0;
    g0 = zeros(dimension, 1);
    m = 1;
    n = 1;
    for poly_i = dimension+2:poly_num
        H = zeros(dimension);
        if m == n % diagonal
            H(m, n) = 2;
        else
            H(m, n) = 1;
            H(n, m) = 1;
        end
        basis(poly_i) = matrices_to_polynomial(c0, g0, H);
        if n < dimension
            n = n + 1;
        else
            m = m + 1;
            n = m;
        end
    end

end