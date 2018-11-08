function [c, g, H] = get_matrices(polynomial)
% Calculates polynomial in the form p(x) = c + g'*x + 0.5*x'*H*x

    dimension = polynomial.dimension;
    coefficients = polynomial.coefficients;

    n_terms =  (dimension + 1)*(dimension + 2)/2;

    [n_coefficients, test_dim] = size(coefficients);
    if (n_coefficients ~= n_terms || test_dim ~= 1)
        error('cmg:wrong_dimension', 'wrong dimension');
    end

    % constant term
    c = coefficients(1);

    % order one term
    idx_coefficients = dimension + 1;
    g = coefficients(2:idx_coefficients);

    % second order term
    H = zeros(dimension);
    for k = 1:dimension
        for m = 1:k
            idx_coefficients = idx_coefficients + 1;
            H(k, m) = coefficients(idx_coefficients);
            H(m, k) = H(k, m);
        end
    end

end
