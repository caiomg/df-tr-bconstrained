function basis = diagonal_basis(dimension)
% Quadratic polynomial basis without crossed products (diagonal
% Hessian)

n_polynomials = 2*dimension + 1;
n_terms =  (dimension + 1)*(dimension + 2)/2;


basis(n_polynomials) = polynomial_zero(dimension);
idx = 1;
idx_sum = 1;
for k = 1:n_polynomials
    coefficients_current = zeros(n_terms, 1);
    coefficients_current(idx) = 1;
    if (idx  < dimension + 1)
        idx = idx+1;
    else
        idx = idx + idx_sum;
        idx_sum = idx_sum + 1;
    end
    basis(k).coefficients = coefficients_current;
    basis(k).dimension = dimension;
end

end


