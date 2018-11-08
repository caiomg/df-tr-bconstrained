function polynomial = polynomial_zero(dimension)
% POLYNOMIAL_ZERO returns a polynomial with zero coefficients

polynomial.dimension = dimension;
n_terms =  (dimension + 1)*(dimension + 2)/2;
polynomial.coefficients = zeros(n_terms, 1);

end

