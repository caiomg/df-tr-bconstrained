function basis = natural_basis(dimension, n_polynomials)
% NATURAL_BASIS calculates the natural basis of monomials

% Number of terms
n_terms =  (dimension + 1)*(dimension + 2)/2;
% Polynomials to be included
if (nargin == 1)
    n_polynomials = n_terms;
end
% Checkint number of polynomials
if(n_polynomials > n_terms)
    error('cmg:npolynomials', 'too high number of polynomials');
end

% Calculating basis of polynomials
basis(n_polynomials) = polynomial_zero(dimension);
for k = 1:n_polynomials
    coefficients_current = zeros(n_terms, 1);
    coefficients_current(k) = 1;
    basis(k).coefficients = coefficients_current;
    basis(k).dimension = dimension;
end

end