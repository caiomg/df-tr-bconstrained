function polynomial = matrices_to_polynomial(c, g, H)

[dim,  ~] = size(g);
n_terms =  (dim + 1)*(dim + 2)/2;
coefficients = zeros(n_terms, 1);

% zero order
coefficients(1) = c;
% first order
ind_coefficients = dim + 1;
coefficients(2:ind_coefficients) = g;
% second order
for k = 1:dim
    for m = 1:k
        ind_coefficients = ind_coefficients + 1;
        coefficients(ind_coefficients) = H(k, m);
        if (H(m, k) ~= H(k, m))
            warning('cmg:h_not_symmetric', 'H not symmetric');
        end
    end
end

polynomial.dimension = dim;
polynomial.coefficients = coefficients;

end

