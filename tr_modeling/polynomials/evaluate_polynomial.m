function value = evaluate_polynomial(polynomial, point)
% Evaluates polynomial in given point

[c, g, H] = coefficients_to_matrices(polynomial.dimension, ...
                                     polynomial.coefficients);

terms = [c, g'*point, 0.5*(point'*H*point)];
terms = sort(terms);

value = (terms(1) + terms(2)) + terms(3);


end

