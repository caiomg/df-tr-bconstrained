function [point, value] = minimize_polynomial_with_bounds(polynomial, bl, ...
                                                          bu, radius)


dimension = polynomial.dimension;
coefficients = polynomial.coefficients;

[c0, g, H] = coefficients_to_matrices(dimension, coefficients);

x0 = zeros(dimension, 1);
s0 = ms_step(H, g, radius);
[point, decr] = pg_path_bounds(H, g, x0, s0, bl, bu, radius);

npoint = norm(point);
if npoint > radius
    point = (point/npoint)*radius;
    value = c0 + g'*point + 0.5*(point'*H*point);
else
    value = c0 + decr;
end

end