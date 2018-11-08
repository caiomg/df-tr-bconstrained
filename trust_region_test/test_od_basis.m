
for dimension = 2:20
    basis = off_diagonal_quadratic(dimension);
    n_polys = length(basis);
    if n_polys ~= (dimension^2 - dimension)/2
        error();
    end
    for k = 1:n_polys
        [c, g, H] = coefficients_to_matrices(basis(k).dimension, basis(k).coefficients);
        if c ~= 0
            error();
        elseif norm(g, inf) ~= 0
            error();
        elseif norm(diag(H), inf) ~= 0
            error();
        elseif norm(H - diag(diag(H)), inf) == 0
            error();
        end
    end
end