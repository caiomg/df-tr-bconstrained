function M = interpolation_matrix(polynomials, points)

    M = zeros(size(points, 2), length(polynomials));
    for m = 1:size(points, 2)
        for n = 1:length(polynomials)
            M(m, n) = evaluate_polynomial(polynomials(n), points(:, m));
        end
    end
end