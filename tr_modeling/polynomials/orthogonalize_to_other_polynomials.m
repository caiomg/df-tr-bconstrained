function polynomial = ...
        orthogonalize_to_other_polynomials(all_polynomials, poly_i, points, last_pt)

    polynomial = all_polynomials(poly_i);
    for n = 1:last_pt
        if n ~= poly_i
            polynomial = zero_at_point(polynomial, all_polynomials(n), points(:, n));
        end
    end
    
end

    