function p = combine_polynomials(polynomials, coefficients)

    terms = length(polynomials);
    if terms == 0 || length(coefficients) ~= terms
        error();
    end

    p = multiply_p(polynomials(1), coefficients(1));
    for k = 2:terms
        p = add_p(p, multiply_p(polynomials(k), coefficients(k)));
    end

end