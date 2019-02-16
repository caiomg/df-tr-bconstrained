function polynomial = normalize_polynomial(polynomial, point)

    val = evaluate_polynomial(polynomial, point);
    for k = 1:3
        if (val - 1) == 0
            break
        elseif val == 0
            warning('cmg:normalize', 'Can not normalize with 0');
            break
        end
        polynomial = multiply_p(polynomial, 1/val);
        val = evaluate_polynomial(polynomial, point);
    end
end