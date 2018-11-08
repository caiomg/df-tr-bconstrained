function polynomial = normalize_polynomial(polynomial, point)

    val = evaluate_polynomial(polynomial, point);
    for k = 1:3
        polynomial = multiply_p(polynomial, 1/val);
        val = evaluate_polynomial(polynomial, point);
        if ((val - 1) == 0)
            break
        end
    end