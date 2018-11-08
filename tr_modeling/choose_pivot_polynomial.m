function [pivot_polynomials, pivot_value, success] = choose_pivot_polynomial(pivot_polynomials, points, initial_i, final_i, tol)

    last_point = initial_i - 1;
    incumbent_point = points(:, initial_i);
    success = false;
    pivot_value = 0;
    for k = initial_i:final_i
        polynomial = orthogonalize_to_other_polynomials(pivot_polynomials, k, points, last_point);
        val = evaluate_polynomial(polynomial, incumbent_point);
        if abs(val) > tol
            % Accept polynomial
            success = true;
            % Swap polynomials
            pivot_polynomials(k) = pivot_polynomials(initial_i);
            pivot_polynomials(initial_i) = polynomial;
            pivot_value = val;
            break
        else
            false; % We won't even save the orthogonalization
        end
    end
    

end