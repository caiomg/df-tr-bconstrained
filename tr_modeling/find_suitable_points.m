function [new_point, new_value] = find_suitable_points(polynomial, bl, bu, attempt)

    if attempt == 1
        % Minimize inside TR
        [new_point, new_value] = ...
           minimize_polynomial_with_bounds(polynomial, bl, bu, smaller_radius);
    elseif attempt == 2
        % Maximize inside TR
        [new_point, new_value] = ...
            minimize_polynomial_with_bounds(multiply_p(polynomial, ...
                                   -1), bl, bu, smaller_radius);
    elseif attempt == 3
        % Minimize in larger region
        [new_point, new_value] = minimize_polynomial_with_bounds(polynomial, bl, bu, 1);
    elseif attempt == 4
        % Maximize in larger region
        [new_point, new_value] = ...
            minimize_polynomial_with_bounds(multiply_p(polynomial, ...
                                   -1), bl, bu, 1);
    elseif attempt == 5
        [cp, gp, Hp] = coefficients_to_matrices(polynomial.dimension, polynomial.coefficients);
        if norm(gp, inf) == 0
            new_point = gp;
            new_value = evaluate_polynomial(polynomial, new_point);
        else
            new_point = round(gp/max(gp));
            new_point = new_point/norm(new_point);
            new_point = min(bu, max(bl, new_point));
            new_value = evaluate_polynomial(polynomial, new_point);
        end
    elseif attempt == 6
        [cp, gp, Hp] = coefficients_to_matrices(polynomial.dimension, polynomial.coefficients);
        if norm(gp, inf) == 0
            new_point = gp;
            new_value = evaluate_polynomial(polynomial, new_point);
        else
            new_point = round(gp/max(gp));
            new_point = -new_point/norm(new_point);
            new_point = min(bu, max(bl, new_point));
            new_value = evaluate_polynomial(polynomial, new_point);
        end
    else
        [cp, gp, Hp] = coefficients_to_matrices(polynomial.dimension, polynomial.coefficients);
        [~, v1] = max(Hp, [], 1);
        [~, v2] = max(Hp, [], 2);
        if v1 ~= v2
            attempt = attempt + 2;
        end
        if attempt == 7
            new_point = zeros(dim, 1);
            new_point(v1) = min(bu(v1), max(bl(v1), 1));
            new_value = evaluate_polynomial(polynomial, new_point);
        elseif attempt == 8
            new_point = zeros(dim, 1);
            new_point(v1) = min(bu(v1), max(bl(v1), -1));
            new_value = evaluate_polynomial(polynomial, new_point);
        elseif attempt == 9
            
        end
        
    elseif attempt == 8
        [cp, gp, Hp] = coefficients_to_matrices(polynomial.dimension, polynomial.coefficients);
        [~, v1] = max(Hp, [], 1);
        [~, v2] = max(Hp, [], 2);

        new_point = zeros(dim, 1);
        new_point(v2) = min(bu(v2), max(bl(v2), 1));
        new_value = evaluate_polynomial(polynomial, new_point);
    elseif attempt == 9
        if v1 == v2
            continue
        end
        new_point = zeros(dim, 1);
        new_point(v1) = 1/sqrt(2);
        new_point(v2) = 1/sqrt(2);
        new_point = min(bu, max(bl, 1));
        new_value = evaluate_polynomial(polynomial, new_point);
    else
        error('cmg:pivot_not_found', 'Pivot not found');
    end
end