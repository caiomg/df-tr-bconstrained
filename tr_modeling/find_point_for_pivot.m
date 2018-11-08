function [status, point_abs, fvalues, f_evaluations] = ...
        find_point_for_pivot(polynomial, bl, bu, radius, ff, tol_pivot, ...
                             unshift_point)
    
    smaller_radius = 0.618*radius;
    n_interpolating_functions = length(ff);
    neg_polynomial = multiply_p(polynomial, -1);
    status = false;
    f_evaluations = 0;
    max_attempts = 3;
    for attempt = 1:max_attempts
        if attempt == 1
            % Minimize inside TR
            [new_point_a, new_value_a] = ...
                minimize_polynomial_with_bounds(polynomial, bl, bu, radius);
            % Maximize inside TR
            [new_point_b, new_value_b] = ...
                minimize_polynomial_with_bounds(neg_polynomial, bl, ...
                                                bu, radius);
            if abs(new_value_a) >= abs(new_value_b)
                new_point = new_point_a;
                new_value = new_value_a;
            else
                new_point = new_point_b;
                new_value = new_value_b;
            end
        elseif attempt == 2
            % Minimize inside TR
            [new_point_a, new_value_a] = ...
                minimize_polynomial_with_bounds(polynomial, bl, bu, ...
                                                smaller_radius);
            % Maximize inside TR
            [new_point_b, new_value_b] = ...
                minimize_polynomial_with_bounds(neg_polynomial, bl, ...
                                                bu, smaller_radius);
            if abs(new_value_a) >= abs(new_value_b)
                new_point = new_point_a;
                new_value = new_value_a;
            else
                new_point = new_point_b;
                new_value = new_value_b;
            end
        elseif attempt == 3
            [new_point, new_value] = ...
                find_other_point_with_bounds(polynomial, bl, bu, radius);
        end
        if abs(new_value) >= tol_pivot
            % Good pivot value, worth evaluating df function
            point_abs = unshift_point(new_point);
            fvalues = zeros(n_interpolating_functions, 1);
            for nf = 1:n_interpolating_functions
                fvalues(nf, 1) = ff{nf}(point_abs);
                f_evaluations = f_evaluations + 1;
                if ~isfinite(fvalues(nf, 1))
                    status = false;
                    break % for(nf)
                end
            end
            status = true;
            break % for(attempt)
        end
    end
    if ~status
        point_abs = [];
        fvalues = [];
    end
