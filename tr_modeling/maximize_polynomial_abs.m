function [new_points, new_pivot_values, new_points_abs] = ...
        maximize_polynomial_abs(polynomial, x_tr_center, radius, bl, ...
                                bu, shift, unshift)
    
    
    [new_point_min, pivot_min, exitflag_min] = minimize_tr(polynomial, ...
                                                      x_tr_center, ...
                                                      radius, bl, ...
                                                      bu);
    new_point_min_abs = unshift(new_point_min);
    new_point_min = shift(new_point_min_abs);
    pivot_min = evaluate_polynomial(polynomial, new_point_min);
    
    polynomial_max = multiply_p(polynomial, -1);
    [new_point_max, pivot_max, exitflag_max] = minimize_tr(polynomial_max, ...
                                                      x_tr_center, ...
                                                      radius, bl, ...
                                                      bu);
    
    new_point_max_abs = unshift(new_point_max);
    new_point_max = shift(new_point_max_abs);
    pivot_min = evaluate_polynomial(polynomial, new_point_max);
    


    if exitflag_min >= 0
        if (exitflag_max >= 0) && (abs(pivot_max)  >= abs(pivot_min))
            new_points = [new_point_max, new_point_min];
            new_pivot_values = [pivot_max, pivot_min];
            new_points_abs = [new_point_max_abs, new_point_min_abs];
        elseif exitflag_max >= 0 && (abs(pivot_max)  < abs(pivot_min))
            new_points = [new_point_min, new_point_max];
            new_pivot_values = [pivot_min, pivot_max];
            new_points_abs = [new_point_min_abs, new_point_max_abs];
        else
            new_points = new_point_min;
            new_pivot_values = pivot_min;
            point_found = true;
        end
    elseif exitflag_max >= 0
        new_points = new_point_max;
        new_pivot_values = pivot_max;
        new_points_abs = new_point_max_abs;
    else
        point_found = false;
        new_points = [];
        new_pivot_values = 0;
    end
        
end
