function [model, succeeded, pt_i] = ...
        exchange_point(model, new_point, new_fvalues, ...
                       relative_pivot_threshold, allow_exchange_center)

    pivot_threshold = min(1, model.radius)*relative_pivot_threshold;

    [dim, last_p] = size(model.points_abs);
    pivot_polynomials = model.pivot_polynomials;
    center_i = model.tr_center;
    
    if last_p < 2
        error('cmg:runtime', 'Error here');
    end
    
    shift_center = model.points_abs(:, 1);
    new_point_shifted = new_point - shift_center;
    points_shifted = model.points_shifted;
    if last_p <= dim + 1
        block_beginning = 2;
        block_end = min(dim+1, last_p);
    else
        block_beginning = dim+2;
        block_end = last_p;
    end
    max_val = 0;
    max_poly_i = 0;
    for poly_i = block_end:-1:block_beginning
        if allow_exchange_center || poly_i ~= center_i
            val = model.pivot_values(poly_i)*...
                  evaluate_polynomial(pivot_polynomials(poly_i), new_point_shifted);
            if abs(max_val) < abs(val)
                max_val = val;
                max_poly_i = poly_i;
            end
        end
    end
    if max_poly_i > 0
        new_pivot_val = model.pivot_values(max_poly_i)*max_val;
        if ~isfinite(new_pivot_val) ...
           && isfinite(max_val) && isfinite(model.pivot_values(max_poly_i))
             % adjustment
             new_pivot_val = sign(new_pivot_val)*realmax;
             warning('cmg:geometry_degenerating', ...
                     'Bad geometry of interpolation set for machine precision');
        end
    else
        new_pivot_val = 0;
    end
    if abs(new_pivot_val) > pivot_threshold
        points_shifted(:, max_poly_i) = new_point_shifted;
        % Normalize polynomial value
        pivot_polynomials(max_poly_i) = ...
            normalize_polynomial(pivot_polynomials(max_poly_i), ...
                                 new_point_shifted);
        % Re-orthogonalize
        pivot_polynomials(max_poly_i) = ...
            orthogonalize_to_other_polynomials(pivot_polynomials, ...
                                               max_poly_i, points_shifted, ...
                                               last_p);

        % Orthogonalize polynomials on present block (until end)
        pivot_polynomials = orthogonalize_block(pivot_polynomials, ...
                                                new_point_shifted, ...
                                                max_poly_i, ...
                                                block_beginning, last_p);

        % Update model
        cache_size = size(model.cached_points, 2);
        cache_size = min(cache_size, model.cache_max - 1);
        model.cached_points = [model.points_abs(:, max_poly_i), ...
                               model.cached_points(:, 1:cache_size)];
        model.cached_fvalues = [model.fvalues(:, max_poly_i), ...
                                model.cached_fvalues(:, 1:cache_size)];
        model.points_abs(:, max_poly_i) = new_point;
        model.fvalues(:, max_poly_i) = new_fvalues;
        model.points_shifted = points_shifted;
        model.pivot_polynomials = pivot_polynomials;
        model.pivot_values(:, max_poly_i) = new_pivot_val;
        model.modeling_polynomials = {};
        succeeded = true;
        pt_i = max_poly_i;
    else
        succeeded = false;
        pt_i = 0;
    end

end
