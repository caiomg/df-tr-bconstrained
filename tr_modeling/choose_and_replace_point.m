function [model, success] = choose_and_replace_point(model, funcs, bl, bu, options)
   
    pivot_threshold = options.exchange_threshold;
    radius = model.radius;
    
    % pivot_threshold = min(1, radius)*rel_pivot_threshold;

    points_shifted = model.points_shifted;
    [dim, p_ini] = size(points_shifted);
    shift_center = model.points_abs(:, 1);
    tr_center = model.tr_center;
    tr_center_x = model.points_shifted(:, tr_center);

    pivot_values = model.pivot_values;
    pivot_polynomials = model.pivot_polynomials;
    
    points_shifted = model.points_shifted;
    [dim, points_num] = size(points_shifted);
    linear_terms = dim+1;
    
    if isempty(bl)
        bl = -inf(dim, 1);
    end
    if isempty(bu)
        bu = inf(dim, 1);
    end
    bl_shifted = bl - shift_center;
    bu_shifted = bu - shift_center;
    unshift_point = @(x) max(min(x + shift_center, bu), bl);
    tol_shift = 10*eps(max(1, norm(shift_center, inf)));

    
    [~, piv_order] = sort(abs(pivot_values(1:points_num)));
    polynomials_num = length(pivot_polynomials);
   
    % Could iterate through all pivots, but will try just dealing
    % with the worst one
    pos = piv_order(1);
    if (pos == 1 || pos == tr_center || ...
        (pos <= linear_terms && points_num > linear_terms))
        % Better to just rebuild model
        success = false;
    else
        current_pivot_value = pivot_values(pos);
        [new_points_shifted, new_pivots, point_found] = ...
            point_new(pivot_polynomials(pos), tr_center_x, radius, ...
                      bl_shifted, bu_shifted, pivot_threshold);
        if point_found
            for found_i = 1:size(new_points_shifted, 2)
                new_point_shifted = new_points_shifted(:, found_i);
                new_pivot_value = new_pivots(found_i);
                new_point_abs = unshift_point(new_point_shifted);
                [new_fvalues, f_succeeded] = ...
                    evaluate_new_fvalues(funcs, new_point_abs);
                if f_succeeded
                    break
                else
                    true;
                end
            end
            if f_succeeded
                % Normalize polynomial value
                pivot_polynomials(pos) = ...
                    normalize_polynomial(pivot_polynomials(pos), ...
                                         new_point_shifted);
                
                % Orthogonalize polynomials on present block (all)
                if pos <= dim + 1
                    block_beginning = 2;
                    block_end = min(points_num, dim + 1);
                else
                    block_beginning = dim + 2;
                    block_end = points_num;
                end
                % Re-orthogonalize
                pivot_polynomials(pos) = ...
                    orthogonalize_to_other_polynomials(pivot_polynomials, ...
                                                       pos, ...
                                                       points_shifted, ...
                                                       block_end);
                % Orthogonalize block
                pivot_polynomials = orthogonalize_block(pivot_polynomials, ...
                                                        new_point_shifted, ...
                                                        pos, block_beginning, ...
                                                        points_num);
                
                % Update model and recompute polynomials
                cache_size = size(model.cached_points, 2);
                cache_size = min(cache_size, model.cache_max - 1);
                model.cached_points = [model.points_abs(:, pos), model.cached_points];
                model.cached_fvalues = [model.fvalues(:, pos), model.cached_fvalues];
                model.points_abs(:, pos) = new_point_abs;
                model.points_shifted(:, pos) = new_point_shifted;
                model.fvalues(:, pos) = new_fvalues;
                model.pivot_polynomials = pivot_polynomials;
                model.pivot_values(:, pos) = new_pivot_value*model.pivot_values(:, pos);
                model.modeling_polynomials = {};
                success = true;
            else
               success = false;
               % Also could try to optimize again
            end
        else
            success = false;
        end
    end
end
