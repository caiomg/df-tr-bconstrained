function model = rebuild_model(model, options)
% REBUILD_MODEL - re-selects the points to include in the model

    % Region considered to use points
    radius_factor = options.radius_factor;
    radius_factor_linear_block = radius_factor*options.radius_factor_extra_tol;

    % Threshold for pivot polynomials
    pivot_threshold_rel = options.pivot_threshold;
    radius = model.radius;
    pivot_threshold = pivot_threshold_rel*min(1, radius);

    % All points we know
    points_abs = [model.points_abs, model.cached_points];
    fvalues = [model.fvalues, model.cached_fvalues];
    [dim, p_ini] = size(points_abs); % dimension, number of points
    
    % Center will be first
    points_abs(:, [1, model.tr_center]) = points_abs(:, [model.tr_center, 1]);
    fvalues(:, [1, model.tr_center]) = fvalues(:, [model.tr_center, 1]);
    
    % Calculate distances
    points_shifted = zeros(dim, p_ini);
    distances = zeros(1, p_ini);
    for n = 2:p_ini
        % Shift all points to TR center
        % Calculate distances
        points_shifted(:, n) = points_abs(:, n) - points_abs(:, 1);
        distances(n) = norm(points_shifted(:, n), inf); % or 2 norm
    end
    % Reorder
    [distances, pt_order] = sort(distances);
    points_shifted = points_shifted(:, pt_order);
    points_abs = points_abs(:, pt_order);
    fvalues = fvalues(:, pt_order);
    % Remove duplicates
    to_remove = false(p_ini, 1);
    for n = 2:p_ini
       if distances(n) == distances(n-1) ...
               && norm(points_abs(:, n) - points_abs(:, n-1), inf) == 0
           to_remove(n) = true;
       end
    end
    if ~isempty(find(to_remove, 1))
       points_shifted(:, to_remove) = [];
       points_abs(:, to_remove) = [];
       fvalues(:, to_remove) = [];
       p_ini = size(points_abs, 2);
    end
    
    % Building model
    pivot_polynomials = nfp_basis(dim); % Basis of Newton
                                        % Fundamental Polynomials
    polynomials_num = length(pivot_polynomials);
    pivot_values = zeros(1, polynomials_num);
    % Constant term
    last_pt_included = 1;
    pivot_values(1) = 1;
    poly_i = 2;
    for iter = 2:polynomials_num
        % Gaussian elimination (using previuos points)
        pivot_polynomials(poly_i) = ...
            orthogonalize_to_other_polynomials(pivot_polynomials, ...
                                                    poly_i, ...
                                                    points_shifted, ...
                                                    last_pt_included);
        if poly_i <= dim+1
            block_beginning = 2;
            block_end = dim+1;
            % Linear block -- we allow more points (*2)
            maxlayer = min(radius_factor_linear_block, distances(end)/radius);
            if iter > dim+1
                % We already tested all linear terms
                % We do not have points to build a FL model
                % How did this happen??? see Comment [1]
                break
            end
        else
            % Quadratic block -- being more carefull
            maxlayer = min(radius_factor, distances(end)/radius);
            block_beginning = dim+2;
            block_end = polynomials_num;
        end
        maxlayer = max(1, maxlayer);
        all_layers = linspace(1, maxlayer, ceil(maxlayer));
        max_absval = 0;
        pt_max = 0;
        for layer = all_layers
            dist_max = layer*radius;
            for n = last_pt_included+1:p_ini
                if distances(n) > dist_max
                    break % for(n)
                end
                val = evaluate_polynomial(pivot_polynomials(poly_i), ...
                                                 points_shifted(:, n));
                val = val/dist_max; % minor adjustment
                if abs(max_absval) < abs(val)
                    max_absval = val;
                    pt_max = n;
                end
            end
            if abs(max_absval) > pivot_threshold
                break % for(layer)
            end
        end
        if abs(max_absval) > pivot_threshold
            % Point accepted
            pt_next = last_pt_included + 1;
            pivot_values(pt_next) = max_absval;
            points_shifted(:, [pt_next:end]) = ...
                points_shifted(:, [pt_max, pt_next:pt_max-1, pt_max+1:end]);
            points_abs(:, [pt_next:end]) = ...
                points_abs(:, [pt_max, pt_next:pt_max-1, pt_max+1:end]);
            fvalues(:, pt_next:end) = ...
                fvalues(:, [pt_max, pt_next:pt_max-1, pt_max+1:end]);
            distances(pt_next:end) = ...
                distances([pt_max, pt_next:pt_max-1, pt_max+1:end]);
            
            % Normalize polynomial value
            pivot_polynomials(poly_i) = ...
                normalize_polynomial(pivot_polynomials(poly_i), ...
                                     points_shifted(:, pt_next));
            % Re-orthogonalize (just to make sure it still assumes
            % 0 in previous points). Unnecessary in infinite precision
            pivot_polynomials(poly_i) = ...
                orthogonalize_to_other_polynomials(pivot_polynomials, ...
                                                        poly_i, ...
                                                        points_shifted, ...
                                                        last_pt_included);
            
            % Orthogonalize polynomials on present block (deffering
            % subsequent ones)
            pivot_polynomials = orthogonalize_block(pivot_polynomials, ...
                                                    points_shifted(:, ...
                                                              poly_i), ...
                                                    poly_i, ...
                                                    block_beginning, poly_i);
            already_included(pt_max) = true;
            last_pt_included = pt_next;
            poly_i = poly_i + 1;
        else
            % These points don't render good pivot value for this
            % specific polynomial
            % Exchange some polynomials and try to advance...
            % Moving this polynomial to the end of the block
            pivot_polynomials(poly_i:block_end) = ...
                pivot_polynomials([poly_i+1:block_end, poly_i]);
            % Comment [1]:
            % If we are on the linear block, this means we won't be
            % able to build a Fully Linear model
        end
    
    end
    
    model.tr_center = 1;
    model.points_abs = points_abs(:, 1:last_pt_included);
    model.points_shifted = points_shifted(:, 1:last_pt_included);
    model.fvalues = fvalues(:, 1:last_pt_included);
    cache_size = min(p_ini - last_pt_included, model.cache_max);
    model.pivot_polynomials = pivot_polynomials;
    model.pivot_values = pivot_values;
    model.modeling_polynomials = {};
    % Points not included
    model.cached_points = points_abs(:, last_pt_included+1:last_pt_included ...
                                     + cache_size);
    model.cached_fvalues = fvalues(:, last_pt_included+1:last_pt_included ...
                                   + cache_size);

end


