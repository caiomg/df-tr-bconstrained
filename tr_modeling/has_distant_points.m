function result = has_distant_points(model, options)
%HAS_DISTANT_POINTS returns wether a model has points distant from center
%    has_distant_points(model, options)

    radius_factor = options.radius_factor;
    radius_factor_extra = radius_factor*options.radius_factor_extra_tol;
    radius = model.radius;
    points_abs = model.points_abs;
    center_x = model.center_point();
    pivot_values = model.pivot_values;
    [dim, point_num] = size(points_abs);
    linear_terms_num = dim + 1;

    % Distance allowed for points in general
    allowed_distance = radius*radius_factor;
    % Extra-distance allowed for completing the linear block
    allowed_distance_extra = radius*radius_factor_extra;

    pt_i = 0;
    for n = 1:length(pivot_values)
        if ~isfinite(pivot_values(n))
            % Just advance n out of this block
            continue
        elseif pivot_values(n) == 0
            % The end. But should've breaked below!
            % If this ever happens, it is a BUG
            error('cmg:pivot_zero', ...
                  'Found pivot zero associated with a point');
        else
            pt_i = pt_i + 1;
            distance = norm(points_abs(:, pt_i) - center_x, inf);
            if distance > allowed_distance_extra ...
                    || (pt_i > linear_terms_num && distance > allowed_distance)
                result = true;
                break
            end
        end
        if pt_i == point_num
            result = false;
            break
        end
    end
end
