function result = has_distant_points(model, options)
%HAS_DISTANT_POINTS returns wether a model has points distant from center
%    has_distant_points(model, options)

    radius_factor = options.radius_factor;
    radius = model.radius;
    points_abs = model.points_abs;
    center_i = model.tr_center;
    points_num = size(points_abs, 2);

    result = false;
    for n = 1:points_num
        distance = norm(points_abs(:, n) - points_abs(:, center_i), inf);
        if distance > radius*radius_factor
            result = true;
            break
        end
    end

end
