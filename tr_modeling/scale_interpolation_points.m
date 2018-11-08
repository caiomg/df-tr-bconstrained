function [points_scaled, scale_factor] = ...
        scale_interpolation_points(points, center, radius)

    center = points(:, 1);
    points_scaled = points;
    max_norm = 0;
    for np = 1:size(points, 2)
        points_scaled(:, np) = points(:, np) - center;
        norm_current = norm(points_scaled(:, np), 2);
        max_norm = max(max_norm, norm_current);
    end
    scale_factor = max(radius, max_norm);
    for np = 1:size(points, 2)
        points_scaled(:, np) = points_scaled(:, np)/scale_factor;
    end
