function keep = discard_points_only(points, discard_radius)
% DISCART_POINTS removes from set of interpolation points that are
% far from center. Does not update polynomial model

[~, n_points] = size(points);
center = points(:, 1);

keep = true(1, n_points); % points to keep
for k = 2:n_points
    if (norm(points(:, k) - center, 2) > discard_radius)
        keep(k) = false;
    end
end

end