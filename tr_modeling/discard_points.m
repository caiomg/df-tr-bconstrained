function model = discard_points(model, discart_radius)
% DISCART_POINTS removes from set of interpolation points that are
% far from center. Does not update polynomial model

center = model.points(:, 1);
points = model.points;
[~, n_points] = size(points);

keep = true(1, n_points); % points to keep
for k = 2:n_points
    if (norm(points(:, k) - center, 2) > discart_radius)
        keep(k) = false;
    end
end
model.points = points(:, keep);
model.fvalues = model.fvalues(:, keep);

end