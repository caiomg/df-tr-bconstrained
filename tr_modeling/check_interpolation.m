function max_diff = check_interpolation(model)

% Tolerance
if model.radius < 1e-3 
    tol_1 = 100*eps;
else
    tol_1 = 10*eps;
end
tol_2 = 10*sqrt(eps);

% Remove shift center from all points
h = model.points_abs;
n_points = size(h, 2);
center_x = model.center_point();
for m = n_points:-1:1
    h(:, m) = h(:, m) - center_x;
end

max_diff = -1;
n_functions = size(model.fvalues, 1);
for k = 1:n_functions
    [c, g, H] = get_model_matrices(model, k-1);
    for m = 1:n_points
        this_value = c + g'*h(:, m) + 0.5*(h(:, m)'*H*h(:, m));
        difference = abs(model.fvalues(k, m) - this_value);
        if difference > max_diff
            max_diff = difference;
        end
        if abs(difference) > max(tol_1*max(abs(model.fvalues(k, :))), tol_2)
            warning('cmg:tr_interpolation_error', 'Interpolation error');
        end
    end
end
