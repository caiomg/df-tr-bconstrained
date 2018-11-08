function result = is_old(model, options)
%IS_OLD Summary of this function goes here
%   Detailed explanation goes here

    radius_factor = options.radius_factor;
    radius = model.radius;
    distance = norm(model.points_abs(:, 1) - model.points_abs(:, model.tr_center), inf);
    
    result = distance > radius*radius_factor;

end

