function result = is_old(model, options)
%IS_OLD Summary of this function goes here
%   Detailed explanation goes here

    radius_factor = options.radius_factor;
    radius = model.radius;
    distance = norm(model.first_point() - model.center_point(), inf);
    
    result = distance > radius*radius_factor;

end

