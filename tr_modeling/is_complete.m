function result = is_complete(model)
%IS_COMPLETE Summary of this function goes here
%   Detailed explanation goes here

    [dim, points_num] = size(model.points_abs);
    
    max_terms = ((dim + 1)*(dim + 2))/2;
    max_terms_unused = length(model.pivot_polynomials);
    
    result = points_num >= max_terms;
    if points_num > max_terms
        warning('cmg:possible_error', 'Too many points');
    end

end

