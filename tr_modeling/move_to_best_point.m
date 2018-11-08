function [model] = move_to_best_point(model, bl, bu, f)
%MOVE_TO_BEST_POINT Summary of this function goes here
%   Detailed explanation goes here


    points = model.points_abs;
    fvalues = model.fvalues;
    [dim, points_num] = size(points);

    if nargin < 2 || isempty(bl)
        bl = -inf(dim);
    end
    if nargin < 3 || isempty(bu)
        bu = inf(dim);
    end
    if nargin < 4 || isempty(f)
        f = @(v) v(1);
    end    
        
    min_f = inf;
    min_i = 0;
    for k = 1:points_num
        if isempty(find(points(:, k) < bl | points(:, k) > bu, 1))
            val = f(fvalues(:, k));
            if val < min_f
                min_f = val;
                min_i = k;
            end
        end
    end
    if min_i ~= model.tr_center
        model.tr_center = min_i;
    end
    % Here should rebuild polynomials!!!
    

end

