function best_i = find_best_point(model, bl, bu, f)
%FIND_BEST_POINT Searches for the best among model interpolation points
%   bl, bu (optional) are lower and upper bounds on variables
%   f (optional) is a function for comparison of points. It receives a
%   vector with the function values of each point in the model and returns
%   a corresponding numeric value

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
        % The first value from the vector
        f = @(v) v(1);
    end    
        
    min_f = inf;
    best_i = 0;
    for k = 1:points_num
        if isempty(find(points(:, k) < bl | points(:, k) > bu, 1))
            val = f(fvalues(:, k));
            if val < min_f
                min_f = val;
                best_i = k;
            end
        end
    end
    
end