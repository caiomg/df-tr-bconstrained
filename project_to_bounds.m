function x = project_to_bounds(x, bl, bu)

    if isempty(bl)
        bl = -inf(size(x));
    end
    if isempty(bu)
        bu = inf(size(x));
    end
    
    x = min(bu, max(bl, x));

end