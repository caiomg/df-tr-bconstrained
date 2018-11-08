classdef tr_model < handle
    
    properties
        points_abs
        fvalues
        points_shifted
        cached_points
        cached_fvalues
        cache_max
        tr_center % Index among points_abs
        radius
        pivot_polynomials % Newton Fundamental Polynomials...
        pivot_values
        modeling_polynomials % Model objective, constraint
    end
    
    methods
        function self = tr_model(points, fvalues, radius)
            self.points_abs = points;
            self.fvalues = fvalues;
            self.radius = radius;
            self.tr_center = 1;
            dim = size(points, 1);
            self.cache_max = 3*dim^2;
        end
    end
    
end

