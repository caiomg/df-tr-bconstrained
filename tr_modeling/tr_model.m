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
        function x = center_point(self)
            center_ind = self.tr_center;
            x = self.points_abs(:, center_ind);
        end
        function x = first_point(self)
            x = self.points_abs(:, 1);
        end
        function fval = center_fvalues(self, fval_ind)
            center_ind = self.tr_center;
            if nargin < 2
                try
                fval = self.fvalues(:, center_ind);
                catch this_error
                    rethrow(this_error);
                end
            else
                fval = self.fvalues(fval_ind, center_ind);
            end
        end
        function n = number_of_points(self)
            n = size(self.points_abs, 2);
        end
    end
    
end

