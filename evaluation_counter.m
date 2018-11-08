classdef evaluation_counter < handle
    %EVALUATION_COUNTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        f
        count
        max_count
    end
    
    methods
        function self = evaluation_counter(f)
            self.f = f;
            self.count = 0;
            self.max_count = inf;
        end
        function reset_count(self)
            self.count = 0;
        end
        function set_max_count(self, n)
            self.max_count = n;
        end
        function n = get_count(self)
           n = self.count;
        end
        function [fx, gx, Hx] = evaluate(self, x)
            if self.count < self.max_count
                self.count = self.count + 1;
                if nargout <= 1
                    fx = self.f(x);
                elseif nargout == 2
                    [fx, gx] = self.f(x);
                elseif nargout == 3
                    [fx, gx, Hx] = self.f(x);
                end
            else
               error('cmg:max_f_count_reached', ...
                     'Exceeded number of function calls');
            end
        end
    end
    
end

