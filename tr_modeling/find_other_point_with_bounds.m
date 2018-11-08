function [point, value] = find_other_point_with_bounds(polynomial, bl, bu, radius)
% FIND_OTHER_POINT tries to find a point in which polynomial
% assumes a sufficiently high value.
%
% Does not perform optimization
% It is assumed one monomial has absolute value of 1
% Tries to find a sufficiently good point by testing possible
% candidates. Based on A.R. Conn's result
    
if nargin < 4 || isempty(radius)
    radius = 1;
end
dimension = polynomial.dimension;
coefficients = polynomial.coefficients;

% Find biggest of monomial coefficients
[~, max_coef] = max(coefficients(2:end));

x0 = zeros(dimension, 1);
x1 = [];
x2 = [];
x3 = [];
x4 = [];

[~, g, H] = coefficients_to_matrices(dimension, coefficients);
if max_coef <= dimension
    % In case biggest coefficient corresponds to linear monomial
    x1 = x0;
    x1(max_coef) = radius;
    x2 = -x1;
else
    % In case biggest coefficient corresponds to quadratic monomial
    [H1, pos_12] = max(H);
    [~, pos_2] = max(H1);
    pos_1 = pos_12(pos_2);
    if pos_1 == pos_2
        x1 = x0;
        x1(pos_1) = radius;
        x2 = -x1;
        x3 = x0;
        x3(pos_1) = -g(pos_1);
        x4 = -x3;
    else
        % If two variables are involved in the offending term
        x1 = x0;
        x1(pos_1) = radius/sqrt(2);
        x1(pos_2) = -radius/sqrt(2);
        x2 = -x1;
        x3 = abs(x1);
        x4 = -x3;
    end
end

% Points to be tested
X = [x0, x1, x2, x3, x4];
% Find point that produces largest absolute value for polynomial
value = 0;
for k = 1:size(X, 2)
    x = project_to_bounds(X(:, k), bl, bu);
    v = evaluate_polynomial(polynomial, x);
    if v > 0
        Hn = -H;
        gn = -(g + H*x);
    else
        Hn = H;
        gn = (g + H*x);
    end
    rradius = radius - norm(x);
    if rradius > 0
        s0 = ms_step(Hn, gn, rradius);
        s = pg_path_bounds(Hn, gn, x, s0, bl, bu, rradius);
        if norm(s) > rradius
            s = (s/norm(s))*rradius;
        end
    else
        s = zeros(dimension, 1);
    end
    x = x + s;
    if norm(x) > radius
        x = (x/norm(x))*radius;
    end
    v = evaluate_polynomial(polynomial, x);
    if abs(v) >= abs(value)
       value = v;
       point = x;
    end
end


end