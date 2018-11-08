function [point, value] = find_other_point_with_bounds(polynomial, bl, bu)
% FIND_OTHER_POINT tries to find a point in which polynomial
% assumes a sufficiently high value.
%
% Does not perform optimization
% It is assumed one monomial has absolute value of 1
% Tries to find a sufficiently good point by testing possible
% candidates. Based on A.R. Conn's result
    

dimension = polynomial.dimension;
coefficients = polynomial.coefficients;

% Find biggest of monomial coefficients
[~, max_coef] = max(coefficients(2:end));

x0 = zeros(dimension, 1);
x1 = [];
x2 = [];
x3 = [];
x4 = [];

if max_coef <= dimension
    % In case biggest coefficient corresponds to linear monomial
    x1 = x0;
    x1(max_coef) = 1;
    x2 = -x1;
else
    % In case biggest coefficient corresponds to quadratic monomial
    [~, ~, H] = coefficients_to_matrices(dimension, coefficients);
    [H1, pos_12] = max(H);
    [~, pos_2] = max(H1);
    pos_1 = pos_12(pos_2);
    if pos_1 == pos_2
        x1 = x0;
        x1(pos_1) = 1;
        x2 = -x1;
    else
        % If two variables are involved in the offending term
        x1 = x0;
        x1(pos_1) = 1/sqrt(2);
        x1(pos_2) = -1/sqrt(2);
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
    if abs(v) > abs(value)
       value = v;
       point = X(:, k);
    end
end


end