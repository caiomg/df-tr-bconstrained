
dim = 2;
nfp_polynomials = natural_basis(dim);
points = [0, 0;
          1, 0;
          0, 1;
          2, 0;
          1, 1;
          0, 2]';

block_beginning = 2;
for position = 2:6
    if position > 3
        block_beginning = 4;
    end
    
    % Orthogonalize
    nfp_polynomials(position) = ...
        orthogonalize_to_other_polynomials(nfp_polynomials, ...
                                                position, points, ...
                                                position - 1);    
    
    % Normalize polynomial value
    nfp_polynomials(position) = ...
        normalize_polynomial(nfp_polynomials(position), ...
                             points(:, position));

    % Re-orthogonalize
    nfp_polynomials(position) = ...
        orthogonalize_to_other_polynomials(nfp_polynomials, ...
                                                position, points, ...
                                                position - 1);
    hopefully_one = evaluate_polynomial(nfp_polynomials(position), points(:, position));
    if abs(hopefully_one - 1) > 10*eps
        1;
    end
    % Orthogonalize polynomials on present block (deffering
    % subsequent ones)
    nfp_polynomials = orthogonalize_block(nfp_polynomials, ...
                                            points(:, position), ...
                                            position, ...
                                            block_beginning, ...
                                            position - 1);
end

M = interpolation_matrix(nfp_polynomials, points);

% Test the interpolation
v1 = [1 1; 2 3];
d = diag([2, 3]);
H = v1*d*inv(v1);
g = [0; 0];
c = 0;

f = @(x) quadratic(H, g, c, x);
for k = 1:6
    fvals(:, k) = f(points(:, k));
end
l_alpha = nfp_finite_differences(points, fvals, nfp_polynomials)
ssol = M\fvals'
