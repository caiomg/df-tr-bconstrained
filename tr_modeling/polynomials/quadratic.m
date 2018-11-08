function [fx, dfx, d2fx] = quadratic(H, g, c, x)

    dim = size(x, 1);

    if (size(x, 2) ~= 1 || length(c) ~= 1 || ...
        size(g, 1) ~= dim || size(H, 1) ~= dim || size(H, 2) ~= dim)
        error('cmg:wrongdim:quadr', 'Wrong dimension');
    end

    fx = c + g'*x + 0.5*(x'*H*x);
    dfx = g + H*x;
    d2fx = H;

end