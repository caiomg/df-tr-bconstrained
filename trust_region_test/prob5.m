
f = @(x) (x(1)-x(2))^2 + (x(2) - x(3))^4;


x0 = [-2.6 ; 2 ; 2];

[x, fval] = trust_region({f}, x0)
