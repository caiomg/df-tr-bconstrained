
f = @(x) (1 - x(1))^2;


x0 = [-1.2;
      2];


[x, fval] = trust_region({f}, x0)