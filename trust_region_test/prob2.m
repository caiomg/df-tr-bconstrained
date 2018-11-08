
f = @(x) log1p(x(1)^2) + x(2)^2;


x0 = [2;
      2];


[x, fval] = trust_region({f}, x0)