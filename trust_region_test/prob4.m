
f = @(x) 0.01*(x(1) - 1)^2 + (x(2) - x(1)^2)^2;

x0 = [2;
      2;
      2];

[x, fval] = trust_region({f}, x0)

