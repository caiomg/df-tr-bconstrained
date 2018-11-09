
f = @(x) (x(1)*x(2)*x(3)*x(4))^2;

x0 = [0.8;
      0.8;
      0.8;
      0.8];


[x, fval] = trust_region({f}, x0)
