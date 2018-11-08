
f = @(x) sin(pi*x(1)/12)*cos(pi*x(2)/16);


x0 = [0;
      0];

[x, fval] = trust_region({f}, x0)
