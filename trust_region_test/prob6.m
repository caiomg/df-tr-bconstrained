
f = @(x) (x(1) + x(2))^2 + (x(2) + x(3))^2;

x0 = [-4;
      1;
      1];
  
  [x, fval] = trust_region({f}, x0)
