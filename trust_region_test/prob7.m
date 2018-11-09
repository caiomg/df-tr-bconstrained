
f = @(x) log1p(x(1)^2) + log1p((x(1) - x(2))^2) + log1p((x(2) - x(3))^2) + log1p((x(3) - x(4))^2);

x0 = [2; 
      2; 
      2; 
      2];
  
  [x, fval] = trust_region({f}, x0)
