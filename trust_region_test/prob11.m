
f = @(x) sum(2*x./(x.*x + 1));

x0 = ones(4, 1);



[x, fval] = trust_region(f, x0)


