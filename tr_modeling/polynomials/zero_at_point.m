function p = zero_at_point(p1, p2, x)
% ZERO_AT_POINT - subtract polynomials so that the result is zero at x
%

p = p1;
px = evaluate_polynomial(p, x);
p2x = evaluate_polynomial(p2, x);

iter = 1;
while px ~= 0
    p = add_p(p, multiply_p(p2, -px/p2x));
    px = evaluate_polynomial(p, x);
    if iter >= 2
        break
    end
    iter = iter + 1;
end

end
