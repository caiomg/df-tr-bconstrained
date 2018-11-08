function [s, decrease] = c_step(H, g, radius)


cauchy_direction = -radius*(g/norm(g));

if g'*H*g <= 0
    s = 1*cauchy_direction;
else
    s = min(1, norm(g)^3/(radius*(g'*H*g)))*cauchy_direction;
end
% trial decrease is the negative of the m(trial) - m(center)
decrease = -(g'*s + 0.5*(s'*H*s));


end