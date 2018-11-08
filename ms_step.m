function [s, fs] = ms_step(H, g, radius)

theta = 0.01;

linopts_l.LT = true;
linopts_u.UT = true;
dim = size(H, 2);

if ~issymmetric(H)
    if norm(H - H') < sqrt(eps)
        H = (H + H')/2;
    end
end
kappa_easy = 0.1;

try
    R = chol(H);
    if min(diag(R)) > dim*sqrt(eps)
        h_pos_def = true;
    else
        h_pos_def = false;
    end
catch exception
   if strcmp(exception.identifier, 'MATLAB:posdef')
       h_pos_def = false;
   else
       rethrow(exception);
   end
end

if h_pos_def
    lambda = 0;
    lambda_l = 0;
else
   % Computing leftmost eigenvalue/vector
%    [~, va] = eigs(H, 1, 'SA');
   [V, Dv] = eig(H);
   [va, ind_va] = min(diag(Dv));
   u = V(:, ind_va);
   % A bit higher
   lambda = abs(va)*1.01 + eps(max(1, norm(H, 1)));
   for k = 1:10
       try
           R = chol(H + lambda*eye(size(H)));
           break
       catch current_exception
           % One more try with other tolerance...
           lambda = lambda*1.05;
       end
   end
   lambda_l = abs(va);
end

y = linsolve(R', -g, linopts_l);
s = linsolve(R, y, linopts_u);

if lambda_l > 0 && norm(s) < radius
    alphas = roots([u'*u, 2*s'*u, s'*s - radius^2]);
    s1 = s + alphas(1)*u;
    s2 = s + alphas(2)*u;
    if 0.5*(s1'*H*s1) + g'*s1 < 0.5*(s2'*H*s2) + g'*s2
        s = s1;
    else
        s = s2;
    end
elseif norm(s) > radius
    lambda_u = norm(g)/radius + norm(H, inf);
    while abs(norm(s) - radius) > kappa_easy*radius
        w = linsolve(R', s, linopts_l);
        lambda = lambda + ((norm(s) - radius)/radius)*((s'*s)/(w'*w));
        if lambda < lambda_l || lambda > lambda_u
           lambda = max(sqrt(lambda_l*lambda_u), lambda_l + theta*(lambda_u - lambda_l)); 
        end
        try
            R = chol(H + lambda*eye(size(H)));
        catch exception
           if strcmp(exception.identifier, 'MATLAB:posdef')
                rethrow(exception);
           else
               rethrow(exception);
           end
        end
        y = linsolve(R', -g, linopts_l);
        s = linsolve(R, y, linopts_u);
        if lambda_u - lambda_l < 10*eps(max(lambda_u, 1))
            break
        end
        if radius - norm(s) > 0
            lambda_u = lambda;
            su = s;
        else
            lambda_l = lambda;
        end
    end
end

fs = g'*s + 0.5*(s'*H*s);

end