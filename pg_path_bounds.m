function [s, fs, hits] = pg_path_bounds(H, g, x0, s0, bl, bu, radius)

if isempty(bl)
    bl = -inf(size(x));
end
if isempty(bu)
    bu = inf(size(x));
end

dim = size(g, 1);
tol = dim*10*eps;

x = x0;

lower_hits = x < bl;
upper_hits = x > bu;
while true
    not_hits = ~(lower_hits | upper_hits);
    x = min(max(x, bl), bu);
    s = x - x0;
    d = s0(not_hits);
    gs = (g + H*s);
    g_red = gs(not_hits);
    H_red = H(not_hits, not_hits);
    d_phony = s0;
    d_phony(lower_hits | upper_hits) = 0;
    if isempty(d) || norm(d) < tol || norm(s) >= radius || norm(g_red) < tol
        break
    end
    
    lower_breakpoints = (bl - x)./d_phony;
    upper_breakpoints = (bu - x)./d_phony;

    tr_breakpoint = roots([d'*d, 2*s(not_hits)'*d, s'*s - radius^2]);
    if ~isreal(tr_breakpoint)
        tr_breakpoint = inf;
    end

    if -g_red'*d >= 0 && (d'*H_red*d) < 0
        t = inf;
    else
        t = -g_red'*d/(d'*H_red*d);
    end
    if t <= 0 || -g_red'*d < 0
        break
    end
    
    min_lower = min(lower_breakpoints(lower_breakpoints > 0 | (lower_breakpoints == 0 & d_phony < 0))); % Not exactly like this
    min_upper = min(upper_breakpoints(upper_breakpoints > 0 | (upper_breakpoints == 0 & d_phony > 0)));
    tr_breakpoint = min(tr_breakpoint(tr_breakpoint > 0));
    
    t = min([t, min_lower, min_upper, tr_breakpoint]);
    
    x(not_hits) = x(not_hits) + t*d;
    if t == tr_breakpoint
        break
    elseif (isempty(min_lower) || t ~= min_lower) && (isempty(min_upper) || t ~= min_upper)
        % Local minimum
        break
    else
        new_l_hits = (t == lower_breakpoints & d_phony < 0);
        x(new_l_hits) = bl(new_l_hits); % small correction
        new_u_hits = (t == upper_breakpoints & d_phony > 0);
        x(new_u_hits) = bu(new_u_hits); % small correction
        lower_hits = (x == bl); % Not exactly like this
        upper_hits = (x == bu);
    end
end

s = x - x0;
fs = g'*s + 0.5*(s'*H*s);
hits = ~not_hits;
if sum(x > bu | x < bl)
    error('cmg:bounds_violated', 'Bounds not satisfied');
end

end