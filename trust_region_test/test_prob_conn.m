
% Objective
Hf = 2*diag([1; 1; 1; 2]);
gf = zeros(4, 1);
cf = 0;
f = @(x) quadratic(Hf, gf, cf, x);


all_con = cell(2, 1);

Hc = zeros(4);
Hc(1, 2) = -1;
Hc(2, 1) = -1;
gc = zeros(4, 1);
cc = 4;
gk = @(x) quadratic(Hc, gc, cc, x);
all_con{1} = gk;

Hc = -2*diag([0; 0; 1; 1]);
gc = zeros(4, 1);
cc = 25;
gk = @(x) quadratic(Hc, gc, cc, x);
all_con{2} = gk;

% Solution
x_sol = [2; 2; 5; 0];

% Initial point
x0 = x_sol + 10*ones(size(x_sol));
x0 = x_sol - 6*ones(size(x_sol));
old_seed = rng('default');
x0 = x_sol + randi(20, 4, 1) - 10;
rng(old_seed);
% x0 = [-1; -1; 5; -1];




% Parameters
mu = 100;
epsilon = 2;
delta = 1e-6;
Lambda = 0.1;

nlcon = @(x) constraints(all_con, {}, x, 1);
fmincon_options = optimoptions(@fmincon, 'Display', 'off', ...
                               'SpecifyObjectiveGradient', true);

%%
[x, hs2] = l1_penalty_article(f, all_con, x0, mu, epsilon, delta, Lambda)


tl1 = @() l1_penalty_article(f, all_con, x0, mu, epsilon, delta, Lambda);

time_exact_penalty = timeit(tl1)

fmincon_options = optimoptions(@fmincon, 'Display', 'off', ...
                               'SpecifyObjectiveGradient', true);
x_fmincon = fmincon(f, x0,[],[],[],[],[],[], nlcon, fmincon_options)
tmlab = @() fmincon(f, x0,[],[],[],[],[],[], nlcon, fmincon_options);
time_fmincon = timeit(tmlab)


