function x = solve_linear_problem(f, varargin)
    
    x = solve_linear_problem_ipopt(f, 0, varargin{:});

end