function polynomial = scale_polynomial(polynomial, scale_factor)
    

    [c, g, H] = get_matrices(polynomial);
    M = diag(scale_factor);
    g = M*g;
    H = M*H*M;
    
    polynomial = matrices_to_polynomial(c, g, H);
    
end
