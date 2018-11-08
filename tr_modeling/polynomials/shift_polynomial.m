function polynomial_shifted = shift_polynomial(polynomial, s)
    
    [c, g, H] = get_matrices(polynomial);

    c_mod = c + g'*s + 0.5*(s'*H*s);
    g_mod = g + H*s;
    
    polynomial_shifted = matrices_to_polynomial(c_mod, g_mod, H);
    
end
