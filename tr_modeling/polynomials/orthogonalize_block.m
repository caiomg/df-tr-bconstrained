function polynomials = orthogonalize_block(polynomials, point, np, ...
                                           orth_beginning, orth_end)
    if nargin < 5
        % All polynomials from this block and higher
        orth_end = length(polynomials);
    end
        
    for p = orth_beginning:orth_end
        if p ~= np
            polynomials(p) = zero_at_point(polynomials(p), ...
                                           polynomials(np), point);
        end
    end
    
end
