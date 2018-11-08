function polynomial = add_p(p1, p2)
% ADD_P adds polynomials: p = p1 + p2

if p1.dimension ~= p2.dimension
    error('cmg:wrongdimension', 'Inconsistent space dimensions');    
end

polynomial.dimension = p1.dimension;
polynomial.coefficients = p1.coefficients + p2.coefficients;

if (max(isnan(polynomial.coefficients)))
    error('cmg:nancoeff', 'NaN coefficient');
end

end

