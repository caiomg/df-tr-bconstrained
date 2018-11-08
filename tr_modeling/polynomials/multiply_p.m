function polynomial = multiply_p(polynomial, factor)
% Multiplies a polynomial by a constant factor

polynomial.coefficients  = factor*polynomial.coefficients;

if (max(isnan(polynomial.coefficients)))
    error('cmg:nancoeff', 'NaN coefficient');
end

end

