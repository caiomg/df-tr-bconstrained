function coeficientes = matrizes_para_coeficientes(polinomio)
%COEFICIENTES_PARA_MATRIZES transforma uma lista de coeficientes dos
% monomios em matrizes c, g e H, tais que o polinomio possa ser
% representado como p(x) = c + g'x + (1/2)*(x'*H*x)

[dimensao_espaco,  ~] = size(polinomio.g);
numero_termos =  (dimensao_espaco + 1)*(dimensao_espaco + 2)/2;
coeficientes = zeros(numero_termos, 1);

% termo constante
coeficientes(1) = polinomio.c;
% termo ordem um
indice_coeficientes = dimensao_espaco + 1;
coeficientes(2:indice_coeficientes) = polinomio.g;
% termo ordem dois
for k = 1:dimensao_espaco
    for m = 1:k
        indice_coeficientes = indice_coeficientes + 1;
        coeficientes(indice_coeficientes) = polinomio.H(k, m);
        if (polinomio.H(m, k) ~= polinomio.H(k, m))
            warning('H nao simetrica');
        end
    end
end

end

