function basis = band_prioritizing_basis(dimension)


    diag_indices = zeros(dimension, 1);
    ind = 0;
    for k = 1:dimension
        ind = ind + k;
        diag_indices(k) = ind;
    end

    all_hessian_indices = [];
    for k = 1:dimension
        all_hessian_indices = [all_hessian_indices;
                               diag_indices];
        diag_indices = diag_indices(2:end) - 1;
    end

    basis = natural_basis(dimension);

    basis_indices = [(1:dimension+1)';
                     all_hessian_indices + dimension+1];

    basis = basis(basis_indices);
