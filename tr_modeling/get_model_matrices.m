function [c, g, H] = get_model_matrices(model, m)
%GET_MODEL_POLYNOMIAL Summary of this function goes here
%   Detailed explanation goes here

    [c, g, H] = get_matrices(model.modeling_polynomials{m + 1});
    

end

