function [fvalues, succeeded] = evaluate_new_fvalues(funcs, point)


    functions_num = length(funcs);
    fvalues = zeros(functions_num, 1);
    succeeded = true;
    for nf = 1:functions_num
        try
            fvalues(nf, 1) = ...
                funcs{nf}(point);
        catch this_error
            fvalues(nf, 1) = nan;
        end
        if ~isfinite(fvalues(nf, 1))
            succeeded = false;
            break
        end
    end

end
