function [model, exitflag] = ensure_improvement(model, funcs, bl, bu, options)

    model_complete = is_complete(model);
    model_fl = is_lambda_poised(model, options);
    model_old = is_old(model, options);
    success = false;
    if ~model_complete && (~model_old || ~model_fl)
        % Calculate a new point to add
        [model, success] = improve_model_nfp(model, funcs, bl, bu, options);
        if success
            exitflag = 1;
        end
    elseif model_complete && ~model_old
        % Replace some point with a new one that improves geometry
        [model, success] = choose_and_replace_point(model, funcs, bl, bu, options);
        if success
            exitflag = 2;
        end
    end
    if ~success
        [model, changed] = rebuild_model(model, options); 
        if ~changed
            if ~model_complete
                % Improve model
                [model, success] = improve_model_nfp(model, funcs, bl, bu, options);
            else
                % Replace point
                [model, success] = choose_and_replace_point(model, funcs, bl, bu, options);
            end
        else
            success = true;
        end
        if model_old
            exitflag = 3;
        else
            exitflag = 4;
        end
    end
end
