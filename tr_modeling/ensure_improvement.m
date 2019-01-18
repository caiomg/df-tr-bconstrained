function [model, exitflag] = ensure_improvement(model, funcs, bl, bu, options)
% ENSURE_IMPROVEMENT - Improves the model by changing set of points
% Calls the appropriate functions to add more points, exchange
% existing points or rebuilding the model from scratch

    STATUS_POINT_ADDED = 1;
    STATUS_POINT_REPLACED = 2;
    STATUS_OLD_MODEL_REBUILT = 3;
    STATUS_MODEL_REBUILT = 4;    
    

    model_complete = is_complete(model);
    model_fl = is_lambda_poised(model, options);
    model_old = is_old(model, options);
    success = false;
    if ~model_complete && (~model_old || ~model_fl)
        % Calculate a new point to add
        [model, success] = improve_model_nfp(model, funcs, bl, bu, options);
        if success
            exitflag = STATUS_POINT_ADDED;
        end
    elseif model_complete && ~model_old
        % Replace some point with a new one that improves geometry
        [model, success] = choose_and_replace_point(model, funcs, bl, bu, options);
        if success
            exitflag = STATUS_POINT_REPLACED;
        end
    end
    if ~success
        model = rebuild_model(model, options); 
        if model_old
            exitflag = STATUS_OLD_MODEL_REBUILT;
        else
            exitflag = STATUS_MODEL_REBUILT;
        end
    end
end
