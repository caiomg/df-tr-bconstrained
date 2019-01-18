function [model, exitflag] = try_to_add_point(model, new_point, ...
                                              new_fvalues, funcs, ...
                                              bl, bu, options)
% TRY_TO_ADD_POINT - tries to add new known point to model
% If it is not possible to add such point to the interpolation set,
% makes a model improvement

    % Step was not accepted. Smaller effort in including point

    STATUS_POINT_ADDED = 1;

    model_is_complete = is_complete(model);
    if ~model_is_complete
        % Add this point
        relative_pivot_threshold = options.add_threshold;
        [model, point_added] = add_point(model, new_point, new_fvalues, ...
                                         relative_pivot_threshold);
        exitflag = STATUS_POINT_ADDED;
    end
    if model_is_complete || ~point_added
        % Save information about this new point, just not to lose
        model.cached_points = [new_point, model.cached_points];
        model.cached_fvalues = [new_fvalues, model.cached_fvalues];
        % Either add a geometry improving point or rebuild model
        [model, exitflag] = ensure_improvement(model, funcs, bl, bu, options);
    end
end