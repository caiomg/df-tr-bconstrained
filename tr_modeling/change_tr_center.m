function [model, exitflag] = change_tr_center(model, new_point, ...
                                              new_fvalues, options)

    point_added = false;
    point_exchanged = false;

    if ~is_complete(model)
        % Add this point
        relative_pivot_threshold = options.pivot_threshold;
        [model, point_added] = add_point(model, new_point, new_fvalues, ...
                                         relative_pivot_threshold);
    end
    if point_added
        % Function add_point adds this as the last.
        % Now we have to set as TR center
        model.tr_center = size(model.points_abs, 2); % Last among
                                                     % the points
        exitflag = 1;
    else
        relative_pivot_threshold = options.pivot_threshold;
        [model, point_exchanged, pt_i] = exchange_point(model, ...
                                                        new_point, ...
                                                        new_fvalues, ...
                                                        relative_pivot_threshold);
        if point_exchanged
            model.tr_center = pt_i;
            point_exchanged = true;
            exitflag = 2;
        else
            % add_point and exchange_point failed
            % We still need to add this new point as TR center
            % Model needs rebuilding
            model.points_abs(:, end+1) = new_point;
            model.fvalues(:, end+1) = new_fvalues;
            model.tr_center = size(model.points_abs, 2); % Last
            model = rebuild_model(model, options);
            exitflag = STATUS_MODEL_REBUILT;
        end
    end

end